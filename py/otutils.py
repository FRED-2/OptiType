import os
import re

from collections import defaultdict, namedtuple
from functools import lru_cache
from itertools import groupby
from operator import attrgetter
from types import SimpleNamespace

import pandas as pd
import numpy as np
import scipy.sparse as sparse
import pysam

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import Locator, MultipleLocator


C1_TO_TYPE = ('A', 'B', 'C')
C2_TO_TYPE = ('DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1', 'DRB345')
LOCI_TO_TYPE = C1_TO_TYPE + C2_TO_TYPE

def locusgroup(locus):
    if locus in LOCI_TO_TYPE:
        return locus
    if locus in ('DRB3', 'DRB4', 'DRB5'):
        return 'DRB345'
    return 'UNTYPED'

basepath = os.path.dirname(os.path.abspath(__file__))

C1_DF = pd.read_csv(os.path.join(basepath, 'hla_ref', 'dom_info_c1.tsv'), sep='\t', index_col=0)
C2_DF = pd.read_csv(os.path.join(basepath, 'hla_ref', 'dom_info_c2.tsv'), sep='\t', index_col=0)

keyfields = ['ref_l', 'allele', 'allele4', 'p_group', 'redist_f', 'imgt_id', 'cumfreq']
LOCUS_INFO = {}
for locus in LOCI_TO_TYPE:
    basedf = C2_DF if locus in C2_TO_TYPE else C1_DF
    LOCUS_INFO[locus] = basedf.loc[(basedf['locus'] == locus), keyfields].set_index('ref_l')


MSA_SIZES = {}
MSA_LOOKUP = {}
with open(os.path.join(basepath, 'hla_ref', 'dna_dom_delpats.tsv'), 'r') as infile:
    for line in infile:
        seqid, delpat = line.strip().split('\t')
        lookup = [i for i, x in enumerate(delpat) if x=='B']
        MSA_LOOKUP[seqid] = lookup
        # redundant below, it would be enough to do it once per locus/feature, doesn't matter
        locus, fttype, _ = seqid.split('_')
        locus = locusgroup(locus)  # DRB345 / untyped
        MSA_SIZES[(locus, fttype)] = len(delpat)


MSA_REPR = {}
with open(os.path.join(basepath, 'hla_ref', 'dna_dom_delpatref.tsv'), 'r') as infile:
    for line in infile:
        allele, allele_rep = line.strip().split('\t')
        MSA_REPR[allele] = allele_rep


MSA_SEQS_FTID = {}
with open(os.path.join(basepath, 'hla_ref', 'dna_dom.msa'), 'r') as infile:
    for line in infile:
        allele, allele_seq = line.strip().split('\t')
        MSA_SEQS_FTID[allele] = allele_seq


# could be list with append instead of [ref_l] since ref_l's are ordered and 0-based
MSA_SEQS = defaultdict(dict)
for _, row in C1_DF.iterrows():
    MSA_SEQS[row['locus']][row['ref_l']] = MSA_SEQS_FTID[row['genseq']]
for _, row in C2_DF.iterrows():
    MSA_SEQS[row['locus']][row['ref_l']] = MSA_SEQS_FTID[row['e2seq']] + MSA_SEQS_FTID[row['e3seq']]
MSA_SEQS = dict(MSA_SEQS)  # avoid defaultdict's silent missing keys


MSA_FTSIZES = {}
with open(os.path.join(basepath, 'hla_ref', 'dna_dom_ftsizes.txt'), 'r') as infile:
    for line in infile:
        locus, feature, *ftsizes = line.strip().split('\t')
        MSA_FTSIZES[(locus, feature)] = [(ftname, int(ftlen)) for ftname, ftlen in zip(ftsizes[::2], ftsizes[1::2])]


MSA_FTBOUNDS = {}
for locus in C1_TO_TYPE:
    ftsizes = MSA_FTSIZES[(locus, 'gen')]
    ftnames = [x[0] for x in ftsizes]
    ftlens = [x[1] for x in ftsizes]
    ftbounds = list(np.cumsum([0] + ftlens))
    MSA_FTBOUNDS[locus] = (ftnames, ftbounds)


for locus in C2_TO_TYPE:
    ftsizes_e2 = MSA_FTSIZES[(locus, 'e2')]
    ftnames_e2 = [x[0] for x in ftsizes_e2]
    ftlens_e2 = [x[1] for x in ftsizes_e2]
    ftsizes_e3 = MSA_FTSIZES[(locus, 'e3')]
    ftnames_e3 = [x[0] for x in ftsizes_e3]
    ftlens_e3 = [x[1] for x in ftsizes_e3]
    ftnames = ftnames_e2 + ftnames_e3  # it's always 1i, 2e, 2i, 2i, 3e, 3i (duplicate 2i)
    ftbounds = list(np.cumsum([0] + ftlens_e2 + ftlens_e3))
    MSA_FTBOUNDS[locus] = (ftnames, ftbounds)


def msa_coord(seq_id, position):
    if seq_id in MSA_REPR:
        return MSA_LOOKUP[MSA_REPR[seq_id]][position]
    return -1


def get_4digit(imgt_id):
    is_null = not imgt_id[-1].isdigit()
    locus, rest = imgt_id.split('*')
    ret_4digit = locus + '*' + ':'.join(rest.split(':')[:2])
    if is_null and ret_4digit[-1].isdigit():
        return ret_4digit + imgt_id[-1]  # A*12:34:56N -> A*12:34N instead of A*12:34
    return ret_4digit


@lru_cache(maxsize=8192)
def length_on_reference(cigar, fullinfo=False):
    cigar_ops = ((opcode, int(length)) for length, opcode in re.findall(r'(\d+)([MIDNSHP=XB])', cigar))
    cue = 0
    length_on_ref = 0
    matches = []  # (seq_slice_start, seq_slice_end)
    insertions = []  # (to_insert_before, seq_slice_start, seq_slice_end)
    deletions = []  # (deletion_after_this, length)
    for operation, length in cigar_ops:
        if operation == 'M':
            length_on_ref += length
            matches.append((cue, cue+length))
            cue += length
        elif operation == 'D':
            deletions.append((length_on_ref, length))
            length_on_ref += length
        elif operation == 'I':
            insertions.append((cue, cue+length))
            cue += length
    if fullinfo:
        return length_on_ref, matches, insertions, deletions
    return length_on_ref


def indel_cigar(cigar):
    return not bool(re.match(r'^\d+M$', cigar))


Hit = namedtuple('Hit', 'read ref locus feature nm cigar start')  # locus, feature, cigar
LocusHit = namedtuple('LocusHit', 'read locus feature bestnm xloci xnm start end')
ReadAln = namedtuple('ReadAln', 'read readlen nm_p ref_p readhash')


def hit_locus_info(hits, ref_ids):
    # Hits of the same read. Calculate how many loci are hit, what is their best NM,
    # and respective MSA start/end coords (for the typed loci at least)
    # Should it even return anything for untyped loci?
    # hits is a list of lists with (read_i, ref_i, nm, cigar, ref_start) elements
    # to be iterated over by itertools.groupby(first element)

    # relies on msa_coord('A_gen_99', 123) function

    # 1. split them based on locus (derived from ref_i)
    # 2. n_loci = # of groups
    # 3. locus_nm = best NM within locus
    # 4. other_nm = best NM from all other loci
    # 5. msa_start, msa_end using the lowest NM hit of the locus and its CIGAR

    sortedhits = sorted(hits, key=lambda x: (x.locus, x.feature, -x.nm))  # big NM better
    locus_bestnm = {}
    locus_hitnum = {}
    locus_startend = {}
    for locus_ft, hitgroup in groupby(sortedhits, attrgetter('locus', 'feature')):
        locushits = list(hitgroup)
        besthit = locushits[0]
        besthit_id = ref_ids[besthit.ref]
        locus_bestnm[locus_ft] = besthit.nm
        locus_hitnum[locus_ft] = len(locushits)
        locus_start = msa_coord(besthit_id, besthit.start)
        # last INCLUDED. Sparse -1 has to be put 1 higher than this
        locus_end = msa_coord(besthit_id, besthit.start + length_on_reference(besthit.cigar) - 1)
        locus_startend[locus_ft] = (locus_start, locus_end)

    xloci = len(locus_bestnm) - 1
    allnms = sorted(locus_bestnm.values(), reverse=True)  # high NM is good
    bestnm = allnms[0]
    secondbestnm = allnms[1] if xloci else 0  # high NM is good

    loci_hitinfo = []

    for locus_ft, (start, end) in locus_startend.items():
        locusbestnm = locus_bestnm[locus_ft]
        locushit = LocusHit(read=sortedhits[0].read, locus=locus_ft[0], feature=locus_ft[1],
                            bestnm=locusbestnm, xloci=xloci,
                            xnm=(secondbestnm if locusbestnm==bestnm else bestnm),
                            start=start, end=end)
        loci_hitinfo.append(locushit)

    return loci_hitinfo


def nmscore(nm):
    # NM=0 -> 5, 1 -> 4, ... 5+ -> 0
    assert nm >= 0, 'unknown-base-deducted NM is negative, how is this possible?'
    # TODO: negligibly rare, but if there is an N inside an insertion, should it be deducted?
    # Depends whether N's can be gaps instead of unknown bases
    return max(5 - nm, 0)


def bam_to_hdf(bam_filename, min_readlength=30):
    def aln_to_hits(aln, i_read, n_deduction=0):
        # uses refname_to_idx, ref_locus, ref_feature from outer scope
        # deducts the number of N's in read from NM values
        hits = []

        # primary alignment
        hits.append(Hit(read=i_read, ref=aln.rname, locus=ref_locus[aln.rname],
                        feature=ref_feature[aln.rname],
                        nm=nmscore(aln.get_tag('NM') - n_deduction),
                        cigar=aln.cigarstring, start=aln.reference_start))

        if aln.has_tag('XA'):  # secondary alignments
            for subtag in aln.get_tag('XA').rstrip(';').split(';'):  # trailing ; always
                refname, start, cigar, nm = subtag.split(',')
                refidx = refname_to_idx[refname]
                hits.append(Hit(read=i_read, ref=refidx, locus=ref_locus[refidx],
                                feature=ref_feature[refidx], nm=nmscore(int(nm) - n_deduction),
                                cigar=cigar, start=abs(int(start))-1))

        # A read should never map to the same reference more than once, only happens with very short reads
        # like TGTGTGTGTG on the end of DRB exon2 right flank. If it happens, it has to be reported and
        # the alternatives removed (or even the whole read? TODO) The only reason it matters is that the sparse
        # mx constructor will ADD repeat entries instead of overwriting.
        if len(set((x.ref for x in hits))) < len(hits):
            print(f'read mapping to the same ref more than once! rl={aln.query_length} ', end='')
            seenrefs = set()
            hits_uniq_ref = []
            for hit in hits:
                if hit.ref in seenrefs:
                    continue
                seenrefs.add(hit.ref)
                hits_uniq_ref.append(hit)
            hits = hits_uniq_ref

        return hits

    bam = pysam.AlignmentFile(bam_filename, 'rb')
    bam_refs = bam.references
    ref_locus = [locusgroup(x.split('_')[0]) for x in bam_refs]  # DRB345 / untyped
    ref_feature = [x.split('_')[1] for x in bam_refs]
    refname_to_idx = {refname: idx for idx, refname in enumerate(bam_refs)}

    readhash = dict()  # -> [read_idx, read_name, duplicate_readname1, duplicate_readname2...]
    read_sequences = []

    all_hits = []
    read_details = []
    locus_hits = []

    i_seq = -1
    print('aln read ', end='')
    for i_aln, aln in enumerate(bam):
        if not (i_aln % 2000):
            print(f'{i_seq+1}/{i_aln} ...', end='')

        if aln.query_length < min_readlength:
            # with sufficiently short read lengths, reads can map to multiple places on
            # the same reference, and are undesirable otherwise as well. Minimum RL of 30
            # is enough to handle issues like the TGTGTGTG exon2 right flank on DRB1
            continue

        if (seq_hash := hash(aln.query_sequence)) in readhash:
            readhash[seq_hash].append(aln.qname)
            continue

        i_seq += 1
        read_sequences.append(aln.query_sequence)
        readhash[seq_hash] = [i_seq, aln.qname]  # eventually, multiplicity: len(x)-1

        read_hits = aln_to_hits(aln, i_seq, n_deduction=aln.query_sequence.count('N'))
        primary = read_hits[0]
        read_details.append(ReadAln(read=i_seq, readlen=aln.query_length, nm_p=primary.nm,
                                    ref_p=primary.ref, readhash=seq_hash))
        all_hits.extend(read_hits)
        locus_hits.extend(hit_locus_info(read_hits, bam_refs))

    refs_df = pd.DataFrame({'ref_id': bam_refs, 'locus': ref_locus, 'feature': ref_feature}
                           ).rename_axis(index='ref')

    dupcounts = pd.Series([len(x)-1 for x in readhash.values()],
                           index=[x[0] for x in readhash.values()],
                           name='dupcount').rename_axis('read')
    duprefs = {}
    for readidx, *readnames in readhash.values():
        for readname in readnames:
            duprefs[readname] = readidx

    locusreads = pd.DataFrame(locus_hits, columns=LocusHit._fields)
    locusreads2 = {locft: subdf.set_index('read')
                   for locft, subdf in locusreads.groupby(['locus', 'feature'])}

    return SimpleNamespace(
        hits=pd.DataFrame(all_hits, columns=Hit._fields),
        reads=pd.DataFrame(read_details, columns=ReadAln._fields).set_index('read').assign(dupcounts=dupcounts),
        refs=refs_df,
        seqs=pd.Series(read_sequences),
        locusreads=locusreads,
        locusreads2=locusreads2,
        dupreads=pd.Series(duprefs).rename_axis('readname'),
        shape=(len(read_details), len(refs_df))
    )


class MinorSymLogLocator(Locator):
    """
    Dynamically find minor tick positions based on the positions of
    major ticks for a symlog scaling.
    """
    def __init__(self, linthresh):
        """
        Ticks will be placed between the major ticks.
        The placement is linear for x between -linthresh and linthresh,
        otherwise its logarithmically
        """
        self.linthresh = linthresh

    def __call__(self):
        'Return the locations of the ticks'
        majorlocs = self.axis.get_majorticklocs()

        # iterate through minor locs
        minorlocs = []

        # handle the lowest part
        for i in range(1, len(majorlocs)):
            majorstep = majorlocs[i] - majorlocs[i-1]
            if abs(majorlocs[i-1] + majorstep/2) < self.linthresh:
                ndivs = 10
            else:
                ndivs = 9
            minorstep = majorstep / ndivs
            locs = np.arange(majorlocs[i-1], majorlocs[i], minorstep)[1:]
            minorlocs.extend(locs)

        return self.raise_if_exceeds(np.array(minorlocs))

    def tick_values(self, vmin, vmax):
        raise NotImplementedError('Cannot get tick locations for a '
                                  '%s type.' % type(self))


def covplot_pair(layers1, layers2, allelename1, allelename2, title='',
                 ftbounds=None, diffpos=None, figsize=(12,6)):

    @ticker.FuncFormatter
    def absolute_labels(yvalue, pos):
        return '%d' % abs(yvalue)

    # layers1 and 2 should be a matrix consisting of these coverage vector layers as rows:
    # nm0_ambig, nm0_xmap_ambig, nm0_uniq, nm0_xmap_uniq, nm1_uniq, nm1_ambig, nm1_xmap_uniq, nm1_xmap_shared
    layers1 =  np.asarray(np.vstack([0*layers1[0], layers1]).cumsum(axis=0))
    layers2 = -np.asarray(np.vstack([0*layers2[0], layers2]).cumsum(axis=0))

    x_axis = np.arange(layers1.shape[1])
    i_colors = ['#51a1cc', '#ffb66d', '#84cc8f', '#f97f7f',
                '#a5d9ad', '#80bad9', '#faa1a1', '#efc581']

    fig = plt.figure(figsize=figsize, dpi=96, facecolor="white")
    ax = fig.add_subplot(111)
    #ax = plt.axes()

    # manually, as stackplot strangely offsets the whole thing by a few pixels up, breaks symmetry
    for i in range(layers1.shape[0] - 1):
        ax.fill_between(x_axis, layers1[i], layers1[i+1], color=i_colors[i], lw=0, zorder=2)
        ax.fill_between(x_axis, layers2[i], layers2[i+1], color=i_colors[i], lw=0, zorder=2)

    if ftbounds is not None:
        ftnames, bound_coords = ftbounds
        for ftbound in bound_coords:
            ax.axvline(ftbound, linewidth=0.6, color='black', alpha=0.15)
        for ftname, ftbound_start, ftbound_end in zip(ftnames, bound_coords, bound_coords[1:]):
            if not 'exon' in ftname:
                ax.axvspan(ftbound_start, ftbound_end, color='grey', alpha=0.1, lw=0, zorder=1)

    if diffpos is not None:
        ax.plot(diffpos, [0]*len(diffpos), '|', zorder=10, color='black')

    ax.grid(color='grey', which='both', axis='both', linestyle='-',
            linewidth=0.6, alpha=0.1, zorder=1)
    ax.axhline(0, color='black', linewidth=0.8)

    ax.set_yscale('symlog', linthresh=10, linscale=0.5)
    max_ylim = max(layers1.max(), -layers2.min()) * 1.05
    max_ylim = 10**np.ceil(np.log10(max_ylim))  # force plot edge to be a power of 10
    ax.set_ylim(-max_ylim, max_ylim)
    ax.set_xlim(0, x_axis[-1])

    ax.yaxis.set_minor_locator(MinorSymLogLocator(10))
    ax.yaxis.set_major_formatter(absolute_labels)

    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.xaxis.set_minor_locator(MultipleLocator(20))
    ax.tick_params(axis='x', labelrotation=45)

    ax.set_title(title)
    ax.text(0.99, 0.99, allelename1, fontsize=10, horizontalalignment='right',
            verticalalignment='top', transform=plt.gca().transAxes)
    ax.text(0.99, 0.01, allelename2, fontsize=10, horizontalalignment='right',
            verticalalignment='bottom', transform=plt.gca().transAxes)
    plt.close()
    return fig


def hits_to_cov(hits, read_endpoints, read_weights=None):
    if read_weights is None:
        res = hits.dot(read_endpoints).todense().cumsum(axis=1).clip(min=0)
    else:
        res = hits.multiply(read_weights).dot(read_endpoints).todense().cumsum(axis=1).clip(min=0)
    return res


def create_cov_layers(a1, a2, endpoint_mx, is_xmap, read_w=None):
    # 8 coverage layers with: (shared,distinguishing) x (perfect,imperfect) x (local,xmapping)

    # a1 and a2 are hit column vectors matching the read endpoint matrix.
    # Hit weights split such that equally good shared hits are split with half weights
    a1w, a2w = splitweight(a1, a2, ternary=True)

    # np.vstack([cov_u0_b, cov_x0_b, cov_u0_1, cov_x0_1, cov_u1_1, cov_u1_b, cov_x1_1, cov_x1_b])
    # u: non-xmapping, x: xmapping
    # 0: perfect match, 1: mismatched
    # _1: unique to allele, _b: shared between alleles

    is_nm0_a1 = (flatten(a1.todense()) == 1)
    is_nm0_a2 = (flatten(a2.todense()) == 1)
    is_shared = (flatten(a1w.todense()) == 0.5)

    a1classes = (4 * is_nm0_a1 + 2 * is_shared + 1 * is_xmap)
    a2classes = (4 * is_nm0_a2 + 2 * is_shared + 1 * is_xmap)

    a1_layers = []
    a2_layers = []

    if True:  # leaving it in for disabling if we want to present NM>0 reads with (0.25^NM) weight
        a1binary = a1.copy()
        a1binary.data = a1binary.data*0 + 1
        a2binary = a2.copy()
        a2binary.data = a2binary.data*0 + 1
    a1 = a1binary if True else a1
    a2 = a2binary if True else a2

    for i in range(8):
        in_layer_a1 = sparse.csr_matrix(a1classes == i).T
        in_layer_a2 = sparse.csr_matrix(a2classes == i).T

        layer_a1 = (a1.multiply(a1w).multiply(in_layer_a1).T.dot(endpoint_mx)
                    .todense().cumsum(axis=1).clip(min=0))

        layer_a2 = (a2.multiply(a2w).multiply(in_layer_a2).T.dot(endpoint_mx)
                    .todense().cumsum(axis=1).clip(min=0))

        a1_layers.append(layer_a1)
        a2_layers.append(layer_a2)

    layer_reorder = [6,7,4,5,0,2,1,3]  # [u0_shared, x0_shared, u0, x0, u1, u1_shared, x1, x1_shared]
    return np.vstack(a1_layers)[layer_reorder], np.vstack(a2_layers)[layer_reorder]


def diff_coords(seq1, seq2):
    return [i for i, (base1, base2) in enumerate(zip(seq1, seq2)) if base1!=base2]


def hits_to_sparse(hits, shape, nm_xform=lambda x: np.power(4, x)):
    return sparse.csr_matrix((nm_xform(hits['nm']), (hits['read'], hits['ref'])), shape=shape)


def flatten(matrix):
    return np.squeeze(np.asarray(matrix))


def splitweight(a1, a2, ternary=True):
    # a1 and a2 are sparse hit vectors for an allele. NM=0 -> 1, NM=1 -> 1/4 etc.
    # we want a1/(a1+a2) where (a1+a2)=0 should still result in 0. If the NMs are even, they
    # will be split half-half. With 1 difference, it's 80%-20%. With 2, 94%-6%.
    # Ternary forces 0 / 0.5 / 1
    a12 = (a1 + a2)
    a12.eliminate_zeros()
    a12.data = 1/(a12.data)  # works because the matrix is sparse and zeros aren't stored
    a1w = a1.multiply(a12)
    a2w = a2.multiply(a12)
    if ternary:
        a1w.data = np.round(2*a1w.data)/2.0
        a2w.data = np.round(2*a2w.data)/2.0
        a1w.eliminate_zeros()
        a2w.eliminate_zeros()
    return a1w, a2w


def splitweight_all(rest, chosen):
    if isinstance(chosen, int):  # either pass the hit vector, or an index
        chosen = rest[:, chosen]
    # chosen: a single column hit vector
    # rest: the entire hit matrix (alleles in columns, reads in rows) with the chosen included
    chosen_expanded = sparse.hstack([chosen] * rest.shape[1])
    total = rest + chosen_expanded
    total.eliminate_zeros()
    total.data = 1/(total.data)
    rest_w = rest.multiply(total)
    chosen_w = chosen_expanded.multiply(total)
    return rest_w, chosen_w, chosen_expanded


def score_pair(a1, a2, endpoint_mx, read_w=None, as_unexplained=None, clip_unexplained=0):
    # a1 and a2 are hit column vectors matching the read endpoint matrix
    # read_w is read weights
    # as_unexplained is best hit score among all alleles.
    #   if NM0 is best in locus but a1,a2>=NM1
    #   then the unexplained penalty will be 1 - 1/4 (assuming NM exponent of 4 in params)
    #   if it's NM=1 vs NM=2 then only a quarter of that.
    a1w, a2w = splitweight(a1, a2)

    a3w = a1w.copy()  # shared hits
    a3w.data *= (a3w.data == 0.5)
    a3w.eliminate_zeros()

    a1all = a1.multiply(a1w)
    a2all = a2.multiply(a2w)
    a3all = a1.multiply(a3w)

    if read_w is not None:
        a1all = a1all.multiply(read_w)
        a2all = a2all.multiply(read_w)
        a3all = a3all.multiply(read_w)

    res = SimpleNamespace() # score1, score2, score0 (unexplained), score

    res.score1 = np.sqrt(a1all.T.dot(endpoint_mx).todense().cumsum(axis=1).clip(min=0)).sum()
    res.score2 = np.sqrt(a2all.T.dot(endpoint_mx).todense().cumsum(axis=1).clip(min=0)).sum()
    res.score3 = np.sqrt(a3all.T.dot(endpoint_mx).todense().cumsum(axis=1).clip(min=0)).sum()
    res.score0 = 0

    if as_unexplained is not None:
        # as_unexplained is always the maximum of the hitmx, so neither a1 or a2 can be larger.
        # a0 is therefore nonnegative
        a0 = as_unexplained - a1.maximum(a2)
        a0.data = a0.data.clip(min=0)
        a0.eliminate_zeros()

        if read_w is not None:
            a0 = a0.multiply(read_w)

        unexplained_coverage = a0.T.dot(endpoint_mx).todense().cumsum(axis=1)
        if clip_unexplained > 0:
            unexplained_coverage -= clip_unexplained

        res.score0 = -np.sqrt(unexplained_coverage.clip(min=0)).sum()

    res.score = res.score1 + res.score2 + res.score0  # total score: allele1 + allele2 + unexplained penalty
    return res


def create_endpoint_mx(locus_reads, msa_width):
    # sparse matrix with +1's at MSA start positions and -1 at MSA end positions
    # created from paired read information DF
    try:
        reads_idx = locus_reads['read_l'].tolist()
    except:
        assert locus_reads.index.name == 'read_l', 'dataframe has to have read_l column or index'
        reads_idx = locus_reads.index.tolist()

    nreads = len(reads_idx)
    sp = sum(
        sparse.csr_matrix(([startorend] * nreads, (reads_idx, locus_reads[colname])),
                            shape=(nreads, msa_width))
        for startorend, colname in [(1, 'start1'), (-1, 'end1'), (1, 'start2'), (-1, 'end2')]
    )
    return sp



def drb_compatibility(df):
    drb_architecture = {'DRB3': ['03', '11', '12', '13', '14'],
                        'DRB4': ['04', '07', '09'],
                        'DRB5': ['15', '16'],
                        'NULL': ['01', '08', '10']}
    drb_haplo = set((drb345, drb1) for drb345, drb1s in drb_architecture.items() for drb1 in drb1s)
    try:
        b1x = df.loc['DRB1', 'allele_x'][5:7]  # DRB1*XX family
        b1y = df.loc['DRB1', 'allele_y'][5:7]
        b3x = df.loc['DRB345', 'allele_x'][:4]  # DRB[345] gene
        b3y = df.loc['DRB345', 'allele_y'][:4]
        b3x_score, b3y_score = df.loc['DRB345', 'score_x'], df.loc['DRB345', 'score_y']

        compatible = []
        for i, (tx, ty) in enumerate([(b3x, b3y), (b3x, 'NULL'), ('NULL', b3y), ('NULL', 'NULL')]):
            if ((tx, b1x) in drb_haplo and (ty, b1y) in drb_haplo) or ((tx, b1y) in drb_haplo and (ty, b1x) in drb_haplo):
                compatible.append(i)
        if compatible == [0]:
            return 0  # both retained
        if compatible == [3]:
            return 3  # both dropped
        if compatible == [1] or (compatible == [1, 2] and b3x_score >= b3y_score):
            return 2  # allele_y to be dropped
        if compatible == [2] or (compatible == [1, 2] and b3x_score < b3y_score):
            return 1  # allele_x to be dropped
        print('Inconsistent DRB haplotypes, DRB345 unreliable')
        return 4
    except:
        # one of the loci (likely DRB345) missing due to no coverage, nothing do verify
        return 0


def swap_x_y(row):
    # swap first and second alleles in a result Series (only really necessary for DRB345)
    def swapper(label):
        if label.endswith('_x'):
            return label[:-2] + '_y'
        if label.endswith('_y'):
            return label[:-2] + '_x'
        return label

    orig_index = row.index
    row.index = [swapper(label) for label in orig_index]
    return row.reindex(orig_index)


def result_table(df, drop_drb345_y=False):
    cols = ['batch', 'sample',
            'allele_x', 'allele4_x', 'p_group_x',
            'allele_y', 'allele4_y', 'p_group_y',
            'score_x', 'score_y', 'score_p', 'score_s', 'score_f', 'imbalance']
    df2 = df[cols].copy()
    if drop_drb345_y:
        df2.loc['DRB345', df2.columns.str.endswith('_y')] = np.nan
        df2.loc['DRB345', 'imbalance'] = np.nan
    return df2.reset_index().set_index(['batch', 'sample', 'locus'])
