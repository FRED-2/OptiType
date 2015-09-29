import pandas as pd
import numpy as np
import re
import pylab
import warnings
from collections import OrderedDict
from datetime import datetime
import sys

try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False

VERBOSE = False

def now(start=datetime.now()):
    # function argument NOT to be ever set! It gets initialized at function definition time
    # so we can calculate difference without further hassle.
    return str(datetime.now()-start)[:-4]


def memoize(f):  # from http://code.activestate.com/recipes/578231-probably-the-fastest-memoization-decorator-in-the-/
    class MemoDict(dict):
        def __missing__(self, key):
            ret = self[key] = f(key)
            return ret
    return MemoDict().__getitem__

# if we used MDINSHP=X chars this could be a full cigar parser, however, we only need M and D to
# determine a read's span on the reference sequence to create a coverage plot.
CIGAR_SLICER = re.compile(r'[0-9]+[MD]')

@memoize
def length_on_reference(cigar_string):
    return sum([int(p[:-1]) for p in re.findall(CIGAR_SLICER, cigar_string)])


def feature_order(feature):
    # to use it in sorted(..., key=feature_order) this returns:
    # 0 for ('UTR', 1)
    # 2 for exon, 1
    # 3 for intron, 1
    # 4 for exon, 2
    # 5 for intron, 2
    # ...
    # 999 for UTR, 2
    feat_type, feat_number = feature
    assert feat_type in ('intron', 'exon', 'UTR'), 'Unknown feature in list (function accepts intron/exon/UTR)'
    assert isinstance(feat_number, int), 'Feature number has to be integer'
    return 0 if feature==('UTR', 1) else 999 if feature==('UTR', 2) else (feat_number*2 + (feat_type=='intron'))


def store_dataframes(out_hdf, **kwargs):
    # DataFrames to serialize have to be passed by keyword arguments. An argument matrix1=DataFrame(...)
    # will be written into table 'matrix1' in the HDF file.

    complevel = kwargs.pop('complevel', 9)   # default complevel & complib values if
    complib = kwargs.pop('complib', 'zlib')  # not explicitly asked for as arguments

    if VERBOSE:
        print now(), 'Storing %d DataFrames in file %s with compression settings %d %s...' % (len(kwargs), out_hdf, complevel, complib)

    store = pd.HDFStore(out_hdf, complevel=complevel, complib=complib)  # TODO: WRITE ONLY? it probably appends now
    for table_name, dataframe in kwargs.iteritems():
        store[table_name] = dataframe
    store.close()

    if VERBOSE:
        print now(), 'DataFrames stored in file.'


def load_hdf(in_hdf, as_dict=False, *args):  # isn't really neccesary, but makes a read-only flag on it to be sure

    store = pd.HDFStore(in_hdf, 'r')
    if len(args):
        if as_dict:
            to_return = {table: store[table] for table in args}
        else:
            to_return = tuple((store[table] for table in args))
        store.close()
        return to_return
    else:
        return store  # store doesn't get closed! Either the user closes it manually or gets closed on exit


def sam_to_hdf(samfile):
    if VERBOSE:
        print now(), 'Loading alleles and read IDs from %s...' % samfile

    # run through the SAM file once to see how many reads and alleles we are dealing with
    # for a one-step DataFrame initialization instead of slow growing

    read_ids, allele_ids = [], []
    first_hit_row = True
    total_hits = 0

    with open(samfile, 'r') as f:
        last_read_id = None
        for line in f:
            if line.startswith('@'):
                if line.startswith('@SQ'):
                    allele_ids.append(line.split('\t')[1][3:]) # SN:HLA:HLA00001
                continue

            total_hits += 1
            read_id = line.split('\t')[0]
            if last_read_id != read_id:
                read_ids.append(read_id)
                last_read_id = read_id

            if first_hit_row:  # analyze SAM file structure and find MD tag column (describes mismatches).
                first_hit_row = False
                columns = line.split()
                try:
                    nm_index = map(lambda x: x.startswith('NM:'), columns).index(True)
                except ValueError:
                    # TODO: we don't really handle the case if NM-tag is not present, code will fail later
                    print '\tNo NM-tag found in SAM file!'
                    nm_index = None

    if VERBOSE:
        print now(), '%d alleles and %d reads found.' % (len(allele_ids), len(read_ids))
        print now(), 'Initializing mapping matrix...'

    # major performance increase if we initialize a numpy zero-matrix and pass that in the constructor
    # than if we just let pandas initialize its default NaN matrix
    matrix_pos = pd.DataFrame(np.zeros((len(read_ids), len(allele_ids)), dtype=np.uint16), columns=allele_ids, index=read_ids)

    # read_details contains NM and read length tuples, calculated from the first encountered CIGAR string.
    read_details = OrderedDict()

    if VERBOSE:
        print now(), '%dx%d mapping matrix initialized. Populating %d hits from SAM file...' % (len(read_ids), len(allele_ids), total_hits)

    milestones = [x * total_hits / 10 for x in range(1, 11)]  # for progress bar

    with open(samfile, 'r') as f:
        counter = 0
        percent = 0
        for line in f:
            if line.startswith('@'):
                continue

            fields = line.strip().split('\t')
            read_id, allele_id, position, cigar, nm = (fields[i] for i in (0, 2, 3, 5, nm_index))

            if read_id not in read_details:
                read_details[read_id] = (int(nm[5:]), length_on_reference(cigar))

            matrix_pos[allele_id][read_id] = int(position)  # SAM indexes from 1, so null elements are not hits

            counter += 1
            if counter in milestones:
                percent += 10
                if VERBOSE:
                    print '\t%d%% completed' % percent
    if VERBOSE:
        print now(), '%d elements filled. Matrix sparsity: 1 in %.2f' % (counter, matrix_pos.shape[0]*matrix_pos.shape[1]/float(counter))

    # convert HLA:HLA00001 identifiers to HLA00001
    matrix_pos.rename(columns=lambda x: x.replace('HLA:', ''), inplace=True)

    details_df = pd.DataFrame.from_dict(read_details, orient='index')
    details_df.columns = ['mismatches', 'read_length']

    return matrix_pos, details_df


def pysam_to_hdf(samfile):
    if not PYSAM_AVAILABLE:
        print "Warning: PySam not available on the system. Falling back to primitive SAM parsing."
        return sam_to_hdf(samfile)

    sam_or_bam = 'rb' if samfile.endswith('.bam') else 'r'
    sam = pysam.AlignmentFile(samfile, sam_or_bam)
    is_yara = (sam.header['PG'][0]['ID'] in ('Yara', 'yara'))
    # If yara produced the sam/bam file, we need to know in what form to expect secondary alignments.
    # If the -os flag was used in the call, they are proper secondary alignments. Otherwise, they are
    # one line per read with a long XA custom tag containing alternative hits.
    xa_tag = is_yara and (' -os ' not in sam.header['PG'][0]['CL'])

    nref = sam.nreferences
    hits = OrderedDict()

    allele_id_to_index = {aa: ii for ii, aa in enumerate(sam.references)}

    # read_details contains NM and read length tuples. We ignore length on reference from now on.
    # Not worth the effort and calculation. Indels are rare and they have to be dealt with differently. The coverage
    # plot wil still be fine and +-1 bp regarding read end isn't concerning.
    read_details = OrderedDict()

    if VERBOSE:
        print now(), 'Loading %s started. Number of HLA reads loaded (updated every thousand):' % samfile

    read_counter = 0
    hit_counter = 0

    for aln in sam:
        # TODO: we could spare on accessing hits[aln.qname] if we were guaranteed alignments of one read come in batches
        # and never mix. Most aligners behave like this but who knows...
        if aln.qname not in hits:  # for primary alignments (first hit of read)
            # not using defaultdict because we have other one-time things to do when new reads come in anyway
            hits[aln.qname] = np.zeros(nref, dtype=np.uint16)  # 16 bits are enough for 65K positions, enough for HLA.
            # TODO: as soon as we don't best-map but let in suboptimal alignments, this below is not good enough,
            # as we need to store NM for every read-allele pair.
            # and then there are cases (usually 1D/1I artefacts at the end) where aligned reference length isn't the same
            # for all alignments of a read. How to handle that? Maybe needless if we re-map reads anyway.
            read_details[aln.qname] = (aln.get_tag('NM'), aln.query_length)  # aln.reference_length it used to be. Soft-trimming is out of question now.
            read_counter += 1
            if VERBOSE and not (read_counter % 1000):
                sys.stdout.write('%dK...' % (len(hits)/1000))
                sys.stdout.flush()
            if xa_tag and aln.has_tag('XA'):
                current_row = hits[aln.qname]  # we may access this hundreds of times, better do it directly
                subtags = aln.get_tag('XA').split(';')[:-1]
                hit_counter += len(subtags)
                for subtag in subtags:
                    allele, pos = subtag.split(',')[:2]  # subtag is like HLA02552,691,792,-,1 (id, start, end, orient, nm)
                    current_row[allele_id_to_index[allele]] = int(pos)  # 1-based positions

        # this runs for primary and secondary alignments as well:
        hits[aln.qname][aln.reference_id] = aln.reference_start + 1  # pysam reports 0-based positions
        hit_counter += 1
        # num_mismatches = aln.get_tag('NM')  # if we ever need suboptimal alignments... doubtful.

    if VERBOSE:
        print '\n', now(), len(hits), 'reads loaded. Creating dataframe...'
    pos_df = pd.DataFrame.from_items(hits.iteritems()).T
    pos_df.columns = sam.references[:]
    details_df = pd.DataFrame.from_dict(read_details, orient='index')
    details_df.columns = ['mismatches', 'read_length']
    if VERBOSE:
        print now(), 'Dataframes created. Shape: %d x %d, hits: %d (%d), sparsity: 1 in %.2f' % (
            pos_df.shape[0], pos_df.shape[1], np.sign(pos_df).sum().sum(), hit_counter, pos_df.shape[0]*pos_df.shape[1]/float(hit_counter)
            )  # TODO: maybe return the binary here right away if we're using it to calculate density anyway.
    return pos_df, details_df


def get_compact_model(hit_df, weak_hit_df=None, weight=None):
# turn a binary hit matrix dataframe into a smaller matrix DF that removes duplicate rows and
# creates the "occurence" vector with the number of rows the representative read represents.
# Note: one can pass "weak" hits (e.g., unpaired reads) and use them with a lower weight.

    hit_df = hit_df.loc[hit_df.any(axis=1)]  # remove all-zero rows
    occurence = {r[0]: len(r) for r in hit_df.groupby(hit_df.columns.tolist()).groups.itervalues()}

    if weak_hit_df is not None:
        weak_hit_df = weak_hit_df.loc[weak_hit_df.any(axis=1)]
        assert 0 < weight <= 1, 'weak hit weight must be in (0, 1]'
        weak_occ = {r[0]: len(r)*weight for r in weak_hit_df.groupby(weak_hit_df.columns.tolist()).groups.itervalues()}
        occurence.update(weak_occ)
        unique_mtx = pd.concat([hit_df.drop_duplicates(), weak_hit_df.drop_duplicates()])
    else:
        unique_mtx = hit_df.drop_duplicates()
    
    return unique_mtx, occurence


def mtx_to_sparse_dict(hit_df):
    # takes a hit matrix and creates a dictionary of (read, allele):1 tuples corresponding to hits
    # (needed by OptiType)
    all_hits = {}
    for read_id, alleles in hit_df.iterrows():
        hit_alleles = alleles[alleles!=0].index
        for hit_allele in hit_alleles:
            all_hits[(read_id, hit_allele)] = 1
    return all_hits
    #return {(read, allele): 1 for read in hit_df.index for allele in hit_df.columns if hit_df[allele][read]>0}  # awfully inefficient


def create_allele_dataframes(imgt_dat, fasta_gen, fasta_nuc):
    from Bio import SeqIO
    if VERBOSE:
        print now(), 'Loading IMGT allele dat file...'

    alleles = OrderedDict()

    with open(imgt_dat, 'rU') as handle:
        for i, record in enumerate(SeqIO.parse(handle, "imgt")):
            # TODO: IMGT has changed the ID system. Now it's HLA00001.1 or HLA12345.2 showing versions. I'll get rid of this now though
            record.id = record.id.split('.')[0]
            alleles[record.id] = record

    if VERBOSE:
        print now(), 'Initializing allele DataFrame...'

    '''
    id      HLA000001
    type    A*01:01:01:01
    locus   A
    class   I --- I, II, TAP, MIC, other
    flags   0 --- problems with allele stored here, see below
    len_gen  3503
    len_nuc  1098
    full_gen  1
    full_nuc  1

    flagging: sum of the below codes
    +1 if HLA type ends in a letter (N, Q, etc.)
    +2 if CDS annotation != exon annotation
    +4 if features don't add up to gen sequence
    +8 if exons don't add up to nuc sequence
    '''

    allele_info = 'id type 4digit locus flags len_dat len_gen len_nuc full_gen full_nuc'

    table = pd.DataFrame(index=alleles.keys(), columns=allele_info.split())
    sequences = []

    if VERBOSE:
        print now(), 'Filling DataFrame with allele data...'

    all_features = []  # contains tuples: (HLA id, feature type, feature number, feature start, feature end)

    for allele in alleles.itervalues():

        allele_type = allele.description.replace('HLA-', '').split(',')[0]

        table.loc[allele.id]['id'] = allele.id
        table.loc[allele.id]['type'] = allele_type
        table.loc[allele.id]['4digit'] = ':'.join(allele_type.split(':')[:2])
        table.loc[allele.id]['locus'] = allele_type.split('*')[0]
        table.loc[allele.id]['flags'] = 0 if allele_type[-1].isdigit() else 1
        table.loc[allele.id]['len_dat'] = len(str(allele.seq))
        table.loc[allele.id]['len_gen'] = 0   # TODO: IT STILL DOESNT SOLVE PERFORMANCEWARNING!
        table.loc[allele.id]['len_nuc'] = 0   # we initialize these nulls so that we don't get a
        table.loc[allele.id]['full_gen'] = 0  # PerformanceWarning + pickling when storing HDF
        table.loc[allele.id]['full_nuc'] = 0  # because of NaNs (they don't map to ctypes)
        sequences.append((allele.id, 'dat', str(allele.seq)))


        # number of features in 2013-04-30 hla.dat:
        # 9296 source (total # of alleles)
        # 9291 CDS, 5 gene (HLA-P pseudogenes)
        # 24493 exons, 3697 introns, 1027 UTRs
        # CDS matches exons in 99+% of cases, some have a single base-pair extension at the end for some
        # weird single-base exons or such (all on unimportant pseudogene loci)
        # so we extract exons, introns and UTRs
        features = [f for f in allele.features if f.type in ('exon', 'intron', 'UTR')]

        for feature in features:
            if feature.type in ('exon', 'intron'):
                feature_num = int(feature.qualifiers['number'][0])
            else:  # UTR
                # UTR either starts at the beginning or ends at the end
                assert feature.location.start == 0 or feature.location.end == len(allele.seq)
                feature_num = 1 if feature.location.start == 0 else 2  # 1 if 5' UTR, 2 if 3' UTR

            all_features.append(
                (allele.id, feature.type, feature_num,
                int(feature.location.start), int(feature.location.end), len(feature), feature_order((feature.type, feature_num)))
                )

        # small sanity check, can be commented out
        cds = [f for f in allele.features if f.type == 'CDS']
        if cds:
            if sum(map(len, [f for f in features if f.type == 'exon'])) != len(cds[0]):
                if VERBOSE:
                    print "\tCDS length doesn't match sum of exons for", allele.id, allele_type
                table.loc[allele.id]['flags'] += 2
        else:
            if VERBOSE:
                print "\tNo CDS found for", allele.id, allele_type
            table.loc[allele.id]['flags'] += 2

    if VERBOSE:
        print now(), 'Loading gen and nuc files...'

    with open(fasta_gen, 'r') as fasta_gen:
        for record in SeqIO.parse(fasta_gen, 'fasta'):
            allele_id = record.id.replace('HLA:', '')
            table.loc[allele_id]['len_gen'] = len(record.seq)
            sequences.append((allele_id, 'gen', str(record.seq)))

    with open(fasta_nuc, 'r') as fasta_nuc:
        for record in SeqIO.parse(fasta_nuc, 'fasta'):
            allele_id = record.id.replace('HLA:', '')
            table.loc[allele_id]['len_nuc'] = len(record.seq)
            sequences.append((allele_id, 'nuc', str(record.seq)))

    # convert list of tuples into DataFrame for features and sequences
    all_features = pd.DataFrame(all_features, columns=['id', 'feature', 'number', 'start', 'end', 'length', 'order'])
    sequences = pd.DataFrame(sequences, columns=['id', 'source', 'sequence'])


    joined = pd.merge(table, all_features, how='inner', on='id')

    exons_for_locus = {}

    for i_locus, i_group in joined.groupby('locus'):
        exons_for_locus[i_locus] = i_group[i_group['feature']=='exon']['number'].max()

    # for i_id, i_group in joined.groupby('id'):
    #     max_exons = exons_for_locus[table.loc[i_id]['locus']]
    #     if len(i_group) >= 2*max_exons-1:
    #         print i_id, 'is fully annotated on all exons and introns and UTRs'

    if VERBOSE:
        print now(), 'Checking dat features vs gen/nuc sequences...'

    for allele, features in joined.groupby('id'):
        row = features.irow(0)  # first row of the features subtable. Contains all allele information because of the join
        sum_features_length = features['length'].sum()
        sum_exons_length = features.loc[features['feature']=='exon']['length'].sum()
        if row['len_gen']>0 and row['len_gen'] != sum_features_length:
            if VERBOSE:
                print "\tFeature lengths don't add up to gen sequence length", allele, row['len_gen'], sum_features_length, row['type']
            table.loc[allele]['flags'] += 4
        if row['len_nuc']>0 and row['len_nuc'] != sum_exons_length:
            if VERBOSE:
                print "\tExon lengths don't add up to nuc sequence length", allele, row['len_nuc'], sum_exons_length, row['type']
            table.loc[allele]['flags'] += 8

    if VERBOSE:
        print now(), 'Sanity check finished. Computing feature sequences...'

    ft_seq_lookup = OrderedDict()
    ft_seq_lookup['---DUMMY---'] = 0  # it will be useful later on if 0 isn't used. lookup_id*boolean operation, etc.
    ft_counter = 1
    all_ft_counter = 0
    all_features['seq_id'] = 0
    for i_id, i_features in all_features.groupby('id'):
        seq = sequences.loc[(sequences['id']==i_id) & (sequences['source']=='dat')].irow(0)['sequence']
        for ft_idx, feature in i_features.iterrows():
            ft_seq = seq[feature['start']:feature['end']]
            all_ft_counter += 1
            if ft_seq not in ft_seq_lookup:
                ft_seq_lookup[ft_seq] = ft_counter
                all_features.loc[ft_idx, 'seq_id'] = ft_counter
                ft_counter += 1
            else:
                all_features.loc[ft_idx, 'seq_id'] = ft_seq_lookup[ft_seq]

    feature_sequences = pd.DataFrame([seq for seq in ft_seq_lookup.keys()], columns=['sequence'])  # , index=ft_seq_lookup.values() but it's 0,1,2,... anyway, as the default

    return table, all_features, sequences, feature_sequences


def prune_identical_alleles(binary_mtx, report_groups=False):
    # return binary_mtx.transpose().drop_duplicates().transpose()
    # # faster:
    hash_columns = binary_mtx.transpose().dot(np.random.rand(binary_mtx.shape[0]))  # dtype np.uint16 okay here because result will be float64
    if report_groups:
        grouper = hash_columns.groupby(hash_columns)
        groups = {g[1].index[0]: g[1].index.tolist() for g in grouper}
    alleles_to_keep = hash_columns.drop_duplicates().index  # relying on keeping the first representative
    # TODO: maybe return sets of alleles that were collapsed into the representative (identical patterned sets)
    return binary_mtx[alleles_to_keep] if not report_groups else (binary_mtx[alleles_to_keep], groups)


def prune_identical_reads(binary_mtx):
    # almost the same as ht.get_compact_model() except it doesn't return an occurence vector.
    # It should only be used to compactify a matrix before passing it to the function finding
    # overshadowed alleles. Final compactifying should be later done on the original binary matrix.
    #
    # return binary_mtx.drop_duplicates()
    # # faster:
    reads_to_keep = binary_mtx.dot(np.random.rand(binary_mtx.shape[1])).drop_duplicates().index  # dtype np.uint16 okay because result will be float64
    return binary_mtx.loc[reads_to_keep]


def prune_overshadowed_alleles(binary_mtx):
    # Calculates B_T*B of the (pruned) binary matrix to determine if certain alleles "overshadow" others, i.e.
    # have the same hits as other alleles plus more. In this case, these "other alleles" can be thrown out early, as
    # they would be never chosen over the one overshadowing them.
    #
    # For a 1000reads x 3600alleles matrix it takes 15 seconds.
    # So you should always prune identical columns and rows before to give it a matrix as small as possible.
    # So ideal usage:
    # prune_overshadowed_alleles(prune_identical_alleles(prune_identical_reads(binary_mtx)))

    # np.dot() is prone to overflow and doesn't expand to int64 like np.sum(). Our binary matrices are np.uint16.
    # If we have less than 65K rows in the binary mtx, we're good with uint16. We're also good if there's no
    # column with 65K+ hits. We check for both with lazy 'or' to avoid calculating the sum if not needed anyway.
    # If we would overflow, we change to uint32.
    if (binary_mtx.shape[0] < np.iinfo(np.uint16).max) or (binary_mtx.sum(axis=0).max() < np.iinfo(np.uint16).max):
        # In case we would reduce the binary mtx to np.uint8 we should to ensure it's at least uint16.
        bb = binary_mtx if all(binary_mtx.dtypes == np.uint16) else binary_mtx.astype(np.uint16)
    else:
        bb = binary_mtx.astype(np.uint32)

    covariance = bb.transpose().dot(bb)
    diagonal = pd.Series([covariance[ii][ii] for ii in covariance.columns], index=covariance.columns)
    new_covariance = covariance[covariance.columns]
    for ii in new_covariance.columns:
        new_covariance[ii][ii] = 0
    overshadowed = []
    for ii in new_covariance.columns:
        potential_superiors = new_covariance[ii][new_covariance[ii]==diagonal[ii]].index
        if any(diagonal[potential_superiors] > diagonal[ii]):
            overshadowed.append(ii)
    non_overshadowed = covariance.columns.difference(overshadowed)
    return non_overshadowed

def create_paired_matrix(binary_1, binary_2, id_cleaning=None):
    # id_cleaning is a function object that turns a read ID that contains pair member information (like XYZ_1234:5678/1) into an
    # ID that identifies the pair itself (like XYZ_1234:5678) so we can match the two DF's IDs. In the above case a suitable
    # id_cleaning function would be lambda x: x[:-2] (gets rid of /1 and /2 from the end). Sometimes it's trickier as pair member
    # information can be somewhere in the middle, and sometimes it's not even required at all as both pair members have the same ID
    # just in two different files (1000 genomes).

    if id_cleaning is not None:
        binary_1.index = map(id_cleaning, binary_1.index)
        binary_2.index = map(id_cleaning, binary_2.index)

    common_read_ids = binary_1.index.intersection(binary_2.index)
    only_1 = binary_1.index.difference(binary_2.index)
    only_2 = binary_2.index.difference(binary_1.index)

    b_1 = binary_1.loc[common_read_ids]
    b_2 = binary_2.loc[common_read_ids]
    b_12 = b_1 * b_2  # elementwise AND
    b_ispaired = b_12.any(axis=1)  # reads with at least one allele w/ paired hits
    b_paired = b_12.loc[b_ispaired]
    b_mispaired = b_1.loc[~b_ispaired] + b_2.loc[~b_ispaired]  # elementwise AND where two ends only hit different alleles
    b_unpaired = pd.concat([binary_1.loc[only_1], binary_2.loc[only_2]])  # concatenation for reads w/ just one end mapping anywhere

    if VERBOSE:
        print now(), ('Alignment pairing completed. %d paired, %d unpaired, %d discordant ' %
            (b_paired.shape[0], b_unpaired.shape[0], b_mispaired.shape[0]))

    return b_paired, b_mispaired, b_unpaired


def get_features(allele_id, features, feature_list):
    # allele_id can be like HLA12345_HLA67890: then we take features from HLA12345 first if possible, then from HLA67890
    # for sequence reconstruction (in that case HLA67890 should be a nearest neighbor of HLA12345)
    if '_' in allele_id:
        partial_allele, complete_allele = allele_id.split('_')
    else:
        complete_allele = allele_id

    feats_complete = {(of['feature'], of['number']): of for _, of in features.loc[features['id']==complete_allele].iterrows()}
    feats_partial = {(of['feature'], of['number']): of for _, of in features.loc[features['id']==partial_allele].iterrows()} if '_' in allele_id else feats_complete

    feats_to_include = []

    for feat in sorted(feature_list, key=feature_order):
        if feat in feats_partial:
            feats_to_include.append(feats_partial[feat])
        elif feat in feats_complete:
            feats_to_include.append(feats_complete[feat])
        else:
            warnings.warn('Feature %s not found for allele %s' % (feat, allele_id))

    return pd.DataFrame(feats_to_include)


def calculate_coverage(alignment, features, alleles_to_plot, features_used):
    assert len(alignment) in (2, 4, 5), ("Alignment tuple either has to have 2, 4 or 5 elements. First four: pos, read_details "
        "once or twice depending on single or paired end, and an optional binary DF at the end for PROPER paired-end plotting")
    has_pairing_info = (len(alignment) == 5)

    if len(alignment) == 2:
        matrix_pos, read_details = alignment[:2]
        pos = matrix_pos[alleles_to_plot]
        pairing = np.sign(pos) * 2  # see explanation on the else branch. These mean unpaired hits.
        hit_counts = np.sign(pos).sum(axis=1)  # how many of the selected alleles does a particular read hit? (unambiguous hit or not)
        # if only a single allele is fed to this function that has no hits at all, we still want to create a null-matrix
        # for it. Its initialization depends on max_ambiguity (it's the size of one of its dimensions) so we set this to 1 at least.
        max_ambiguity = max(1, hit_counts.max())
        to_process = [(pos, read_details, hit_counts, pairing)]
    elif len(alignment) == 4:  # TODO: deprecated. We don't use this, it should probably be taken out.
        matrix_pos1, read_details1, matrix_pos2, read_details2 = alignment
        pos1 = matrix_pos1[alleles_to_plot]
        pos2 = matrix_pos2[alleles_to_plot]
        pairing1 = np.sign(pos1) * 2  # every read is considered unpaired
        pairing2 = np.sign(pos2) * 2
        hit_counts1 = np.sign(pos1).sum(axis=1)
        hit_counts2 = np.sign(pos2).sum(axis=1)
        max_ambiguity = max(1, hit_counts1.max(), hit_counts2.max())
        to_process = [(pos1, read_details1, hit_counts1, pairing1), (pos2, read_details2, hit_counts2, pairing2)]
    else:
        matrix_pos1, read_details1, matrix_pos2, read_details2, pairing_binaries = alignment
        bin_p, bin_u, bin_m = pairing_binaries
        pos1 = matrix_pos1[alleles_to_plot]
        pos2 = matrix_pos2[alleles_to_plot]
        # pairing's values - 0: no hit, 1: paired hit, 2: unpaired hit, 3: discordantly paired hit
        pairing = pd.concat([bin_p[alleles_to_plot], bin_u[alleles_to_plot]*2, bin_m[alleles_to_plot]*3])
        pairing1 = pairing.loc[pos1.index]  # pairing information for first ends' hit matrix
        pairing2 = pairing.loc[pos2.index]  # pairing information for second ends' hit matrix
        hit_counts1 = np.sign(pairing1).sum(axis=1)
        hit_counts2 = np.sign(pairing2).sum(axis=1)
        max_ambiguity = max(1, hit_counts1.max(), hit_counts2.max())
        to_process = [(pos1, read_details1, hit_counts1, pairing1), (pos2, read_details2, hit_counts2, pairing2)]

    coverage_matrices = []

    for allele in alleles_to_plot:
        allele_features = get_features(allele, features, features_used)
        allele_length = allele_features['length'].sum()

        # Dimensions:
        #   1st - 0: perfect hit, 1: hit w/ mismatch(es)
        #   2nd - 0: paired, 1: unpaired, 2: discordantly paired  # maybe introduce unpaired despite being paired on other alleles
        #   3rd - 0: unique hit, 1: ambiguous hit between two alleles, ... n: ambiguous bw/ n+1 alleles
        #   4th - 0 ... [sequence length]
        coverage = np.zeros((2, 3, max_ambiguity, allele_length), dtype=int)

        # to_process is a list containing tuples with mapping/position/hit property dataframes to superimpose.
        # For single-end plotting to_process has just one element. For paired end, to_process has two elements that
        # were pre-processed to only plot reads where both pairs map, etc. but the point is that superimposing the two
        # will provide the correct result.

        for pos, read_details, hit_counts, pairing_info in to_process:
            reads = pos[pos[allele]!=0].index  # reads hitting allele: coverage plot will be built using these
            for i_pos, i_read_length, i_mismatches, i_hitcount, i_pairing in zip(
                    pos.loc[reads][allele],
                    read_details.loc[reads]['read_length'],
                    read_details.loc[reads]['mismatches'],
                    hit_counts[reads],
                    pairing_info.loc[reads][allele]):
                if not i_pairing:
                    continue  # or i_pairing = 4. Happens if one end maps to the allele but a different allele has a full paired hit.
                coverage[bool(i_mismatches)][i_pairing-1][i_hitcount-1][i_pos-1:i_pos-1+i_read_length] += 1

        coverage_matrices.append((allele, coverage))
    return coverage_matrices


def plot_coverage(outfile, coverage_matrices, allele_data, features, features_used, columns=2):

    def start_end_zeros(cov_array):
        # the areaplot polygon gets screwed up if the series don't end/start with zero
        # this adds a zero to the sides to circumvent that
        return np.append(np.append([0], cov_array), [0])

    def allele_sorter(allele_cov_mtx_tuple):
        allele, _ = allele_cov_mtx_tuple
        return allele_data.loc[allele.split('_')[0]]['type']

    def get_allele_locus(allele):
        return allele_data.loc[allele.split('_')[0]]['locus']

    number_of_loci = len(set((get_allele_locus(allele) for allele, _ in coverage_matrices)))

    dpi = 50
    box_size = (7, 1)

    # subplot_rows = len(coverage_matrices)/columns + bool(len(coverage_matrices) % columns)  # instead of float division + rounding up
    # subplot_rows = 3 # TODO
    subplot_rows = 3 * number_of_loci + 1  # rowspan=3 for each plot, legend at the bottom with third height and colspan=2

    area_colors = [  # stacked plot colors
        (0.26, 0.76, 0.26), # perfect, paired, unique
        (0.40, 0.84, 0.40), # perfect, paired, shared

        (0.99, 0.75, 0.20), # perfect, unpaired, unique
        (0.99, 0.75, 0.20), # perfect, mispaired, unique
        (0.99, 0.85, 0.35), # perfect, unpaired, shared
        (0.99, 0.85, 0.35), # perfect, mispaired, shared

        (0.99, 0.23, 0.23), # mismatch, paired, unique
        (0.99, 0.49, 0.49), # mismatch, paired, shared
        
        (0.14, 0.55, 0.72), # mismatch, unpaired, unique
        (0.14, 0.55, 0.72), # mismatch, mispaired, unique
        (0.33, 0.70, 0.88), # mismatch, unpaired, shared
        (0.33, 0.70, 0.88)] # mismatch, mispaired, shared

    figure = pylab.figure(figsize=(box_size[0]*columns, box_size[1]*subplot_rows), dpi=dpi)  # TODO: dpi doesn't seem to do shit. Is it stuck in 100?

    coverage_matrices = sorted(coverage_matrices, key=allele_sorter)  # so that the A alleles come first, then B, and so on.
    prev_locus = ''
    i_locus = -1

    for allele, coverage in coverage_matrices:

        if '_' in allele:
            partial, complete = allele.split('_')
            plot_title = '%s (introns from %s)' % (allele_data.loc[partial]['type'], allele_data.loc[complete]['type']) # , allele for debugging (original ID)
        else:
            plot_title = allele_data.loc[allele]['type'] # + allele for debugging original ID

        if prev_locus != get_allele_locus(allele):  # new locus, start new row
            i_locus += 1
            i_allele_in_locus = 0
        else:
            i_allele_in_locus = 1

        prev_locus = get_allele_locus(allele)

        plot = pylab.subplot2grid((subplot_rows, columns), (3*i_locus, i_allele_in_locus), rowspan=3, adjustable='box')

        _, _, max_ambig, seq_length = coverage.shape  # first two dimensions known (mismatched[2], pairing[3])

        shared_weighting = np.reciprocal(np.arange(max_ambig)+1.0)  # --> 1, 1/2, 1/3...
        shared_weighting[0] = 0  # --> 0, 1/2, 1/3, so the unique part doesn't get mixed in

        perfect_paired_unique = start_end_zeros(coverage[0][0][0])
        mismatch_paired_unique = start_end_zeros(coverage[1][0][0])
        perfect_unpaired_unique = start_end_zeros(coverage[0][1][0])
        mismatch_unpaired_unique = start_end_zeros(coverage[1][1][0])
        perfect_mispaired_unique = start_end_zeros(coverage[0][2][0])
        mismatch_mispaired_unique = start_end_zeros(coverage[1][2][0])

        perfect_paired_shared = start_end_zeros(shared_weighting.dot(coverage[0][0]))
        mismatch_paired_shared = start_end_zeros(shared_weighting.dot(coverage[1][0]))
        perfect_unpaired_shared = start_end_zeros(shared_weighting.dot(coverage[0][1]))
        mismatch_unpaired_shared = start_end_zeros(shared_weighting.dot(coverage[1][1]))
        perfect_mispaired_shared = start_end_zeros(shared_weighting.dot(coverage[0][2]))
        mismatch_mispaired_shared = start_end_zeros(shared_weighting.dot(coverage[1][2]))

        # Exon annotation
        i_start = 1  # position of last feature's end. It's one because we padded with zeros above
        for _, ft in get_features(allele, features, features_used).iterrows():
            if ft['feature'] == 'exon':
                plot.axvspan(i_start, i_start + ft['length'], facecolor='black', alpha=0.1, linewidth=0, zorder=1)
            i_start += ft['length']


        areas = plot.stackplot(np.arange(seq_length+2), # seq_length+2 because of the zero padding at either end
            perfect_paired_unique + 0.001,  # so that 0 isn't -inf on the logplot, but still below cutoff
            perfect_paired_shared,

            perfect_unpaired_unique,
            perfect_mispaired_unique,
            perfect_unpaired_shared,
            perfect_mispaired_shared,

            mismatch_paired_unique,
            mismatch_paired_shared,
            
            mismatch_unpaired_unique,
            mismatch_mispaired_unique,
            mismatch_unpaired_shared,
            mismatch_mispaired_shared,

            linewidth=0, colors=area_colors, zorder=5)

        for aa in areas:
            # if you output to pdf it strangely doesn't respect linewidth=0 perfectly, you'll have
            # to set line colors identical to area color to avoid having a black line between them
            aa.set_edgecolor(aa.get_facecolor())

        plot.tick_params(axis='both', labelsize=10, direction='out', which='both', top=False)

        plot.text(.015, 0.97, plot_title, horizontalalignment='left', verticalalignment='top', transform=plot.transAxes, fontsize=10, zorder=6)
        _, _, _, y2 = plot.axis()
        plot.axis((0, seq_length, 1, y2))  # enforce y axis minimum at 10^0. This corresponds to zero coverage because of the +1 above
        plot.set_yscale('log')
        plot.set_ylim(bottom=0.5)

    legend = pylab.subplot2grid((subplot_rows, columns), (subplot_rows-1, 0), colspan=2, adjustable='box')
    ppp = pylab.matplotlib.patches
    legend.add_patch(ppp.Rectangle((0, 2), 2, 2, color=area_colors[0]))
    legend.add_patch(ppp.Rectangle((0, 0), 2, 2, color=area_colors[1]))
    legend.add_patch(ppp.Rectangle((25, 2), 2, 2, color=area_colors[2]))
    legend.add_patch(ppp.Rectangle((25, 0), 2, 2, color=area_colors[4]))
    legend.add_patch(ppp.Rectangle((50, 2), 2, 2, color=area_colors[6]))
    legend.add_patch(ppp.Rectangle((50, 0), 2, 2, color=area_colors[7]))
    legend.add_patch(ppp.Rectangle((75, 2), 2, 2, color=area_colors[8]))
    legend.add_patch(ppp.Rectangle((75, 0), 2, 2, color=area_colors[10]))
    legend.text( 2.5, 3, 'paired, no mismatches, unique', va='center', size='smaller')
    legend.text( 2.5, 1, 'paired, no mismatches, ambiguous', va='center', size='smaller')
    legend.text(27.5, 3, 'unpaired, no mismatches, unique', va='center', size='smaller')
    legend.text(27.5, 1, 'unpaired, no mismatches, ambiguous', va='center', size='smaller')
    legend.text(52.5, 3, 'paired, mismatched, unique', va='center', size='smaller')
    legend.text(52.5, 1, 'paired, mismatched, ambiguous', va='center', size='smaller')
    legend.text(77.5, 3, 'unpaired, mismatched, unique', va='center', size='smaller')
    legend.text(77.5, 1, 'unpaired, mismatched, ambiguous', va='center', size='smaller')
    legend.set_xlim(0, 100)
    legend.set_ylim(0, 4)
    legend.axison = False

    figure.tight_layout()
    figure.savefig(outfile)
