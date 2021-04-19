from datetime import datetime
def now(start=datetime.now()):
    return '>>>' + str(datetime.now()-start)[:-4]

import os
from argparse import ArgumentParser
from types import SimpleNamespace

import pandas as pd
import numpy as np
import scipy.sparse as sparse

from otutils import (bam_to_hdf, hits_to_sparse, flatten, create_endpoint_mx, score_pair,
    drb_compatibility, swap_x_y, diff_coords, result_table, covplot_pair, create_cov_layers,
    LOCI_TO_TYPE, C2_TO_TYPE, C2_DF, LOCUS_INFO, MSA_SEQS, MSA_FTBOUNDS, MSA_SIZES)
    

print(now(), 'imports ready')

parser = ArgumentParser(description="OptiType HLA Ligand atlas version")
parser.add_argument("-b1", dest="bam1")
parser.add_argument("-b2", dest="bam2")
parser.add_argument("-o", dest="outdir")
parser.add_argument("-s", dest="sample_id", type=str, default="")
parser.add_argument("-r", dest="batch_id", type=str, default="")
parser.add_argument("--indiv-plots", dest="indiv", action="store_true")
parser.add_argument("--png-plots", dest="png", action="store_true")
parser.add_argument("--dpi", dest="dpi", type=int, default=144)
parser.add_argument("--pickle-out", dest="pickle", action="store_true")

args = parser.parse_args()

params = SimpleNamespace(
            clip_unexplained = 1.0,
            nm_exponent = 4.0,
            unpaired_weight = 0.25,
            xmap_none = 1.0,
            xmap_sup = 0.8,
            xmap_eq = 0.5,
            bound_leniency = 0.95,
            report_leniency = 0.9,
    )


print(now(), 'start bam loading')
b1 = bam_to_hdf(args.bam1, min_readlength=30)
print(now(), 'bam1 loaded')
b2 = bam_to_hdf(args.bam2, min_readlength=30)
print(now(), 'bam2 loaded')

b1.hitmx = hits_to_sparse(b1.hits, b1.shape, nm_xform=lambda x: x)
b2.hitmx = hits_to_sparse(b2.hits, b2.shape, nm_xform=lambda x: x)


refs_l = {}
for (locus, feature), subdf in b1.refs.groupby(['locus', 'feature']):
    refs_l[(locus, feature)] = subdf.reset_index().rename_axis('ref_l')
    # TODO: For Class II ref_l should correspond to alleles and not feature references.
    # It isn't misused right now, but it's misleading and should be renamed to avoid inconsistency.


# TODO: the stuff below could be wrapped into a pair_ends(b1, b2) function and create a
# new namespace object with: pairing info, new hitmatrices, locus-wise MSA endpoint matrices
# (with e2+e3 concatenated for Class II)

# Hitmx construction for paired-up reads
pairs = (b1.dupreads.to_frame('r1').merge(b2.dupreads.rename('r2'), on='readname', how='inner')
           .reset_index(drop=True).rename_axis('read_p'))  # drop=False would keep read names

p1 = b1.hitmx[pairs['r1']] # hitmx of 1st ends of pairs
p2 = b2.hitmx[pairs['r2']] # hitmx of 2nd ends of pairs
h12 = p1.minimum(p2)

# Hitmx construction for lone-end reads
idx_offset = len(pairs)
single_1 = b1.reads.loc[b1.reads.index.difference(pairs['r1'])].reset_index().rename_axis('read_p')
single_1.index = single_1.index + idx_offset

idx_offset += len(single_1)
single_2 = b2.reads.loc[b2.reads.index.difference(pairs['r2'])].reset_index().rename_axis('read_p')
single_2.index = single_2.index + idx_offset

h1 = b1.hitmx[single_1['read']]  # hitmx of lone 1st ends
h2 = b2.hitmx[single_2['read']]  # hitmx of lone 2nd ends

h = sparse.vstack([h12, h1, h2])  # hitmx of paired, lone 1st, lone 2nd ends

# read index tracker for the individual ends of pairs
read_idmapper = pd.concat([pairs,
                           single_1['read'].to_frame().rename(columns={'read': 'r1'}).assign(r2=-1),
                           single_2['read'].to_frame().rename(columns={'read': 'r2'}).assign(r1=-1)],
                           axis=0)

# Paired-up sparse hitmx slices for every feature. Rows contain all reads,
# whether they hit any reference in the group or not
hitmx_feature_slices = {locft: h[:, subdf['ref']] for locft, subdf in refs_l.items()}

# rows: paired-up read IDs, columns: locft, values: best NM for read for that locft
bestnm_per_feature = pd.DataFrame({locft: flatten(subhits.max(axis=1).todense())
                                    for locft, subhits in hitmx_feature_slices.items()}
                                  ).rename_axis(index='read_p')

all_locfts = bestnm_per_feature.columns

# This block creates the:
# * smaller, feature-wise hit matrices (with only the reads that actually hit them)
# * corresponding read information (best NM for feature, best NM for xmapped features,
#   original r1/r2 ids, MSA start/ends)
# * sparse endpoint matrices (rows: reads, columns: msa positions. +1 is start, -1 is end.
#   Paired reads have two +1s and -1s. TODO: clean middle ones if overlap?)

reads_lp = {}
hits_lp = {}
endpoints_lp = {}

for locft, i_refs in refs_l.items():
    if locft[0] not in LOCI_TO_TYPE:
        continue

    needed_reads = bestnm_per_feature.loc[:, locft] > 0
    needed_refs = i_refs['ref'].tolist()

    # read information columns
    # best NM for this locus/feature
    nm_locus = bestnm_per_feature.loc[needed_reads, locft].rename('bestnm')

    # best NM for any other locus/feature
    nm_xmap = (bestnm_per_feature.loc[needed_reads, all_locfts.difference([locft])]
                .max(axis=1).rename('xmapnm'))

    # reads_locft DataFrame:
    # index: read_p
    # columns: bestnm (for this feature), xmapnm (best nm for any other feature),
    #          r1 (read integer idx in b1), r2 (read integer idx in b2)
    reads_locft = pd.concat([nm_locus, nm_xmap], axis=1).join(read_idmapper)

    # Extend it further with MSA start/end coordinates (start1/end1 for r1 reads,
    # start2/end2 for r2 reads. Zeros if end isn't available)

    # Endpoint information for all reads of locus.
    # Initialize empty DF if there's an absence of hits on either end
    r1_ends = (b1.locusreads2.get(locft, pd.DataFrame(columns=['start', 'end']))[['start', 'end']]
                .rename(columns={'start': 'start1', 'end': 'end1'})
                .reindex(reads_locft['r1'], fill_value=0))

    r2_ends = (b2.locusreads2.get(locft, pd.DataFrame(columns=['start', 'end']))[['start', 'end']]
                .rename(columns={'start': 'start2', 'end': 'end2'})
                .reindex(reads_locft['r2'], fill_value=0))

    r1_ends.index = reads_locft.index
    r2_ends.index = reads_locft.index

    reads_locft = (pd.concat([reads_locft, r1_ends, r2_ends], axis=1)
                    .reset_index().rename_axis(index='read_l'))

    # Finalize three blocks for read pairs: read information, hit matrix, endpoint matrix
    reads_lp[locft] = reads_locft
    hits_lp[locft] = h[needed_reads][:, needed_refs]
    endpoints_lp[locft] = create_endpoint_mx(reads_locft, MSA_SIZES[locft])



# Combine the Class II split e2-e3-wise hitmatrices, endpoint matrices and read info DFs into
# locus-wise hitmatrices, endpoint matrices, and read info DFs.
# C2_DF information DataFrame with reused e2 and e3 reference combinations for alleles
for locus in C2_TO_TYPE:
    alleles = C2_DF.loc[C2_DF['locus'] == locus]

    e2hits = hits_lp[(locus, 'e2')]
    e3hits = hits_lp[(locus, 'e3')]

    e2hits_reidx = e2hits[:, alleles['ref_l2']]
    e3hits_reidx = e3hits[:, alleles['ref_l3']]

    hits_e23 = sparse.vstack([e2hits_reidx, e3hits_reidx])

    e2endpoints = endpoints_lp[(locus, 'e2')]
    e3endpoints = endpoints_lp[(locus, 'e3')]

    endpoints_e23 = sparse.block_diag([e2endpoints, e3endpoints]) # could it fail if e3 was empty?
    reads_e23 = pd.concat([reads_lp[(locus, 'e2')], reads_lp[(locus, 'e3')]], axis=0)

    # new Class II feature 'e23' to behave quasi-identical to Class I gen
    reads_lp[(locus, 'e23')] = reads_e23
    hits_lp[(locus, 'e23')] = hits_e23
    endpoints_lp[(locus, 'e23')] = endpoints_e23

print(now(), 'preprocessing ready, typing starts')


# used for typing:
# reads_lp      dictionary with ('A', 'gen') keys for class1 and ('DRB1', 'e23') keys for class II
# hits_lp       same
# endpoints_lp  same
# locus_info     dictionary with locus as key, containing allele names, frequencies and G/P groups

locus_typinginfo = {}
sample_results = []

for locus in LOCI_TO_TYPE:
    print(now(), 'typing locus', locus)

    locft = (locus, ('e23' if locus in C2_TO_TYPE else 'gen'))

    # shorthands for vectorized ops
    # rr: read info (best NM, end indices, best NM on other loci)
    # ee: sparse read endpoint mx
    # hh: sparse hit matrix
    # ww: read weights (paired-ness, xmapping status)
    # uu: read weight when unexplained
    # cc: coverage vectors
    # vv: coverage scores

    rr = reads_lp[locft]
    ee = endpoints_lp[locft]
    hh = hits_lp[locft].copy()  # NM to weight conversion to be applied directly on the .data attr
    hh.data = np.sign(hh.data) * np.power(params.nm_exponent, -(5-hh.data))  # sign * just in case there are un-eliminated zeros

    # cross-mapping read status
    no_xmap = rr['bestnm'] > (rr['xmapnm'] + 1)  # 2+ better than xmap, if any
    good_xmap = rr['bestnm'] == (rr['xmapnm'] + 1)  # 1 better than xmap
    equal_xmap = rr['bestnm'] == rr['xmapnm']

    is_paired = ~((rr['r1'] == -1) | (rr['r2'] == -1))
    is_paired = (1 * is_paired).clip(lower=params.unpaired_weight)  # paired: 1, unpaired: 0.25

    weight_xmap = (params.xmap_none * no_xmap
                  + params.xmap_sup * good_xmap
                  + params.xmap_eq  * equal_xmap)

    read_weight = weight_xmap * is_paired

    # weights for the penalty layer based on the best possible match on this locus
    as_unexplained_weight = pd.Series(flatten(hh.max(axis=1).todense()), index=read_weight.index)
    # equal xmaps don't count as unexplained
    as_unexplained_weight = as_unexplained_weight * (weight_xmap > params.xmap_eq)

    # read weight vectors for the penalty and regular layers
    uu = sparse.csr_matrix(as_unexplained_weight).T
    ww = sparse.csr_matrix(read_weight).T

    # coverage vector from hitmx, weight vector and endpoint matrix
    cc = hh.multiply(ww).T.dot(ee).todense().cumsum(axis=1)
    cc = pd.DataFrame(cc, index=LOCUS_INFO[locus].index).clip(lower=0)  # index is ref_l
    cc = np.sqrt(cc)
    vv = cc.sum(axis=1).rename('covscore')

    res1 = LOCUS_INFO[locus].join(vv)

    #########################################################################################

    best_single = res1['covscore'].idxmax()  # index is ref_l (0-based in this scope)
    singlescores = res1['covscore'].tolist()

    # maximum possible score of any pair (no shared reads, zero unexplained reads)
    theomax = np.add.outer(singlescores, singlescores)

    # evaluate pairs within this ratio of the momentary lower bound for output with secondary calls
    bound_leniency = params.bound_leniency
    # report evaluated pairs above this range compared to the top result
    report_leniency = params.report_leniency

    # lower_bound = max(singlescores) * np.sqrt(2)  # best HOM solution w/o unexplained penalty
    lower_bound = 0  # cheap enough to start from 0, goes up quick anyway

    candidates = sparse.coo_matrix(theomax * (theomax >= lower_bound * bound_leniency))
    candidate_pairs = list((iix, iiy) for iix, iiy in zip(*candidates.nonzero()) if iix <= iiy)

    # bring the pairs containing the best_single ahead to get a good lower bound quickly
    candidate_pairs = ([iipair for iipair in candidate_pairs if best_single in iipair] +
                       [iipair for iipair in candidate_pairs if best_single not in iipair])

    evaluated = 0
    scored_pairs = []

    for i, (cand1, cand2) in enumerate(candidate_pairs):
        if (singlescores[cand1] + singlescores[cand2]) < lower_bound:
            continue
        evaluated += 1
        sp = score_pair(hh[:, cand1], hh[:, cand2], endpoint_mx=ee, read_w=ww, as_unexplained=uu,
                        clip_unexplained=params.clip_unexplained)

        scored_pairs.append((evaluated, cand1, cand2,
                             sp.score1, sp.score2, sp.score0, sp.score3, sp.score))

        # we don't push the lower bound beyond 95% so decent pairs still make it into the report
        lower_bound = max(lower_bound, sp.score * bound_leniency)

    scored_pairs_df = pd.DataFrame(scored_pairs, columns=['i_try', 'ref_l_x', 'ref_l_y', 'score_x',
                                                        'score_y', 'score_p', 'score_s', 'score'])
    # TODO: ref_l2 name misleading in Class II (also used to refer to the e2 reference elsewhere)

    scored_pairs_df.insert(0, 'batch', args.batch_id)
    scored_pairs_df.insert(1, 'sample', args.sample_id)
    scored_pairs_df.insert(2, 'locus', locus)

    res = (scored_pairs_df
            .merge(LOCUS_INFO[locus], left_on='ref_l_x', right_index=True)
            .merge(LOCUS_INFO[locus], left_on='ref_l_y', right_index=True))

    res = res.loc[res['score'] >= res['score'].max() * report_leniency]  # clean early poor winners

    # Additional descriptors. Imbalance of non-shared coverage
    # (0: evenly split, 1: one allele hogs all unique hits)
    res['imbalance'] = 2 * np.abs((res['score_x'] - res['score_s'])
                         / (res['score_x'] + res['score_y'] - 2*res['score_s']) - 0.5).fillna(0)

    # prior-based tiebreaker (<0.1% bias towards the most frequent alleles)
    res['score_f'] = (res['score_x'] * (1 + 0.001*res['cumfreq_x']) +  # allele x
                      res['score_y'] * (1 + 0.001*res['cumfreq_y']) +  # allele y
                      res['score_p']  # unexplained reads penalty
                      )

    sample_results.append(res)
    locus_typinginfo[locus] = SimpleNamespace(
        locft=locft, res=res, rr=rr, ee=ee, hh=hh, uu=uu, ww=ww, vv=vv,
        no_xmap=no_xmap, good_xmap=good_xmap, equal_xmap=equal_xmap, is_paired=is_paired,
        weight_xmap=weight_xmap, read_weight=read_weight,
        as_unexplained_weight=as_unexplained_weight)

print(now(), 'finished typing')


sample_result_df = pd.concat(sample_results, ignore_index=True)
locus_typinginfo['all'] = sample_result_df

winners = (sample_result_df.loc[sample_result_df.groupby('locus')['score_f'].idxmax()]
            .set_index('locus'))



# DRB345 consistency check and removal
drbx_compat = drb_compatibility(winners)
if drbx_compat == 3:
    winners = winners.drop('DRB345')
elif drbx_compat == 1:
    winners.loc['DRB345'] = swap_x_y(winners.loc['DRB345'])  # swap retained allele into 1st place

drop_drb345_y = drbx_compat in (1, 2)  # second call to be dropped and not displayed on plot/CSV

final_out = result_table(winners, drop_drb345_y)
locus_typinginfo['winners'] = final_out


# Outputs: tsv, pickle, coverage plots
print(now(), 'writing outputs')

os.makedirs(args.outdir, exist_ok=True)
final_out.to_csv(os.path.join(args.outdir, 'result.tsv'), sep='\t', float_format='%.2f')

if args.pickle:
    import pickle
    with open(os.path.join(args.outdir, 'run.pkl'), 'wb') as outfile:
        pickle.dump(locus_typinginfo, outfile)
        print(now(), 'pickled results written')


from matplotlib.backends.backend_pdf import PdfPages
with PdfPages(os.path.join(args.outdir, 'covplots.pdf')) as pdf:

    for locus, row in winners.iterrows():
        l = locus_typinginfo[locus]
        hh, ee, weight_xmap, res = l.hh, l.ee, l.weight_xmap, l.res

        rl1 = row['ref_l_x']
        rl2 = row['ref_l_y']

        allelename1 = LOCUS_INFO[locus].loc[rl1]['allele']
        allelename2 = LOCUS_INFO[locus].loc[rl2]['allele']
        diffpos = diff_coords(MSA_SEQS[locus][rl1], MSA_SEQS[locus][rl2])
        hits1 = hh[:, rl1]
        hits2 = hh[:, rl2]

        if locus == 'DRB345' and drop_drb345_y:
            allelename2 = ''
            diffpos = []
            hits2 *= 0

        layers1, layers2 = create_cov_layers(hits1, hits2, endpoint_mx=ee,
                                            is_xmap=(flatten(weight_xmap) <= 0.7))#params.xmap_eq))
        fig = covplot_pair(layers1, layers2,
                    allelename1=allelename1,
                    allelename2=allelename2,
                    ftbounds=MSA_FTBOUNDS[locus],
                    diffpos=diffpos,
                    title=f"{args.batch_id} {args.sample_id}")
        if args.png:
            fig.savefig(os.path.join(args.outdir, f'covplot_{locus}.png'), dpi=args.dpi)
        if args.indiv:
            fig.savefig(os.path.join(args.outdir, f'covplot_{locus}.pdf'), dpi=args.dpi)
        pdf.savefig(fig, dpi=args.dpi)
    print(now(), 'finished coverage plots')
