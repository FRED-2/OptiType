import pandas as pd
import numpy as np
import re
import pylab
import warnings
from collections import OrderedDict
from Bio import SeqIO
from datetime import datetime


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

    if kwargs.get("verbosity", 0):
        print now(), 'Storing %d DataFrames in file %s with compression settings %d %s...' % (len(kwargs), out_hdf, complevel, complib)

    store = pd.HDFStore(out_hdf, complevel=complevel, complib=complib)  # TODO: WRITE ONLY? it probably appends now
    for table_name, dataframe in kwargs.iteritems():
        store[table_name] = dataframe
    store.close()

    if kwargs.get("verbosity", 0):
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


def sam_to_hdf(samfile, **kwargs):
    if kwargs.get("verbosity", 0):
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
                    md_index = map(lambda x: x.startswith('MD:'), columns).index(True)
                except ValueError:
                    # TODO: we don't really handle the case if MD-tag is not present, code will fail later
                    print '\tNo MD-tag found in SAM file!'
                    md_index = None
                try:
                    nm_index = map(lambda x: x.startswith('NM:'), columns).index(True)
                except ValueError:
                    # TODO: we don't really handle the case if NM-tag is not present, code will fail later
                    print '\tNo NM-tag found in SAM file!'
                    nm_index = None

    if kwargs.get("verbosity", 0):
        print now(), '%d alleles and %d reads found.' % (len(allele_ids), len(read_ids))
        print now(), 'Initializing mapping matrix...'

    # major performance increase if we initialize a numpy zero-matrix and pass that in the constructor
    # than if we just let pandas initialize its default NaN matrix
    matrix_pos = pd.DataFrame(np.zeros((len(read_ids), len(allele_ids)), dtype=int), columns=allele_ids, index=read_ids)
    matrix_etc = pd.DataFrame(np.zeros((len(read_ids), len(allele_ids)), dtype=int), columns=allele_ids, index=read_ids)

    if kwargs.get("verbosity", 0):
        print now(), '%dx%d mapping matrix initialized. Populating %d hits from SAM file...' % (len(read_ids), len(allele_ids), total_hits)

    milestones = [x * total_hits / 10 for x in range(1, 11)]  # for progress bar

    with open(samfile, 'r') as f:

        # NM, CIGAR and MD fields are extremely redundant. We'll build an ordered set of their triplets
        # and create a second matrix identical to the first one where fields reference positions in this
        # ordered set of hit descriptors. Reason to use OrderedDict: each key triplet's value is their
        # position index in the set. So if I do a lookup on a triplet, I can tell its index, and later on
        # I can retreive the triplet from its index by nm_cigar_mismatch_set.keys()[index]
        # A minor hack: to be able to initialize a fast numpy zero matrix and not a slow pandas NaN matrix
        # I add a first dummy element to this set (index 0) that correspond to unset matrix fields
        # and the actual triplets will start from 1.
        nm_cigar_mismatch_set = OrderedDict()
        nm_cigar_mismatch_set[(0, '', '')] = 0

        counter = 0
        percent = 0
        for line in f:
            if line.startswith('@'):
                continue

            fields = line.strip().split('\t')
            read_id, allele_id, position, cigar, nm, mismatches = (fields[i] for i in (0, 2, 3, 5, nm_index, md_index))

            position = int(position)
            nm = int(nm[5:])  # NM:i:2 --> 2
            mismatches = mismatches[5:]  # MD:Z:1T32A67 --> 1T32A67

            matrix_pos[allele_id][read_id] = position  # SAM indexes from 1, so null elements are not hits

            if (nm, cigar, mismatches) not in nm_cigar_mismatch_set:
                descriptor_index = len(nm_cigar_mismatch_set)
                nm_cigar_mismatch_set[(nm, cigar, mismatches)] = descriptor_index
            else:
                descriptor_index = nm_cigar_mismatch_set[(nm, cigar, mismatches)]

            matrix_etc[allele_id][read_id] = descriptor_index

            counter += 1
            if counter in milestones:
                percent += 10
                if kwargs.get("verbosity", 0):
                    print '\t%d%% completed' % percent
    if kwargs.get("verbosity", 0):
        print now(), '%d elements filled. Matrix sparsity: 1 in %.2f' % (counter, matrix_pos.shape[0]*matrix_pos.shape[1]/float(counter))

    # convert HLA:HLA00001 identifiers to HLA00001
    matrix_pos.rename(columns=lambda x: x.replace('HLA:', ''), inplace=True)
    matrix_etc.rename(columns=lambda x: x.replace('HLA:', ''), inplace=True)

    # make a DataFrame from descriptor triplet set
    hit_descriptors = pd.DataFrame(nm_cigar_mismatch_set.keys(), columns=['NM', 'CIGAR', 'MD'])

    return matrix_pos, matrix_etc, hit_descriptors


def get_compact_model(hit_df, to_bool=False):
# turn a hit matrix dataframe (can be mapping position matrix if used with to_bool=True)
# into a smaller matrix DF that removes duplicate rows and creates the "occurence" vector
# with the number of rows the representative read represents

    hit_df = hit_df.ix[hit_df.any(axis=1)]  # remove all-zero rows
    if to_bool:
        hit_df = hit_df.applymap(bool)
    occurence = {r[0]: len(r) for r in hit_df.groupby(hit_df.columns.tolist()).groups.itervalues()}
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


def create_allele_dataframes(imgt_dat, fasta_gen, fasta_nuc, **kwargs):
    if kwargs.get("verbosity", 0):
        print now(), 'Loading IMGT allele dat file...'

    alleles = OrderedDict()

    with open(imgt_dat, 'rU') as handle:
        for i, record in enumerate(SeqIO.parse(handle, "imgt")):
            alleles[record.id] = record

    if kwargs.get("verbosity", 0):
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

    if kwargs.get("verbosity", 0):
        print now(), 'Filling DataFrame with allele data...'

    all_features = []  # contains tuples: (HLA id, feature type, feature number, feature start, feature end)

    for allele in alleles.itervalues():

        allele_type = allele.description.replace('HLA-', '').split(',')[0]

        table.ix[allele.id]['id'] = allele.id
        table.ix[allele.id]['type'] = allele_type
        table.ix[allele.id]['4digit'] = ':'.join(allele_type.split(':')[:2])
        table.ix[allele.id]['locus'] = allele_type.split('*')[0]
        table.ix[allele.id]['flags'] = 0 if allele_type[-1].isdigit() else 1
        table.ix[allele.id]['len_dat'] = len(str(allele.seq))
        table.ix[allele.id]['len_gen'] = 0   # TODO: IT STILL DOESNT SOLVE PERFORMANCEWARNING!
        table.ix[allele.id]['len_nuc'] = 0   # we initialize these nulls so that we don't get a
        table.ix[allele.id]['full_gen'] = 0  # PerformanceWarning + pickling when storing HDF
        table.ix[allele.id]['full_nuc'] = 0  # because of NaNs (they don't map to ctypes)
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
                int(feature.location.start), int(feature.location.end), len(feature))
                )

        # small sanity check, can be commented out
        cds = [f for f in allele.features if f.type == 'CDS']
        if cds:
            if sum(map(len, [f for f in features if f.type == 'exon'])) != len(cds[0]):
                print "\tCDS length doesn't match sum of exons for", allele.id, allele_type
                table.ix[allele.id]['flags'] += 2
        else:
            print "\tNo CDS found for", allele.id, allele_type
            table.ix[allele.id]['flags'] += 2

    if kwargs.get("verbosity", 0):
        print now(), 'Loading gen and nuc files...'

    with open(fasta_gen, 'r') as fasta_gen:
        for record in SeqIO.parse(fasta_gen, 'fasta'):
            allele_id = record.id.replace('HLA:', '')
            table.ix[allele_id]['len_gen'] = len(record.seq)
            sequences.append((allele_id, 'gen', str(record.seq)))

    with open(fasta_nuc, 'r') as fasta_nuc:
        for record in SeqIO.parse(fasta_nuc, 'fasta'):
            allele_id = record.id.replace('HLA:', '')
            table.ix[allele_id]['len_nuc'] = len(record.seq)
            sequences.append((allele_id, 'nuc', str(record.seq)))

    # convert list of tuples into DataFrame for features and sequences
    all_features = pd.DataFrame(all_features, columns=['id', 'feature', 'number', 'start', 'end', 'length'])
    sequences = pd.DataFrame(sequences, columns=['id', 'source', 'sequence'])


    joined = pd.merge(table, all_features, how='inner', on='id')

    exons_for_locus = {}

    for i_locus, i_group in joined.groupby('locus'):
        exons_for_locus[i_locus] = i_group[i_group['feature']=='exon']['number'].max()

    # for i_id, i_group in joined.groupby('id'):
    #     max_exons = exons_for_locus[table.ix[i_id]['locus']]
    #     if len(i_group) >= 2*max_exons-1:
    #         print i_id, 'is fully annotated on all exons and introns and UTRs'

    if kwargs.get("verbosity", 0):
        print now(), 'Checking dat features vs gen/nuc sequences...'

    for allele, features in joined.groupby('id'):
        row = features.irow(0)  # first row of the features subtable. Contains all allele information because of the join
        sum_features_length = features['length'].sum()
        sum_exons_length = features.ix[features['feature']=='exon']['length'].sum()
        if row['len_gen']>0 and row['len_gen'] != sum_features_length:
            print "\tFeature lengths don't add up to gen sequence length", allele, row['len_gen'], sum_features_length, row['type']
            table.ix[allele]['flags'] += 4
        if row['len_nuc']>0 and row['len_nuc'] != sum_exons_length:
            print "\tExon lengths don't add up to nuc sequence length", allele, row['len_nuc'], sum_exons_length, row['type']
            table.ix[allele]['flags'] += 8

    if kwargs.get("verbosity", 0):
        print now(), 'Sanity check finished...'

    return table, all_features, sequences


def prune_identical_alleles(binary_mtx, report_groups=False):
    # return binary_mtx.transpose().drop_duplicates().transpose()
    # # faster:
    hash_columns = binary_mtx.transpose().dot(np.random.rand(binary_mtx.shape[0]))
    if report_groups:
        grouper = hash_columns.groupby(hash_columns)
        groups = {g[1].index[0]: g[1].index.tolist() for g in grouper}
    alleles_to_keep = hash_columns.drop_duplicates().index  # relying on keeping the first representative
    # TODO: maybe return sets of alleles that were collapsed into the representative (identical patterned sets)
    return binary_mtx[alleles_to_keep] if not report_groups else (binary_mtx[alleles_to_keep], groups)


def prune_identical_reads(binary_mtx):
    # almost the same as ht.get_compact_model() except it doesn't return an occurence vector.
    # It should only be used to compactify a matrix before passing it to the function finding
    # overshadowed alleles as it uses a super expensive square product operation.
    # Final compactifying should be later done on the original binary matrix.
    #
    # return binary_mtx.drop_duplicates()
    # # faster:
    reads_to_keep = binary_mtx.dot(np.random.rand(binary_mtx.shape[1])).drop_duplicates().index
    return binary_mtx.ix[reads_to_keep]


def prune_overshadowed_alleles(binary_mtx):
    # USES COVARIANCE which can be slow for huge matrices.
    # For a 1000reads x 3600alleles matrix it takes 15 seconds.
    # So you should always prune identical columns and rows before to give it a matrix as small as possible.
    # So ideal usage:
    # prune_overshadowed_alleles(prune_identical_alleles(prune_identical_reads(binary_mtx)))

    bb = binary_mtx.applymap(int)  # make sure dot products are not boolean but numbers
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
    non_overshadowed = covariance.columns.diff(overshadowed)
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

    b_1 = binary_1.ix[common_read_ids]
    b_2 = binary_2.ix[common_read_ids]
    return b_1 * b_2  # elementwise AND


def get_features(allele_id, features, feature_list):
    # allele_id can be like HLA12345_HLA67890: then we take features from HLA12345 first if possible, then from HLA67890
    # for sequence reconstruction (in that case HLA67890 should be a nearest neighbor of HLA12345)
    if '_' in allele_id:
        partial_allele, complete_allele = allele_id.split('_')
    else:
        complete_allele = allele_id

    feats_complete = {(of['feature'], of['number']): of for _, of in features.ix[features['id']==complete_allele].iterrows()}
    feats_partial = {(of['feature'], of['number']): of for _, of in features.ix[features['id']==partial_allele].iterrows()} if '_' in allele_id else feats_complete

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
    assert len(alignment) in (3, 6, 7), ("Alignment tuple either has to have 3, 6 or 7 elements. First six: pos, etc, desc "
        "once or twice depending on single or paired end, and an optional binary DF at the end for PROPER paired-end plotting")

    if len(alignment) == 3:
        matrix_pos, matrix_etc, hit_descriptors = alignment[:3]
        pos = matrix_pos[alleles_to_plot]
        etc = matrix_etc[alleles_to_plot]
        hit = hit_descriptors  # TODO: just want a shorter name
        hit_counts = pos.applymap(bool).sum(axis=1)  # how many of the selected alleles does a particular read hit? (unambiguous hit or not)
        # if only a single allele is fed to this function that has no hits at all, we still want to create a null-matrix
        # for it. Its initialization depends on max_ambiguity (it's the size of one of its dimensions) so we set this to 1 at least.
        max_ambiguity = max(1, hit_counts.max())
        to_process = [(pos, etc, hit, hit_counts)]
    elif len(alignment) == 6:
        matrix_pos1, matrix_etc1, hit_descriptors1, matrix_pos2, matrix_etc2, hit_descriptors2 = alignment
        pos1 = matrix_pos1[alleles_to_plot]
        etc1 = matrix_etc1[alleles_to_plot]
        hit1 = hit_descriptors1
        pos2 = matrix_pos2[alleles_to_plot]
        etc2 = matrix_etc2[alleles_to_plot]
        hit2 = hit_descriptors2
        hit_counts1 = pos1.applymap(bool).sum(axis=1)
        hit_counts2 = pos2.applymap(bool).sum(axis=1)
        max_ambiguity = max(1, hit_counts1.max(), hit_counts2.max())
        to_process = [(pos1, etc1, hit1, hit_counts1), (pos2, etc2, hit2, hit_counts2)]
    else:
        matrix_pos1, matrix_etc1, hit_descriptors1, matrix_pos2, matrix_etc2, hit_descriptors2, binary = alignment
        binx = binary[alleles_to_plot]
        pos1 = matrix_pos1[alleles_to_plot].ix[binx.index] * binx
        pos2 = matrix_pos2[alleles_to_plot].ix[binx.index] * binx
        etc1 = matrix_etc1[alleles_to_plot].ix[binx.index] * binx
        etc2 = matrix_etc2[alleles_to_plot].ix[binx.index] * binx
        hit1 = hit_descriptors1
        hit2 = hit_descriptors2
        hit_counts = binx.sum(axis=1)
        max_ambiguity = max(1, hit_counts.max())
        to_process = [(pos1, etc1, hit1, hit_counts), (pos2, etc2, hit2, hit_counts)]

    coverage_matrices = []

    for allele in alleles_to_plot:
        allele_features = get_features(allele, features, features_used)
        allele_length = allele_features['length'].sum()

        # typically 1-3 rows x sequence_length columns, plus a dimension for perfect hits vs hits with mismatches.
        # Unique hits go in first row, reads hitting 2 alleles in second, etc.
        # Nth row stands for reads hitting N alleles (ie. read might have come from one of N alleles)
        coverage = np.zeros((2, max_ambiguity, allele_length), dtype=int)

        # to_process is a list containing tuples with mapping/position/hit property dataframes to superimpose.
        # For single-end plotting to_process has just one element. For paired end, to_process has two elements that
        # were pre-processed to only plot reads where both pairs map, etc. but the point is that superimposing the two
        # will provide the correct result.

        for pos, etc, hit, hit_counts in to_process:
            reads = pos[pos[allele]!=0].index  # reads hitting allele: coverage plot will be built using these
            for i_pos, i_cigar, i_mismatches, i_hitcount in zip(
                    pos.ix[reads][allele],
                    etc.ix[reads][allele].map(hit['CIGAR']),
                    etc.ix[reads][allele].map(hit['NM']),
                    hit_counts[reads]):
                read_length = length_on_reference(i_cigar)  # this does not neccessarily equal actual read length because of indel mismatches
                coverage[bool(i_mismatches)][i_hitcount-1][i_pos-1:i_pos-1+read_length] += 1

        coverage_matrices.append((allele, coverage))
    return coverage_matrices


def plot_coverage(outfile, coverage_matrices, allele_data, features, features_used, columns=2):

    def start_end_zeros(cov_array):
        # the areaplot polygon gets screwed up if the series don't end/start with zero
        # this adds a zero to the sides to circumvent that
        return np.append(np.append([0], cov_array), [0])

    def allele_sorter(allele_cov_mtx_tuple):
        allele, _ = allele_cov_mtx_tuple
        return allele_data.ix[allele.split('_')[0]]['type']

    def get_allele_locus(allele):
        return allele_data.ix[allele.split('_')[0]]['locus']

    number_of_loci = len(set((get_allele_locus(allele) for allele, _ in coverage_matrices)))

    dpi = 50
    box_size = (7, 3)

    # subplot_rows = len(coverage_matrices)/columns + bool(len(coverage_matrices) % columns)  # instead of float division + rounding up
    # subplot_rows = 3 # TODO
    subplot_rows = number_of_loci

    area_colors = [(.4, .7, .4), (.6, .9, .6), (.9, .3, .3), (.9, .6, .6)]  # colors: unique perfect hits, shared perfect hits, unique mismatch hits, shared mismatch hits
    figure = pylab.figure(figsize=(box_size[0]*columns, box_size[1]*subplot_rows), dpi=dpi)  # TODO: dpi doesn't seem to do shit. Is it stuck in 100?

    coverage_matrices = sorted(coverage_matrices, key=allele_sorter)  # so that the A alleles come first, then B, and so on.
    prev_locus = ''
    i_subplot = 0

    for allele, coverage in coverage_matrices:

        if '_' in allele:
            partial, complete = allele.split('_')
            plot_title = '%s (introns from %s)' % (allele_data.ix[partial]['type'], allele_data.ix[complete]['type']) # , allele for debugging (original ID)
        else:
            plot_title = allele_data.ix[allele]['type'] # + allele for debugging original ID

        if prev_locus != get_allele_locus(allele):  # new locus, start new row
            i_subplot = i_subplot + columns - ((i_subplot-1) % columns)
        else:
            i_subplot += 1

        #print i_subplot

        prev_locus = get_allele_locus(allele)

        plot = figure.add_subplot(subplot_rows, columns, i_subplot, adjustable='box')

        _, max_ambig, seq_length = coverage.shape  # first dimension size is 2 (perfect vs mismatch)

        shared_weighting = np.reciprocal(np.arange(max_ambig)+1.0)  # --> 1, 1/2, 1/3...
        shared_weighting[0] = 0  # --> 0, 1/2, 1/3, so the unique part doesn't get mixed in

        unique_perfect  = start_end_zeros(coverage[0][0])
        unique_mismatch = start_end_zeros(coverage[1][0])
        shared_perfect  = start_end_zeros(shared_weighting.dot(coverage[0]))
        shared_mismatch = start_end_zeros(shared_weighting.dot(coverage[1]))

        # Exon annotation
        i_start = 1  # position of last feature's end. It's one because we padded with zeros above
        for _, ft in get_features(allele, features, features_used).iterrows():
            if ft['feature'] == 'exon':
                plot.axvspan(i_start, i_start + ft['length'], facecolor='blue', alpha=0.1, linewidth=0, zorder=1)
            i_start += ft['length']


        areas = plot.stackplot(np.arange(seq_length+2), unique_perfect+0.001, shared_perfect, unique_mismatch,  # seq_length+2 because of the zero padding at either end
            shared_mismatch, linewidth=0, colors=area_colors, zorder=5)  # +0.001 so zeros are not -inf on the logplot, but somewhere well below the cutoff line. Makes everything nicer

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

    figure.tight_layout()
    figure.savefig(outfile)
