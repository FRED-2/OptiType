"""SAM/BAM file parsing utilities."""

import logging
import re
import sys
from collections import OrderedDict
from functools import cache

import numpy as np
import pandas as pd

from optitype._logging import elapsed

try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False

logger = logging.getLogger(__name__)

# CIGAR pattern for calculating read span on reference
CIGAR_SLICER = re.compile(r"[0-9]+[MD]")


@cache
def _length_on_reference(cigar_string: str) -> int:
    """Calculate the length a read spans on the reference from its CIGAR string."""
    return sum(int(p[:-1]) for p in CIGAR_SLICER.findall(cigar_string))


def sam_to_dataframe(samfile: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Parse a SAM file into position and read detail DataFrames.

    This is a fallback parser when pysam is not available.

    Args:
        samfile: Path to the SAM file.

    Returns:
        Tuple of (matrix_pos, read_details) DataFrames.
    """
    logger.debug("%s Loading alleles and read IDs from %s...", elapsed(), samfile)

    read_ids, allele_ids = [], []
    first_hit_row = True
    total_hits = 0
    nm_index = None

    with open(samfile, "rb") as f:
        last_read_id = None
        for line in f:
            line = line.decode("utf-8")
            if line.startswith("@"):
                if line.startswith("@SQ"):
                    allele_ids.append(line.split("\t")[1][3:])  # SN:HLA:HLA00001
                continue

            total_hits += 1
            read_id = line.split("\t")[0]
            if last_read_id != read_id:
                read_ids.append(read_id)
                last_read_id = read_id

            if first_hit_row:
                first_hit_row = False
                columns = line.split()
                try:
                    nm_index = [x.startswith("NM:") for x in columns].index(True)
                except ValueError:
                    logger.warning("No NM-tag found in SAM file!")
                    nm_index = None

    logger.debug("%s %d alleles and %d reads found.", elapsed(), len(allele_ids), len(read_ids))
    logger.debug("%s Initializing mapping matrix...", elapsed())

    matrix_pos = pd.DataFrame(
        np.zeros((len(read_ids), len(allele_ids)), dtype=np.uint16),
        columns=allele_ids,
        index=read_ids,
    )

    read_details = OrderedDict()

    logger.debug(
        "%s %dx%d mapping matrix initialized. Populating %d hits from SAM file...",
        elapsed(), len(read_ids), len(allele_ids), total_hits,
    )

    milestones = [x * total_hits // 10 for x in range(1, 11)]

    with open(samfile, "rb") as f:
        counter = 0
        percent = 0
        for line in f:
            line = line.decode("utf-8")
            if line.startswith("@"):
                continue

            fields = line.strip().split("\t")
            if int(fields[1]) & 4:  # unmapped
                continue
            read_id, allele_id, position, cigar = fields[0], fields[2], fields[3], fields[5]
            nm = fields[nm_index] if nm_index is not None else "NM:i:0"

            if read_id not in read_details:
                read_details[read_id] = (int(nm[5:]), _length_on_reference(cigar))

            matrix_pos.at[read_id, allele_id] = int(position)

            counter += 1
            if counter in milestones:
                percent += 10
                logger.debug("\t%d%% completed", percent)

    if counter > 0:
        logger.debug(
            "%s %d elements filled. Matrix sparsity: 1 in %.2f",
            elapsed(), counter, matrix_pos.shape[0] * matrix_pos.shape[1] / float(counter),
        )
    else:
        logger.debug("%s No alignments found in SAM file.", elapsed())

    matrix_pos.rename(columns=lambda x: x.replace("HLA:", ""), inplace=True)

    details_df = pd.DataFrame.from_dict(read_details, orient="index")
    details_df.columns = ["mismatches", "read_length"]

    return matrix_pos, details_df


def pysam_to_dataframe(samfile: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Parse a SAM/BAM file into position and read detail DataFrames using pysam.

    Args:
        samfile: Path to the SAM or BAM file.

    Returns:
        Tuple of (matrix_pos, read_details) DataFrames.
    """
    if not PYSAM_AVAILABLE:
        logger.warning("PySam not available. Falling back to primitive SAM parsing.")
        return sam_to_dataframe(samfile)

    sam_or_bam = "rb" if samfile.endswith(".bam") else "r"
    sam = pysam.AlignmentFile(samfile, sam_or_bam)

    # Check if yara produced the file (for XA tag handling)
    is_yara = sam.header.get("PG", [{}])[0].get("ID", "").lower() == "yara"
    xa_tag = is_yara and (" -os " not in sam.header.get("PG", [{}])[0].get("CL", ""))

    nref = sam.nreferences
    hits = OrderedDict()
    allele_id_to_index = {aa: ii for ii, aa in enumerate(sam.references)}
    read_details = OrderedDict()

    logger.debug(
        "%s Loading %s started. Number of HLA reads loaded (updated every thousand):",
        elapsed(), samfile,
    )

    read_counter = 0
    hit_counter = 0

    for aln in sam:
        if aln.is_unmapped:
            continue
        if aln.qname not in hits:
            hits[aln.qname] = np.zeros(nref, dtype=np.uint16)
            nm = aln.get_tag("NM") if aln.has_tag("NM") else 0
            read_details[aln.qname] = (nm, aln.query_length)
            read_counter += 1

            if logger.isEnabledFor(logging.DEBUG) and not (read_counter % 1000):
                sys.stdout.write(f"{len(hits) // 1000}K...")
                sys.stdout.flush()

            if xa_tag and aln.has_tag("XA"):
                current_row = hits[aln.qname]
                subtags = aln.get_tag("XA").split(";")[:-1]
                hit_counter += len(subtags)
                for subtag in subtags:
                    allele, pos = subtag.split(",")[:2]
                    current_row[allele_id_to_index[allele]] = int(pos)

        hits[aln.qname][aln.reference_id] = aln.reference_start + 1
        hit_counter += 1

    logger.debug("\n %s %d reads loaded. Creating dataframe...", elapsed(), len(hits))

    pos_df = pd.DataFrame.from_dict(hits, orient="index", columns=sam.references[:])

    details_df = pd.DataFrame.from_dict(
        read_details,
        orient="index",
        columns=["mismatches", "read_length"],
    )

    if hit_counter > 0:
        logger.debug(
            "%s Dataframes created. Shape: %d x %d, hits: %s (%d), sparsity: 1 in %.2f",
            elapsed(), pos_df.shape[0], pos_df.shape[1],
            np.sign(pos_df).sum().sum(), hit_counter,
            pos_df.shape[0] * pos_df.shape[1] / float(hit_counter),
        )
    else:
        logger.debug("%s No alignments found in %s.", elapsed(), samfile)

    return pos_df, details_df
