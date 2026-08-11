"""
Core pipeline orchestration for OptiType.

This module contains the main HLA typing pipeline logic, coordinating
read mapping, matrix construction, and ILP solving.
"""

import datetime
import logging
import multiprocessing
import os
import shlex
import shutil
import subprocess
import tempfile
import time
from collections import defaultdict
from configparser import RawConfigParser
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from optitype import hlatyper as ht
from optitype._logging import configure_logging, elapsed
from optitype.io.data import get_reference_fasta, load_reference_data
from optitype.io.readers import PYSAM_AVAILABLE, pysam_to_dataframe
from optitype.model import OptiType

logger = logging.getLogger(__name__)


# Frequent alleles list for filtering
FREQ_ALLELES = """A*01:01 A*01:02 A*01:03 A*01:06 A*01:09 A*01:23 A*01:38 A*01:44 A*02:01 A*02:02 A*02:03 A*02:04 A*02:05 A*02:06 A*02:07 A*02:08 A*02:09 A*02:10 A*02:11 A*02:12 A*02:13 A*02:133 A*02:14 A*02:141 A*02:146 A*02:15N A*02:16 A*02:17 A*02:18 A*02:19 A*02:20 A*02:21 A*02:22 A*02:226N A*02:24 A*02:25 A*02:26 A*02:27 A*02:28 A*02:29 A*02:30 A*02:33 A*02:34 A*02:35 A*02:36 A*02:37 A*02:38 A*02:40 A*02:42 A*02:44 A*02:45 A*02:46 A*02:48 A*02:49 A*02:51 A*02:53N A*02:54 A*02:55 A*02:57 A*02:58 A*02:60 A*02:64 A*02:67 A*02:74 A*02:85 A*02:90 A*02:93 A*03:01 A*03:02 A*03:05 A*03:07 A*03:08 A*03:10 A*03:12 A*03:22 A*03:25 A*03:65 A*03:69N A*03:97 A*11:01 A*11:02 A*11:03 A*11:04 A*11:05 A*11:06 A*11:08 A*11:10 A*11:12 A*11:13 A*11:18 A*11:19 A*11:20 A*11:29 A*11:40 A*23:01 A*23:02 A*23:03 A*23:04 A*23:05 A*23:09 A*24:02 A*24:03 A*24:04 A*24:05 A*24:06 A*24:07 A*24:08 A*24:09N A*24:10 A*24:13 A*24:14 A*24:15 A*24:17 A*24:18 A*24:20 A*24:21 A*24:22 A*24:23 A*24:24 A*24:25 A*24:26 A*24:27 A*24:28 A*24:29 A*24:31 A*24:35 A*24:46 A*24:51 A*24:56 A*24:63 A*24:93 A*25:01 A*25:02 A*25:04 A*26:01 A*26:02 A*26:03 A*26:04 A*26:05 A*26:06 A*26:07 A*26:08 A*26:09 A*26:10 A*26:11N A*26:12 A*26:14 A*26:15 A*26:16 A*26:17 A*26:18 A*26:20 A*26:49 A*29:01 A*29:02 A*29:03 A*29:04 A*29:10 A*29:12 A*30:01 A*30:02 A*30:03 A*30:04 A*30:06 A*30:08 A*30:09 A*30:10 A*30:11 A*30:12 A*30:16 A*31:01 A*31:02 A*31:03 A*31:04 A*31:05 A*31:06 A*31:08 A*31:09 A*31:12 A*32:01 A*32:02 A*32:03 A*32:04 A*32:05 A*32:06 A*32:08 A*32:13 A*32:20 A*32:22 A*33:01 A*33:03 A*33:04 A*33:05 A*33:10 A*33:26 A*34:01 A*34:02 A*34:03 A*34:05 A*36:01 A*36:03 A*43:01 A*66:01 A*66:02 A*66:03 A*68:01 A*68:02 A*68:03 A*68:04 A*68:05 A*68:06 A*68:07 A*68:08 A*68:12 A*68:13 A*68:15 A*68:16 A*68:17 A*68:18N A*68:23 A*68:24 A*68:38 A*69:01 A*74:01 A*74:02 A*74:03 A*74:04 A*74:06 A*74:09 A*74:11 A*80:01 A*80:02 B*07:02 B*07:03 B*07:04 B*07:05 B*07:06 B*07:07 B*07:08 B*07:09 B*07:10 B*07:12 B*07:13 B*07:14 B*07:15 B*07:17 B*07:20 B*07:22 B*07:26 B*07:33 B*07:36 B*07:47 B*07:53 B*07:85 B*08:01 B*08:02 B*08:03 B*08:04 B*08:05 B*08:09 B*08:12 B*08:18 B*08:23 B*13:01 B*13:02 B*13:03 B*13:04 B*13:07N B*13:09 B*13:11 B*13:13 B*14:01 B*14:02 B*14:03 B*14:04 B*14:05 B*14:06 B*15:01 B*15:02 B*15:03 B*15:04 B*15:05 B*15:06 B*15:07 B*15:08 B*15:09 B*15:10 B*15:108 B*15:11 B*15:12 B*15:123 B*15:125 B*15:13 B*15:135 B*15:15 B*15:153 B*15:16 B*15:17 B*15:18 B*15:20 B*15:21 B*15:23 B*15:24 B*15:25 B*15:27 B*15:28 B*15:29 B*15:30 B*15:31 B*15:32 B*15:33 B*15:34 B*15:35 B*15:36 B*15:37 B*15:38 B*15:39 B*15:40 B*15:42 B*15:45 B*15:46 B*15:47 B*15:48 B*15:50 B*15:52 B*15:53 B*15:54 B*15:55 B*15:56 B*15:58 B*15:61 B*15:63 B*15:67 B*15:68 B*15:70 B*15:71 B*15:73 B*15:82 B*15:86 B*18:01 B*18:02 B*18:03 B*18:04 B*18:05 B*18:06 B*18:07 B*18:08 B*18:09 B*18:11 B*18:13 B*18:14 B*18:18 B*18:19 B*18:20 B*18:28 B*18:33 B*27:01 B*27:02 B*27:03 B*27:04 B*27:05 B*27:06 B*27:07 B*27:08 B*27:09 B*27:10 B*27:11 B*27:12 B*27:13 B*27:14 B*27:19 B*27:20 B*27:21 B*27:30 B*27:39 B*35:01 B*35:02 B*35:03 B*35:04 B*35:05 B*35:06 B*35:08 B*35:09 B*35:10 B*35:11 B*35:12 B*35:13 B*35:14 B*35:15 B*35:16 B*35:17 B*35:18 B*35:19 B*35:20 B*35:21 B*35:22 B*35:23 B*35:24 B*35:25 B*35:27 B*35:28 B*35:29 B*35:30 B*35:31 B*35:32 B*35:33 B*35:34 B*35:36 B*35:43 B*35:46 B*35:51 B*35:77 B*35:89 B*37:01 B*37:02 B*37:04 B*37:05 B*38:01 B*38:02 B*38:04 B*38:05 B*38:06 B*38:15 B*39:01 B*39:02 B*39:03 B*39:04 B*39:05 B*39:06 B*39:07 B*39:08 B*39:09 B*39:10 B*39:11 B*39:12 B*39:13 B*39:14 B*39:15 B*39:23 B*39:24 B*39:31 B*39:34 B*40:01 B*40:02 B*40:03 B*40:04 B*40:05 B*40:06 B*40:07 B*40:08 B*40:09 B*40:10 B*40:11 B*40:12 B*40:14 B*40:15 B*40:16 B*40:18 B*40:19 B*40:20 B*40:21 B*40:23 B*40:27 B*40:31 B*40:35 B*40:36 B*40:37 B*40:38 B*40:39 B*40:40 B*40:42 B*40:44 B*40:49 B*40:50 B*40:52 B*40:64 B*40:80 B*41:01 B*41:02 B*41:03 B*42:01 B*42:02 B*44:02 B*44:03 B*44:04 B*44:05 B*44:06 B*44:07 B*44:08 B*44:09 B*44:10 B*44:12 B*44:13 B*44:15 B*44:18 B*44:20 B*44:21 B*44:22 B*44:26 B*44:27 B*44:29 B*44:31 B*44:59 B*45:01 B*45:02 B*45:04 B*45:06 B*46:01 B*46:02 B*46:13 B*47:01 B*47:02 B*47:03 B*48:01 B*48:02 B*48:03 B*48:04 B*48:05 B*48:06 B*48:07 B*48:08 B*49:01 B*49:02 B*49:03 B*50:01 B*50:02 B*50:04 B*50:05 B*51:01 B*51:02 B*51:03 B*51:04 B*51:05 B*51:06 B*51:07 B*51:08 B*51:09 B*51:10 B*51:12 B*51:13 B*51:14 B*51:15 B*51:18 B*51:21 B*51:22 B*51:27N B*51:29 B*51:31 B*51:32 B*51:33 B*51:34 B*51:36 B*51:37 B*51:63 B*51:65 B*52:01 B*52:02 B*52:06 B*53:01 B*53:02 B*53:03 B*53:04 B*53:05 B*53:07 B*53:08 B*54:01 B*54:02 B*55:01 B*55:02 B*55:03 B*55:04 B*55:07 B*55:08 B*55:10 B*55:12 B*55:16 B*55:46 B*56:01 B*56:02 B*56:03 B*56:04 B*56:05 B*56:06 B*56:07 B*56:09 B*56:11 B*57:01 B*57:02 B*57:03 B*57:04 B*57:05 B*57:06 B*57:10 B*58:01 B*58:02 B*58:06 B*59:01 B*67:01 B*67:02 B*73:01 B*78:01 B*78:02 B*78:03 B*78:05 B*81:01 B*81:02 B*82:01 B*82:02 B*83:01 C*01:02 C*01:03 C*01:04 C*01:05 C*01:06 C*01:08 C*01:14 C*01:17 C*01:30 C*01:32 C*02:02 C*02:03 C*02:04 C*02:06 C*02:08 C*02:09 C*02:10 C*02:19 C*02:20 C*02:27 C*03:02 C*03:03 C*03:04 C*03:05 C*03:06 C*03:07 C*03:08 C*03:09 C*03:10 C*03:13 C*03:14 C*03:15 C*03:16 C*03:17 C*03:19 C*03:21 C*03:32 C*03:42 C*03:43 C*03:56 C*03:67 C*03:81 C*04:01 C*04:03 C*04:04 C*04:05 C*04:06 C*04:07 C*04:10 C*04:11 C*04:14 C*04:15 C*04:24 C*04:29 C*04:33 C*04:37 C*05:01 C*05:03 C*05:04 C*05:07N C*05:09 C*06:02 C*06:03 C*06:04 C*06:06 C*06:08 C*06:09 C*06:17 C*06:24 C*06:53 C*07:01 C*07:02 C*07:03 C*07:04 C*07:05 C*07:06 C*07:07 C*07:08 C*07:09 C*07:10 C*07:109 C*07:123 C*07:13 C*07:14 C*07:17 C*07:18 C*07:19 C*07:21 C*07:22 C*07:27 C*07:29 C*07:32N C*07:37 C*07:43 C*07:46 C*07:49 C*07:56 C*07:66 C*07:67 C*07:68 C*07:80 C*07:95 C*08:01 C*08:02 C*08:03 C*08:04 C*08:05 C*08:06 C*08:13 C*08:15 C*08:20 C*08:21 C*08:27 C*12:02 C*12:03 C*12:04 C*12:05 C*12:07 C*12:12 C*12:15 C*12:16 C*14:02 C*14:03 C*14:04 C*15:02 C*15:03 C*15:04 C*15:05 C*15:06 C*15:07 C*15:08 C*15:09 C*15:11 C*15:12 C*15:13 C*15:17 C*16:01 C*16:02 C*16:04 C*16:08 C*16:09 C*17:01 C*17:02 C*17:03 C*17:04 C*18:01 C*18:02 A*30:07 B*15:64 B*18:12""".split()


# Default mapper arguments
# RazerS3: -i identity%, -m max_matches, --distance-range, -pa purge_ambiguous
RAZERS3_DEFAULT_ARGS = ["-i", "97", "-m", "99999", "--distance-range", "0", "-pa"]
# YARA (nf-core/hlatyping params): -e error_rate%, -f output_format
YARA_DEFAULT_ARGS = ["-e", "3", "-f", "bam"]


@dataclass
class PipelineConfig:
    """Configuration for the OptiType pipeline."""

    razers3_path: str = "razers3"
    mapper: str = "razers3"
    razers3_args: list[str] = field(default_factory=lambda: RAZERS3_DEFAULT_ARGS[:])
    yara_args: list[str] = field(default_factory=lambda: YARA_DEFAULT_ARGS[:])
    mapping_threads: int = 4
    solver: str = "glpk"
    ilp_threads: int = 1
    delete_bam: bool = True
    unpaired_weight: float = 0.0
    use_discordant: bool = False


@dataclass
class PipelineResult:
    """Results from the OptiType pipeline."""

    result_df: pd.DataFrame
    result_4digit: pd.DataFrame
    output_csv: str
    output_plot: str


class LowCoverageError(RuntimeError):
    """Raised when HLA typing cannot proceed due to insufficient mapped HLA reads."""


def _write_low_coverage_result(out_csv: str, message: str) -> None:
    """Write a structured failure TSV so callers get an explicit output artifact."""
    failed = pd.DataFrame([
        {
            "A1": None,
            "A2": None,
            "B1": None,
            "B2": None,
            "C1": None,
            "C2": None,
            "Reads": 0,
            "Objective": None,
            "Status": "FAILED",
            "Message": message,
        }
    ])
    failed.to_csv(out_csv, sep="\t", index=False)


def get_num_threads(configured_threads: int) -> int:
    """Get the number of threads to use, capped by available CPUs."""
    try:
        cpu_count = multiprocessing.cpu_count()
    except (ImportError, NotImplementedError):
        return 2

    return min(cpu_count, configured_threads)


def _is_frequent(allele_id: str, table: pd.DataFrame) -> bool:
    """Check if an allele is in the frequent alleles list."""
    allele_id = allele_id.split("_")[0]
    return (
        table.loc[allele_id]["4digit"] in FREQ_ALLELES
        and table.loc[allele_id]["flags"] == 0
    ) or (table.loc[allele_id]["locus"] in "HGJ")


def _get_4digit(allele_id: str, table: pd.DataFrame) -> str:
    """Get the 4-digit type for an allele."""
    allele_id = allele_id.split("_")[0]
    return table.loc[allele_id]["4digit"]


def _get_types(allele_id: str | float, table: pd.DataFrame) -> str | float:
    """Get the 4-digit type, handling non-string inputs."""
    if not isinstance(allele_id, str):
        return allele_id
    aa = allele_id.split("_")
    return table.loc[aa[0]]["4digit"]


def _parse_mapper_args(parser: RawConfigParser) -> tuple[list[str], list[str]]:
    """Extract mapper args from a parsed config, falling back to defaults."""
    razers3_args = shlex.split(
        parser.get("razers3", "args", fallback=" ".join(RAZERS3_DEFAULT_ARGS))
    )
    yara_args = shlex.split(
        parser.get("yara", "args", fallback=" ".join(YARA_DEFAULT_ARGS))
    )
    return razers3_args, yara_args


def load_mapper_config(config_path: str | Path | None = None) -> tuple[list[str], list[str]]:
    """
    Load mapper arguments from a config file.

    The config file uses INI format with [razers3] and [yara] sections,
    each containing an 'args' key with command-line arguments.

    Args:
        config_path: Path to mapper config file. If None, returns defaults.

    Returns:
        Tuple of (razers3_args, yara_args) as lists of strings.
    """
    parser = RawConfigParser()
    if config_path:
        parser.read(config_path)
    return _parse_mapper_args(parser)


def load_config(config_path: str | Path | None = None) -> PipelineConfig:
    """
    Load configuration from an INI file.

    Args:
        config_path: Path to config file. If None, uses environment or defaults.

    Returns:
        PipelineConfig with loaded settings.
    """
    parser = RawConfigParser()

    if config_path and Path(config_path).exists():
        parser.read(config_path)

    razers3_args, yara_args = _parse_mapper_args(parser)

    return PipelineConfig(
        razers3_path=parser.get("mapping", "razers3", fallback="razers3"),
        mapper=parser.get("mapping", "mapper", fallback="razers3"),
        razers3_args=razers3_args,
        yara_args=yara_args,
        mapping_threads=parser.getint("mapping", "threads", fallback=4),
        solver=parser.get("ilp", "solver", fallback="glpk"),
        ilp_threads=parser.getint("ilp", "threads", fallback=1),
        delete_bam=parser.getboolean("behavior", "deletebam", fallback=True),
        unpaired_weight=parser.getfloat("behavior", "unpaired_weight", fallback=0.0),
        use_discordant=parser.getboolean("behavior", "use_discordant", fallback=False),
    )


def _check_binary(name: str, install_hint: str) -> None:
    """Raise FileNotFoundError if binary is not in PATH."""
    if not shutil.which(name):
        raise FileNotFoundError(f"{name} not found in PATH. {install_hint}")


def _run_cmd(cmd: list[str], label: str) -> None:
    """Run a subprocess command, raising RuntimeError with stderr on failure."""
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"{label} failed (exit {proc.returncode}): {proc.stderr.strip()}")


def _run_razers3(
    config: PipelineConfig,
    mapping_ref: Path,
    input_files: list[str],
    bam_paths: list[str],
    threads: int,
    ref_type: str,
) -> None:
    """Map reads to HLA reference using RazerS3."""
    _check_binary(config.razers3_path, "Install with: conda install -c bioconda razers3")
    logger.debug("Mapping with RazerS3 using %d threads...", threads)

    for sample, outbam in zip(input_files, bam_paths, strict=True):
        logger.debug("%s Mapping %s to %s reference...", elapsed(), os.path.basename(sample), ref_type.upper())
        _run_cmd([
            config.razers3_path, *config.razers3_args,
            "-tc", str(threads), "-o", str(outbam), str(mapping_ref), str(sample),
        ], "RazerS3")


def _run_yara(
    config: PipelineConfig,
    mapping_ref: Path,
    input_files: list[str],
    bam_paths: list[str],
    threads: int,
    ref_type: str,
) -> None:
    """Map reads to HLA reference using YARA (index + map)."""
    _check_binary("yara_mapper", "Install with: conda install -c bioconda yara")
    _check_binary("yara_indexer", "Install with: conda install -c bioconda yara")
    logger.debug("Mapping with YARA using %d threads...", threads)

    with tempfile.TemporaryDirectory(prefix="optitype_yara_") as tmpdir:
        index_prefix = os.path.join(tmpdir, "ref")
        logger.debug("%s Building YARA index for %s reference...", elapsed(), ref_type.upper())
        _run_cmd(["yara_indexer", "-o", index_prefix, str(mapping_ref)], "yara_indexer")

        for sample, outbam in zip(input_files, bam_paths, strict=True):
            logger.debug("%s Mapping %s to %s reference...", elapsed(), os.path.basename(sample), ref_type.upper())
            _run_cmd([
                "yara_mapper", *config.yara_args,
                "-t", str(threads), "-o", str(outbam), index_prefix, str(sample),
            ], "yara_mapper")


def run_pipeline(
    input_files: list[str],
    seq_type: str,
    output_dir: str,
    prefix: str | None = None,
    beta: float = 0.009,
    enumerate_count: int = 1,
    config: PipelineConfig | None = None,
    verbose: bool = False,
) -> PipelineResult:
    """
    Run the OptiType HLA typing pipeline.

    Args:
        input_files: List of input FASTQ or BAM files (1 or 2 files).
        seq_type: Sequencing type ('dna' or 'rna').
        output_dir: Output directory path.
        prefix: Output file prefix (defaults to timestamp).
        beta: Homozygosity detection parameter (0.0 to 0.1).
        enumerate_count: Number of solutions to enumerate.
        config: Pipeline configuration.
        verbose: Enable verbose output.

    Returns:
        PipelineResult with typing results.
    """
    if config is None:
        config = PipelineConfig()

    # Set verbosity
    configure_logging(verbose)

    # Validate inputs
    if len(input_files) not in (1, 2):
        raise ValueError("Number of input files must be 1 (single-end) or 2 (paired-end)")

    if seq_type not in ("dna", "rna"):
        raise ValueError("seq_type must be 'dna' or 'rna'")

    if not 0.0 <= beta < 0.1:
        raise ValueError("beta must be between 0.0 and 0.1")

    if enumerate_count <= 0:
        raise ValueError("enumerate_count must be > 0")

    # Check input file extensions
    input_extension = input_files[0].split(".")[-1].lower()
    bam_input = input_extension in ("sam", "bam")

    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Set up paths
    if prefix is None:
        prefix = datetime.datetime.fromtimestamp(time.time()).strftime("%Y_%m_%d_%H_%M_%S")
        out_dir = os.path.join(output_dir, prefix)
    else:
        out_dir = output_dir

    Path(out_dir).mkdir(parents=True, exist_ok=True)

    # YARA always produces BAM; RazerS3 can output SAM when pysam is unavailable
    if config.mapper == "yara" or PYSAM_AVAILABLE:
        extension = "bam"
    else:
        extension = "sam"
    bam_paths = (
        input_files
        if bam_input
        else [os.path.join(out_dir, f"{prefix}_{i + 1}.{extension}") for i in range(len(input_files))]
    )

    ref_type = "nuc" if seq_type == "rna" else "gen"
    is_paired = len(input_files) > 1

    out_csv = os.path.join(out_dir, f"{prefix}_result.tsv")
    out_plot = os.path.join(out_dir, f"{prefix}_coverage_plot.pdf")

    # Run mapping if needed
    if not bam_input:
        threads = get_num_threads(config.mapping_threads)
        mapping_ref = get_reference_fasta(seq_type)

        if config.mapper == "yara":
            _run_yara(config, mapping_ref, input_files, bam_paths, threads, ref_type)
        elif config.mapper == "razers3":
            _run_razers3(config, mapping_ref, input_files, bam_paths, threads, ref_type)
        else:
            raise ValueError(f"Unknown mapper: {config.mapper!r}. Must be 'razers3' or 'yara'.")

    # Load reference data
    table, features = load_reference_data()

    logger.debug("%s Generating binary hit matrix.", elapsed())

    # Process alignments
    if is_paired:
        pos, read_details = pysam_to_dataframe(bam_paths[0])
        binary1 = np.sign(pos)

        pos2, read_details2 = pysam_to_dataframe(bam_paths[1])
        binary2 = np.sign(pos2)

        if pos.empty and pos2.empty:
            msg = (
                "HLA typing was not possible due to low coverage: no mapped HLA reads "
                "were found in either input read file."
            )
            logger.error(msg)
            _write_low_coverage_result(out_csv, msg)
            raise LowCoverageError(f"{msg} See {out_csv}.")

        if not bam_input and config.delete_bam:
            os.remove(bam_paths[0])
            os.remove(bam_paths[1])

        id1 = set(binary1.index)
        id2 = set(binary2.index)

        # Handle read ID suffixes
        if len(set(r[-1] for r in id1)) == 1 and len(set(r[-1] for r in id2)) == 1:
            def cut_last_char(x): return x[:-1]
            binary1.index = list(map(cut_last_char, binary1.index))
            binary2.index = list(map(cut_last_char, binary2.index))
            pos.index = list(map(cut_last_char, pos.index))
            pos2.index = list(map(cut_last_char, pos2.index))
            read_details.index = list(map(cut_last_char, read_details.index))
            read_details2.index = list(map(cut_last_char, read_details2.index))

        binary_p, binary_mis, binary_un = ht.create_paired_matrix(binary1, binary2)

        if binary_p.shape[0] < len(id1) * 0.1:
            logger.warning(
                "Less than 10%% of reads could be paired. Consider an appropriate "
                "unpaired_weight setting in your config file (currently %.3f), "
                "because you may need to resort to using unpaired reads.",
                config.unpaired_weight,
            )

        if config.unpaired_weight > 0:
            if config.use_discordant:
                binary = pd.concat([binary_p, binary_un, binary_mis])
            else:
                binary = pd.concat([binary_p, binary_un])
        else:
            binary = binary_p
    else:
        pos, read_details = pysam_to_dataframe(bam_paths[0])

        if pos.empty:
            msg = (
                "HLA typing was not possible due to low coverage: no mapped HLA reads "
                "were found in the input file."
            )
            logger.error(msg)
            _write_low_coverage_result(out_csv, msg)
            raise LowCoverageError(f"{msg} See {out_csv}.")

        if not bam_input and config.delete_bam:
            os.remove(bam_paths[0])

        binary = np.sign(pos)

    if binary.shape[0] == 0:
        msg = (
            "HLA typing was not possible due to low coverage: no usable mapped HLA reads "
            "remain after read pairing/selection."
        )
        logger.error(msg)
        _write_low_coverage_result(out_csv, msg)
        raise LowCoverageError(f"{msg} See {out_csv}.")

    # Filter to frequent alleles
    alleles_to_keep = [col for col in binary.columns if _is_frequent(col, table)]
    binary = binary[alleles_to_keep]

    if binary.shape[1] == 0:
        msg = (
            "HLA typing was not possible due to low coverage: no informative HLA allele "
            "hits remained after filtering."
        )
        logger.error(msg)
        _write_low_coverage_result(out_csv, msg)
        raise LowCoverageError(f"{msg} See {out_csv}.")

    logger.debug("%s Temporary pruning of identical rows and columns", elapsed())

    unique_col, _representing = ht.prune_identical_alleles(binary, report_groups=True)

    temp_pruned = ht.prune_identical_reads(unique_col)

    logger.debug("%s Size of mtx with unique rows and columns: %s", elapsed(), temp_pruned.shape)
    logger.debug("%s Determining minimal set of non-overshadowed alleles", elapsed())

    minimal_alleles = ht.prune_overshadowed_alleles(temp_pruned)

    if len(minimal_alleles) == 0:
        msg = (
            "HLA typing was not possible due to low coverage: no minimal informative "
            "alleles could be retained."
        )
        logger.error(msg)
        _write_low_coverage_result(out_csv, msg)
        raise LowCoverageError(f"{msg} See {out_csv}.")

    logger.debug("%s Keeping only the minimal number of required alleles %s", elapsed(), minimal_alleles.shape)

    binary = binary[minimal_alleles]

    logger.debug("%s Creating compact model...", elapsed())

    if is_paired and config.unpaired_weight > 0:
        if config.use_discordant:
            compact_mtx, compact_occ = ht.get_compact_model(
                binary_p[minimal_alleles],
                pd.concat([binary_un, binary_mis])[minimal_alleles],
                weight=config.unpaired_weight,
            )
        else:
            compact_mtx, compact_occ = ht.get_compact_model(
                binary_p[minimal_alleles],
                binary_un[minimal_alleles],
                weight=config.unpaired_weight,
            )
    else:
        compact_mtx, compact_occ = ht.get_compact_model(binary)

    allele_ids = binary.columns

    groups_4digit = defaultdict(list)
    for allele in allele_ids:
        type_4digit = _get_4digit(allele, table)
        groups_4digit[type_4digit].append(allele)

    sparse_dict = ht.mtx_to_sparse_dict(compact_mtx)
    threads = get_num_threads(config.ilp_threads)

    logger.debug("Starting ILP solver with %d threads...", threads)
    logger.debug("%s Initializing OptiType model...", elapsed())

    op = OptiType(
        sparse_dict,
        compact_occ,
        groups_4digit,
        table,
        beta,
        2,
        config.solver,
        threads,
        verbosity=verbose,
    )
    result = op.solve(enumerate_count)

    logger.debug("%s Result dataframe has been constructed...", elapsed())

    result_4digit = result.map(lambda x: _get_types(x, table))
    for col in ["A1", "A2", "B1", "B2", "C1", "C2"]:
        if col not in result_4digit:
            result_4digit[col] = None

    r = result_4digit[["A1", "A2", "B1", "B2", "C1", "C2", "nof_reads", "obj"]]

    # Write results
    r.to_csv(
        out_csv,
        sep="\t",
        columns=["A1", "A2", "B1", "B2", "C1", "C2", "nof_reads", "obj"],
        header=["A1", "A2", "B1", "B2", "C1", "C2", "Reads", "Objective"],
    )

    # Generate coverage plot
    hlatype = result.iloc[0].reindex(["A1", "A2", "B1", "B2", "C1", "C2"]).drop_duplicates().dropna()
    features_used = (
        [("intron", 1), ("exon", 2), ("intron", 2), ("exon", 3), ("intron", 3)]
        if seq_type == "dna"
        else [("exon", 2), ("exon", 3)]
    )

    plot_variables = (
        [pos, read_details, pos2, read_details2, (binary_p, binary_un, binary_mis)]
        if is_paired
        else [pos, read_details]
    )

    coverage_mat = ht.calculate_coverage(plot_variables, features, hlatype, features_used)
    ht.plot_coverage(out_plot, coverage_mat, table, features, features_used)

    return PipelineResult(
        result_df=result,
        result_4digit=r,
        output_csv=out_csv,
        output_plot=out_plot,
    )
