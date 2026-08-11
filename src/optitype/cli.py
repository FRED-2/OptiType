"""
Click-based command-line interface for OptiType.

Usage:
    optitype run -i sample_1.fq -i sample_2.fq --dna -o output/
    optitype check-deps
    optitype init-config
"""

import shutil
import sys
from pathlib import Path

import click

from optitype import __version__
from optitype.io.data import get_data_path, get_reference_fasta


@click.group()
@click.version_option(version=__version__, prog_name="optitype")
def main():
    """OptiType: Precision HLA typing from next-generation sequencing data.

    OptiType is a novel HLA genotyping algorithm based on integer linear
    programming, capable of producing accurate 4-digit HLA genotyping
    predictions from NGS data.

    \b
    External dependencies:
      - RazerS3: For read mapping (conda install -c bioconda razers3)
      - YARA: Optional alternative mapper (conda install -c bioconda yara)
      - ILP solver: GLPK, CBC, or CPLEX

    \b
    Example usage:
      optitype run -i reads_1.fq -i reads_2.fq --dna -o results/
      optitype run -i sample.bam --rna -o results/
    """
    pass


@main.command()
@click.option(
    "-i", "--input",
    "input_files",
    multiple=True,
    required=True,
    type=click.Path(exists=True),
    help="Input FASTQ or BAM files. One file for single-end, two for paired-end.",
)
@click.option(
    "-r", "--rna",
    "seq_type",
    flag_value="rna",
    help="Input data is RNA sequencing.",
)
@click.option(
    "-d", "--dna",
    "seq_type",
    flag_value="dna",
    default=True,
    help="Input data is DNA sequencing (default).",
)
@click.option(
    "-o", "--outdir",
    required=True,
    type=click.Path(),
    help="Output directory for results.",
)
@click.option(
    "-p", "--prefix",
    default=None,
    help="Output filename prefix. Default: timestamp.",
)
@click.option(
    "-b", "--beta",
    type=float,
    default=0.009,
    help="Homozygosity detection parameter (0.0-0.1). Default: 0.009.",
)
@click.option(
    "-e", "--enumerate",
    "enumerate_count",
    type=int,
    default=1,
    help="Number of solutions to enumerate. Default: 1.",
)
@click.option(
    "--solver",
    type=click.Choice(["glpk", "cbc", "cplex"]),
    default="glpk",
    envvar="OPTITYPE_SOLVER",
    help="ILP solver to use. Default: glpk.",
)
@click.option(
    "--mapper",
    type=click.Choice(["razers3", "yara"]),
    default="razers3",
    help="Read mapper to use. Default: razers3.",
)
@click.option(
    "--razers3",
    type=click.Path(),
    default=None,
    envvar="OPTITYPE_RAZERS3",
    help="Path to RazerS3 binary. Default: search PATH.",
)
@click.option(
    "--threads",
    type=int,
    default=4,
    help="Number of threads for mapping. Default: 4.",
)
@click.option(
    "--ilp-threads",
    type=int,
    default=1,
    help="Number of threads for ILP solver. Default: 1.",
)
@click.option(
    "--keep-bam/--delete-bam",
    default=False,
    help="Keep or delete intermediate BAM files. Default: delete.",
)
@click.option(
    "-c", "--config",
    type=click.Path(exists=True),
    default=None,
    help="Path to config.ini file for additional settings.",
)
@click.option(
    "--mapper-config",
    type=click.Path(exists=True),
    default=None,
    help="Path to mapper config file with custom alignment parameters.",
)
@click.option(
    "-v", "--verbose",
    is_flag=True,
    help="Enable verbose output.",
)
def run(
    input_files,
    seq_type,
    outdir,
    prefix,
    beta,
    enumerate_count,
    solver,
    mapper,
    razers3,
    threads,
    ilp_threads,
    keep_bam,
    config,
    mapper_config,
    verbose,
):
    """Run HLA typing analysis.

    \b
    Examples:
      # Paired-end DNA analysis
      optitype run -i reads_1.fq -i reads_2.fq --dna -o results/

      # Single-end RNA analysis
      optitype run -i sample.fastq --rna -o results/

      # Re-analyze from BAM file
      optitype run -i mapped.bam --dna -o results/

      # With custom settings
      optitype run -i reads.fq --dna -o results/ --solver cbc --threads 8
    """
    from optitype.pipeline import LowCoverageError, PipelineConfig, load_config, load_mapper_config, run_pipeline

    # Resolve razers3 path (--razers3 option or search PATH)
    if razers3 is None:
        razers3 = shutil.which("razers3") or "razers3"

    # Load mapper-specific arguments
    razers3_args, yara_args = load_mapper_config(mapper_config)

    # Build config: start from file config (or defaults), override with CLI args
    pipeline_config = load_config(config) if config else PipelineConfig()
    pipeline_config.razers3_path = razers3
    pipeline_config.mapper = mapper
    pipeline_config.razers3_args = razers3_args
    pipeline_config.yara_args = yara_args
    pipeline_config.mapping_threads = threads
    pipeline_config.solver = solver
    pipeline_config.ilp_threads = ilp_threads
    pipeline_config.delete_bam = not keep_bam

    # Run the pipeline
    if verbose:
        click.echo(f"OptiType {__version__}")
        click.echo(f"Input files: {list(input_files)}")
        click.echo(f"Sequence type: {seq_type}")
        click.echo(f"Output directory: {outdir}")
        click.echo(f"Mapper: {mapper}")
        click.echo(f"Solver: {solver}")
        click.echo()

    try:
        result = run_pipeline(
            input_files=list(input_files),
            seq_type=seq_type,
            output_dir=outdir,
            prefix=prefix,
            beta=beta,
            enumerate_count=enumerate_count,
            config=pipeline_config,
            verbose=verbose,
        )

        # Print results summary
        click.echo()
        click.echo("=" * 60)
        click.echo("HLA Typing Results")
        click.echo("=" * 60)

        best = result.result_4digit.iloc[0]
        for locus in ["A", "B", "C"]:
            a1 = best.get(f"{locus}1", "-")
            a2 = best.get(f"{locus}2", "-")
            if a1 == a2:
                click.echo(f"  HLA-{locus}: {a1} (homozygous)")
            else:
                click.echo(f"  HLA-{locus}: {a1}, {a2}")

        click.echo()
        click.echo(f"Reads: {best['nof_reads']:.0f}")
        click.echo(f"Objective: {best['obj']:.2f}")
        click.echo()
        click.echo(f"Results written to: {result.output_csv}")
        click.echo(f"Coverage plot: {result.output_plot}")

    except LowCoverageError as e:
        click.secho(f"WARNING: {e}", fg="yellow", err=True)
        return
    except Exception as e:
        raise click.ClickException(str(e)) from None


@main.command("check-deps")
def check_deps():
    """Check that all external dependencies are available.

    Verifies that RazerS3 and an ILP solver are installed and accessible.
    """
    click.echo("Checking OptiType dependencies...")
    click.echo()

    all_ok = True

    # Check RazerS3
    razers3 = shutil.which("razers3")
    if razers3:
        click.echo(click.style("  [OK]", fg="green") + f" RazerS3: {razers3}")
    else:
        click.echo(click.style("  [MISSING]", fg="red") + " RazerS3")
        click.echo("    Install with: conda install -c bioconda razers3")
        all_ok = False

    # Check YARA (optional alternative mapper)
    yara_mapper = shutil.which("yara_mapper")
    yara_indexer = shutil.which("yara_indexer")
    if yara_mapper and yara_indexer:
        click.echo(click.style("  [OK]", fg="green") + f" YARA mapper: {yara_mapper}")
        click.echo(click.style("  [OK]", fg="green") + f" YARA indexer: {yara_indexer}")
    else:
        click.echo(click.style("  [N/A]", fg="cyan") + " YARA (optional alternative mapper: conda install -c bioconda yara)")

    # Check GLPK
    glpsol = shutil.which("glpsol")
    if glpsol:
        click.echo(click.style("  [OK]", fg="green") + f" GLPK: {glpsol}")
    else:
        click.echo(click.style("  [MISSING]", fg="yellow") + " GLPK")
        click.echo("    Install with: apt install glpk-utils  (or conda install -c conda-forge glpk)")

    # Check CBC
    cbc = shutil.which("cbc")
    if cbc:
        click.echo(click.style("  [OK]", fg="green") + f" CBC: {cbc}")
    else:
        click.echo(click.style("  [MISSING]", fg="yellow") + " CBC")
        click.echo("    Install with: apt install coinor-cbc  (or conda install -c conda-forge coincbc)")

    # Check CPLEX (optional)
    cplex = shutil.which("cplex")
    if cplex:
        click.echo(click.style("  [OK]", fg="green") + f" CPLEX: {cplex}")
    else:
        click.echo(click.style("  [N/A]", fg="cyan") + " CPLEX (optional, commercial solver)")

    # Check reference data
    click.echo()
    try:
        data_path = get_data_path()
        click.echo(click.style("  [OK]", fg="green") + f" Reference data: {data_path}")
    except FileNotFoundError:
        click.echo(click.style("  [MISSING]", fg="red") + " Reference data")
        click.echo("    Set OPTITYPE_DATA environment variable to data directory")
        all_ok = False

    click.echo()
    if all_ok and (glpsol or cbc):
        click.echo(click.style("All required dependencies are available!", fg="green"))
    else:
        click.echo(click.style("Some dependencies are missing.", fg="red"))
        if not (glpsol or cbc):
            click.echo("At least one ILP solver (GLPK or CBC) is required.")
        raise SystemExit(1)


@main.command("init-config")
@click.option(
    "-o", "--output",
    type=click.Path(),
    default="config.ini",
    help="Output path for config file. Default: config.ini",
)
@click.option(
    "--force",
    is_flag=True,
    help="Overwrite existing config file.",
)
def init_config(output, force):
    """Generate a default configuration file.

    Creates a config.ini file with default settings that can be customized.
    """
    config_content = """[mapping]

# Read mapper to use: razers3 or yara
# YARA is faster and more memory-efficient for large WGS datasets.
# mapper=razers3

# Path to RazerS3 binary. If not specified, searched in PATH.
# razers3=/path/to/razers3

# Number of threads for read mapping
threads=4

[ilp]

# ILP solver to use: glpk, cbc, or cplex
solver=glpk

# Number of threads for ILP solver
threads=1

[behavior]

# Delete intermediate BAM files after processing
deletebam=true

# Weight for unpaired reads in paired-end mode (0-1)
# 0 means unpaired reads are ignored
# 1 means unpaired reads are weighted equally to paired reads
unpaired_weight=0

# Use discordant read pairs (where ends map to different alleles)
use_discordant=false

[razers3]

# Arguments passed to razers3 (excluding threads and output path).
# -i: percent identity threshold (0-100)
# -m: max number of matches per read
# --distance-range: distance range for best mapping
# -pa: purge ambiguous alignments
args = -i 97 -m 99999 --distance-range 0 -pa

[yara]

# Arguments passed to yara_mapper (excluding threads and output path).
# Uses nf-core/hlatyping defaults.
# -e: error rate in percent (0-10)
# -f: output format (bam or sam)
args = -e 3 -f bam
"""

    output_path = Path(output)
    if output_path.exists() and not force:
        raise click.ClickException(
            f"Config file already exists: {output}\n"
            "Use --force to overwrite."
        )

    output_path.write_text(config_content)
    click.echo(f"Config file created: {output}")
    click.echo("Edit the file to customize settings.")


@main.command("info")
def info():
    """Show information about the OptiType installation."""
    click.echo(f"OptiType version: {__version__}")
    click.echo()

    # Python info
    click.echo(f"Python: {sys.version}")
    click.echo()

    # Check data path
    try:
        data_path = get_data_path()
        click.echo(f"Data directory: {data_path}")

        # Check for reference files
        for ref_type in ["dna", "rna"]:
            try:
                ref_path = get_reference_fasta(ref_type)
                click.echo(f"  {ref_type.upper()} reference: {ref_path}")
            except FileNotFoundError:
                click.echo(f"  {ref_type.upper()} reference: NOT FOUND")

        # Check for allele data
        alleles_path = data_path / "alleles"
        if alleles_path.exists():
            parquet_files = list(alleles_path.glob("*.parquet"))
            click.echo(f"  Allele data: {len(parquet_files)} parquet files")
        else:
            click.echo("  Allele data: NOT FOUND")
    except FileNotFoundError as e:
        click.echo(f"Data directory: NOT FOUND ({e})")

    click.echo()

    # Package dependencies
    click.echo("Dependencies:")
    for pkg in ["click", "numpy", "pandas", "pyarrow", "pyomo", "pysam", "matplotlib"]:
        try:
            mod = __import__(pkg)
            version = getattr(mod, "__version__", "unknown")
            click.echo(f"  {pkg}: {version}")
        except ImportError:
            click.echo(f"  {pkg}: NOT INSTALLED")


if __name__ == "__main__":
    main()
