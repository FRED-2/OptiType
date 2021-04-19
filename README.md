## Installation

1\. Set up an environment with all dependencies using the `conda` package manager

    conda create -n ot2env -c conda-forge -c bioconda python=3.8 numpy pandas scipy matplotlib snakemake pysam samtools yara
    conda activate ot2env


2\. Clone the repository

    git clone --single-branch --branch hla-ligand-atlas https://github.com/FRED-2/OptiType.git ot2-snake
    cd ot2-snake

## Usage

The HLA typing pipeline's input is a tsv file residing in `datasets/<batch_id>.tsv` containing sample ID's and the corresponding paired-end fastq files. Output will be organized into a folder structure derived from `<batch_id>` and the sample IDs inside the data sheet.

Example files and folder structure:

`datasets/project1.tsv`

    sample  reads1  reads2
    ABC123  ABC123_R1.fq.gz ABC123_R1.fq.gz
    DEF456  DEF456_R1.fq.gz DEF456_R1.fq.gz
    ZZZ001  ZZZ001_R1.fq.gz ZZZ001_R1.fq.gz

`datasets/project2.tsv`

    sample  reads1  reads2
    sample1 linked_dir/SMP001_1.fastq   linked_dir/SMP001_2.fastq
    sample2 linked_dir/SMP002_1.fastq   linked_dir/SMP002_2.fastq
    sample3 linked_dir/SMP003_1.fastq   linked_dir/SMP003_2.fastq

Directory structure

    ┌───datasets/                     ┆
    │   ├───project1.tsv              ┆
    │   └───project2.tsv              ┆
    │                                 ┆
    ├───fastq/                        ┆
    │   ├──linked_dir/                ┆
    │   │   ├───SMP001_1.fastq        ┆
    │   │   ├───SMP001_2.fastq        ┆
    │   │   ├───SMP002_1.fastq        ┆ input
    │   │   ├───SMP002_2.fastq        ┆
    │   │   ├───SMP003_1.fastq        ┆
    │   │   └───SMP003_2.fastq        ┆
    │   ├───ABC123_R1.fq.gz           ┆
    │   ├───ABC123_R2.fq.gz           ┆
    │   ├───DEF456_R1.fq.gz           ┆
    │   ├───DEF456_R2.fq.gz           ┆
    │   ├───ZZZ001_R1.fq.gz           ┆
    │   └───ZZZ001_R2.fq.gz           ┆
    │
    ├───mapped/                              ┆ intermediate
    │
    ├───results/                                    ┆
    │   ├──project1/                                ┆
    │   │   ├───ABC123/  result.tsv  covplots.pdf   ┆
    │   │   ├───DEF456/  result.tsv  covplots.pdf   ┆
    │   │   ├───ZZZ001/  result.tsv  covplots.pdf   ┆
    │   │   └───all.tsv                             ┆ output
    │   └──project2/                                ┆
    │       ├───sample1/ result.tsv  covplots.pdf   ┆
    │       ├───sample2/ result.tsv  covplots.pdf   ┆
    │       ├───sample3/ result.tsv  covplots.pdf   ┆
    │       └───all.tsv                             ┆
    │
    ├───ref_idx/
    ├───py/
    └───Snakefile


The command `snakemake -k -j 8 results/project1/all.tsv` builds the target file `results/project1/all.tsv` by reading `datasets/project1.tsv`, processing all listed samples in `project1.tsv` and gathering their outputs in a single file, using up to 8 threads for the task. Individual samples can be processed by building their specific target files with `snakemake -j 8 results/project2/sample2/result.tsv`.

The fastq paths in the dataset tsv's shall be relative to the `fastq` directory. They can be hard- or symbolic links as well. Absolute paths pointing outside the directory structure altogether are also accepted in the dataset tsv's.

The code is bundled with an example sheet and corresponding fastq files which can be run with `snakemake -k -j 8 results/example/all.tsv`

## Notes

Given a conda environment with `snakemake` available, but not the other dependencies, the necessary environment with all dependencies can be created on the fly with the `--use-conda` flag. The first execution of the pipeline will initialize a new environment with all dependencies, and subsequent executions will re-use this environment.

Other powerful features of snakemake, such as [deploying the pipeline on a cluster](https://snakemake.readthedocs.io/en/stable/executing/cluster.html) can also be leveraged.

---

Coverage plots of individual samples are placed in their respective output folders. The coverage depth chart of each locus consists of the individual coverage plots of the two alleles mirrored about the x-axis, which also marks their distinguishing polymorphisms as ticks.

Reads contribute to the colored bands based on the following scheme:

* Reads exclusive to locus
  * <span style="background-color: #84cc8f">green</span>: exclusive to one of the two alleles
  * <span style="background-color: #51a1cc">blue</span>: attributed to both alleles (due to allele sequence identity)
* Reads mapping to multiple loci
  * <span style="background-color: #f97f7f">red</span>: exclusive to one allele
  * <span style="background-color: #ffb66d">orange</span>: attributed to both alleles

Bright shades stand for perfect alignments, light shades for mismatched alignments. <span style="background-color: #f2f2f2">Grey</span> bands in the background correspond to intron regions. HLA Class II predictions are limited to the exon 2+3 regions of the loci. Full allele names are displayed nonetheless, but they shall be seen as representatives of their ambiguity group.

---

The output tsv files also contain the [P-group](http://hla.alleles.org/alleles/p_groups.html) antigen binding domain designation of the predicted alleles. In certain contexts this representation may be preferable to four-digit genotypes.

---

The present release is limited to the detection of non-rare HLA alleles based on their [CWD 3.0](https://doi.org/10.1111/tan.13811) status, covering ~99% of known allelic diversity.