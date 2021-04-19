# snakemake --snakefile Snakefile-drb345 --use-conda --conda-frontend mamba --conda-prefix=/home/szolek/condaenvs/snakeenv -j8 -n
# /home/szolek/condaenvs/snakeenv/6f765271/conda-meta/
# snakemake --use-conda --conda-prefix=/home/szolek/condaenvs/snakeenv --config fqprefix=highcov -j8 results/highcov/all.tsv
import os
import pandas as pd

batch_sheets = {}

# Allowing a hackish way to prefix the fastq paths using --config fqprefix=/path/to/fqs in the
# snakemake call. Note: this must not be the last flag argument in the call, a further flag
# such as -j8 should come afterwards, otherwise the parser get confused.
FQ_PATH_PREFIX = config.get('fqprefix', '')

def read_sheet(batch_name):
    batch_sheet = pd.read_csv(f"datasets/{batch_name}.tsv", sep="\t", comment="#")
    if batch_sheet.columns.tolist() != ['sample', 'reads1', 'reads2']:
        # missing or malformed header in {batch_name}.tsv, assuming headerless sample/read1/read2 structure
        batch_sheet = batch_sheet.T.reset_index().T
        batch_sheet.columns = ['sample', 'reads1', 'reads2']
    batch_sheet.insert(0, "batch", batch_name)
    batch_sheet.set_index(["batch", "sample"], inplace=True)
    return batch_sheet

def batch_outputs(wc):
    batch_sheet = read_sheet(wc.batch)
    batch_sheets[wc.batch] = batch_sheet
    target_files = [f"results/{batch}/{sample}/result.tsv" for batch, sample in batch_sheet.index]
    return target_files

def fetch_fq(wc):
    # if run locally, batch_sheets is loaded once and available forever, but when deployed to
    # a cluster, the Snakefile is run from scratch before each job, and the already loaded sheet
    # lost. Workaround, load the entire sheet if not present:
    batch_sheet = batch_sheets.get(wc.batch, read_sheet(wc.batch))
    fq_path = batch_sheet.loc[(wc.batch, wc.sample), f"reads{wc.readend}"]
    # for relative paths, assume fastq as the base directory. os.path.join leaves absolutes alone
    return os.path.join('fastq', FQ_PATH_PREFIX, fq_path) 


rule all:
    input:
        'datasets/{batch}.tsv',  # for timestamp comparison whether it needs to be re-run
        batch_outputs
    output:
        out="results/{batch}/all.tsv"
    run:
        pd.concat([pd.read_csv(infile, sep="\t") for infile in input[1:]], ignore_index=True).to_csv(output.out, sep="\t", index=False)


rule ref_index:
    input:
        "py/hla_ref/dna_dom.fasta"
    output:
        touch("ref_idx/ok")
    params:
        target="ref_idx/dna_dom"
    threads:
        1
    resources:
        mem_mb="1G",
    conda:
        "py/environment.yaml"    
    shell:
        "yara_indexer -o {params.target} {input}"


rule map_dom:
    input:
        fq=fetch_fq,
        ref="ref_idx/ok"
    output:
        out="mapped/{batch}/{sample}/r{readend}.bam",
        ok=touch("mapped/{batch}/{sample}/r{readend}.ok")  # samtools doesn't remove output if broken
    params:
        errate=5,
        sensitivity="full",
        stratacount=2,
        reference="ref_idx/dna_dom"
    threads:
        4
    resources:
        mem_mb="2G"
    log:
        "mapped/{batch}/{sample}/r{readend}.log"
    benchmark:
        "mapped/{batch}/{sample}/r{readend}.benchmark"
    conda:
        "py/environment.yaml"
    shell:
        "yara_mapper -t {threads} -e {params.errate} -as -sc {params.stratacount} -y {params.sensitivity} -vv {params.reference} {input.fq} 2> {log} | samtools view -b -F 4 - > {output.out}"


rule runtype:
    input:
        r1="mapped/{batch}/{sample}/r1.bam",
        r2="mapped/{batch}/{sample}/r2.bam",
        ok1="mapped/{batch}/{sample}/r1.ok",
        ok2="mapped/{batch}/{sample}/r2.ok"
    output:
        "results/{batch}/{sample}/result.tsv",
    log:
        "results/{batch}/{sample}/run.log"
    benchmark:
        "results/{batch}/{sample}/run.benchmark"
    params:
        outdir=lambda w: f"results/{w.batch}/{w.sample}",
    threads:
        1
    resources:
        mem_mb="2G"
    conda:
        "py/environment.yaml"
    shell:
        "python py/enumtyper.py -b1 {input.r1} -b2 {input.r2} -o {params.outdir} -s {wildcards.sample} -r {wildcards.batch} --png-plots --dpi 200 &> {log}"