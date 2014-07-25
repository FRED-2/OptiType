OptiType
========

Precision HLA typing from next-generation sequencing data

Authors: András Szolek, Benjamin Schubert, Christopher Mohr  
Date: April 2014  
Version: 1.0  
License: OptiType is released under a three-clause BSD license


Introduction:
-------------
OptiType, is a novel HLA genotyping algorithm based on integer linear
programming, capable of producing accurate 4-digit HLA genotyping predictions
from NGS data by simultaneously selecting all minor and major HLA-I alleles.


Requirement:
-------------
OptiType uses the following software and libraries:  

1. Python 2.7
2. Biopython 1.63
3. Coopr 3.3
4. Matplotlib 1.3.1
5. Pandas 0.12 (with HDF5 support)
6. HDF5 1.8.11
7. RazerS 3.1
8. Cplex 12.5 (or other ILP solver suppored by Coopr)

Please make sure you have installed said software/libraries
and their dependencies.


Installation:
-------------
First install all required software and libraries and register the static path
in the configuration file for RazerS 3.1. CPLEX should be globally executable
via command line. Alternative ILP solver supported by Cooper can also be used 
by changing the config file accordingly. CPLEX is free for academic use. For more 
details see IBMs Academic Initiative (http://www-304.ibm.com/ibm/university/academic/pub/page/academic_initiative)
Please do not change the folder structure or make sure you changed the necessary
entries in the config file.


Step-by-Step Installation Guide:
-------------
1) Install Python 2.7 on your system  
Either download it from https://www.python.org/download/ or just install it by easy_install e.g. depending on your system  

***The following installation steps require a working pip installation on your system.***

2) Install the Biopython package (https://pypi.python.org/pypi/biopython)
```
pip install biopython
```

3) Install the Coopr python package (https://pypi.python.org/pypi/Coopr)
```
pip install Coopr
```

4) Install the Matplotlib python package (https://pypi.python.org/pypi/matplotlib)
```
pip install matplotlib
```

5) Install the Pandas python package (https://pypi.python.org/pypi/pandas/)
```
pip install pandas
```
5.1) Install HDF5 (http://www.hdfgroup.org/HDF5/)  
Either install the pre-built binaries or build from source as following
```
cd <top HDF5 source code directory>
./configure --prefix=<location for HDF5 software> 
make >& make.out
make check >& check.out
make install 
```
5.2) Install pyTables (http://www.pytables.org)  
Add HDF5 to LD_LIBRARY PATH.
```
pip install tables
```

6) Install RazerS (http://www.seqan.de/projects/razers/)  
Go to the project webpage and download a precompiled binary of RazerS 3 
or follow the instruction in the Compilation from Source Code section.

7) Install an ILP solver supported by Coopr (e.g. CPLEX or GLPK)


Usage:
-------------
1) First filter the read files with the following settings:
```
>razers3 --percent-identity 90 --max-hits 1 --distance-range 0 --output-format sam --output sample_fished.sam
         ./data/hla_reference.fasta sample.fastq
```
where reference.fasta is either nuc_reference.fasta or gen_reference.fasta
depending on the type of NGS data. The references can be found in the ./data
sub-folder or in the supplementary material. To use the results as input
for OptiType the sam-files have to be converted into fastq format. On Unix-
based operating system you can convert from sam to fastq with the following
command:
```
>cat sample_fished.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > sample_fished.fastq
```
For paired-end data pre-process each file individually.

2) After pre-filtering, OptiType can be called as follows:
```
>python OptiTypePipeline.py -i sample_fished_1.fastq [sample_fished_2.fastq]
                    (--rna | --dna) [--beta BETA] [--enumerate ENUMERATE]
                    --o ./out_dir/
```
This will produce a CSV with the optimal typing and possible sub-optimal
typings if specified, as well as a coverage plot of the genotype for
diagnostic purposes and a HTML file containing a summary of the results.
```
>python OptiTypePipeline.py --help  
usage: OptiType [-h] --input INPUT [INPUT ...] (--rna | --dna) [--beta BETA]  
                [--enumerate ENUMERATE] --outdir OUTDIR [--verbose]

OptiType: 4-digit HLA typer

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT [INPUT ...], -i INPUT [INPUT ...]
                        Fastq files with fished HLA reads. Max two files (for
                        paired-end)
  --rna, -r             Specifiying the mapped data as RNA.
  --dna, -d             Specifiying the mapped data as DNA.
  --beta BETA, -b BETA  The beta value for for homozygosity detection.
  --enumerate ENUMERATE, -e ENUMERATE
                        The number of enumerations.
  --outdir OUTDIR, -o   OUTDIR
                        Specifies the out directory to which all files should
                        be written
  --verbose, -v         Set verbose mode on.
```
Examples:
-------------
DNA data (paired-end):
```
python OptiTypePipeline.py -i ./test/exome/NA11995_SRR766010_1_fished.fastq ./test/exome/NA11995_SRR766010_2_fished.fastq -d -v -o ./test/exome/
```
RNA data (paired-end):
```
python OptiTypePipeline.py -i ./test/rna/CRC_81_N_2_fished.fastq ./test/rna/CRC_81_N_2_fished.fastq -r -v -o ./test/rna/
```
Contacts:
-------------
András Szolek  
szolek@informatik.uni-tuebingen.de  
University of Tübingen, Applied Bioinformatics,  
Center for Bioinformatics, Quantitative Biology Center,  
and Dept. of Computer Science,  
Sand 14, 72076 Tübingen, Germany


Downloads:
-------------

1. http://python.org/download/
2. http://biopython.org/
3. http://software.sandia.gov/trac/coopr
4. http://matplotlib.org/
5. http://pandas.pydata.org/
6. http://www.hdfgroup.org/HDF5/
7. https://www.seqan.de/projects/razers/
8. http://www-01.ibm.com/software/info/ilog/


Reference:
-------------
 Szolek, A., Schubert, B., Mohr, C., Feldhahn, M., Strum, M.,
 & Kohlbacher, O. (2014). OptiType: precision HLA typing from next-generation
 sequencing data. (in review)
