#!/usr/bin/env python
# coding=utf-8
"""
###################################################################

OptiType: precision HLA typing from next-generation sequencing data

###################################################################

Authors: András Szolek, Benjamin Schubert, Christopher Mohr
Date: August 2017
Version: 1.3.1
License: OptiType is released under a three-clause BSD license


Introduction:
-------------
OptiType, is a novel HLA genotyping algorithm based on integer linear
programming, capable of producing accurate 4-digit HLA genotyping predictions
from NGS data by simultaneously selecting all minor and major HLA-I alleles.


Requirements:
-------------
OptiType uses the following software and libraries:
1) Python 2.7
2) Biopython 1.63
3) Pyomo 4.1
4) Matplotlib 1.3.1
5) Pandas 0.12 (with HDF5 support)
6) HDF5 1.8.11
7) RazerS 3.1
8) Cplex 12.5

Please make sure you have installed said software/libraries
and their dependencies.


Installation:
-------------
First install all required software and libraries and register the static path
in the configuration file for RazerS 3.1. CPLEX should be globally executable
via command line. Alternative ILP solver supported by Cooper are also usable.
Please do not change the folder structure or make sure you changed the necessary
entries in the config file.


Usage:
-------------
1) First filter the read files with the following settings:

>razers3 --percent-identity 90 --max-hits 1 --distance-range 0
         --output-format sam --output sample_fished.sam
         ./data/hla_reference.fasta sample.fastq

where reference.fasta is either nuc_reference.fasta or gen_reference.fasta
depending on the type of NGS data. The references can be found in the ./data
sub-folder or in the supplementary material. To use the results as input
for OptiType the sam-files have to be converted into fastq format. On Unix-
based operating system you can convert from sam to fastq with the following
command:

>cat sample_fished.sam | grep -v ^@
	| awk '{print "@"$1"\n"$10"\n+\n"$11}' > sample_fished.fastq

For paired-end data pre-process each file individually.

2) After pre-filtering, OptiType can be called as follows:

>python OptiTypePipeline.py -i sample_fished_1.fastq [sample_fished_2.fastq]
                    (--rna | --dna) [--beta BETA] [--enumerate ENUMERATE]
                    --o ./out_dir/

This will produce a CSV with the optimal typing and possible sub-optimal
typings if specified, as well as a coverage plot of the genotype for
diagnostic purposes and a HTML file containing a summary of the results.

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
"""
from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import zip
from builtins import filter
from builtins import map
from builtins import range

## eliminate dependency on an X11 server
import matplotlib
matplotlib.use('Agg')

import sys
import subprocess
import os
import argparse
if sys.version_info > (3,0):
    import configparser
else:
    import ConfigParser as configparser
import time
import datetime
import pandas as pd
import hlatyper as ht

import numpy as np
from model import OptiType
from collections import defaultdict

# Try to import pysam, important for RazerS 3 output format
try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False


freq_alleles = '''A*01:01 A*01:02 A*01:03 A*01:06 A*01:09 A*01:23 A*01:38 A*01:44 A*02:01 A*02:02 A*02:03 A*02:04 A*02:05 A*02:06 A*02:07 A*02:08 A*02:09 A*02:10 A*02:11 A*02:12 A*02:13 A*02:133 A*02:14 A*02:141 A*02:146 A*02:15N A*02:16 A*02:17 A*02:18 A*02:19 A*02:20 A*02:21 A*02:22 A*02:226N A*02:24 A*02:25 A*02:26 A*02:27 A*02:28 A*02:29 A*02:30 A*02:33 A*02:34 A*02:35 A*02:36 A*02:37 A*02:38 A*02:40 A*02:42 A*02:44 A*02:45 A*02:46 A*02:48 A*02:49 A*02:51 A*02:53N A*02:54 A*02:55 A*02:57 A*02:58 A*02:60 A*02:64 A*02:67 A*02:74 A*02:85 A*02:90 A*02:93 A*03:01 A*03:02 A*03:05 A*03:07 A*03:08 A*03:10 A*03:12 A*03:22 A*03:25 A*03:65 A*03:69N A*03:97 A*11:01 A*11:02 A*11:03 A*11:04 A*11:05 A*11:06 A*11:08 A*11:10 A*11:12 A*11:13 A*11:18 A*11:19 A*11:20 A*11:29 A*11:40 A*23:01 A*23:02 A*23:03 A*23:04 A*23:05 A*23:09 A*24:02 A*24:03 A*24:04 A*24:05 A*24:06 A*24:07 A*24:08 A*24:09N A*24:10 A*24:13 A*24:14 A*24:15 A*24:17 A*24:18 A*24:20 A*24:21 A*24:22 A*24:23 A*24:24 A*24:25 A*24:26 A*24:27 A*24:28 A*24:29 A*24:31 A*24:35 A*24:46 A*24:51 A*24:56 A*24:63 A*24:93 A*25:01 A*25:02 A*25:04 A*26:01 A*26:02 A*26:03 A*26:04 A*26:05 A*26:06 A*26:07 A*26:08 A*26:09 A*26:10 A*26:11N A*26:12 A*26:14 A*26:15 A*26:16 A*26:17 A*26:18 A*26:20 A*26:49 A*29:01 A*29:02 A*29:03 A*29:04 A*29:10 A*29:12 A*30:01 A*30:02 A*30:03 A*30:04 A*30:06 A*30:08 A*30:09 A*30:10 A*30:11 A*30:12 A*30:16 A*31:01 A*31:02 A*31:03 A*31:04 A*31:05 A*31:06 A*31:08 A*31:09 A*31:12 A*32:01 A*32:02 A*32:03 A*32:04 A*32:05 A*32:06 A*32:08 A*32:13 A*32:20 A*32:22 A*33:01 A*33:03 A*33:04 A*33:05 A*33:10 A*33:26 A*34:01 A*34:02 A*34:03 A*34:05 A*36:01 A*36:03 A*43:01 A*66:01 A*66:02 A*66:03 A*68:01 A*68:02 A*68:03 A*68:04 A*68:05 A*68:06 A*68:07 A*68:08 A*68:12 A*68:13 A*68:15 A*68:16 A*68:17 A*68:18N A*68:23 A*68:24 A*68:38 A*69:01 A*74:01 A*74:02 A*74:03 A*74:04 A*74:06 A*74:09 A*74:11 A*80:01 A*80:02 B*07:02 B*07:03 B*07:04 B*07:05 B*07:06 B*07:07 B*07:08 B*07:09 B*07:10 B*07:12 B*07:13 B*07:14 B*07:15 B*07:17 B*07:20 B*07:22 B*07:26 B*07:33 B*07:36 B*07:47 B*07:53 B*07:85 B*08:01 B*08:02 B*08:03 B*08:04 B*08:05 B*08:09 B*08:12 B*08:18 B*08:23 B*13:01 B*13:02 B*13:03 B*13:04 B*13:07N B*13:09 B*13:11 B*13:13 B*14:01 B*14:02 B*14:03 B*14:04 B*14:05 B*14:06 B*15:01 B*15:02 B*15:03 B*15:04 B*15:05 B*15:06 B*15:07 B*15:08 B*15:09 B*15:10 B*15:108 B*15:11 B*15:12 B*15:123 B*15:125 B*15:13 B*15:135 B*15:15 B*15:153 B*15:16 B*15:17 B*15:18 B*15:20 B*15:21 B*15:23 B*15:24 B*15:25 B*15:27 B*15:28 B*15:29 B*15:30 B*15:31 B*15:32 B*15:33 B*15:34 B*15:35 B*15:36 B*15:37 B*15:38 B*15:39 B*15:40 B*15:42 B*15:45 B*15:46 B*15:47 B*15:48 B*15:50 B*15:52 B*15:53 B*15:54 B*15:55 B*15:56 B*15:58 B*15:61 B*15:63 B*15:67 B*15:68 B*15:70 B*15:71 B*15:73 B*15:82 B*15:86 B*18:01 B*18:02 B*18:03 B*18:04 B*18:05 B*18:06 B*18:07 B*18:08 B*18:09 B*18:11 B*18:13 B*18:14 B*18:18 B*18:19 B*18:20 B*18:28 B*18:33 B*27:01 B*27:02 B*27:03 B*27:04 B*27:05 B*27:06 B*27:07 B*27:08 B*27:09 B*27:10 B*27:11 B*27:12 B*27:13 B*27:14 B*27:19 B*27:20 B*27:21 B*27:30 B*27:39 B*35:01 B*35:02 B*35:03 B*35:04 B*35:05 B*35:06 B*35:08 B*35:09 B*35:10 B*35:11 B*35:12 B*35:13 B*35:14 B*35:15 B*35:16 B*35:17 B*35:18 B*35:19 B*35:20 B*35:21 B*35:22 B*35:23 B*35:24 B*35:25 B*35:27 B*35:28 B*35:29 B*35:30 B*35:31 B*35:32 B*35:33 B*35:34 B*35:36 B*35:43 B*35:46 B*35:51 B*35:77 B*35:89 B*37:01 B*37:02 B*37:04 B*37:05 B*38:01 B*38:02 B*38:04 B*38:05 B*38:06 B*38:15 B*39:01 B*39:02 B*39:03 B*39:04 B*39:05 B*39:06 B*39:07 B*39:08 B*39:09 B*39:10 B*39:11 B*39:12 B*39:13 B*39:14 B*39:15 B*39:23 B*39:24 B*39:31 B*39:34 B*40:01 B*40:02 B*40:03 B*40:04 B*40:05 B*40:06 B*40:07 B*40:08 B*40:09 B*40:10 B*40:11 B*40:12 B*40:14 B*40:15 B*40:16 B*40:18 B*40:19 B*40:20 B*40:21 B*40:23 B*40:27 B*40:31 B*40:35 B*40:36 B*40:37 B*40:38 B*40:39 B*40:40 B*40:42 B*40:44 B*40:49 B*40:50 B*40:52 B*40:64 B*40:80 B*41:01 B*41:02 B*41:03 B*42:01 B*42:02 B*44:02 B*44:03 B*44:04 B*44:05 B*44:06 B*44:07 B*44:08 B*44:09 B*44:10 B*44:12 B*44:13 B*44:15 B*44:18 B*44:20 B*44:21 B*44:22 B*44:26 B*44:27 B*44:29 B*44:31 B*44:59 B*45:01 B*45:02 B*45:04 B*45:06 B*46:01 B*46:02 B*46:13 B*47:01 B*47:02 B*47:03 B*48:01 B*48:02 B*48:03 B*48:04 B*48:05 B*48:06 B*48:07 B*48:08 B*49:01 B*49:02 B*49:03 B*50:01 B*50:02 B*50:04 B*50:05 B*51:01 B*51:02 B*51:03 B*51:04 B*51:05 B*51:06 B*51:07 B*51:08 B*51:09 B*51:10 B*51:12 B*51:13 B*51:14 B*51:15 B*51:18 B*51:21 B*51:22 B*51:27N B*51:29 B*51:31 B*51:32 B*51:33 B*51:34 B*51:36 B*51:37 B*51:63 B*51:65 B*52:01 B*52:02 B*52:06 B*53:01 B*53:02 B*53:03 B*53:04 B*53:05 B*53:07 B*53:08 B*54:01 B*54:02 B*55:01 B*55:02 B*55:03 B*55:04 B*55:07 B*55:08 B*55:10 B*55:12 B*55:16 B*55:46 B*56:01 B*56:02 B*56:03 B*56:04 B*56:05 B*56:06 B*56:07 B*56:09 B*56:11 B*57:01 B*57:02 B*57:03 B*57:04 B*57:05 B*57:06 B*57:10 B*58:01 B*58:02 B*58:06 B*59:01 B*67:01 B*67:02 B*73:01 B*78:01 B*78:02 B*78:03 B*78:05 B*81:01 B*81:02 B*82:01 B*82:02 B*83:01 C*01:02 C*01:03 C*01:04 C*01:05 C*01:06 C*01:08 C*01:14 C*01:17 C*01:30 C*01:32 C*02:02 C*02:03 C*02:04 C*02:06 C*02:08 C*02:09 C*02:10 C*02:19 C*02:20 C*02:27 C*03:02 C*03:03 C*03:04 C*03:05 C*03:06 C*03:07 C*03:08 C*03:09 C*03:10 C*03:13 C*03:14 C*03:15 C*03:16 C*03:17 C*03:19 C*03:21 C*03:32 C*03:42 C*03:43 C*03:56 C*03:67 C*03:81 C*04:01 C*04:03 C*04:04 C*04:05 C*04:06 C*04:07 C*04:10 C*04:11 C*04:14 C*04:15 C*04:24 C*04:29 C*04:33 C*04:37 C*05:01 C*05:03 C*05:04 C*05:07N C*05:09 C*06:02 C*06:03 C*06:04 C*06:06 C*06:08 C*06:09 C*06:17 C*06:24 C*06:53 C*07:01 C*07:02 C*07:03 C*07:04 C*07:05 C*07:06 C*07:07 C*07:08 C*07:09 C*07:10 C*07:109 C*07:123 C*07:13 C*07:14 C*07:17 C*07:18 C*07:19 C*07:21 C*07:22 C*07:27 C*07:29 C*07:32N C*07:37 C*07:43 C*07:46 C*07:49 C*07:56 C*07:66 C*07:67 C*07:68 C*07:80 C*07:95 C*08:01 C*08:02 C*08:03 C*08:04 C*08:05 C*08:06 C*08:13 C*08:15 C*08:20 C*08:21 C*08:27 C*12:02 C*12:03 C*12:04 C*12:05 C*12:07 C*12:12 C*12:15 C*12:16 C*14:02 C*14:03 C*14:04 C*15:02 C*15:03 C*15:04 C*15:05 C*15:06 C*15:07 C*15:08 C*15:09 C*15:11 C*15:12 C*15:13 C*15:17 C*16:01 C*16:02 C*16:04 C*16:08 C*16:09 C*17:01 C*17:02 C*17:03 C*17:04 C*18:01 C*18:02 A*30:07 B*15:64 B*18:12'''.split(' ')


def is_frequent(allele_id):
    # prepare for HLA12345_HLA67890 and use the first part
    allele_id = allele_id.split('_')[0]
    return table.loc[allele_id]['4digit'] in freq_alleles and table.loc[allele_id]['flags'] == 0 or (table.loc[allele_id]['locus'] in 'HGJ')


def get_4digit(allele_id):
    allele_id = allele_id.split('_')[0]  # for reconstructed IDs like HLA12345_HLA67890 return HLA12345's 4-digit type
    return table.loc[allele_id]['4digit']


def get_types(allele_id):
    if not isinstance(allele_id, str):
        return allele_id
    else:
        aa = allele_id.split('_')
        if len(aa) == 1:
            return table.loc[aa[0]]['4digit']
        else:
            return table.loc[aa[0]]['4digit']  #+ '/' + table.loc[aa[1]]['4digit']

def get_num_threads(configured_threads):
  try:
    import multiprocessing
  except (ImportError, NotImplementedError):
    return 2
  if(multiprocessing.cpu_count() < configured_threads):
    return multiprocessing.cpu_count()
  return configured_threads


if __name__ == '__main__':

    this_dir = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    config_default = os.path.join(this_dir, 'config.ini')

    parser = argparse.ArgumentParser(description=' OptiType: 4-digit HLA typer', prog='OptiType')
    parser.add_argument('--input','-i',
                      nargs='+',
                      required=True,
                      metavar='FQ',
                      help=(".fastq file(s) (fished or raw) or .bam files stored for re-use, generated by "
                          "an earlier OptiType run. One file: single-end mode, two files: paired-end mode.")
                      )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--rna','-r',
                      action="store_true",
                      help="Use with RNA sequencing data."
                      )
    group.add_argument('--dna','-d',
                      action="store_true",
                      help="Use with DNA sequencing data."
                      )
    parser.add_argument('--beta','-b',
                      type=float,
                      metavar='B',
                      default=0.009,
                      help="The beta value for for homozygosity detection (see paper). Default: 0.009. Handle with care."
                      )
    parser.add_argument('--enumerate','-e',
                      type=int,
                      default=1,
                      metavar='N',
                      help=("Number of enumerations. OptiType will output the optimal solution and "
                          "the top N-1 suboptimal solutions in the results CSV. Default: 1")
                      )
    parser.add_argument('--outdir','-o',
                      required=True,
                      help="Specifies the out directory to which all files should be written."
                      )
    parser.add_argument('--prefix', '-p',
                      default=None, dest="prefix", type=str,
                      help="Specifies prefix of output files"
                      )
    parser.add_argument('--verbose','-v',
                      required=False,
                      action="store_true",
                      help="Set verbose mode on."
                      )
    parser.add_argument('--config', '-c',
                      type=argparse.FileType('r'),
                      default=config_default,
                      help="Path to config file. Default: config.ini in the same directory as this script"
                      )

    args = parser.parse_args()

    if not os.path.isfile(args.config.name):
        print ("Config file not found. Place config.ini either alongside this script or use the -c option. "
            "See config.ini.example and note that its fields have changed recently.")
        sys.exit(-1)

    config = configparser.ConfigParser(os.environ)
    config.read(args.config.name)

    unpaired_weight = config.getfloat('behavior', 'unpaired_weight')
    use_discordant = config.getboolean('behavior', 'use_discordant')

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)        

    # test if inputs are legit:
    if args.beta < 0.0 or args.beta >= 0.1:
        print("Beta value is not correctly chosen. Please choose another beta value between [0,0.1]")
        sys.exit(-1)

    if args.enumerate <= 0:
        print("The specified number of enumerations must be bigger than %i"%args.enumeration)
        sys.exit(-1)

    if len(args.input) not in (1, 2):
        print("Number of input files can only be 1 (single-end) or 2 (paired-end)")
        sys.exit(-1)

    input_extension = args.input[0].split('.')[-1]
    assert all(ii.endswith('.' + input_extension) for ii in args.input), 'Mixed input file extensions'

    bam_input = (input_extension in ('sam', 'bam', 'SAM', 'BAM'))  # otherwise treated as fastq

    # Constants
    VERBOSE = ht.VERBOSE = bool(args.verbose)  # set verbosity setting in hlatyper too
    COMMAND = "-i 97 -m 99999 --distance-range 0 -pa -tc %d -o %s %s %s"
    ALLELE_HDF = os.path.join(this_dir, 'data/alleles.h5')
    MAPPING_REF = {'gen': os.path.join(this_dir, 'data/hla_reference_dna.fasta'),
                   'nuc': os.path.join(this_dir, 'data/hla_reference_rna.fasta')}
    MAPPING_CMD = config.get("mapping", "razers3") + " " + COMMAND
    date = datetime.datetime.fromtimestamp(time.time()).strftime('%Y_%m_%d_%H_%M_%S')
    if args.prefix == None:
        prefix = date
        out_dir = os.path.join(args.outdir, date)
    else:
        prefix = args.prefix
        out_dir = args.outdir

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if PYSAM_AVAILABLE:
        extension = 'bam'
    else:
        extension = 'sam'

    bam_paths = args.input if bam_input else [os.path.join(out_dir, ("%s_%i.%s" % (date, i+1, extension))) for i in range(len(args.input))]

    # SETUP variables and OUTPUT samples
    ref_type = "nuc" if args.rna else "gen"
    is_paired = len(args.input) > 1
    
    out_csv = os.path.join(out_dir, ("%s_result.tsv" % prefix))
    out_plot = os.path.join(out_dir, ("%s_coverage_plot.pdf" % prefix))

    # mapping fished file to reference
    if not bam_input:
        threads = get_num_threads(config.getint("mapping", "threads"))
        if VERBOSE:
          print("\nmapping with %s threads..." % threads)
        for (i, sample), outbam in zip(enumerate(args.input), bam_paths):          
            if VERBOSE:
                print("\n", ht.now(), "Mapping %s to %s reference..." % (os.path.basename(sample), ref_type.upper()))

            subprocess.call(MAPPING_CMD % (threads, outbam,
                                           MAPPING_REF[ref_type], sample), shell=True)

    # sam-to-hdf5
    table, features = ht.load_hdf(ALLELE_HDF, False, 'table', 'features')
    if VERBOSE:
        print("\n", ht.now(), "Generating binary hit matrix.")

    if is_paired:
        # combine matrices for paired-end mapping
        pos, read_details = ht.pysam_to_hdf(bam_paths[0])
        binary1 = np.sign(pos)  # dtype=np.uint16
        
        pos2, read_details2 = ht.pysam_to_hdf(bam_paths[1])
        binary2 = np.sign(pos2)  # dtype=np.uint16

        if not bam_input and config.getboolean('behavior', 'deletebam'):
            os.remove(bam_paths[0])
            os.remove(bam_paths[1])
        
        id1 = set(binary1.index)
        id2 = set(binary2.index)
        
        '''
            test if we actually can do paired-end mapping
            1) look at the last character 2-1 character if they are always the same if so -> proon them away and do
               paired-end
            2) if not test if the intersection of ID-binary1 and ID-binary2 has at least 10% of the former read
               number -> do paired-end
            3) if nothing worked ( perhaps pair-end ID was in the middle or something) raise flag and do single-end
               mapping on first input
        '''
        if len(set([r[-1] for r in id1])) == 1 and len(set([r[-1] for r in id2])) == 1:
            # if this case is true you have to edit also all pos,etc,desc indices such that the plotting works correctly
            # again .. maybe it is also neccessary to test for the last two characters
            cut_last_char = lambda x: x[:-1]
            binary1.index = list(map(cut_last_char, binary1.index))
            binary2.index = list(map(cut_last_char, binary2.index))
            pos.index = list(map(cut_last_char, pos.index))
            pos2.index = list(map(cut_last_char, pos2.index))
            read_details.index = list(map(cut_last_char, read_details.index))
            read_details2.index = list(map(cut_last_char, read_details2.index))
            
        binary_p, binary_mis, binary_un =  ht.create_paired_matrix(binary1, binary2)

        if binary_p.shape[0] < len(id1) * 0.1:
            print(("\nWARNING: Less than 10%% of reads could be paired. Consider an appropriate unpaired_weight setting "
             "in your config file (currently %.3f), because you may need to resort to using unpaired reads.") % unpaired_weight)

        if unpaired_weight > 0:
            if use_discordant:
                binary = pd.concat([binary_p, binary_un, binary_mis])
            else:
                binary = pd.concat([binary_p, binary_un])
        else:
            binary = binary_p

    else:
        pos, read_details = ht.pysam_to_hdf(bam_paths[0])

        if not bam_input and config.getboolean('behavior', 'deletebam'):
            os.remove(bam_paths[0])

        binary = np.sign(pos)  # dtype=np.uint16

    # dimensionality reduction and typing

    alleles_to_keep = list(filter(is_frequent, binary.columns))
    binary = binary[alleles_to_keep]

    if VERBOSE:
        print("\n", ht.now(), 'temporary pruning of identical rows and columns')
    unique_col, representing = ht.prune_identical_alleles(binary, report_groups=True)
    representing_df = pd.DataFrame([[a1, a2] for a1, a_l in representing.items() for a2 in a_l],
                                   columns=['representative', 'represented'])

    temp_pruned = ht.prune_identical_reads(unique_col)

    if VERBOSE:
        print("\n", ht.now(), 'Size of mtx with unique rows and columns:', temp_pruned.shape)
        print(ht.now(), 'determining minimal set of non-overshadowed alleles')

    minimal_alleles = ht.prune_overshadowed_alleles(temp_pruned)

    if VERBOSE:
        print("\n", ht.now(), 'Keeping only the minimal number of required alleles', minimal_alleles.shape)

    binary = binary[minimal_alleles]

    if VERBOSE:
        print("\n", ht.now(), 'Creating compact model...')

    if is_paired and unpaired_weight > 0:
        if use_discordant:
            compact_mtx, compact_occ = ht.get_compact_model(binary_p[minimal_alleles],
                pd.concat([binary_un, binary_mis])[minimal_alleles], weight=unpaired_weight)
        else:
            compact_mtx, compact_occ = ht.get_compact_model(binary_p[minimal_alleles],
                binary_un[minimal_alleles], weight=unpaired_weight)
    else:
        compact_mtx, compact_occ = ht.get_compact_model(binary)

    allele_ids = binary.columns

    groups_4digit = defaultdict(list)
    for allele in allele_ids:
        type_4digit = get_4digit(allele)
        groups_4digit[type_4digit].append(allele)

    sparse_dict = ht.mtx_to_sparse_dict(compact_mtx)
    threads = get_num_threads(config.getint("ilp", "threads"))               
    if VERBOSE:
        print("\nstarting ilp solver with %s threads..." % threads)
        print("\n", ht.now(), 'Initializing OptiType model...')

    op = OptiType(sparse_dict, compact_occ, groups_4digit, table, args.beta, 2,
                  config.get("ilp", "solver"), threads, verbosity=VERBOSE)
    result = op.solve(args.enumerate)

    if VERBOSE:
        print("\n", ht.now(), 'Result dataframe has been constructed...')

    result_4digit = result.applymap(get_types)
    for iii in ["A1", "A2", "B1", "B2", "C1", "C2"]:
        if not iii in result_4digit:
            result_4digit[iii] = None
    r = result_4digit[["A1", "A2", "B1", "B2", "C1", "C2", "nof_reads", "obj"]]

    # write CSV to out. And generate plots
    r.to_csv(out_csv, sep="\t",
                         columns=["A1", "A2", "B1", "B2", "C1", "C2", "nof_reads", "obj"],
                         header=["A1", "A2", "B1", "B2", "C1", "C2", "Reads", "Objective"])
    
    hlatype = result.iloc[0].reindex(["A1", "A2", "B1", "B2", "C1", "C2"]).drop_duplicates().dropna()
    features_used = [('intron', 1), ('exon', 2), ('intron', 2), ('exon', 3), ('intron', 3)] \
                     if not args.rna else [('exon',2),('exon',3)]
    plot_variables = [pos, read_details, pos2, read_details2, (binary_p, binary_un, binary_mis)] if is_paired else [pos, read_details]
    coverage_mat = ht.calculate_coverage(plot_variables, features, hlatype, features_used)
    ht.plot_coverage(out_plot, coverage_mat, table, features, features_used)
