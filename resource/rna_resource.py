#!/usr/bin/python

#
# Job scripts for rna_pipeline.py
#
tophat2 = """
#!/bin/bash
#
#  Copyright Human Genome Center, Institute of Medical Science, the University of Tokyo
#  @since 2012
#
# Set SGE
#
#$ -S /bin/bash         # set shell in UGE
#$ -cwd                 # execute at the submitted dir
#$ -e {log}             # log file directory
#$ -o {log}             # log file directory
pwd                     # print current working directory
hostname                # print hostname
date                    # print date
set -xv

{tophat2} -p 8 \
          -o {output_file} \
          -G {gtf_file} \
         {bowtie2_dataabes} \
         {fastq_file}          
"""


cufflinks = """
#!/bin/bash
#
#  Copyright Human Genome Center, Institute of Medical Science, the University of Tokyo
#  @since 2012
#
# Set SGE
#
#$ -S /bin/bash         # set shell in UGE
#$ -cwd                 # execute at the submitted dir
#$ -e {log}             # log file directory
#$ -o {log}             # log file directory
pwd                     # print current working directory
hostname                # print hostname
date                    # print date
set -xv

cufflinks -p 8 \
          -o {output_dir} \
          -G {gtf_file} \
          {bam_file}
"""

cuffdiff = """
#!/bin/bash
#
#  Copyright Human Genome Center, Institute of Medical Science, the University of Tokyo
#  @since 2012
#
# Set SGE
#
#$ -S /bin/bash         # set shell in UGE
#$ -cwd                 # execute at the submitted dir
#$ -e {log}             # log file directory
#$ -o {log}             # log file directory
pwd                     # print current working directory
hostname                # print hostname
date                    # print date
set -xv

cuffdiff -p $THREAD_NUM \
         {gtf_file} \
         -o {output_dir} \
         -L {data_label} \
         {bam_files}
"""

