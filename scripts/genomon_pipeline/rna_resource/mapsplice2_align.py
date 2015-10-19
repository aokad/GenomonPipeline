#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Mapsplice2_align(Stage_task):

    task_name = "mapsplice2_align"

    script_template = """
#!/bin/bash
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

{python} {mapsplice2} --fusion-non-canonical -p 16 --bam --gene-gtf {ref_gtf} -c {ref_fasta} -x {bow_ind_mp2} -1 {fastq1} -2 {fastq2} -o {out_dir} 
"""

    def __init__(self, qsub_option, script_dir):
        super(Mapsplice2_align, self).__init__(qsub_option, script_dir)
