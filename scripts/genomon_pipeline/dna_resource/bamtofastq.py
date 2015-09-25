#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Bam2Fastq(Stage_task):

    task_name = "bam2fastq"

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

{biobambam}/bamtofastq collate=1 exclude=QCFAIL,SECONDARY,SUPPLEMENTARY filename={input_bam} F={f1_name} F2={f2_name} T={t} S={s} O={o1_name} 2={o2_name}

"""

    def __init__(self, qsub_option, script_dir):
        super(Bam2Fastq, self).__init__(qsub_option, script_dir)


