#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Bwa_align(Stage_task):

    task_name = "bwa_align"

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

{bwa} mem {bwa_params} {ref_fa} {fastq1} {fastq2}  > {sam} || exit $?

{biobambam}/bamsort index=1 level=1 inputthreads=2 outputthreads=2 calmdnm=1 calmdnmrecompindentonly=1 calmdnmreference={ref_fa} tmpfile={bam}.tmp inputformat=sam indexfilename={bam}.bai I={sam} O={bam}

"""

    def __init__(self, qsub_option, script_dir):
        super(Bwa_align, self).__init__(qsub_option, script_dir)


