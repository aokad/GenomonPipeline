#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class TopHat2_align(Stage_task):

    task_name = "tophat2_align"

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

export PATH={bowtie_path}:$PATH
export PATH={samtools_path}:$PATH

{tophat2} -o {output_dir} \
          --GTF {ref_gtf} \
         {additional_params} \
         {bowtie2_database} \
         {fastq1} \
         {fastq2}     
"""

    def __init__(self, qsub_option, script_dir):
        super(TopHat2_align, self).__init__(qsub_option, script_dir)
