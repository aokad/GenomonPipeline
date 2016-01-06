#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Bam_index(Stage_task):

    task_name = "bam_index"

    script_template = """
#!/bin/bash
#
# Set SGE
#
#$ -S /bin/bash         # set shell in UGE
#$ -cwd                 # execute at the submitted dir
pwd                     # print current working directory
hostname                # print hostname
date                    # print date
set -xv

{samtools} index {input_bam}
"""

    def __init__(self, qsub_option, script_dir):
        super(Bam_index, self).__init__(qsub_option, script_dir)
