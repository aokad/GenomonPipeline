#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Markduplicates(Stage_task):

    task_name = "markduplicates"

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

{biobambam}/bammarkduplicates M={out_prefix}.metrics tmpfile={out_prefix}.tmp markthreads=2 rewritebam=1 rewritebamlevel=1 index=1 md5=1 {input_bam_files} O={out_bam}


"""

    def __init__(self, qsub_option, script_dir):
        super(Markduplicates, self).__init__(qsub_option, script_dir)


