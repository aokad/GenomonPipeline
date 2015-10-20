#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Star_fusion(Stage_task):

    task_name = "star_fusion"

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

export PERL5LIB={environment_variables}
{star_fusion} \
    --chimeric_out_sam {chimeric_sam} \
    --chimeric_junction {chimeric_junction} \
    --ref_GTF {gtf_file} \
    --out_prefix {out_prefix} \
    {additional_params}
"""

    def __init__(self, qsub_option, script_dir):
        super(Star_fusion, self).__init__(qsub_option, script_dir)
