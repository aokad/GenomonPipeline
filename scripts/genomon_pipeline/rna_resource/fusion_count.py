#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Fusion_count(Stage_task):

    task_name = "fusion_count"

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

# set python environment
export PYTHONHOME={pythonhome}
export PATH=$PYTHONHOME/bin:$PATH
export LD_LIBRARY_PATH={ld_library_path}
export PYTHONPATH={pythonpath}

{chimera_utils} count {additional_params} {chimeric_sam} {output}
"""

    def __init__(self, qsub_option, script_dir):
        super(Fusion_count, self).__init__(qsub_option, script_dir)
