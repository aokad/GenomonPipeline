#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class SV_parse(Stage_task):

    task_name = "sv_parse"

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

# set python environment
export PYTHONHOME={pythonhome}
export PATH=$PYTHONHOME/bin:$PATH
export LD_LIBRARY_PATH={ld_library_path}
export PYTHONPATH={pythonpath}

{genomon_sv} parse {sample_conf} {param_conf}

"""

    def __init__(self, qsub_option, script_dir):
        super(SV_parse, self).__init__(qsub_option, script_dir)


