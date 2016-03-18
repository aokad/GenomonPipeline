#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Res_PostAnalysis(Stage_task):

    task_name = "post_analysis"

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
export LD_LIBRARY_PATH={ld_library_path}:$LD_LIBRARY_PATH
export PYTHONPATH={pythonpath}

{genomon_pa} run {mode} {output_dir} {genomon_root} --config_file {config_file} --input_file "{input_file}"
"""

    def __init__(self, qsub_option, script_dir):
        super(Res_PostAnalysis, self).__init__(qsub_option, script_dir)
        

    def list_to_string(self, li):
        txt = ""
        for i in li:
            if len(txt) > 0:
                txt += ","
            txt += i
            
        txt += ""
        return txt
