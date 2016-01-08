#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Res_PA_Merge(Stage_task):

    task_name = "post_analysis_merge"

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

{merge_result}
{merge_igv}
"""

    script_template_result = """
genomon_pa merge {mode} {input_files} {output_file} --config_text "{config}"
"""

    script_template_igv = """
genomon_pa merge igv {input_files} {output_file} --config_text "{config}"
"""

    script_template_config = "[result_format_{mode}] sept={sept} header={header} suffix={suffix} [merge_format_{mode}] filters={filters}"
    
    def __init__(self, qsub_option, script_dir):
        super(Res_PA_Merge, self).__init__(qsub_option, script_dir)
        
