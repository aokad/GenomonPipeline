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

{genomon_pa} run {mode} {output_dir} {genomon_root} {sample_sheet} \
--config_file {config_file} \
--samtools {samtools} --bedtools {bedtools} \
--input_file_case1 "{input_file_case1}" \
--input_file_case2 "{input_file_case2}" \
--input_file_case3 "{input_file_case3}" \
--input_file_case4 "{input_file_case4}"
"""

    def __init__(self, qsub_option, script_dir):
        super(Res_PostAnalysis, self).__init__(qsub_option, script_dir)
        
    def sample_split_case(self, li, unpair=True, unpanel=True):
        tmr_nrml_list = []
        tmr_nrml_none = []
        tmr_none_list = []
        tmr_none_none = []
    
        for item in li:
            if item[1]== None:
                if unpair == True:
                    if item[2] == None:
                        if unpanel == True: tmr_none_none.append(item[0])
                    else:
                        tmr_none_list.append(item[0])
            else:
                if item[2] == None:
                    if unpanel == True: tmr_nrml_none.append(item[0])
                else:
                    tmr_nrml_list.append(item[0])
        
        return {"case1": tmr_nrml_list, "case2": tmr_nrml_none, "case3": tmr_none_list, "case4": tmr_none_none}

            
