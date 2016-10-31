#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Res_PrePmsignature(Stage_task):

    task_name = "pre_pmsignature"

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
export PATH={r_path}:$PATH
export LD_LIBRARY_PATH={r_ld_library_path}:$LD_LIBRARY_PATH
export R_LIBS={r_libs}
export R_PATH={r_path}

if [ -e {output_file}.tmp ]; then
    rm {output_file}.tmp || exit $?
fi

for input_file in {input_files}
do
    $R_PATH/R --vanilla --slave --args $input_file {output_file}.tmp < {script_path}/pmsignature/cut_mutation.R  || exit $?
done

mv {output_file}.tmp {output_file}
"""

    def __init__(self, qsub_option, script_dir):
        super(Res_PrePmsignature, self).__init__(qsub_option, script_dir)

