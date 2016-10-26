#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Res_Pmsignature(Stage_task):

    task_name = "pmsignature"

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

# export PATH=/usr/local/package/r/current3/bin:${PATH}
# export LD_LIBRARY_PATH=/usr/local/package/r/current3/lib64/R/lib:${LD_LIBRARY_PATH}
# export R_LIBS=/home/aiokada/.R_3_2_5
# export R_PATH=/usr/local/package/r/current3/bin

export PATH={R_PATH}:$PATH
export LD_LIBRARY_PATH={R_LD_LIBRARY_PATH}:$LD_LIBRARY_PATH
export R_LIBS={R_LIBS}
export R_PATH={R_PATH}

sig_list=({sig_list})
sig_num=${{sig_list[$SGE_TASK_ID]}}

{command} 
"""

    ind_template = """
${R_PATH}/R --vanilla --slave --args ${INPUTFILE} ${OUTPUTFILE} ${SIGNUM} ${TRDIRFLAG} ${TRIALNUM} < perform_signature.R
{R_PATH}/R --vanilla --slave --args {INPUTFILE} {OUTPUTFILE} $sig_num {TRDIRFLAG} {TRIALNUM} < {SCRIPT_PATH}/pmsignature/perform_ind.R
"""

    full_template = """
{R_PATH}/R --vanilla --slave --args {INPUTFILE} {OUTPUTFILE} {SIGNUM} {TRDIRFLAG} {TRIALNUM} < {SCRIPT_PATH}/pmsignature/perform_full.R
"""

    def __init__(self, qsub_option, script_dir):
        super(Res_Pmsignature, self).__init__(qsub_option, script_dir)

