#! /usr/bin/env python

from genomon_qc.stage_task import *

class Res_QC_Bamstats(Stage_task):

    task_name = "qc_bamstats"

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

export PERL5LIB={PERL5LIB}

{bamstats} -i {input} -o {output}.tmp || exit $?
mv {output}.tmp {output}
"""

    def __init__(self, qsub_option, script_dir):
        super(Res_QC_Bamstats, self).__init__(qsub_option, script_dir)
