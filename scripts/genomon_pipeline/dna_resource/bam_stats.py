#! /usr/bin/env python

from genomon_summary.stage_task import *

class Res_Bamstats(Stage_task):

    task_name = "bam_stats"

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

export PERL5LIB={PERL5LIB}

{PCAP}/bin/bam_stats.pl -i {input} -o {output}.tmp
mv {output}.tmp {output}
"""

    def __init__(self, qsub_option, script_dir):
        super(Res_Bamstats, self).__init__(qsub_option, script_dir)
