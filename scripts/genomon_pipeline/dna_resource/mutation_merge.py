#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Mutation_merge(Stage_task):

    task_name = "mutation_merge"

    script_template = """
#!/bin/bash
#
#$ -S /bin/bash         # set shell in UGE
#$ -cwd                 # execute at the submitted dir
#$ -e {log}             # log file directory
#$ -o {log}             # log file directory
pwd                     # print current working directory
hostname                # print hostname
date                    # print date
set -xv


echo -n > {out_prefix}_genomon_mutations.result.txt || exit $?
for i in `seq 1 1 {filecount}`
do
    if [ _{active_annovar_flag} = _True ]
    then
        awk 'NR>1 {{print}}' {out_prefix}_mutations_candidate.${{i}}.hg19_multianno.txt >> {out_prefix}_genomon_mutations.result.txt || exit $?
    else
        cat {out_prefix}_mutations_candidate.${{i}}.hg19_multianno.txt >> {out_prefix}_genomon_mutations.result.txt || exit $?
    fi
done

"""
    def __init__(self, qsub_option, script_dir):
        super(Mutation_merge, self).__init__(qsub_option, script_dir)


