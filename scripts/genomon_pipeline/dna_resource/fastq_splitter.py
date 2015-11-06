#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Fastq_splitter(Stage_task):

    task_name = "fastq_splitter"

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

if [ -f {target_dir}/1_${{SGE_TASK_ID}}.gz ]; then
    zcat {target_dir}/1_${{SGE_TASK_ID}}.gz | split -a 4 -d -l {lines} - {target_dir}/${{SGE_TASK_ID}}_
    status=("${{PIPESTATUS[@]}}")
    [ ${{PIPESTATUS[0]}} -ne 0 ] || echo ${{PIPESTATUS[0]}}
else
    split -a 4 -d -l {lines} {target_dir}/1_${{SGE_TASK_ID}}.fastq {target_dir}/${{SGE_TASK_ID}}_ || exit $?
fi

ls -1 {target_dir}/${{SGE_TASK_ID}}_[0-9][0-9][0-9][0-9] | while read filename; do
    mv $filename $filename.fastq_split || exit $?
done


"""

    def __init__(self, qsub_option, script_dir):
        super(Fastq_splitter, self).__init__(qsub_option, script_dir)


