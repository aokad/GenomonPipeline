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

input_files=`ls {target_dir}/*_${{SGE_TASK_ID}}{ext}`

if [ "_{fastq_filter}" = "_True" ]; then

    if [ "_{ext}" = "_.gz" ]; then
        gzip -dc $input_files | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v '^--$' | split -a 4 -d -l {lines} - {target_dir}/${{SGE_TASK_ID}}_ || exit $?
    else
        cat $input_files | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v '^--$' | split -a 4 -d -l {lines} - {target_dir}/${{SGE_TASK_ID}}_ || exit $?
    fi
else
    if [ "_{ext}" = "_.gz" ]; then
        gzip -dc $input_files | split -a 4 -d -l {lines} - {target_dir}/${{SGE_TASK_ID}}_ || exit $?
    else
        cat $input_files | split -a 4 -d -l {lines} - {target_dir}/${{SGE_TASK_ID}}_ || exit $?
    fi
fi

ls -1 {target_dir}/${{SGE_TASK_ID}}_[0-9][0-9][0-9][0-9] | while read filename; do
    mv $filename $filename.fastq_split || exit $?
done


"""

    def __init__(self, qsub_option, script_dir):
        super(Fastq_splitter, self).__init__(qsub_option, script_dir)


