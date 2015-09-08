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

split -a 3 -d -l {lines} {input_file} {out_name_stem}/{pair_id}_ || exit $?

ls -1 {out_name_stem}/{pair_id}_[0-9][0-9][0-9] | while read filename; do
    mv $filename $filename.fastq_split || exit $?
done

"""

    def __init__(self, qsub_option, script_dir):
        super(Fastq_splitter, self).__init__(qsub_option, script_dir)


