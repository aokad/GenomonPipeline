#! /usr/bin/env python

class Star_align(Stage_task):

    self.task_name = "bwa_mem"
    self.script_resource = """
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

{star} \
    --genomeDir {star_genome} \
    --readFilesIn {fastq1} {fastq2} \
    --outFileNamePrefix {out_prefix} \
    {additional_params}
"""

    def __init__(self, qsub_option, script_dir):
        super(Star_align, self).__init__(qsub_option, script_dir)

