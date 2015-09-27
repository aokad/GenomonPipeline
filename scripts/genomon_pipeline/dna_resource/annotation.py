#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Annotation(Stage_task):

    task_name = "annotation"

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


"""

    def __init__(self, qsub_option, script_dir):
        super(Annotation, self).__init__(qsub_option, script_dir)


