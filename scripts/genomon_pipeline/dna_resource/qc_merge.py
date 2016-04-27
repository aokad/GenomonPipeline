#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Res_QC_Merge(Stage_task):
    
    task_name = "qc_merge"

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

# cat {input1} {input2} {input3} {input4} > {output}
"""

    def __init__(self, qsub_option, script_dir):
        super(Res_QC_Merge, self).__init__(qsub_option, script_dir)
        
    def write_qc(self, files, output, meta):

        # read .bamstat
        f = open(files[0])
        text = f.read()
        f.close()
        bamstat = text.split("\n")
        
        if len(bamstat) < 2:
            bamstat = ["",""]
            
        # read .coverage
        f = open(files[1])
        text = f.read()
        f.close()
        coverage = text.split("\n")
        
        if len(coverage) < 2:
            coverage = ["",""]
            
        f = open(output, "w")
        f.write(meta + "\n")
        f.write(bamstat[0] + "\t" + coverage[0] + "\n")    # header
        f.write(bamstat[1] + "\t" + coverage[1])    # data
        f.close()
   
