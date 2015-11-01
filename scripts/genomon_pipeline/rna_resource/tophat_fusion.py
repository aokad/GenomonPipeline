#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class TopHat_fusion(Stage_task):

    task_name = "tophat_fusion"

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

export PATH={bowtie_path}:$PATH
export PATH={samtools_path}:$PATH
export PATH={blast_path}:$PATH

cd {output_dir}

ln -s {tophatfusion_db}/refGene.txt ./
ln -s {tophatfusion_db}/ensGene.txt ./
ln -s {tophatfusion_db}/blast ./

sample_dir_name=`basename {tophat_result}`
ln -s {tophat_result} tophat_$sample_dir_name

{tophat2_dir}/tophat-fusion-post {additional_params} {bowtie1_db} 
"""

    def __init__(self, qsub_option, script_dir):
        super(TopHat_fusion, self).__init__(qsub_option, script_dir)
