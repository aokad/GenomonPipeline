#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Res_PA_Plot(Stage_task):

    task_name = "paplot"

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
# set python environment
export PYTHONHOME={pythonhome}
export PATH=$PYTHONHOME/bin:$PATH
export PYTHONPATH={pythonpath}

{command} 
"""

    qc_template = """
{paplot} qc '{inputs}' {output_dir} {title} --config_file {config_file}
"""
    sv_template = """
{paplot} ca '{inputs}' {output_dir} {title} --config_file {config_file}
"""
    mutation_template = """
if test '{annovar}' == 'True'; then
  if test '{inputs_mutation}' != ''; then
    {paplot} mutation '{inputs}' {output_dir} {title} --config_file {config_file}
  fi
else
  echo 'paplot: [annotation] active_annovar_flag = False in genomon_conf_file, skip mutation-matrix.'
fi
"""
    full_template = """
{paplot} signature {input} {output_dir} {title} {signature_num} --config_file {config_file}
"""
    ind_template = """
{paplot} pmsignature {input} {output_dir} {title} {pmsignature_num} --config_file {config_file}
"""

    index_template = """
{paplot} index {output_dir} --config_file {config_file} --remarks '{remarks}'
"""

    def __init__(self, qsub_option, script_dir):
        super(Res_PA_Plot, self).__init__(qsub_option, script_dir)

        
