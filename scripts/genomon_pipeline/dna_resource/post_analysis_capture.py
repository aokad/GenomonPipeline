#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Res_PA_Capture(Stage_task):

    task_name = "post_analysis_capture"

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
export LD_LIBRARY_PATH={ld_library_path}:$LD_LIBRARY_PATH
export PYTHONPATH={pythonpath}

{bam_script}
{igv_script}
"""

    script_template_bam = """
# bam trimming
genomon_pa bam {mode} {input_file} {id} {output_file} {output_bam_dir} {output_log_dir} {genomon_root} --config_text "{config}"
bash {output_file}
"""

    script_template_igv = """
# igv capture
genomon_pa igv {mode} {input_file} {id} {output_file} {output_dir} {genomon_root} --config_text "{config}"
"""

    script_template_config = "[capture] capture_max={capture_max} capture_width={capture_width} [pickup] pickup_width={pickup_width} markdup_bam_suffix={markdup_bam_suffix} pickup_bam_suffix={pickup_bam_suffix} [result_format_{mode}] sept={sept} header={header} col_pos_chr1={col_pos_chr1} col_pos_start={col_pos_start} col_pos_chr2={col_pos_chr2} col_pos_end={col_pos_end} [SOFTWARE] samtools={samtools} bedtools={bedtools}"

    def __init__(self, qsub_option, script_dir):
        super(Res_PA_Capture, self).__init__(qsub_option, script_dir)
        
