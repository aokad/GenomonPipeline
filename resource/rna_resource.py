#!/usr/bin/python

#
# Job scripts for rna_pipeline.py
#
tophat2 = """
#!/bin/bash
#
#  Copyright Human Genome Center, Institute of Medical Science, the University of Tokyo
#  @since 2012
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

PATH=$PATH:{bowtie_path}

{tophat2} -p 8 \
          -o {output_dir} \
          -G {ref_gtf} \
         {bowtie2_database} \
         {input_fastq}          
"""


cufflinks = """
#!/bin/bash
#
#  Copyright Human Genome Center, Institute of Medical Science, the University of Tokyo
#  @since 2012
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

{cufflinks} -p 8 \
          -o {output_dir} \
          -G {ref_gtf} \
          {bam_file}
"""

cuffdiff_merge = """
#!/bin/bash
#
#  Copyright Human Genome Center, Institute of Medical Science, the University of Tokyo
#  @since 2012
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
#
# <parameters>
# log
# env_variables
# array_data
# merge_bam_flag
# input_bam_files
# merge_bam_flag
# samtools
# biobambam
#
echo SGE_TASK_ID:$SGE_TASK_ID
echo SGE_TASK_FIRST:$SGE_TASK_FIRST
echo SGE_TASK_LAST:$SGE_TASK_LAST
echo SGE_TASK_STEPSIZE:$SGE_TASK_STEPSIZE

{env_variables}

{array_data}

OUTPUT_BAM_PREFIX=`echo {merged_bam_file} | sed 's/\.[^\.]\+$//'`
mkdir -p {out_dir}
if [ "{merge_bam_flag}" = "True" ]
then
    if [ "{use_biobambam}" = "True" ]
    then
        {biobambam}/bammerge  \
                        {input_bam_files} \
                        md5filename="$OUT_BAM_PREFIX".metrics \
                        tmpfile="$OUT_BAM_PREFIX".tmp \
                        indexfilename={merged_bam_file}.bai \
                        md5=1 \
                        index=1 \
                        > {merged_bam_file}
    else
        {samtools} merge \
                    "$OUTPUT_BAM_PREFIX".bam \
                    {input_bam_files};
        {samtools} index {merged_bam_file};
    fi
else
    cp {input_bam_files} {merged_bam_file}
fi
"""

cuffdiff = """
#!/bin/bash
#
#  Copyright Human Genome Center, Institute of Medical Science, the University of Tokyo
#  @since 2012
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

{env_variables}

{cuffdiff} -p 4\
         {ref_gtf} \
         -o {output_dir} \
         -L {data_labels} \
         {merged_control_bam_file} {merged_disease_bam_file}

"""

cummeRbund = """
#!/bin/bash
#
#  Copyright Human Genome Center, Institute of Medical Science, the University of Tokyo
#  @since 2012
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

{env_variables}

{R} -q --vanilla --args {sample_name} {input_dir} {output_dir} < {script_dir}/cummeRbund.R

"""

