#!/usr/bin/python

#
# Job scripts for rna_pipeline.py
#
star_genome = """
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

{star} \
    --runMode genomeGenerate \
    --genomeDir {out_dir} \
    --genomeFastaFiles {ref_fasta} \
    --sjdbGTFfile {ref_gtf} \
    {additional_params} \

"""

extract_fastq = """
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

echo SGE_TASK_ID:$SGE_TASK_ID
echo SGE_TASK_FIRST:$SGE_TASK_FIRST
echo SGE_TASK_LAST:$SGE_TASK_LAST
echo SGE_TASK_STEPSIZE:$SGE_TASK_STEPSIZE

{array_data}

case {input_file} in
*\.gz)
    gzip -dc {input_file} >> {output_file}
    ;;

*\.bz2)

    bzip2 -dc  {input_file} >> {output_file}
    ;;
esac

"""


star_map = """
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

{star} \
    --genomeDir {star_genome} \
    --readFilesIn {fastq1} {fastq2} \
    --outFileNamePrefix {out_prefix} \
    {additional_params}

"""


star_fusion = """
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

export PERL5LIB={environment_variables}

{star_fusion} \
    --chimeric_out_sam {chimeric_sam} \
    --chimeric_junction {chimeric_junction} \
    --ref_GTF {gtf_file} \
    --out_prefix {out_prefix} \
    {additional_params}

"""

