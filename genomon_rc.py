#!/usr/bin/python
"""

Genome Analysis Pipeline
resource file


"""
#
# General
#
date_format = "{year:0>4d}{month:0>2d}{day:0>2d}"

#
# Commands
#
qsub_cmd = "qsub -sync yes -now no -l {job_type},s_vmem={s_vmem},mem_req={mem_req} {cmd}"

shell_script_format = "{name}_{year:0>4d}{month:0>2d}{day:0>2d}_{hour:0>2d}{min:0>2d}_{msecond:0>6d}"

splitfile = \
"""
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

split -a {suffix_len} -d -l {lines_per_file} {input_file} {output_prefix}_
for FILE in `ls {output_prefix}_*`
do
    mv $FILE $FILE{output_suffix}
done

"""

bwa_mem = \
"""
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

{bwa} mem \
    -t 1 \
    {hg19_fa} \
    {fastq1} \
    {fastq2} | \
{samtools} view -Sb - \
    > {bam}
"""

merge_bam = \
"""
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

{samtools} merge \
    -nr \
    {output_bam_file} \
    {input_bam_files}
"""


#
# Misc
#
dir_tree_resource = \
"""
project_directory:
    data:
        data_date:
            sample_name
    results:
        run_date:
            data_date_sample_name:
                - config
                - script
                - log
                - out:
                    - fastq
                    - bam
                    - annovar
                    - fisher
                    - cnv
                    - fusion
                    - sv
"""

end_dir_list = ( 'config', 'script', 'log', 'fastq', 'bam', 'annovar', 'fisher', 'cnv', 'fusion', 'sv' )

