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

file_timestamp_format = "{name}_{year:0>4d}{month:0>2d}{day:0>2d}_{hour:0>2d}{min:0>2d}_{msecond:0>6d}"

splitfile = """
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

split -a {suffix_len} -d -l {lines_per_file} {input_file} {output_prefix}_
for FILE in `ls {output_prefix}_*`
do
    mv $FILE $FILE{output_suffix}
done

"""

cutadapt = """
#!/bin/bash
#
#  Copyright Human Genome Center, Institute of Medical Science, the University of Tokyo
#  @since 2012
#
#$ -S /bin/bash
#$ -cwd
#$ -e {log}             # log file directory
#$ -o {log}             # log file directory
pwd                     # print current working directory
hostname                # print hostname
date                    # print date
set -xv

#infastq=$1
#outfastq=$2
#tmpoutfastq=$3
#casavacoDE=$4
#adapters=$5
#cutadapt=$6
#scriptdir=$7


source {scriptdir}/utility.sh

# sleep 
sh {scriptdir}/sleep.sh

ADAPTERS=`echo "{adapters}" | sed -e 's/,/ /g'`
for adapter in ${adapters}
do
    optadapters="${optadapters}"" ""-a ${adapter}"
done

echo "{cutadapt} ${optadapters} {infastq} > {tmpoutfastq}"
{cutadapt} ${optadapters} {infastq} > {tmpoutfastq}
check_error $?

echo "{scriptdir}/fastqNPadding.pl {casavacode} {tmpoutfastq} > {outfastq}"
{scriptdir}/fastqNPadding.pl {casavacode} {tmpoutfastq} > {outfastq}
check_error $?

"""



bwa_mem = """
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

{samtools} merge \
    -nr \
    {output_bam_file} \
    {input_bam_files}
"""

fisher_mutation_call = """
#!/bin/bash
#
#  Copyright Human Genome Center, Institute of Medical Science, the University of Tokyo
#  @since 2012
#
#$ -S /bin/bash
#$ -cwd
#$ -e {log}             # log file directory
#$ -o {log}             # log file directory
pwd                     # print current working directory
hostname                # print hostname
date                    # print date
set -xv


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

