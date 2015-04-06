#!/usr/bin/python

#
# Job scripts for wgs_pipeline.py
#
bamtofastq_s = """
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

{bamtofastq}  exclude=QCFAIL,SECONDARY,SUPPLEMENTARY\
            T={tmpfastq}\
            F={outfastq1}\
            collate=1\
            filename={bamfile}
"""

bamtofastq_p = """
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

{bamtofastq}  exclude=QCFAIL,SECONDARY,SUPPLEMENTARY\
            T={tmpfastq}\
            F={outfastq1}\
            F2={outfastq2}\
            collate=1\
            filename={bamfile}
"""
extract_gz = """
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

gzip -dc {input_file} >> {output_file}

"""


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

echo SGE_TASK_ID:$SGE_TASK_ID
echo SGE_TASK_FIRST:$SGE_TASK_FIRST
echo SGE_TASK_LAST:$SGE_TASK_LAST
echo SGE_TASK_STEPSIZE:$SGE_TASK_STEPSIZE

{array_data}

if [ "{file_ext}" = ".gz" ]
then
    zcat {input_file} | split -a {suffix_len} -d -l {lines_per_file} - {output_prefix}
elif [ "{file_ext}" = ".bz2" ]
then
    bzip2 -dc  {input_file} | split -a {suffix_len} -d -l {lines_per_file} - {output_prefix}
else
    split -a {suffix_len} -d -l {lines_per_file} {input_file} {output_prefix}
fi

for FILE in `ls {output_prefix}*`
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

echo SGE_TASK_ID:$SGE_TASK_ID
echo SGE_TASK_FIRST:$SGE_TASK_FIRST
echo SGE_TASK_LAST:$SGE_TASK_LAST
echo SGE_TASK_STEPSIZE:$SGE_TASK_STEPSIZE

#infastq
#outfastq
#tmpoutfastq
#casavacode
#optadapters
#cutadapt
#scriptdir

source {scriptdir}/utility.sh

# sleep 
sh {scriptdir}/sleep.sh

{array_data}

echo "{cutadapt} {optadapters} {infastq} > {tmpoutfastq}"
{cutadapt} {optadapters}  {infastq} > {tmpoutfastq}
check_error $?

echo "{scriptdir}/fastqNPadding.pl {casavacode} {tmpoutfastq} > {outfastq}"
{scriptdir}/fastqNPadding.pl {casavacode} {tmpoutfastq} > {outfastq}
check_error $?

"""



bwa_mem = r"""
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

TMP_BAM=`echo {bam} | sed 's/\.bam/_unsorted.bam/'`
BAM_WITHOUT_SUFFIX=`echo {bam} | sed 's/\.bam//'`

{bwa} mem \
    -t 2 \
    -T {min_score} \
    -R '@RG\tID:{rg_id}\tSM:{sample_desc}\tLB:{library}\tPL:{platform}\tPU:{platform_unit}\tCN:{seq_center}\tPI:{pred_med_insert}' \
    {hg19_fa} \
    {fastq1} \
    {fastq2} | \
{samtools} view -Sb - \
    > $TMP_BAM;

{samtools} sort $TMP_BAM $BAM_WITHOUT_SUFFIX;

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

OUTPUT_BAM_PREFIX=`echo {output_bam_file} | cut -d '.' -f1`

{samtools} merge \
    -nr \
    "$OUTPUT_BAM_PREFIX".bam \
    {input_bam_files};

{samtools} index {output_bam_file};
"""

fisher_mutation_call = \
"""
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

CONTROL_FILTERED_BAM=`echo {control_input_bam} | cut -d '.' -f1`_filtered.bam
DISEASE_FILTERED_BAM=`echo {disease_input_bam} | cut -d '.' -f1`_filtered.bam

if [ -f {control_input_bam} -a ! -f $CONTROL_FILTERED_BAM ]
then
    python {script_dir}/bamfilter.py \
               --ref_fa {ref_fa} \
               --max_indel {max_indel} \
               --max_distance {max_distance} \
               --map_quality {map_quality} \
            {control_input_bam} \
            $CONTROL_FILTERED_BAM
fi

if [ -f {disease_input_bam} -a ! -f $DISEASE_FILTERED_BAM ]
then
    python {script_dir}/bamfilter.py \
               --ref_fa {ref_fa} \
               --max_indel {max_indel} \
               --max_distance {max_distance} \
               --map_quality {map_quality} \
            {disease_input_bam} \
            $DISEASE_FILTERED_BAM
fi

if [ -f $CONTROL_FILTERED_BAM -a -f $DISEASE_FILTERED_BAM -a ! -f {output_txt} ]
then
    samtools mpileup \
                -BQ0 \
                -d 10000000 \
                -f {ref_fa} \
                $CONTROL_FILTERED_BAM \
                $DISEASE_FILTERED_BAM |\
    python {script_dir}/fisher.py \
               --output {output_txt}\
               --ref_fa {ref_fa}\
               --base_quality {base_quality} \
               --mismatch_rate {mismatch_rate}\
               --min_depth {min_depth}

elif [ -f $CONTROL_FILTERED_BAM -a ! -f {output_txt} ]
then
    samtools mpileup \
                -BQ0 \
                -d 10000000 \
                -f {ref_fa} \
                $CONTROL_FILTERED_BAM |\
    python {script_dir}/fisher.py \
               --output {output_txt} \
               --ref_fa {ref_fa} \
               --base_quality {base_quality} \
               --mismatch_rate {mismatch_rate} \
               --min_depth {min_depth}

else           
    exit 1
fi

"""

