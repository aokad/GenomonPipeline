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

case {input_file} in
*\.gz)
    if [ "{fastq_filter}" = "True" ]
    then
        zcat {input_file} | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v '^--$' | split -a {suffix_len} -d -l {lines_per_file} - {output_prefix}
    else
        zcat {input_file} | split -a {suffix_len} -d -l {lines_per_file} - {output_prefix}
    fi
    ;;

*\.bz2)
    if [ "{fastq_filter}" = "True" ]
    then
        bzip2 -dc  {input_file} | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v '^--$' | split -a {suffix_len} -d -l {lines_per_file} - {output_prefix}
    else
        bzip2 -dc  {input_file} | split -a {suffix_len} -d -l {lines_per_file} - {output_prefix}
    fi
    ;;

*)
    if [ "{fastq_filter}" = "True" ]
    then
        cat {input_file} | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v '^--$' | split -a {suffix_len} -d -l {lines_per_file} - {output_prefix}
    else
        split -a {suffix_len} -d -l {lines_per_file} {input_file} {output_prefix}
    fi
    ;;
esac

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

UNSORTED_BAM=`echo {bam} | sed 's/\.bam/_unsorted.bam/'`
BAM_WITHOUT_SUFFIX=`echo {bam} | sed 's/\.bam//'`

if [ ! -e {ref_fa}.fai ]
then
    {samtools} faidx {ref_fa}
    if [ $1 -ne 0 ]
    then
        exit $1
    fi
fi

{bwa} mem \
    -t 2 \
    -T {min_score} \
    -R '{read_group}' \
    {ref_fa} \
    {fastq1} \
    {fastq2} | \
{samtools} view -Sb - \
    > $UNSORTED_BAM;

{samtools} sort $UNSORTED_BAM "$BAM_WITHOUT_SUFFIX"_sorted;

"""

bwa_mem_biobambam = r"""
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

echo SGE_TASK_ID:$SGE_TASK_ID
echo SGE_TASK_FIRST:$SGE_TASK_FIRST
echo SGE_TASK_LAST:$SGE_TASK_LAST
echo SGE_TASK_STEPSIZE:$SGE_TASK_STEPSIZE

{array_data}

BAM_WITHOUT_SUFFIX=`echo {bam} | sed 's/\.bam//'`

{bwa} mem \
    -t 2 \
    -T {min_score} \
    -R '{read_group}' \
    {ref_fa} \
    {fastq1} \
    {fastq2} | \
{samtools} view -Sb - \
    > "$BAM_WITHOUT_SUFFIX"_unsorted.bam;

{biobambam}/bamsort index=1 \
                    level=1 \
                    inputthreads=2 \
                    outputthreads=2 \
                    calmdnm=1 \
                    calmdnmrecompindentonly=1 \
                    calmdnmreference={ref_fa} \
                    tmpfile="$BAM_WITHOUT_SUFFIX".tmp \
                    inputformat=bam\
                    indexfilename="$BAM_WITHOUT_SUFFIX"_bamsorted.bam.bai \
                    I="$BAM_WITHOUT_SUFFIX"_unsorted.bam \
                    O="$BAM_WITHOUT_SUFFIX"_bamsorted.bam

"""

samtools_merge_bam = \
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

NUM_FILES=`ls -1 {input_bam_files} | wc -l `

if [ $NUM_FILES  -ge 2 ]
then
    {samtools} merge \
        {output_bam_file}\
        {input_bam_files};
else
    cp {input_bam_files} {output_bam_file}
fi
{samtools} index {output_bam_file};

"""

biobambam_merge_bam = \
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

{env_variables}

OUT_BAM_PREFIX=`echo {output_bam_file} | sed 's/\.[^\.]\+$//'`

{biobambam}/bammerge  \
        {input_bam_files} \
        md5filename="$OUT_BAM_PREFIX".metrics \
        tmpfile="$OUT_BAM_PREFIX".tmp \
        indexfilename={output_bam_file}.bai \
        md5=1 \
        index=1 \
        > {output_bam_file}

# older version of pysam does not take care of PP in @PG
#samtools view -H {output_bam_file} | sed 's/PP:[^ 	]\+//' | samtools reheader - {output_bam_file}

"""

markduplicates = \
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

java -Xmx{memory} -Xms1G -jar {picard}/MarkDuplicates.jar \
         ASSUME_SORTED=true \
         I={input_bam} \
         O={output_bam} \
         M={output_bam}.met

{samtools} index {output_bam}

"""

biobambam_markduplicates = \
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

{env_variables}

OUT_BAM_PREFIX=`echo {output_bam} | sed 's/\.[^\.]\+$//'`

{biobambam}/bammarkduplicates  \
    M="$OUT_BAM_PREFIX".metrics \
    tmpfile="$OUT_BAM_PREFIX".tmp \
    markthreads=8 \
    rewritebam=1 \
    rewritebamlevel=1 \
    index=1 \
    md5=1\
    {input_bam_files} \
    O={output_bam}

# older version of pysam does not take care of PP in @PG
#samtools view -H {output_bam} | sed 's/PP:[^ 	]\+//' | samtools reheader - {output_bam}

"""
bam_stats_calc = \
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

# -t 1-13:1
echo SGE_TASK_ID:$SGE_TASK_ID
echo SGE_TASK_FIRST:$SGE_TASK_FIRST
echo SGE_TASK_LAST:$SGE_TASK_LAST
echo SGE_TASK_STEPSIZE:$SGE_TASK_STEPSIZE

if [ "$SGE_TASK_ID" = "1" ]
then
    {pcap}/bin/bam_stats.pl \
            --threads 4 \
            -i {bam_file} \
            -o {output_txt}.$SGE_TASK_ID

elif [ "$SGE_TASK_ID" = "2" ]
then
    {samtools} depth {bam_file} |\
    awk '{{sum+=$3; sumsq+=$3*$3}} END {{ print "Average:\t",sum/NR; print "Stdev:\t",sqrt(sumsq/NR - (sum/NR)**2)}}' \
         > {output_txt}.$SGE_TASK_ID
elif [ "$SGE_TASK_ID" = "3" ]
then
    #
    # metrics files are generated only by biobambam bammarkduplicates
    #
    MET_FILE=`echo {bam_file} | sed 's/\.bam/.metrics/'`
    if [ -e $MET_FILE ]
    then
        grep -A1 LIBRARY $MET_FILE > {output_txt}.$SGE_TASK_ID
    fi
elif [ "$SGE_TASK_ID" = "4" ]
then
    echo "samtools_flagstat" > {output_txt}.$SGE_TASK_ID
    {samtools} flagstat {bam_file} >> {output_txt}.$SGE_TASK_ID
else
    #SGE_TASK_ID = 5~14
    {python} {script_dir}/coverage.py \
            -i {bam_file} \
            -t {output_txt}.$SGE_TASK_ID.tmp \
            -f {ref_fa} \
            -g {genome_size} \
            -e {bed_file} \
            -b 1000 \
            -n 1000 \
            -p 10 \
            {chr_str_in_fa} \
            -c {coverage} \
            -s {samtools} \
            > {output_txt}.tmp.$SGE_TASK_ID
fi
"""


bam_stats_merge = \
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

OUT_TSV=`echo {output_txt} | sed 's/\.txt/.tsv/'`
OUT_XLS=`echo {output_txt} | sed 's/\.txt/.xls/'`

awk '/ratio/ {{ print }}' {output_txt}.tmp.5 > {output_txt}.5-1
for TMP_FILE in `ls {output_txt}.tmp.*`
do
    tail -n1 $TMP_FILE >> {output_txt}.5-1
    rm $TMP_FILE
done

{python} {script_dir}/merge_cov.py -i {output_txt}.5-1 -o {output_txt}.5
cat {output_txt}.1 {output_txt}.2 {output_txt}.3 {output_txt}.4 {output_txt}.5 > {output_txt}
{python} {script_dir}/mkxls.py -i {output_txt} -x $OUT_XLS
{python} {script_dir}/xl2tsv.py -t $OUT_TSV -x $OUT_XLS

"""

fisher_merge_bams = """
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

OUT_BAM_PREFIX=`echo {merged_bam_file} | sed 's/\.[^\.]\+$//'`
mkdir -p {out_dir}
if [ "{merge_bam_flag}" = "True" ]
then
    if [ "{use_biobambam}" = "True" ]
    then
        {biobambam}/bammerge  \
                        {bambam_input_bam_files} \
                        md5filename="$OUT_BAM_PREFIX".metrics \
                        tmpfile="$OUT_BAM_PREFIX".tmp \
                        indexfilename={merged_bam_file}.bai \
                        md5=1 \
                        index=1 \
                        > {merged_bam_file}
    else
        {samtools} merge \
                    {merged_bam_file} \
                    {samtools_input_bam_files};
        {samtools} index {merged_bam_file};
    fi
else
    cp {input_bam_file} {merged_bam_file}
    cp {input_bam_file}.bai {merged_bam_file}.bai
fi
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

echo SGE_TASK_ID:$SGE_TASK_ID
echo SGE_TASK_FIRST:$SGE_TASK_FIRST
echo SGE_TASK_LAST:$SGE_TASK_LAST
echo SGE_TASK_STEPSIZE:$SGE_TASK_STEPSIZE

echo {control_input_bam}
echo {disease_input_bam}
echo {output_txt}

source {script_dir}/interval.sh
source {script_dir}/interval_list.sh

CONTROL_INTERVAL_BAM=`echo {control_input_bam} | sed "s/\.bam/_${{INTERVAL[$SGE_TASK_ID]}}.bam/"`
DISEASE_INTERVAL_BAM=`echo {disease_input_bam} | sed "s/\.bam/_${{INTERVAL[$SGE_TASK_ID]}}.bam/"`

CONTROL_FILTERED_BAM=`echo {control_input_bam} | sed "s/\.bam/_${{INTERVAL[$SGE_TASK_ID]}}_filtered.bam/"`
DISEASE_FILTERED_BAM=`echo {disease_input_bam} | sed "s/\.bam/_${{INTERVAL[$SGE_TASK_ID]}}_filtered.bam/"`

if [ -e "{control_input_bam}" -a ! -e "$CONTROL_FILTERED_BAM" ]
then

    {samtools} view -bh \
                  -q {map_quality} \
                  {control_input_bam} \
                  ${{INTERVAL[$SGE_TASK_ID]}} \
              > $CONTROL_INTERVAL_BAM

    {python} {script_dir}/bamfilter.py \
               --ref_fa {ref_fa} \
               --max_indel {max_indel} \
               --max_distance {max_distance} \
               --map_quality {map_quality} \
            $CONTROL_INTERVAL_BAM \
            $CONTROL_FILTERED_BAM

    if [ "{remove_intermediate}" = "True" ]
    then
        rm $CONTROL_INTERVAL_BAM
    fi
fi

if [ -e "{disease_input_bam}" -a ! -e "$DISEASE_FILTERED_BAM" ]
then
    {samtools} view -bh \
                  -q {map_quality} \
                  {disease_input_bam} \
                  ${{INTERVAL[$SGE_TASK_ID]}} \
              > $DISEASE_INTERVAL_BAM

    {python} {script_dir}/bamfilter.py \
               --ref_fa {ref_fa} \
               --max_indel {max_indel} \
               --max_distance {max_distance} \
               --map_quality {map_quality} \
            $DISEASE_INTERVAL_BAM \
            $DISEASE_FILTERED_BAM

    if [ "{remove_intermediate}" = "True" ]
    then
        rm $DISEASE_INTERVAL_BAM
    fi
fi

INTERVAL_OUT=`echo {output_txt} | sed "s/\.txt/_${{INTERVAL[$SGE_TASK_ID]}}.txt/"`

if [ -e "$CONTROL_FILTERED_BAM"  ]
then
    SIZE_OF_CONTROL_BAM=`du -b $CONTROL_FILTERED_BAM | sed 's/\([0-9\.]\+\).\+/\\1/'`
else
    SIZE_OF_CONTROL_BAM=0
fi
if [ -e "$DISEASE_FILTERED_BAM"  ]
then
    SIZE_OF_DISEASE_BAM=`du -b $DISEASE_FILTERED_BAM | sed 's/\([0-9\.]\+\).\+/\\1/'`
else
    SIZE_OF_DISEASE_BAM=0
fi

if [ $SIZE_OF_CONTROL_BAM -gt 300 -a $SIZE_OF_DISEASE_BAM -gt 300 -a ! -e $INTERVAL_OUT  ]
then
    {samtools} mpileup \
                -BQ0 \
                -d 10000000 \
                -f {ref_fa} \
                $CONTROL_FILTERED_BAM \
                $DISEASE_FILTERED_BAM |\
    {python} {script_dir}/fisher.py \
               --output $INTERVAL_OUT \
               --ref_fa {ref_fa} \
               --base_quality {base_quality} \
               --mismatch_rate {mismatch_rate} \
               --min_depth {min_depth}\
               --compare
fi

if [ $SIZE_OF_CONTROL_BAM -gt 300 -a ! -e $INTERVAL_OUT ]
then
    {samtools} mpileup \
                -BQ0 \
                -d 10000000 \
                -f {ref_fa} \
                $CONTROL_FILTERED_BAM |\
    {python} {script_dir}/fisher.py \
               --output $INTERVAL_OUT \
               --ref_fa {ref_fa} \
               --base_quality {base_quality} \
               --mismatch_rate {mismatch_rate} \
               --min_depth {min_depth}

fi

if [ $SIZE_OF_DISEASE_BAM -gt 300 -a ! -e $INTERVAL_OUT ]
then
    {samtools} mpileup \
                -BQ0 \
                -d 10000000 \
                -f {ref_fa} \
                $DISEASE_FILTERED_BAM |\
    {python} {script_dir}/fisher.py \
               --output $INTERVAL_OUT \
               --ref_fa {ref_fa} \
               --base_quality {base_quality} \
               --mismatch_rate {mismatch_rate} \
               --min_depth {min_depth}

fi

if [ "{remove_intermediate}" = "True" ]
then
    rm $CONTROL_FILTERED_BAM
    rm $DISEASE_FILTERED_BAM
fi

"""

merge_fisher_result= \
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

source {script_dir}/interval_list.sh

FIRST_INTERVAL=`echo $INTERVAL_LIST | head -1 | sed 's/^\([^ ]\+\) .\+$/\\1/g'`
INTERVAL_OUT=`echo {output_txt} | sed "s/\.txt/_$FIRST_INTERVAL.txt/"`

head -1 $INTERVAL_OUT > {output_txt}

for INTERVAL in $INTERVAL_LIST
do
    INTERVAL_OUT=`echo {output_txt} | sed "s/\.txt/_${{INTERVAL[$SGE_TASK_ID]}}.txt/"`
    if [ -f $INTERVAL_OUT ]
    then
        tail -n+2 $INTERVAL_OUT >> {output_txt}
        rm $INTERVAL_OUT
    fi
done

"""

interval_num = 274

itd_detection = \
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

echo SGE_TASK_ID:$SGE_TASK_ID
echo SGE_TASK_FIRST:$SGE_TASK_FIRST
echo SGE_TASK_LAST:$SGE_TASK_LAST
echo SGE_TASK_STEPSIZE:$SGE_TASK_STEPSIZE

# {itd_inhouse_files}

{array}

bash {itd_detector}/detectITD.sh \
        {bam_file} \
        {output_file} \
        {name} \
        {itd_inhouse_dir}

"""

annotation = \
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

if [ "{output_in_vcf}" = "True" ]
then
    {python} {scriptdir}/ToVCF.py \
        --ref_fasta {ref_fa} \
        --output_vcf {output_vcf} \
        --annovar_input {input_file}

    {annovar}/table_annovar.pl \
        --vcfinput \
        --outfile {output_prefix} \
        {output_vcf} \
        {annovar}/humandb \
        {table_annovar_params}

else

    if [ "{use_table_annovar}" = "True" ]
    then
        {annovar}/table_annovar.pl \
            {input_file} \
            --outfile {output_prefix} \
            {annovar}/humandb \
            {table_annovar_params}
    else

        {annovar}/summarize_annovar.pl \
            {summarize_annovar_params} \
            {input_file} \
            --outfile {output_prefix} \
            {annovar}/humandb

        {python} {scriptdir}/ToVCF.py \
            --ref_fasta {ref_fa} \
            --output_vcf {output_vcf} \
            --tsv_annovar {input_file}
    fi
fi

"""

