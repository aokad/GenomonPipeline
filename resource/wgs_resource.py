#!/usr/bin/python

#
# Job scripts for wgs_pipeline.py
#
bamtofastq_s = """
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

source {scriptdir}/utility.sh
{biobambam}/bamtofastq  exclude=QCFAIL,SECONDARY,SUPPLEMENTARY\
            T={tmpfastq}\
            F={outfastq1}\
            collate=1\
            filename={bamfile}
check_error $?

"""

bamtofastq_p = """
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

source {scriptdir}/utility.sh

{biobambam}/bamtofastq  exclude=QCFAIL,SECONDARY,SUPPLEMENTARY\
            T={tmpfastq}\
            F={outfastq1}\
            F2={outfastq2}\
            collate=1\
            filename={bamfile}
check_error $?

"""
extract_gz = """
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

echo SGE_TASK_ID:$SGE_TASK_ID
echo SGE_TASK_FIRST:$SGE_TASK_FIRST
echo SGE_TASK_LAST:$SGE_TASK_LAST
echo SGE_TASK_STEPSIZE:$SGE_TASK_STEPSIZE

source {scriptdir}/utility.sh

{array_data}

gzip -dc {input_file} >> {output_file}

check_error $?
"""


splitfile = """
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

echo SGE_TASK_ID:$SGE_TASK_ID
echo SGE_TASK_FIRST:$SGE_TASK_FIRST
echo SGE_TASK_LAST:$SGE_TASK_LAST
echo SGE_TASK_STEPSIZE:$SGE_TASK_STEPSIZE

source {scriptdir}/utility.sh

{array_data}

case {input_file} in
*\.gz)
    if [ "{fastq_filter}" = "True" ]
    then
        zcat {input_file} | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v '^--$' | split -a {suffix_len} -d -l {lines_per_file} - {output_prefix}
        check_error $?
    else
        zcat {input_file} | split -a {suffix_len} -d -l {lines_per_file} - {output_prefix}
        check_error $?
    fi
    ;;

*\.bz2)
    if [ "{fastq_filter}" = "True" ]
    then
        bzip2 -dc  {input_file} | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v '^--$' | split -a {suffix_len} -d -l {lines_per_file} - {output_prefix}
        check_error $?
    else
        bzip2 -dc  {input_file} | split -a {suffix_len} -d -l {lines_per_file} - {output_prefix}
        check_error $?
    fi
    ;;

*)
    if [ "{fastq_filter}" = "True" ]
    then
        cat {input_file} | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v '^--$' | split -a {suffix_len} -d -l {lines_per_file} - {output_prefix}
        check_error $?
    else
        split -a {suffix_len} -d -l {lines_per_file} {input_file} {output_prefix}
        check_error $?
    fi
    ;;
esac

for FILE in `ls {output_prefix}*`
do
    mv $FILE $FILE{output_suffix}
    check_error $?
done

"""

cutadapt = """
#!/bin/bash
#
# Set SGE
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

source {scriptdir}/utility.sh

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
    check_error $?
fi

{bwa} mem \
    -T {min_score} \
    -R '{read_group}' \
    {additional_params} \
    {ref_fa} \
    {fastq1} \
    {fastq2} | \
{samtools} view -Sb - \
    > $UNSORTED_BAM;
check_error $?

{samtools} sort $UNSORTED_BAM "$BAM_WITHOUT_SUFFIX"_sorted;
check_error $?

"""

bwa_mem_biobambam = r"""
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

{env_variables}

echo SGE_TASK_ID:$SGE_TASK_ID
echo SGE_TASK_FIRST:$SGE_TASK_FIRST
echo SGE_TASK_LAST:$SGE_TASK_LAST
echo SGE_TASK_STEPSIZE:$SGE_TASK_STEPSIZE

source {scriptdir}/utility.sh

{array_data}

BAM_WITHOUT_SUFFIX=`echo {bam} | sed 's/\.bam//'`

{bwa} mem \
    -T {min_score} \
    -R '{read_group}' \
    {additional_params} \
    {ref_fa} \
    {fastq1} \
    {fastq2} | \
{samtools} view -Sb - \
    > "$BAM_WITHOUT_SUFFIX"_unsorted.bam;
check_error $?

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
check_error $?

"""

samtools_merge_bam = \
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

source {scriptdir}/utility.sh

NUM_FILES=`ls -1 {input_bam_files} | wc -l `

if [ $NUM_FILES  -ge 2 ]
then
    {samtools} merge \
        {output_bam_file}\
        {input_bam_files};
    check_error $?

else
    cp {input_bam_files} {output_bam_file}
    check_error $?
fi

{samtools} index {output_bam_file};
check_error $?

"""

biobambam_merge_bam = \
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

{env_variables}

source {scriptdir}/utility.sh

OUT_BAM_PREFIX=`echo {output_bam_file} | sed 's/\.[^\.]\+$//'`

{biobambam}/bammerge  \
        {input_bam_files} \
        md5filename="$OUT_BAM_PREFIX".metrics \
        tmpfile="$OUT_BAM_PREFIX".tmp \
        indexfilename={output_bam_file}.bai \
        md5=1 \
        index=1 \
        > {output_bam_file}

check_error $?

# older version of pysam does not take care of PP in @PG
#samtools view -H {output_bam_file} | sed 's/PP:[^ 	]\+//' | samtools reheader - {output_bam_file}

"""

markduplicates = \
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

source {scriptdir}/utility.sh

java -Xmx{memory} -Xms1G -jar {picard}/MarkDuplicates.jar \
         ASSUME_SORTED=true \
         I={input_bam} \
         O={output_bam} \
         M={output_bam}.met
check_error $?

{samtools} index {output_bam}
check_error $?

"""

biobambam_markduplicates = \
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

{env_variables}

source {scriptdir}/utility.sh

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
check_error $?

# older version of pysam does not take care of PP in @PG
#samtools view -H {output_bam} | sed 's/PP:[^ 	]\+//' | samtools reheader - {output_bam}

"""
bam_stats_calc = \
"""
#!/bin/bash
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

{env_variables}
source {scriptdir}/utility.sh

if [ "$SGE_TASK_ID" = "1" ]
then
    {pcap}/bin/bam_stats.pl \
            --threads 4 \
            -i {bam_file} \
            -o {output_txt}.$SGE_TASK_ID
    check_error $?

elif [ "$SGE_TASK_ID" = "2" ]
then
    {samtools} depth {bam_file} |\
    awk '{{sum+=$3; sumsq+=$3*$3}} END {{ print "Average:\t",sum/NR; print "Stdev:\t",sqrt(sumsq/NR - (sum/NR)**2)}}' \
         > {output_txt}.$SGE_TASK_ID
    check_error $?

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
    check_error $?

else
    #SGE_TASK_ID = 5~14
    {python} {scriptdir}/coverage.py \
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
    check_error $?
fi
"""


bam_stats_merge = \
"""
#!/bin/bash
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

source {scriptdir}/utility.sh

awk '/ratio/ {{ print }}' {output_txt}.tmp.5 > {output_txt}.5-1
for TMP_FILE in `ls {output_txt}.tmp.*`
do
    tail -n1 $TMP_FILE >> {output_txt}.5-1
    rm -f $TMP_FILE
done

{python} {scriptdir}/merge_cov.py -i {output_txt}.5-1 -o {output_txt}.5
check_error $?

cat {output_txt}.1 {output_txt}.2 {output_txt}.3 {output_txt}.4 {output_txt}.5 > {output_txt}
{python} {scriptdir}/mkxls.py -i {output_txt} -x $OUT_XLS
check_error $?

{python} {scriptdir}/xl2tsv.py -t $OUT_TSV -x $OUT_XLS
check_error $?

"""

fisher_merge_bams = """
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

source {scriptdir}/utility.sh

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
        check_error $?
    else
        {samtools} merge \
                    -f \
                    {merged_bam_file} \
                    {samtools_input_bam_files};
        check_error $?

        {samtools} index {merged_bam_file};
        check_error $?
    fi
else
    cp {input_bam_file} {merged_bam_file}
    check_error $?

    cp {input_bam_file}.bai {merged_bam_file}.bai
    check_error $?
fi
"""

fisher_mutation_call = \
"""
#!/bin/bash
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

source {scriptdir}/utility.sh
source {scriptdir}/interval.sh
source {scriptdir}/interval_list.sh

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
    check_error $?

    {python} {scriptdir}/bamfilter.py \
               --ref_fa {ref_fa} \
               --max_indel {max_indel} \
               --max_distance {max_distance} \
               --map_quality {map_quality} \
            $CONTROL_INTERVAL_BAM \
            $CONTROL_FILTERED_BAM
    check_error $?

    if [ "{remove_intermediate}" = "True" ]
    then
        rm -f $CONTROL_INTERVAL_BAM
    fi
fi

if [ -e "{disease_input_bam}" -a ! -e "$DISEASE_FILTERED_BAM" ]
then
    {samtools} view -bh \
                  -q {map_quality} \
                  {disease_input_bam} \
                  ${{INTERVAL[$SGE_TASK_ID]}} \
              > $DISEASE_INTERVAL_BAM
    check_error $?

    {python} {scriptdir}/bamfilter.py \
               --ref_fa {ref_fa} \
               --max_indel {max_indel} \
               --max_distance {max_distance} \
               --map_quality {map_quality} \
            $DISEASE_INTERVAL_BAM \
            $DISEASE_FILTERED_BAM
    check_error $?

    if [ "{remove_intermediate}" = "True" ]
    then
        rm -f $DISEASE_INTERVAL_BAM
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
    {python} {scriptdir}/fisher.py \
               --output $INTERVAL_OUT \
               --ref_fa {ref_fa} \
               --base_quality {base_quality} \
               --mismatch_rate {mismatch_rate} \
               --min_depth {min_depth}\
               --compare
    check_error $?
fi

if [ $SIZE_OF_CONTROL_BAM -gt 300 -a ! -e $INTERVAL_OUT ]
then
    {samtools} mpileup \
                -BQ0 \
                -d 10000000 \
                -f {ref_fa} \
                $CONTROL_FILTERED_BAM |\
    {python} {scriptdir}/fisher.py \
               --output $INTERVAL_OUT \
               --ref_fa {ref_fa} \
               --base_quality {base_quality} \
               --mismatch_rate {mismatch_rate} \
               --min_depth {min_depth}
    check_error $?

fi

if [ $SIZE_OF_DISEASE_BAM -gt 300 -a ! -e $INTERVAL_OUT ]
then
    {samtools} mpileup \
                -BQ0 \
                -d 10000000 \
                -f {ref_fa} \
                $DISEASE_FILTERED_BAM |\
    {python} {scriptdir}/fisher.py \
               --output $INTERVAL_OUT \
               --ref_fa {ref_fa} \
               --base_quality {base_quality} \
               --mismatch_rate {mismatch_rate} \
               --min_depth {min_depth}
    check_error $?

fi

if [ "{remove_intermediate}" = "True" ]
then
    rm -f $CONTROL_FILTERED_BAM
    rm -f ${{CONTROL_FILTERED_BAM}}.bai
    rm -f $DISEASE_FILTERED_BAM
    rm -f ${{DISEASE_FILTERED_BAM}}.bai
fi

"""

merge_fisher_result= \
"""
#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -e {log}             # log file directory
#$ -o {log}             # log file directory
pwd                     # print current working directory
hostname                # print hostname
date                    # print date
set -xv

source {scriptdir}/utility.sh
source {scriptdir}/interval_list.sh

FIRST_INTERVAL=`echo $INTERVAL_LIST | head -1 | sed 's/^\([^ ]\+\) .\+$/\\1/g'`
INTERVAL_OUT=`echo {output_txt} | sed "s/\.txt/_$FIRST_INTERVAL.txt/"`

head -1 $INTERVAL_OUT > {output_txt}

for INTERVAL in $INTERVAL_LIST
do
    INTERVAL_OUT=`echo {output_txt} | sed "s/\.txt/_${{INTERVAL[$SGE_TASK_ID]}}.txt/"`
    if [ -f $INTERVAL_OUT ]
    then
        tail -n+2 $INTERVAL_OUT >> {output_txt}
        rm -f $INTERVAL_OUT
    fi
done

"""

interval_num = 274

itd_detection = \
"""
#!/bin/bash
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
source {scriptdir}/utility.sh

{array}

bash {itd_detector}/detectITD.sh \
        {bam_file} \
        {output_file} \
        {name} \
        {itd_inhouse_dir}
check_error $?

"""

sv_parse_filt = \
"""
#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -e {log}             # log file directory
#$ -o {log}             # log file directory
pwd                     # print current working directory
hostname                # print hostname
date                    # print date
set -xv

# set python environment
export PYTHONHOME={pythonhome}
export PATH=$PYTHONHOME/bin:$PATH
export LD_LIBRARY_PATH={ld_library_path}
export PYTHONPATH={pythonpath}


source {scriptdir}/utility.sh

{genomon_sv} {method} {sample_conf} {param_conf}

check_error $?

"""

annotation = \
"""
#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -e {log}             # log file directory
#$ -o {log}             # log file directory
pwd                     # print current working directory
hostname                # print hostname
date                    # print date
set -xv

source {scriptdir}/utility.sh

if [ "{output_in_vcf}" = "True" ]
then
    {python} {scriptdir}/ToVCF.py \
        --ref_fasta {ref_fa} \
        --output_vcf {output_vcf} \
        --annovar_input {input_file}
    check_error $?

    {annovar}/table_annovar.pl \
        --vcfinput \
        --outfile {output_prefix} \
        {output_vcf} \
        {annovar}/humandb \
        {table_annovar_params}
    check_error $?

else

    if [ "{use_table_annovar}" = "True" ]
    then
        {annovar}/table_annovar.pl \
            {input_file} \
            --outfile {output_prefix} \
            {annovar}/humandb \
            {table_annovar_params}
        check_error $?
    else

        {annovar}/summarize_annovar.pl \
            {summarize_annovar_params} \
            {input_file} \
            --outfile {output_prefix} \
            {annovar}/humandb
        check_error $?

        {python} {scriptdir}/ToVCF.py \
            --ref_fasta {ref_fa} \
            --output_vcf {output_vcf} \
            --tsv_annovar {input_file}
        check_error $?
    fi
fi

"""

