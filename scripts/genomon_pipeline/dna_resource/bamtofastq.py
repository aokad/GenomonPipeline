#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Bam2Fastq(Stage_task):

    task_name = "bam2fastq"

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

echo -n > {f1_name}.tmp
echo -n > {f2_name}.tmp
echo -n > {t}.tmp
echo -n > {s}.tmp
echo -n > {o1_name}.tmp
echo -n > {o2_name}.tmp
echo "{input_bam}" | awk -F ";" '{{n=split($0,f); for(i=1;i<=n;i++){{print f[i];}}}}' | while read bam; do
    echo $bam
    {biobambam}/bamtofastq collate=1 exclude=QCFAIL,SECONDARY,SUPPLEMENTARY filename=${{bam}} F={f1_name} F2={f2_name} T={t} S={s} O={o1_name} O2={o2_name} || exit $?
    cat {f1_name} >> {f1_name}.tmp || exit $?
    cat {f2_name} >> {f2_name}.tmp || exit $?
    if [ -s {t} ]; then
        cat {t} >> {t}.tmp || exit $?
    fi
    if [ -s {s} ]; then
        cat {s} >> {s}.tmp || exit $?
    fi
    if [ -s {o1_name} ]; then
        cat {o1_name} >> {o1_name}.tmp || exit $?
    fi
    if [ -s {o2_name} ]; then
        cat {o2_name} >> {o2_name}.tmp || exit $?
    fi
done
mv {f1_name}.tmp {f1_name} || exit $?
mv {f2_name}.tmp {f2_name} || exit $?
mv {t}.tmp {t} || exit $?
mv {s}.tmp {s} || exit $?
mv {o1_name}.tmp {o1_name} || exit $?
mv {o2_name}.tmp {o2_name} || exit $?

"""

    def __init__(self, qsub_option, script_dir):
        super(Bam2Fastq, self).__init__(qsub_option, script_dir)


