#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Mutation_merge(Stage_task):

    task_name = "mutation_merge"

    script_template = """
#!/bin/bash
#
#$ -S /bin/bash         # set shell in UGE
#$ -cwd                 # execute at the submitted dir
#$ -e {log}             # log file directory
#$ -o {log}             # log file directory
pwd                     # print current working directory
hostname                # print hostname
date                    # print date
set -xv

mut_header=""

if [ _{control_bam} = "_None" ]
then
    mut_header="Chr,Start,End,Ref,Alt,depth,variantNum,bases,A_C_G_T,misRate,strandRatio,10%_posterior_quantile,posterior_mean,90%_posterior_quantile,readPairNum,variantPairNum,otherPairNum,10%_posterior_quantile(realignment),posterior_mean(realignment),90%_posterior_quantile(realignment),simple_repeat_pos,simple_repeat_seq,P-value(EBCall)"
else
    mut_header="Chr,Start,End,Ref,Alt,depth_tumor,variantNum_tumor,depth_normal,variantNum_normal,bases_tumor,bases_normal,A_C_G_T_tumor,A_C_G_T_normal,misRate_tumor,strandRatio_tumor,misRate_normal,strandRatio_normal,P-value(fisher),readPairNum_tumor,variantPairNum_tumor,otherPairNum_tumor,readPairNum_normal,variantPairNum_normal,otherPairNum_normal,P-value(fhsher_realignment),indel_mismatch_rate,indel_mismatch_rate,bp_mismatch_count,distance_from_breakpoint,simple_repeat_pos,simple_repeat_seq,P-value(EBCall)"
fi

if [ _{control_bam_list} = "_None" ]
then
    mut_header=`echo ${{mut_header}} | awk -F"," -v OFS="," '{{$NF=""; sub(/.$/,""); print $0}}'` || exit $?
fi

print_header=""

if [ _{active_annovar_flag} = _True ]
then
    tmp_header=`head -n 1 {out_prefix}_mutations_candidate.1.hg19_multianno.txt | awk -F"\t" -v OFS="\t" '{{$NF=""; sub(/.$/,""); print $0}}'` || exit $?
    print_header=${{tmp_header}}

    tmp_header=`echo $mut_header | cut -d "," -f 6- | tr "," "\t"` || exit $?
    print_header=${{print_header}}'\t'${{tmp_header}}
else
    tmp_header=`echo $mut_header | tr "," "\t"` || exit $?
    print_header=${{tmp_header}}
fi

if [ _{active_inhouse_normal_flag} = "_True" ]
then
    inhouse_normal_header="inhouseSNV"
    print_header=${{print_header}}'\t'${{inhouse_normal_header}} || exit $?
fi

if [ _{active_inhouse_tumor_flag} = "_True" ]
then
    inhouse_tumor_header="inhouseCandidateSNV"
    print_header=${{print_header}}'\t'${{inhouse_tumor_header}} || exit $?
fi

if [ _{active_HGVD_flag} = "_True" ]
then
    HGVD_header="HGVD_20131010:rsID,HGVD_20131010:AF,HGVD_20131010:JPN,HGVD_20131010:JPN_NS,HGVD_20131010:ASN,HGVD_20131010:EUR,HGVD_20131010:AMR,HGVD_20131010:AFR,HGVD_20131010:NumberOfSamples,HGVD_20131010:Filter,HGVD_20131010:Mean_depth,HGVD_20131010:SD_depth,HGVD_20131010:NR,HGVD_20131010:NA,HGVD_20131010:Gene"
    tmp_header=`echo $HGVD_header | tr "," "\t"` || exit $?
    print_header=${{print_header}}'\t'${{tmp_header}} || exit $?
fi

echo "$print_header" > {out_prefix}_genomon_mutations.result.txt || exit $?

for i in `seq 1 1 {filecount}`
do
    if [ _{active_annovar_flag} = "_True" ]
    then
        awk 'NR>1 {{print}}' {out_prefix}_mutations_candidate.${{i}}.hg19_multianno.txt >> {out_prefix}_genomon_mutations.result.txt || exit $?
    else
        cat {out_prefix}_mutations_candidate.${{i}}.hg19_multianno.txt >> {out_prefix}_genomon_mutations.result.txt || exit $?
    fi
done

"""
    def __init__(self, qsub_option, script_dir):
        super(Mutation_merge, self).__init__(qsub_option, script_dir)


