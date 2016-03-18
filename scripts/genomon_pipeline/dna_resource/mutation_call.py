#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Mutation_call(Stage_task):

    task_name = "mutation_call"

    script_template = """
#!/bin/bash
#
#$ -S /bin/bash         # set shell in UGE
#$ -cwd                 # execute at the submitted dir
pwd                     # print current working directory
hostname                # print hostname
date                    # print date
set -xv

# set python environment
export PYTHONHOME={pythonhome}
samtools_home={samtools}
export PATH=${{samtools_home%/*}}:$PYTHONHOME/bin:$PATH
export LD_LIBRARY_PATH={ld_library_path}
export PYTHONPATH={pythonpath}

REGION=`head -n $SGE_TASK_ID {interval_list} | tail -n 1`

if [ _{control_bam} = "_None" ]; then 
    {fisher} single -R ${{REGION}} -o {out_prefix}.fisher_mutations.${{SGE_TASK_ID}}.txt --ref_fa {ref_fa} --mapping_quality {map_quality} --base_quality {base_quality}  --min_allele_freq {min_allele_freq} --post_10_q {post_10_q} --min_depth {min_depth} -1 {disease_bam} --samtools_path {samtools} || exit $?

    {mutfilter} realignment --tumor_min_mismatch {realign_min_mismatch} --normal_max_mismatch {realign_max_mismatch} --score_difference {realign_score_diff} --window_size {realign_window_size} --max_depth {realign_max_depth} --target_mutation_file {out_prefix}.fisher_mutations.${{SGE_TASK_ID}}.txt -1 {disease_bam} --output {out_prefix}.realignment_mutations.${{SGE_TASK_ID}}.txt --ref_genome {ref_fa} --blat_path {blat} || exit $?

    {mutfilter} simplerepeat --target_mutation_file {out_prefix}.realignment_mutations.${{SGE_TASK_ID}}.txt --output {out_prefix}.simplerepeat_mutations.${{SGE_TASK_ID}}.txt --simple_repeat_db {simple_repeat_db} || exit $?

else
    {fisher} comparison -R ${{REGION}} -o {out_prefix}.fisher_mutations.${{SGE_TASK_ID}}.txt --ref_fa {ref_fa} --mapping_quality {map_quality} --base_quality {base_quality}  --min_allele_freq {min_allele_freq} --max_allele_freq {max_allele_freq} --min_depth {min_depth} --fisher_value {fisher_thres} -2 {control_bam} -1 {disease_bam} --samtools_path {samtools} || exit $?

    {mutfilter} realignment --tumor_min_mismatch {realign_min_mismatch} --normal_max_mismatch {realign_max_mismatch} --score_difference {realign_score_diff} --window_size {realign_window_size} --max_depth {realign_max_depth} --target_mutation_file {out_prefix}.fisher_mutations.${{SGE_TASK_ID}}.txt -1 {disease_bam} -2 {control_bam} --output {out_prefix}.realignment_mutations.${{SGE_TASK_ID}}.txt --ref_genome {ref_fa} --blat_path {blat} || exit $?

    {mutfilter} indel --search_length {indel_search_length} --neighbor {indel_neighbor} --base_qual {indel_base_quality} --min_depth {indel_min_depth} --min_mismatch {indel_min_mismatch} --af_thres {indel_min_allele_freq} --target_mutation_file {out_prefix}.realignment_mutations.${{SGE_TASK_ID}}.txt -2 {control_bam} --output {out_prefix}.indel_mutations.${{SGE_TASK_ID}}.txt || exit $?

    {mutfilter} breakpoint --max_depth {bp_max_depth} --min_clip_size {bp_min_clip_size} --junc_num_thres {bp_junc_num_thres} --mapq_thres {bp_map_quality} --target_mutation_file {out_prefix}.indel_mutations.${{SGE_TASK_ID}}.txt -2 {control_bam} --output {out_prefix}.breakpoint_mutations.${{SGE_TASK_ID}}.txt || exit $?

    {mutfilter} simplerepeat --target_mutation_file {out_prefix}.breakpoint_mutations.${{SGE_TASK_ID}}.txt --output {out_prefix}.simplerepeat_mutations.${{SGE_TASK_ID}}.txt --simple_repeat_db {simple_repeat_db} || exit $?
fi

if [ _{control_bam_list} != "_None" ]; then 
    {EBFilter} --loption --region ${{REGION}} -f anno -q {eb_map_quality} -Q {eb_base_quality} {out_prefix}.simplerepeat_mutations.${{SGE_TASK_ID}}.txt {disease_bam} {control_bam_list} {out_prefix}.ebfilter_mutations.${{SGE_TASK_ID}}.txt || exit $?
else
    cp {out_prefix}.simplerepeat_mutations.${{SGE_TASK_ID}}.txt {out_prefix}.ebfilter_mutations.${{SGE_TASK_ID}}.txt
fi

if [ _{active_inhouse_normal_flag} = "_True" ]; then 
    {mutanno} mutation -t {out_prefix}.ebfilter_mutations.${{SGE_TASK_ID}}.txt -o {out_prefix}.inhouse_normal.${{SGE_TASK_ID}}.txt -d {inhouse_normal_database} || exit $?
else
    cp {out_prefix}.ebfilter_mutations.${{SGE_TASK_ID}}.txt {out_prefix}.inhouse_normal.${{SGE_TASK_ID}}.txt
fi

if [ _{active_inhouse_tumor_flag} = "_True" ]; then 
    {mutanno} mutation -t {out_prefix}.inhouse_normal.${{SGE_TASK_ID}}.txt -o {out_prefix}.inhouse_tumor.${{SGE_TASK_ID}}.txt -d {inhouse_tumor_database} || exit $?
else
    cp {out_prefix}.inhouse_normal.${{SGE_TASK_ID}}.txt {out_prefix}.inhouse_tumor.${{SGE_TASK_ID}}.txt
fi

if [ _{active_HGVD_flag} = "_True" ]; then 
    {mutanno} mutation -t {out_prefix}.inhouse_tumor.${{SGE_TASK_ID}}.txt -o {out_prefix}.HGVD.${{SGE_TASK_ID}}.txt -d {HGVD_database} -c 15 || exit $?
else
    cp {out_prefix}.inhouse_tumor.${{SGE_TASK_ID}}.txt {out_prefix}.HGVD.${{SGE_TASK_ID}}.txt
fi

if [ _{active_HGMD_flag} = "_True" ]; then 
    {mutanno} mutation -t {out_prefix}.HGVD.${{SGE_TASK_ID}}.txt -o {out_prefix}.HGMD.${{SGE_TASK_ID}}.txt -d {HGMD_database} -c 7 || exit $?
else
    cp {out_prefix}.HGVD.${{SGE_TASK_ID}}.txt {out_prefix}.HGMD.${{SGE_TASK_ID}}.txt
fi

if [ _{active_annovar_flag} = "_True" ];then
    {annovar}/table_annovar.pl --outfile {out_prefix}_mutations_candidate.${{SGE_TASK_ID}} {table_annovar_params} {out_prefix}.HGMD.${{SGE_TASK_ID}}.txt {annovar}/humandb || exit $?
else
    cp {out_prefix}.HGMD.${{SGE_TASK_ID}}.txt {out_prefix}_mutations_candidate.${{SGE_TASK_ID}}.hg19_multianno.txt
fi

"""
    def __init__(self, qsub_option, script_dir):
        super(Mutation_call, self).__init__(qsub_option, script_dir)


