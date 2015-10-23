#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Mutation_call(Stage_task):

    task_name = "mutation_call"

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

# set python environment
export PYTHONHOME={pythonhome}
export PATH=$PYTHONHOME/bin:$PATH
export LD_LIBRARY_PATH={ld_library_path}
export PYTHONPATH={pythonpath}

if [ _{control_bam} = "_None" ]; then 
    {fisher} single -R {region} -o {out_prefix}.fisher_mutations.{task_id}.txt --ref_fa {ref_fa} --mapping_quality {map_quality} --base_quality {base_quality}  --min_allele_freq {min_allele_freq} --post_10_q {post_10_q} --min_depth {min_depth} -1 {disease_bam} --samtools_path {samtools} || exit $?

    {mutfilter} realignment --tumor_min_mismatch {realign_min_mismatch} --normal_max_mismatch {realign_max_mismatch} --score_difference {realign_score_diff} --window_size {realign_window_size} --max_depth {realign_max_depth} --target_mutation_file {out_prefix}.fisher_mutations.{task_id}.txt -1 {disease_bam} --output {out_prefix}.realignment_mutations.{task_id}.txt --ref_genome {ref_fa} --blat_path {blat} || exit $?

    {mutfilter} simplerepeat --target_mutation_file {out_prefix}.realignment_mutations.{task_id}.txt --output {out_prefix}.simplerepeat_mutations.{task_id}.txt --simple_repeat_db {simple_repeat_db} || exit $?

else
    {fisher} comparison -R {region} -o {out_prefix}.fisher_mutations.{task_id}.txt --ref_fa {ref_fa} --mapping_quality {map_quality} --base_quality {base_quality}  --min_allele_freq {min_allele_freq} --max_allele_freq {max_allele_freq} --min_depth {min_depth} --fisher_value {fisher_thres} -2 {control_bam} -1 {disease_bam} --samtools_path {samtools} || exit $?

    {mutfilter} realignment --tumor_min_mismatch {realign_min_mismatch} --normal_max_mismatch {realign_max_mismatch} --score_difference {realign_score_diff} --window_size {realign_window_size} --max_depth {realign_max_depth} --target_mutation_file {out_prefix}.fisher_mutations.{task_id}.txt -1 {disease_bam} -2 {control_bam} --output {out_prefix}.realignment_mutations.{task_id}.txt --ref_genome {ref_fa} --blat_path {blat} || exit $?

    {mutfilter} indel --search_length {indel_search_length} --neighbor {indel_neighbor} --base_qual {indel_base_quality} --min_depth {indel_min_depth} --min_mismatch {indel_min_mismatch} --af_thres {indel_min_allele_freq} --target_mutation_file {out_prefix}.realignment_mutations.{task_id}.txt -2 {control_bam} --output {out_prefix}.indel_mutations.{task_id}.txt || exit $?

    {mutfilter} breakpoint --max_depth {bp_max_depth} --min_clip_size {bp_min_clip_size} --junc_num_thres {bp_junc_num_thres} --mapq_thres {bp_map_quality} --target_mutation_file {out_prefix}.indel_mutations.{task_id}.txt -2 {control_bam} --output {out_prefix}.simplerepeat_mutations.{task_id}.txt || exit $?

    {mutfilter} simplerepeat --target_mutation_file {out_prefix}.simplerepeat_mutations.{task_id}.txt --output {out_prefix}.simplerepeat_mutations.{task_id}.txt --simple_repeat_db {simple_repeat_db} || exit $?
fi

if [ _{control_bam_list} != "_None" ]; then 
    {EBFilter} -f anno -q {eb_map_quality} -Q {eb_base_quality} {out_prefix}.simplerepeat_mutations.{task_id}.txt {disease_bam} {control_bam_list} {out_prefix}.ebfilter_mutations.{task_id}.txt || exit $?
else
    cp {out_prefix}.simplerepeat_mutations.{task_id}.txt {out_prefix}.ebfilter_mutations.{task_id}.txt
fi

if [ _{active_annovar_flag} = "_True" ];then
    {annovar}/table_annovar.pl --outfile {out_prefix}_mutations_candidate.{task_id} {table_annovar_params} {out_prefix}.ebfilter_mutations.{task_id}.txt {annovar}/humandb || exit $?
else
    cp {out_prefix}.ebfilter_mutations.{task_id}.txt {out_prefix}_mutations_candidate.{task_id}.hg19_multation.txt
fi

"""
    def __init__(self, qsub_option, script_dir):
        super(Mutation_call, self).__init__(qsub_option, script_dir)


