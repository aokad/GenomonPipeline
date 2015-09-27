import os
import shutil
from ruffus import *
from genomon_pipeline.config.run_conf import *
from genomon_pipeline.config.genomon_conf import *
from genomon_pipeline.config.task_conf import *
from genomon_pipeline.config.sample_conf import *
from genomon_pipeline.dna_resource.fastq_splitter import *
from genomon_pipeline.dna_resource.bwa_align import *
from genomon_pipeline.dna_resource.markduplicates import *
from genomon_pipeline.dna_resource.mutation_call import *
from genomon_pipeline.dna_resource.bamtofastq import *

# set task classes
fastq_splitter = Fastq_splitter(task_conf.get("split_fast", "qsub_option"), run_conf.project_root + '/script')
bwa_align = Bwa_align(task_conf.get("bwa_mem", "qsub_option"), run_conf.project_root + '/script')
markduplicates = Markduplicates(task_conf.get("markduplicates", "qsub_option"), run_conf.project_root + '/script')
mutation_call = Mutation_call(task_conf.get("mutation_call", "qsub_option"), run_conf.project_root + '/script')
bamtofastq = Bam2Fastq(task_conf.get("bam2fastq", "qsub_option"), run_conf.project_root + '/script')

# generate list of linked_fastq file path
linked_fastq_list = []
for sample in sample_conf.fastq:
    fastq_prefix, ext = os.path.splitext(sample_conf.fastq[sample][0][0])
    linked_fastq_list.append([run_conf.project_root + '/fastq/' + sample + '/1_1' + ext,
                              run_conf.project_root + '/fastq/' + sample + '/1_2' + ext])

markdup_bam_list = []
for complist in sample_conf.compare:
    markdup_bam_list.append([run_conf.project_root + '/bam/' + complist[0] + '/' + complist[0] + '.markdup.bam',
                             run_conf.project_root + '/bam/' + complist[1] + '/' + complist[1] + '.markdup.bam',
                             run_conf.project_root + '/control_panel/' + complist[2] + ".control_panel.txt"])


bam2fastq_output_list = []
for sample in sample_conf.bam_tofastq:
    bam2fastq_output_list.append([run_conf.project_root + '/fastq/' + sample + '/1_1.fastq',
                              run_conf.project_root + '/fastq/' + sample + '/1_2.fastq'])
sample_list_fastq = sample_conf.fastq
control_panel_keys = sample_conf.control_panel.keys()

# prepare output directories
if not os.path.isdir(run_conf.project_root): os.mkdir(run_conf.project_root)
if not os.path.isdir(run_conf.project_root + '/script'): os.mkdir(run_conf.project_root + '/script')
if not os.path.isdir(run_conf.project_root + '/log'): os.mkdir(run_conf.project_root + '/log')
if not os.path.isdir(run_conf.project_root + '/fastq'): os.mkdir(run_conf.project_root + '/fastq')
if not os.path.isdir(run_conf.project_root + '/bam'): os.mkdir(run_conf.project_root + '/bam')
if not os.path.isdir(run_conf.project_root + '/mutation'): os.mkdir(run_conf.project_root + '/mutation')
if not os.path.isdir(run_conf.project_root + '/control_panel'): os.mkdir(run_conf.project_root + '/control_panel')

# bamtofastq
# @originate(bam2fastq_output_list, bam2fastq_input_list)
@originate(bam2fastq_output_list)
def bam2fastq(outputfiles):
    
    sample = os.path.basename(os.path.dirname(outputfiles[0]))
    output_dir = run_conf.project_root + '/fastq/' + sample
    bam_dir = run_conf.project_root + '/bam/' + sample
    if not os.path.isdir(output_dir): os.mkdir(output_dir)
    if not os.path.isdir(bam_dir): os.mkdir(bam_dir)
            
    arguments = {"biobambam": genomon_conf.get("SOFTWARE", "biobambam"),
                 "input_bam": sample_conf.bam_tofastq[sample],
                 "f1_name": outputfiles[0],
                 "f2_name": outputfiles[1],
                 "o1_name": output_dir + '/unmatched_first_output.txt',
                 "o2_name": output_dir + '/unmatched_second_output.txt',
                 "t": output_dir + '/temp.txt',
                 "s": output_dir + '/single_end_output.txt',
                 "log": run_conf.project_root + '/log'}

    bamtofastq.task_exec(arguments)

# link the input fastq files
@originate(linked_fastq_list, sample_list_fastq)
def link_input_fastq(output_file, sample_list_fastq):
    sample = os.path.basename(os.path.dirname(output_file[0]))
    link_dir = run_conf.project_root + '/fastq/' + sample
    bam_dir = run_conf.project_root + '/bam/' + sample
   
    if not os.path.isdir(link_dir): os.mkdir(link_dir)
    if not os.path.isdir(bam_dir): os.mkdir(bam_dir)

    fastq_prefix, ext = os.path.splitext(sample_list_fastq[sample][0][0])
    # Todo
    # 1. should compare the timestamps between input and linked file
    # 2. check md5sum ?
    if not os.path.exists(link_dir + '/1_1' + ext): os.symlink(sample_list_fastq[sample][0][0], link_dir + '/1_1' + ext)
    if not os.path.exists(link_dir + '/1_2' + ext): os.symlink(sample_list_fastq[sample][1][0], link_dir + '/1_2' + ext)


##################
#  split stage
@subdivide([bam2fastq, link_input_fastq], formatter(), "{path[0]}/*_*.fastq_split", "{path[0]}")
def split_files(input_files, output_files, output_name_stem):
  
    for oo in output_files:
        os.unlink(oo)

    pair_id = 0
    for input_file in input_files:
        pair_id += 1

        arguments = {"lines": task_conf.get("split_fast", "split_fastq_line_number"),
                     "input_file": input_files[(pair_id - 1)],
                     "out_name_stem": output_name_stem,
                     "pair_id": pair_id,
                     "log": run_conf.project_root + '/log'}

        fastq_splitter.task_exec(arguments)


###################
# mapping stage
# @transform(split_files, formatter(".+/1_(?P<NAME>[0-9]+).fastq_split"), add_inputs("{path[0]}/2_{NAME[0]}.fastq_split"), "{subpath[0][2]}/bam/{subdir[0][0]}/{NAME[0]}.bamsorted.bam", {subdir[0][0]})
@transform(split_files, formatter(".+/1_(?P<NAME>[0-9]+).fastq_split"), add_inputs("{path[0]}/2_{NAME[0]}.fastq_split"), "{subpath[0][2]}/bam/{subdir[0][0]}/{NAME[0]}.sorted.bam")
def map_dna_sequence(input_files, output_file):
   
#    if os.path.exists(run_conf.project_root + '/bam/' + key + '/' + key + '.bammarkdup.bam'): continue

    dir_name = os.path.dirname(output_file)
    output_bwa_sam = output_file.replace('.sorted.bam', '.bwa.sam')
         
    arguments = {"bwa": genomon_conf.get("SOFTWARE", "bwa"),
                 "bwa_params": task_conf.get("bwa_mem", "bwa_params"),
                 "ref_fa":genomon_conf.get("REFERENCE", "ref_fasta"),
                 "fastq1": input_files[0],
                 "fastq2": input_files[1],
                 "sam": output_bwa_sam,
                 "biobambam": genomon_conf.get("SOFTWARE", "biobambam"),
                 "bam": output_file,
                 "log": run_conf.project_root + '/log'}

    bwa_align.task_exec(arguments) 
 

###################
# merge stage
@collate(map_dna_sequence, formatter(".+/(?P<SAMPLE>.+)/([0-9]+).sorted.bam"), "{path[0]}/{SAMPLE[0]}.markdup.bam")
def markdup(input_files, output_file):

    output_prefix, ext = os.path.splitext(output_file)

    input_bam_files = ""
    for input_file in input_files:
        input_bam_files = input_bam_files + " I=" + input_file

    arguments = {"biobambam": genomon_conf.get("SOFTWARE", "biobambam"),
                 "out_prefix": output_prefix,
                 "input_bam_files": input_bam_files,
                 "out_bam": output_file,
                 "log": run_conf.project_root + '/log'}

    markduplicates.task_exec(arguments)


# link the input fastq files
@originate(control_panel_keys)
def make_control_panel_file(control_panel_name):
    control_panel_file = run_conf.project_root + '/control_panel/' + control_panel_name + ".control_panel.txt"
    with open(control_panel_file,  "w") as out_handle:
        for bam in sample_conf.get_control_panel_list(control_panel_name):
            out_handle.write(bam + "\n")

##################
#  split stage
###################
# mutation calling stage
# @transform(markdup_bam_list, formatter(".+/(?P<SAMPLE>.+).markdup.bam"), "{subpath[0][2]}/mutation/{subdir[0][0]}/{SAMPLE[0]}.candidate_mutations.txt")
# def identify_mutations(input_files, output_file):
@follows( make_control_panel_file )
@follows( markdup )
@subdivide(markdup_bam_list, formatter(".+/(?P<SAMPLE>.+).markdup.bam"), "{subpath[0][2]}/mutation/{subdir[0][0]}/{SAMPLE[0]}_mutations_candidate.*.hg19_multianno.txt", "{subpath[0][2]}/mutation/{subdir[0][0]}")
def identify_mutations(input_files, output_files, output_dir):
    print "        mutation %s -> %s" % (input_files, output_files)

    # for oo in output_files:
    #    os.unlink(oo)

    sample_name = os.path.basename(output_dir)

    arguments = {
        # fisher mutation
        "fisher": genomon_conf.get("SOFTWARE", "fisher"),
        "map_quality": task_conf.get("fisher_mutation_call", "map_quality"),
        "base_quality": task_conf.get("fisher_mutation_call", "base_quality"),
        "min_allele_freq": task_conf.get("fisher_mutation_call", "disease_min_allele_frequency"),
        "max_allele_freq": task_conf.get("fisher_mutation_call", "control_max_allele_frequency"),
        "min_depth": task_conf.get("fisher_mutation_call", "min_depth"),
        # realignment filter
        "mutfilter": genomon_conf.get("SOFTWARE", "mutfilter"),
        "realign_min_mismatch": task_conf.get("realignment_filter","disease_min_mismatch"),
        "realign_max_mismatch": task_conf.get("realignment_filter","control_max_mismatch"),
        "realign_score_diff": task_conf.get("realignment_filter","score_diff"),
        "realign_window_size": task_conf.get("realignment_filter","window_size"),
        # indel filter
        "indel_search_length": task_conf.get("indel_filter","search_length"),
        "indel_neighbor": task_conf.get("indel_filter","neighbor"),
        "indel_base_quality": task_conf.get("indel_filter","base_quality"),
        "indel_min_depth": task_conf.get("indel_filter","min_depth"),
        "indel_min_mismatch": task_conf.get("indel_filter","max_mismatch"),
        "indel_min_allele_freq": task_conf.get("indel_filter","max_allele_freq"),
        # breakpoint filter
        "bp_max_depth": task_conf.get("breakpoint_filter","max_depth"),
        "bp_min_clip_size": task_conf.get("breakpoint_filter","min_clip_size"),
        "bp_junc_num_thres": task_conf.get("breakpoint_filter","junc_num_thres"),
        "bp_map_quality": task_conf.get("breakpoint_filter","map_quality"),
        # eb_filter
        "EBFilter": genomon_conf.get("SOFTWARE", "ebfilter"),
        "eb_map_quality": task_conf.get("eb_filter","map_quality"),
        "eb_base_quality": task_conf.get("eb_filter","base_quality"),
        "control_bam_list": input_files[2],
        # annovar
        "active_annovar_flag": "True",
        "annovar": genomon_conf.get("SOFTWARE", "annovar"),
        "table_annovar_params": task_conf.get("annotation", "table_annovar_params"),
        # commmon
        "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
        "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),   
        "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
        "interval_list": genomon_conf.get("REFERENCE", "interval_list"),
        "ref_fa":genomon_conf.get("REFERENCE", "ref_fasta"),
        "simple_repeat_db":genomon_conf.get("REFERENCE", "simple_repeat_tabix_db"),
        "disease_bam": input_files[0],
        "control_bam": input_files[1],
        "out_prefix": output_dir + '/' + sample_name,
        "samtools": genomon_conf.get("SOFTWARE", "samtools"),
        "blat": genomon_conf.get("SOFTWARE", "blat"),
        "log": run_conf.project_root + '/log'}

    interval_list = genomon_conf.get("REFERENCE", "interval_list")
    num_lines = sum(1 for line in open(interval_list))

    if not os.path.isdir(output_dir): os.mkdir(output_dir)
    mutation_call.task_exec(arguments, num_lines)


@collate(identify_mutations, formatter(".+/(?P<SAMPLE>.+)_mutations_candidate.(?P<NAME>[0-9]+).hg19_multianno.txt"), "{path[0]}/{SAMPLE[0]}_mutations.txt")
def merge_mutation(input_files, output_file):
    with open(output_file,  "w") as out_handle:
        for input_file in input_files:
            with open(input_file) as in_handle:
                for line in in_handle:
                   out_handle.write(line)
            
        

# @transform(identify_mutations. formatter(".+/(?P<SAMPLE>.+).candidate_mutations.(?P<NAME>[0-9]+).txt"), "{SAMPLE[0]}.annotated_mutations.{NAME[0]}.txt")
# @transform(map_dna_sequence,   formatter(".+/(?P<SAMPLE>.+)/([0-9]+).sorted.bam"), "{path[0]}/{SAMPLE[0]}.markdup.bam")
# @transform(identify_mutations, formatter(".+/(?P<SAMPLE>.+)_mutations_candidate.(?P<NAME>[0-9]+).txt"), "{path[0]}/{SAMPLE[0]}_mutation_annotation.{NAME[0]}.txt", "{path[0]}")
# def annotator(input_file, output_file, output_dir):
#     print "        annotator %s -> %s" % (input_file, output_file)
#     print "        annotator %s -> %s" % (input_file, output_dir)
# 
#     # if (run_conf.annovar_path == ""):
#         # shutil.copy(input_file, output_file)
#     # else:
#     sample_name = os.path.basename(output_dir)
#     arguments = {
#                 "annovar": "/home/kchiba/work_genomonprj/annovar",
#                 "table_annovar_params": task_conf.get("annotation", "table_annovar_params"),
#                 "input_file": input_file,
#                 "out_prefix": output_dir + '/' + sample_name,
#                 "log": run_conf.project_root + '/log'}
# 
#     annotation.task_exec(arguments)
 