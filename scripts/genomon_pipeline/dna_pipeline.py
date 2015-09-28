import os
import shutil
import yaml
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
from genomon_pipeline.dna_resource.sv_parse import *
from genomon_pipeline.dna_resource.sv_merge import *
from genomon_pipeline.dna_resource.sv_filt import *


# set task classes
fastq_splitter = Fastq_splitter(task_conf.get("split_fast", "qsub_option"), run_conf.project_root + '/script')
bwa_align = Bwa_align(task_conf.get("bwa_mem", "qsub_option"), run_conf.project_root + '/script')
markduplicates = Markduplicates(task_conf.get("markduplicates", "qsub_option"), run_conf.project_root + '/script')
mutation_call = Mutation_call(task_conf.get("mutation_call", "qsub_option"), run_conf.project_root + '/script')
bamtofastq = Bam2Fastq(task_conf.get("bam2fastq", "qsub_option"), run_conf.project_root + '/script')
sv_parse = SV_parse(task_conf.get("sv_parse", "qsub_option"), run_conf.project_root + '/script')
sv_merge = SV_merge(task_conf.get("sv_merge", "qsub_option"), run_conf.project_root + '/script')
sv_filt = SV_filt(task_conf.get("sv_filt", "qsub_option"), run_conf.project_root + '/script')

# generate list of linked_fastq file path
linked_fastq_list = []
for sample in sample_conf.fastq:
    fastq_prefix, ext = os.path.splitext(sample_conf.fastq[sample][0][0])
    linked_fastq_list.append([run_conf.project_root + '/fastq/' + sample + '/1_1' + ext,
                              run_conf.project_root + '/fastq/' + sample + '/1_2' + ext])

bam2fastq_output_list = []
for sample in sample_conf.bam_tofastq:
    bam2fastq_output_list.append([run_conf.project_root + '/fastq/' + sample + '/1_1.fastq',
                                  run_conf.project_root + '/fastq/' + sample + '/1_2.fastq'])

markdup_bam_list = []
for complist in sample_conf.compare:
    markdup_bam_list.append([run_conf.project_root + '/bam/' + complist[0] + '/' + complist[0] + '.markdup.bam',
                             run_conf.project_root + '/bam/' + complist[1] + '/' + complist[1] + '.markdup.bam',
                             run_conf.project_root + '/mutation/control_panel/' + complist[2] + ".control_panel.txt"])

parse_bedpe_list = []
for complist in sample_conf.compare:
    parse_bedpe_list.append(run_conf.project_root+ "/sv/"+ complist[0] +"/"+ complist[0] +".junction.clustered.bedpe.gz")

control_bedpe_list = []
for control_panel_name in sample_conf.control_panel.keys():
    tmp_list = []
    tmp_list.append(run_conf.project_root + '/sv/config/' + control_panel_name + ".control.yaml")
    for sample in sample_conf.control_panel[control_panel_name]:
       tmp_list.append(run_conf.project_root+ "/sv/"+ sample +"/"+ sample +".junction.clustered.bedpe.gz")
    control_bedpe_list.append(tmp_list)


# prepare output directories
if not os.path.isdir(run_conf.project_root): os.makedirs(run_conf.project_root)
if not os.path.isdir(run_conf.project_root + '/script'): os.mkdir(run_conf.project_root + '/script')
if not os.path.isdir(run_conf.project_root + '/log'): os.mkdir(run_conf.project_root + '/log')
if not os.path.isdir(run_conf.project_root + '/fastq'): os.mkdir(run_conf.project_root + '/fastq')
if not os.path.isdir(run_conf.project_root + '/bam'): os.mkdir(run_conf.project_root + '/bam')
if not os.path.isdir(run_conf.project_root + '/mutation'): os.mkdir(run_conf.project_root + '/mutation')
if not os.path.isdir(run_conf.project_root + '/mutation/control_panel'): os.mkdir(run_conf.project_root + '/mutation/control_panel')
if not os.path.isdir(run_conf.project_root + '/sv'): os.mkdir(run_conf.project_root + '/sv')
if not os.path.isdir(run_conf.project_root + '/sv/non_matched_control_panel'): os.mkdir(run_conf.project_root + '/sv/non_matched_control_panel')
if not os.path.isdir(run_conf.project_root + '/sv/config'): os.mkdir(run_conf.project_root + '/sv/config')


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


#############################
# link the input fastq files
@originate(linked_fastq_list, sample_conf.fastq)
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


############################
# link the import bam files
@originate(sample_conf.bam_import.keys())
def link_import_bam(sample):
    bam = sample_conf.bam_import[sample]
    link_dir = run_conf.project_root + '/bam/' + sample
    bam_prefix, ext = os.path.splitext(bam)
    
    if not os.path.isdir(link_dir): os.mkdir(link_dir)

    if not os.path.exists(link_dir +'/'+ sample +'.bam'): os.symlink(bam, link_dir +'/'+ sample +'.bam')
    if not os.path.exists(link_dir +'/'+ sample +'.bam.bai'): 
        if (os.path.exists(bam +'.bai')):
            os.symlink(bam +'.bai', link_dir +'/'+ sample +'.bam.bai')
        elif (os.path.exists(bam_prefix +'.bai')):
            os.symlink(bam_prefix +'.bai', link_dir +'/'+ sample +'.bam.bai')


##########################
# make control panel txt
@originate(sample_conf.control_panel.keys())
def make_control_panel_file(control_panel_name):
    control_panel_file = run_conf.project_root + '/mutation/control_panel/' + control_panel_name + ".control_panel.txt"
    with open(control_panel_file,  "w") as out_handle:
        for bam in sample_conf.get_control_panel_list(control_panel_name):
            out_handle.write(bam + "\n")


##########################
# make sv control panel txt
@originate(sample_conf.control_panel.keys())
def make_sv_control_panel_file(control_panel_name):
    control_conf = run_conf.project_root + '/sv/config/' + control_panel_name + ".control.yaml"
    with open(control_conf,  "w") as out_handle:
        for sample in sample_conf.control_panel[control_panel_name]:
            out_handle.write(sample +": "+run_conf.project_root+ "/sv/"+ sample +"/"+ sample +".junction.clustered.bedpe.gz\n")


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
@transform(split_files, formatter(".+/1_(?P<NAME>[0-9]+).fastq_split"), add_inputs("{path[0]}/2_{NAME[0]}.fastq_split"), "{subpath[0][2]}/bam/{subdir[0][0]}/{NAME[0]}.sorted.bam")
def map_dna_sequence(input_files, output_file):
   
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



###################
# mutation calling stage
@follows( link_import_bam )
@follows( make_control_panel_file )
@follows( markdup )
@subdivide(markdup_bam_list, formatter(".+/(?P<SAMPLE>.+).markdup.bam"), "{subpath[0][2]}/mutation/{subdir[0][0]}/{subdir[0][0]}_mutations_candidate.*.hg19_multianno.txt", "{subpath[0][2]}/mutation/{subdir[0][0]}")
def identify_mutations(input_files, output_files, output_dir):
    print "        mutation %s -> %s" % (input_files, output_files)

    for oo in output_files:
       os.unlink(oo)

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
        # eb filter
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


###################
# merge mutation result
@collate(identify_mutations, formatter(".+/(?P<SAMPLE>.+)_mutations_candidate.(?P<NAME>[0-9]+).hg19_multianno.txt"), "{path[0]}/{SAMPLE[0]}_genomon_mutations.result.txt")
def merge_mutation(input_files, output_file):
    with open(output_file,  "w") as out_handle:
        for input_file in input_files:
            with open(input_file) as in_handle:
                next(in_handle)
                for line in in_handle:
                   out_handle.write(line)


###################
# SV parse
@follows( link_import_bam )
@follows( markdup )
@transform(sample_conf.get_disease_and_control_panel_bam(), formatter(), "{subpath[0][2]}/sv/{subdir[0][0]}/{subdir[0][0]}.junction.clustered.bedpe.gz")
def parse_sv(input_file, output_file):
    print "        parse %s -> %s" % (input_file, output_file)

    dir_name = os.path.dirname(output_file)
    sample_name = os.path.basename(dir_name)

    if not os.path.isdir(dir_name): os.mkdir(dir_name)
    idx = -1
    sv_sampleConf = {"target": {}, "matched_control": {}, "non_matched_control_panel": {}}
    for complist in sample_conf.compare:
        if sample_name in complist:
            idx = complist.index(sample_name)

        if idx == 0: # is disease
            sv_sampleConf["target"]["label"] = sample_name
            sv_sampleConf["target"]["path_to_bam"] = input_file
            sv_sampleConf["target"]["path_to_output_dir"] = dir_name
            sv_sampleConf["matched_control"]["use"] = True
            sv_sampleConf["matched_control"]["path_to_bam"] = sample_conf.sample2bam(complist[1])
            sv_sampleConf["non_matched_control_panel"]["use"] = True
            sv_sampleConf["non_matched_control_panel"]["matched_control_label"] = complist[1]
            sv_sampleConf["non_matched_control_panel"]["data_path"] = run_conf.project_root +"/sv/non_matched_control_panel/"+ complist[2] +".merged.junction.control.bedpe.gz"
            break

    if idx != 0: # are control or non matched control
        sv_sampleConf["target"]["label"] = sample_name
        sv_sampleConf["target"]["path_to_bam"] = input_file
        sv_sampleConf["target"]["path_to_output_dir"] = dir_name
        sv_sampleConf["matched_control"]["use"] = False
        sv_sampleConf["non_matched_control_panel"]["use"] = False

    sample_yaml = run_conf.project_root + "/sv/config/" + sample_name + ".yaml"
    hOUT = open(sample_yaml, "w")
    print >> hOUT, yaml.dump(sv_sampleConf, default_flow_style = False)
    hOUT.close()

    arguments = {"genomon_sv": genomon_conf.get("SOFTWARE", "genomon_sv"),
                 "sample_conf": sample_yaml,
                 "param_conf": task_conf.get("genomon_sv", "param_file"),
                 "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),   
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "log": run_conf.project_root + '/log'}
    sv_parse.task_exec(arguments)


###################
# SV merge
@follows( make_sv_control_panel_file )
@follows( parse_sv )
@transform(control_bedpe_list, formatter(".+/(?P<NAME>.+).control.yaml",".+/(.+).junction.clustered.bedpe.gz"), "{subpath[0][2]}/sv/non_matched_control_panel/{NAME[0]}.merged.junction.control.bedpe.gz")
def merge_sv(input_files,  output_file):
    print "        merge %s -> %s" % (input_files, output_file)

    arguments = {"genomon_sv": genomon_conf.get("SOFTWARE", "genomon_sv"),
                 "control_conf": input_files[0],
                 "bedpe": output_file,
                 "param_conf": task_conf.get("genomon_sv", "param_file"),
                 "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),   
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "log": run_conf.project_root + '/log'}
    sv_merge.task_exec(arguments)


@follows( merge_sv )
@transform(parse_bedpe_list, formatter(), "{subpath[0][2]}/sv/{subdir[0][0]}/{subdir[0][0]}.genomonSV.result.txt")
def filt_sv(input_files,  output_file):
    print "        filt %s -> %s" % (input_files, output_file)

    dir_name = os.path.dirname(output_file)
    sample_name = os.path.basename(dir_name)

    sample_yaml = run_conf.project_root + "/sv/config/" + sample_name + ".yaml"

    arguments = {"genomon_sv": genomon_conf.get("SOFTWARE", "genomon_sv"),
                 "sample_conf": sample_yaml,
                 "param_conf": task_conf.get("genomon_sv", "param_file"),
                 "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),   
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "log": run_conf.project_root + '/log'}
    sv_filt.task_exec(arguments)


