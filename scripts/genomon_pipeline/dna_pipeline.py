import os
import shutil
import yaml
import linecache
import glob
from ruffus import *
from genomon_pipeline.config.run_conf import *
from genomon_pipeline.config.genomon_conf import *
from genomon_pipeline.config.task_conf import *
from genomon_pipeline.config.sample_conf import *
from genomon_pipeline.dna_resource.bamtofastq import *
from genomon_pipeline.dna_resource.fastq_splitter import *
from genomon_pipeline.dna_resource.bwa_align import *
from genomon_pipeline.dna_resource.markduplicates import *
from genomon_pipeline.dna_resource.mutation_call import *
from genomon_pipeline.dna_resource.mutation_merge import *
from genomon_pipeline.dna_resource.sv_parse import *
from genomon_pipeline.dna_resource.sv_merge import *
from genomon_pipeline.dna_resource.sv_filt import *
from genomon_pipeline.dna_resource.bam_stats import *
from genomon_pipeline.dna_resource.coverage import *
from genomon_pipeline.dna_resource.merge import *

# set task classes
bamtofastq = Bam2Fastq(task_conf.get("bam2fastq", "qsub_option"), run_conf.drmaa)
fastq_splitter = Fastq_splitter(task_conf.get("split_fastq", "qsub_option"), run_conf.drmaa)
bwa_align = Bwa_align(task_conf.get("bwa_mem", "qsub_option"), run_conf.drmaa)
markduplicates = Markduplicates(task_conf.get("markduplicates", "qsub_option"), run_conf.drmaa)
mutation_call = Mutation_call(task_conf.get("mutation_call", "qsub_option"), run_conf.drmaa)
mutation_merge = Mutation_merge(task_conf.get("mutation_merge", "qsub_option"), run_conf.drmaa)
sv_parse = SV_parse(task_conf.get("sv_parse", "qsub_option"), run_conf.drmaa)
sv_merge = SV_merge(task_conf.get("sv_merge", "qsub_option"), run_conf.drmaa)
sv_filt = SV_filt(task_conf.get("sv_filt", "qsub_option"), run_conf.drmaa)
r_bamstats = Res_Bamstats(task_conf.get("bam_stats", "qsub_option"), run_conf.drmaa)
r_coverage = Res_Coverage(task_conf.get("coverage", "qsub_option"), run_conf.drmaa)
r_merge = Res_Merge(task_conf.get("merge", "qsub_option"), run_conf.drmaa)
r_pa_capture = Res_PA_Capture(task_conf.get("pa_capture", "qsub_option"), run_conf.drmaa)
r_pa_merge = Res_PA_Merge(task_conf.get("pa_merge", "qsub_option"), run_conf.drmaa)

# generate output list of 'linked fastq'
linked_fastq_list = []
for sample in sample_conf.fastq:
    if os.path.exists(run_conf.project_root + '/bam/' + sample + '/1.sorted.bam'): continue
    if os.path.exists(run_conf.project_root + '/bam/' + sample + '/' + sample + '.markdup.bam'): continue

    link_fastq_arr1 = []
    link_fastq_arr2 = []
    for (count, fastq_file) in enumerate(sample_conf.fastq[sample][0]):
        fastq_prefix, ext = os.path.splitext(fastq_file)
        link_fastq_arr1.append(run_conf.project_root + '/fastq/' + sample + '/' + str(count+1) + '_1' + ext)
        link_fastq_arr2.append(run_conf.project_root + '/fastq/' + sample + '/' + str(count+1) + '_2' + ext)
    linked_fastq_list.append([link_fastq_arr1,link_fastq_arr2])

# generate output list of 'bam2fastq'
bam2fastq_output_list = []
for sample in sample_conf.bam_tofastq:
    if os.path.exists(run_conf.project_root + '/bam/' + sample + '/1.sorted.bam'): continue
    if os.path.exists(run_conf.project_root + '/bam/' + sample + '/' + sample + '.markdup.bam'): continue
    bam2fastq_arr1 = []
    bam2fastq_arr2 = []
    bam2fastq_arr1.append(run_conf.project_root + '/fastq/' + sample + '/1_1.fastq')
    bam2fastq_arr2.append(run_conf.project_root + '/fastq/' + sample + '/1_2.fastq')
    bam2fastq_output_list.append([bam2fastq_arr1,bam2fastq_arr2])

# generate input list of 'mutation call'
markdup_bam_list = []
merge_mutation_list = []
for complist in sample_conf.mutation_call:
     if os.path.exists(run_conf.project_root + '/mutation/' + complist[0] + '/' + complist[0] + '_genomon_mutations.result.txt'): continue
     tumor_bam  = run_conf.project_root + '/bam/' + complist[0] + '/' + complist[0] + '.markdup.bam'
     normal_bam = run_conf.project_root + '/bam/' + complist[1] + '/' + complist[1] + '.markdup.bam' if complist[1] != None else None
     panel = run_conf.project_root + '/mutation/control_panel/' + complist[2] + ".control_panel.txt" if complist[2] != None else None
     markdup_bam_list.append([tumor_bam, normal_bam, panel])


# generate input list of 'SV parse'
parse_sv_bam_list = []
all_target_bams = []
unique_bams = []
for complist in sample_conf.sv_detection:
    tumor_sample = complist[0]
    if tumor_sample != None:
        all_target_bams.append(run_conf.project_root + '/bam/' + tumor_sample + '/' + tumor_sample + '.markdup.bam')
    normal_sample = complist[1]
    if normal_sample != None:
        all_target_bams.append(run_conf.project_root + '/bam/' + normal_sample + '/' + normal_sample + '.markdup.bam')
    panel_name = complist[2]
    if panel_name != None:
        for panel_sample in sample_conf.control_panel[panel_name]:
            all_target_bams.append(run_conf.project_root + '/bam/' + panel_sample + '/' + panel_sample + '.markdup.bam')
    unique_bams = list(set(all_target_bams))       
    
for bam in unique_bams:
    dir_name = os.path.dirname(bam)
    sample_name = os.path.basename(dir_name)
    if os.path.exists(run_conf.project_root + '/sv/' + sample_name + '/' + sample_name + '.junction.clustered.bedpe.gz'): continue
    parse_sv_bam_list.append(bam)

# generate input list of 'SV merge'
unique_complist = []
merge_bedpe_list = []
for complist in sample_conf.sv_detection:
    control_panel_name = complist[2]
    if control_panel_name != None and control_panel_name not in unique_complist:
        unique_complist.append(control_panel_name)

for control_panel_name in unique_complist:
    if os.path.exists(run_conf.project_root + '/sv/non_matched_control_panel/' + control_panel_name + '.merged.junction.control.bedpe.gz'): continue
    tmp_list = []
    tmp_list.append(run_conf.project_root + '/sv/config/' + control_panel_name + ".control.yaml")
    for sample in sample_conf.control_panel[control_panel_name]:
        tmp_list.append(run_conf.project_root+ "/sv/"+ sample +"/"+ sample +".junction.clustered.bedpe.gz")
    merge_bedpe_list.append(tmp_list)

# generate input list of 'SV filt'
filt_bedpe_list = []
for complist in sample_conf.sv_detection:
    if os.path.exists(run_conf.project_root + '/sv/' + complist[0] +'/'+ complist[0] +'.genomonSV.result.txt'): continue
    filt_bedpe_list.append(run_conf.project_root+ "/sv/"+ complist[0] +"/"+ complist[0] +".junction.clustered.bedpe.gz")

# generate input list of 'summary'
summary_bamstats_list = []
summary_coverage_list = []
summary_merge_list = []
for sample in sample_conf.summary:
    if os.path.exists(run_conf.project_root + '/summary/' + sample + '/' + sample + '.tsv'): continue
    summary_merge_list.append(
        [run_conf.project_root + '/summary/' + sample + '/' + sample + '.bamstats',
         run_conf.project_root + '/summary/' + sample + '/' + sample + '.coverage'])
    if not os.path.exists(run_conf.project_root + '/summary/' + sample + '/' + sample + '.bamstats'):
        summary_bamstats_list.append(run_conf.project_root + '/bam/' + sample +'/'+ sample +'.markdup.bam')
    if not os.path.exists(run_conf.project_root + '/summary/' + sample + '/' + sample + '.coverage'):
        summary_coverage_list.append(run_conf.project_root + '/bam/' + sample +'/'+ sample +'.markdup.bam')

# generate input list of 'post analysis for mutation'
pa_dir = run_conf.project_root + '/post_analysis/mutation'
pa_capt_list_mutation = []
for complist in sample_conf.mutation_call:
    if (complist[1] == None): continue
    if (not os.path.exists(pa_dir + '/capture_script/' + complist[0] + '.bat')) \
       or (not os.path.exists(pa_dir + '/bam/' + complist[0] + '/' + complist[0] + '.markdup.pickup.bam')):
        pa_capt_list_mutation.append(run_conf.project_root + '/mutation/' + complist[0] + '/' + complist[0] + '_genomon_mutations.result.txt')

pa_merge_list_mutation = []
if (not os.path.exists(pa_dir + '/capture_script/capture.bat')) \
  or (not os.path.exists(run_conf.project_root + '/post_analysis/merge.mutation.csv')):
    for complist in sample_conf.mutation_call:
        if (complist[1] == None): continue
        pa_merge_list_mutation.append(run_conf.project_root + '/mutation/' + complist[0] + '/' + complist[0] + '_genomon_mutations.result.txt')

# generate input list of 'post analysis for SV'
pa_dir = run_conf.project_root + '/post_analysis/sv'
pa_capt_list_sv = []
for complist in sample_conf.sv_detection:
    if (complist[1] == None): continue
    if (not os.path.exists(pa_dir + '/capture_script/' + complist[0] + '.bat')) \
      or (not os.path.exists(pa_dir + '/bam/' + complist[0] + '/' + complist[0] + '.markdup.pickup.bam')):
        pa_capt_list_sv.append(run_conf.project_root + '/sv/' + complist[0] +'/'+ complist[0] +'.genomonSV.result.txt')

pa_merge_list_sv = []
if (not os.path.exists(pa_dir + '/capture_script/capture.bat')) \
  or (not os.path.exists(run_conf.project_root + '/post_analysis/merge.sv.csv')):
    for complist in sample_conf.sv_detection:
        if (complist[1] == None): continue
        pa_merge_list_sv.append(run_conf.project_root + '/sv/' + complist[0] +'/'+ complist[0] +'.genomonSV.result.txt')

# generate input list of 'post analysis for summary'
pa_summary_list = []
if not os.path.exists(run_conf.project_root + '/post_analysis/merge.summary.csv'):
    for sample in sample_conf.summary:
        pa_summary_list.append(run_conf.project_root + '/summary/' + sample + '/' + sample + '.tsv')

# prepare output directories
if not os.path.isdir(run_conf.project_root): os.mkdir(run_conf.project_root)
if not os.path.isdir(run_conf.project_root + '/script'): os.mkdir(run_conf.project_root + '/script')
if not os.path.isdir(run_conf.project_root + '/log'): os.mkdir(run_conf.project_root + '/log')
if not os.path.isdir(run_conf.project_root + '/fastq'): os.mkdir(run_conf.project_root + '/fastq')
if not os.path.isdir(run_conf.project_root + '/bam'): os.mkdir(run_conf.project_root + '/bam')
if not os.path.isdir(run_conf.project_root + '/mutation'): os.mkdir(run_conf.project_root + '/mutation')
if not os.path.isdir(run_conf.project_root + '/mutation/control_panel'): os.mkdir(run_conf.project_root + '/mutation/control_panel')
if not os.path.isdir(run_conf.project_root + '/sv'): os.mkdir(run_conf.project_root + '/sv')
if not os.path.isdir(run_conf.project_root + '/sv/non_matched_control_panel'): os.mkdir(run_conf.project_root + '/sv/non_matched_control_panel')
if not os.path.isdir(run_conf.project_root + '/sv/config'): os.mkdir(run_conf.project_root + '/sv/config')
if not os.path.isdir(run_conf.project_root + '/summary'): os.mkdir(run_conf.project_root + '/summary')
if (task_conf.getboolean("post_analysis", "enable") == True):
    if not os.path.isdir(run_conf.project_root + '/post_analysis/mutation/capture'):
        os.makedirs(run_conf.project_root + '/post_analysis/mutation/capture')
    if not os.path.isdir(run_conf.project_root + '/post_analysis/mutation/capture_script'):
        os.makedirs(run_conf.project_root + '/post_analysis/mutation/capture_script')
    if not os.path.isdir(run_conf.project_root + '/post_analysis/mutation/bam'):
        os.mkdir(run_conf.project_root + '/post_analysis/mutation/bam')
    if not os.path.isdir(run_conf.project_root + '/post_analysis/mutation/bam_script'):
        os.mkdir(run_conf.project_root + '/post_analysis/mutation/bam_script')
    if not os.path.isdir(run_conf.project_root + '/post_analysis/sv/capture'):
        os.makedirs(run_conf.project_root + '/post_analysis/sv/capture')
    if not os.path.isdir(run_conf.project_root + '/post_analysis/sv/capture_script'):
        os.makedirs(run_conf.project_root + '/post_analysis/sv/capture_script')
    if not os.path.isdir(run_conf.project_root + '/post_analysis/sv/bam'):
        os.mkdir(run_conf.project_root + '/post_analysis/sv/bam')
    if not os.path.isdir(run_conf.project_root + '/post_analysis/sv/bam_script'):
        os.mkdir(run_conf.project_root + '/post_analysis/sv/bam_script')
    
for outputfiles in (bam2fastq_output_list, linked_fastq_list):
    for outputfile in outputfiles:
        sample = os.path.basename(os.path.dirname(outputfile[0][0]))
        fastq_dir = run_conf.project_root + '/fastq/' + sample
        bam_dir = run_conf.project_root + '/bam/' + sample
        if not os.path.isdir(fastq_dir): os.mkdir(fastq_dir)
        if not os.path.isdir(bam_dir): os.mkdir(bam_dir)


# prepare output directory for each sample and make mutation control panel file
for complist in sample_conf.mutation_call:
    # make dir
    mutation_dir = run_conf.project_root + '/mutation/' + complist[0]
    if not os.path.isdir(mutation_dir): os.mkdir(mutation_dir)
    # make the control panel text 
    control_panel_name = complist[2]
    if control_panel_name != None:
        control_panel_file = run_conf.project_root + '/mutation/control_panel/' + control_panel_name + ".control_panel.txt"
        with open(control_panel_file,  "w") as out_handle:
            for panel_sample in sample_conf.control_panel[control_panel_name]:
                out_handle.write(run_conf.project_root + '/bam/' + panel_sample + '/' + panel_sample + '.markdup.bam' + "\n")

# make SV configuration file
for complist in sample_conf.sv_detection:
    # make the control yaml file
    control_panel_name = complist[2]
    if control_panel_name != None:
        control_conf = run_conf.project_root + '/sv/config/' + control_panel_name + ".control.yaml"
        with open(control_conf,  "w") as out_handle:
            for sample in sample_conf.control_panel[control_panel_name]:
                out_handle.write(sample +": "+run_conf.project_root+ "/sv/"+ sample +"/"+ sample +".junction.clustered.bedpe.gz\n")

# link the import bam to project directory
@originate(sample_conf.bam_import.keys())
def link_import_bam(sample):
    bam = sample_conf.bam_import[sample]
    link_dir = run_conf.project_root + '/bam/' + sample
    bam_prefix, ext = os.path.splitext(bam)
    
    if not os.path.isdir(link_dir): os.mkdir(link_dir)
    if (not os.path.exists(link_dir +'/'+ sample +'.markdup.bam')) and (not os.path.exists(link_dir +'/'+ sample +'.markdup.bam.bai')): 
        os.symlink(bam, link_dir +'/'+ sample +'.markdup.bam')
        if (os.path.exists(bam +'.bai')):
            os.symlink(bam +'.bai', link_dir +'/'+ sample +'.markdup.bam.bai')
        elif (os.path.exists(bam_prefix +'.bai')):
            os.symlink(bam_prefix +'.bai', link_dir +'/'+ sample +'.markdup.bam.bai')

# convert bam to fastq
@originate(bam2fastq_output_list)
def bam2fastq(outputfiles):
    sample = os.path.basename(os.path.dirname(outputfiles[0][0]))
    output_dir = run_conf.project_root + '/fastq/' + sample
            
    arguments = {"biobambam": genomon_conf.get("SOFTWARE", "biobambam"),
                 "input_bam": sample_conf.bam_tofastq[sample],
                 "f1_name": outputfiles[0][0],
                 "f2_name": outputfiles[1][0],
                 "o1_name": output_dir + '/unmatched_first_output.txt',
                 "o2_name": output_dir + '/unmatched_second_output.txt',
                 "t": output_dir + '/temp.txt',
                 "s": output_dir + '/single_end_output.txt'}
    bamtofastq.task_exec(arguments, run_conf.project_root + '/log', run_conf.project_root + '/script')


# link the input fastq to project directory
@originate(linked_fastq_list)
def link_input_fastq(output_file):
    sample = os.path.basename(os.path.dirname(output_file[0][0]))
    fastq_dir = run_conf.project_root + '/fastq/' + sample
    fastq_prefix, ext = os.path.splitext(sample_conf.fastq[sample][0][0])
    # Todo
    # 1. should compare the timestamps between input and linked file
    # 2. check md5sum ?
    for (count, fastq_files) in enumerate(sample_conf.fastq[sample][0]):
        fastq_prefix, ext = os.path.splitext(fastq_files)
        if not os.path.exists(fastq_dir + '/'+str(count+1)+'_1'+ ext): os.symlink(sample_conf.fastq[sample][0][count], fastq_dir + '/'+str(count+1)+'_1'+ ext)
        if not os.path.exists(fastq_dir + '/'+str(count+1)+'_2'+ ext): os.symlink(sample_conf.fastq[sample][1][count], fastq_dir + '/'+str(count+1)+'_2'+ ext)


# split fastq
@subdivide([bam2fastq, link_input_fastq], formatter(), "{path[0]}/*_*.fastq_split", "{path[0]}")
def split_files(input_files, output_files, target_dir):

    for oo in output_files:
        os.unlink(oo)

    input_prefix, ext = os.path.splitext(input_files[0][0])
    arguments = {"lines": task_conf.get("split_fastq", "split_fastq_line_number"),
                 "fastq_filter": task_conf.get("split_fastq", "fastq_filter"),
                 "target_dir": target_dir,
                 "ext": ext}
    
    fastq_splitter.task_exec(arguments, run_conf.project_root + '/log', run_conf.project_root + '/script', 2)
   
    all_line_num = 0
    for fastq in glob.glob(target_dir + '/1_*.fastq_split'):
        all_line_num += sum(1 for line in open(fastq))
    with open(target_dir + "/fastq_line_num.txt",  "w") as out_handle:
        out_handle.write(str(all_line_num)+"\n")
    
    for input_fastq in input_files[0]:
        os.unlink(input_fastq)
    for input_fastq in input_files[1]:
        os.unlink(input_fastq)


#bwa
@subdivide(split_files, formatter(".+/(.+)/1_0000.fastq_split"), add_inputs("{subpath[0][2]}/fastq/{subdir[0][0]}/2_0000.fastq_split"), "{subpath[0][2]}/bam/{subdir[0][0]}/{subdir[0][0]}_*.sorted.bam", "{subpath[0][2]}/fastq/{subdir[0][0]}", "{subpath[0][2]}/bam/{subdir[0][0]}")
def map_dna_sequence(input_files, output_files, input_dir, output_dir):

    sample_name = os.path.basename(output_dir)

    all_line_num = 0
    with open(input_dir + "/fastq_line_num.txt") as in_handle:
        tmp_num = in_handle.read()
        all_line_num = int(tmp_num)
    split_lines = task_conf.get("split_fastq", "split_fastq_line_number")

    ans_quotient = all_line_num / int(split_lines)
    ans_remainder = all_line_num % int(split_lines)
    max_task_id = ans_quotient if ans_remainder == 0 else ans_quotient + 1
    
    arguments = {"input_dir": input_dir,
                 "output_dir": output_dir,
                 "sample_name": sample_name,
                 "bwa": genomon_conf.get("SOFTWARE", "bwa"),
                 "bwa_params": task_conf.get("bwa_mem", "bwa_params"),
                 "ref_fa":genomon_conf.get("REFERENCE", "ref_fasta"),
                 "biobambam": genomon_conf.get("SOFTWARE", "biobambam")}

    bwa_align.task_exec(arguments, run_conf.project_root + '/log', run_conf.project_root + '/script', max_task_id) 

    for task_id in range(max_task_id):
        num = str(task_id).zfill(4)
        os.unlink(input_dir +'/1_'+str(num)+'.fastq_split')
        os.unlink(input_dir +'/2_'+str(num)+'.fastq_split')
        os.unlink(output_dir+'/'+sample_name+'_'+str(num)+'.bwa.sam')


# merge sorted bams into one and mark duplicate reads with biobambam
@collate(map_dna_sequence, formatter(), "{subpath[0][2]}/bam/{subdir[0][0]}/{subdir[0][0]}.markdup.bam")
def markdup(input_files, output_file):

    output_prefix, ext = os.path.splitext(output_file)

    input_bam_files = ""
    for input_file in input_files:
        input_bam_files = input_bam_files + " I=" + input_file

    arguments = {"biobambam": genomon_conf.get("SOFTWARE", "biobambam"),
                 "out_prefix": output_prefix,
                 "input_bam_files": input_bam_files,
                 "out_bam": output_file}

    markduplicates.task_exec(arguments, run_conf.project_root + '/log', run_conf.project_root + '/script')

    for input_file in input_files:
        os.unlink(input_file)
        os.unlink(input_file + ".bai")


# identify mutations
@follows( markdup )
@follows( link_import_bam )
@transform(markdup_bam_list, formatter(), "{subpath[0][2]}/mutation/{subdir[0][0]}/{subdir[0][0]}_genomon_mutations.result.txt", "{subpath[0][2]}/mutation/{subdir[0][0]}")
def identify_mutations(input_file, output_file, output_dir):

    sample_name = os.path.basename(output_dir)

    active_inhouse_normal_flag = False
    inhouse_normal_tabix_db = ""
    try:
        active_inhouse_normal_flag = task_conf.get("annotation", "active_inhouse_normal_flag")
        inhouse_normal_tabix_db = genomon_conf.get("REFERENCE", "inhouse_normal_tabix_db")
    except ConfigParser.NoOptionError:
        print "assumed cause: [inhouse normal] is not defined"

    active_inhouse_tumor_flag = False
    inhouse_tumor_tabix_db = ""
    try:
        active_inhouse_tumor_flag = task_conf.get("annotation", "active_inhouse_tumor_flag")
        inhouse_tumor_tabix_db = genomon_conf.get("REFERENCE", "inhouse_tumor_tabix_db")
    except ConfigParser.NoOptionError:
        print "assumed cause: [inhouse tumor] is not defined"

    active_HGMD_flag = False
    HGMD_tabix_db = ""
    try:
        active_HGMD_flag = task_conf.get("annotation", "active_HGMD_flag")
        HGMD_tabix_db = genomon_conf.get("REFERENCE", "HGMD_tabix_db")
    except ConfigParser.NoOptionError:
        print "assumed cause: [HGMD] is not defined"

    arguments = {
        # fisher mutation
        "fisher": genomon_conf.get("SOFTWARE", "fisher"),
        "map_quality": task_conf.get("fisher_mutation_call", "map_quality"),
        "base_quality": task_conf.get("fisher_mutation_call", "base_quality"),
        "min_allele_freq": task_conf.get("fisher_mutation_call", "disease_min_allele_frequency"),
        "max_allele_freq": task_conf.get("fisher_mutation_call", "control_max_allele_frequency"),
        "min_depth": task_conf.get("fisher_mutation_call", "min_depth"),
        "fisher_thres": task_conf.get("fisher_mutation_call", "fisher_thres_hold"),
        "post_10_q": task_conf.get("fisher_mutation_call", "post_10_q"),
        # realignment filter
        "mutfilter": genomon_conf.get("SOFTWARE", "mutfilter"),
        "realign_min_mismatch": task_conf.get("realignment_filter","disease_min_mismatch"),
        "realign_max_mismatch": task_conf.get("realignment_filter","control_max_mismatch"),
        "realign_score_diff": task_conf.get("realignment_filter","score_diff"),
        "realign_window_size": task_conf.get("realignment_filter","window_size"),
        "realign_max_depth": task_conf.get("realignment_filter","max_depth"),
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
        # simplerepeat filter
        "simple_repeat_db":genomon_conf.get("REFERENCE", "simple_repeat_tabix_db"),
        # EB filter
        "EBFilter": genomon_conf.get("SOFTWARE", "ebfilter"),
        "eb_map_quality": task_conf.get("eb_filter","map_quality"),
        "eb_base_quality": task_conf.get("eb_filter","base_quality"),
        "control_bam_list": input_file[2],
        # original_annotations
        "mutanno": genomon_conf.get("SOFTWARE", "mutanno"),
        "active_inhouse_normal_flag": active_inhouse_normal_flag,
        "inhouse_normal_database":inhouse_normal_tabix_db,
        "active_inhouse_tumor_flag": active_inhouse_tumor_flag,
        "inhouse_tumor_database":inhouse_tumor_tabix_db,
        "active_HGVD_flag": task_conf.get("annotation", "active_HGVD_flag"),
        "HGVD_database":genomon_conf.get("REFERENCE", "HGVD_tabix_db"),
        "active_HGMD_flag": active_HGMD_flag,
        "HGMD_database": HGMD_tabix_db,
        # annovar
        "active_annovar_flag": task_conf.get("annotation", "active_annovar_flag"),
        "annovar": genomon_conf.get("SOFTWARE", "annovar"),
        "table_annovar_params": task_conf.get("annotation", "table_annovar_params"),
        # commmon
        "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
        "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),   
        "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
        "ref_fa":genomon_conf.get("REFERENCE", "ref_fasta"),
        "interval_list": genomon_conf.get("REFERENCE", "interval_list"),
        "disease_bam": input_file[0],
        "control_bam": input_file[1],
        "out_prefix": output_dir + '/' + sample_name,
        "samtools": genomon_conf.get("SOFTWARE", "samtools"),
        "blat": genomon_conf.get("SOFTWARE", "blat")}

    interval_list = genomon_conf.get("REFERENCE", "interval_list")
    max_task_id = sum(1 for line in open(interval_list))

    mutation_call.task_exec(arguments, run_conf.project_root + '/log', run_conf.project_root + '/script', max_task_id)
    
    arguments = {
        "control_bam": input_file[1],
        "control_bam_list": input_file[2],
        "active_annovar_flag": task_conf.get("annotation", "active_annovar_flag"),
        "active_HGVD_flag": task_conf.get("annotation", "active_HGVD_flag"),
        "active_HGMD_flag": active_HGMD_flag,
        "active_inhouse_normal_flag": active_inhouse_normal_flag,
        "active_inhouse_tumor_flag": active_inhouse_tumor_flag,
        "filecount": max_task_id,
        "out_prefix": output_dir + '/' + sample_name}

    mutation_merge.task_exec(arguments, run_conf.project_root + '/log', run_conf.project_root + '/script')
     
    for task_id in range(1,(max_task_id + 1)):
        input_file = output_dir+'/'+sample_name+'_mutations_candidate.'+str(task_id)+'.hg19_multianno.txt'
        os.unlink(input_file)

    for task_id in range(1,(max_task_id + 1)):
        if os.path.exists(output_dir+'/'+sample_name+'.fisher_mutations.'+str(task_id)+'.txt'):
            os.unlink(output_dir+'/'+sample_name+'.fisher_mutations.'+str(task_id)+'.txt')
        if os.path.exists(output_dir+'/'+sample_name+'.realignment_mutations.'+str(task_id)+'.txt'):
            os.unlink(output_dir+'/'+sample_name+'.realignment_mutations.'+str(task_id)+'.txt')
        if os.path.exists(output_dir+'/'+sample_name+'.indel_mutations.'+str(task_id)+'.txt'):
            os.unlink(output_dir+'/'+sample_name+'.indel_mutations.'+str(task_id)+'.txt')
        if os.path.exists(output_dir+'/'+sample_name+'.breakpoint_mutations.'+str(task_id)+'.txt'):
            os.unlink(output_dir+'/'+sample_name+'.breakpoint_mutations.'+str(task_id)+'.txt')
        if os.path.exists(output_dir+'/'+sample_name+'.simplerepeat_mutations.'+str(task_id)+'.txt'):
            os.unlink(output_dir+'/'+sample_name+'.simplerepeat_mutations.'+str(task_id)+'.txt')
        if os.path.exists(output_dir+'/'+sample_name+'.ebfilter_mutations.'+str(task_id)+'.txt'):
            os.unlink(output_dir+'/'+sample_name+'.ebfilter_mutations.'+str(task_id)+'.txt')
        if os.path.exists(output_dir+'/'+sample_name+'.inhouse_normal.'+str(task_id)+'.txt'):
           os.unlink(output_dir+'/'+sample_name+'.inhouse_normal.'+str(task_id)+'.txt')
        if os.path.exists(output_dir+'/'+sample_name+'.inhouse_tumor.'+str(task_id)+'.txt'):
           os.unlink(output_dir+'/'+sample_name+'.inhouse_tumor.'+str(task_id)+'.txt')
        if os.path.exists(output_dir+'/'+sample_name+'.HGVD.'+str(task_id)+'.txt'):
           os.unlink(output_dir+'/'+sample_name+'.HGVD.'+str(task_id)+'.txt')
        if os.path.exists(output_dir+'/'+sample_name+'.HGMD.'+str(task_id)+'.txt'):
           os.unlink(output_dir+'/'+sample_name+'.HGMD.'+str(task_id)+'.txt')


# parse SV 
@follows( link_import_bam )
@follows( markdup )
@transform(parse_sv_bam_list, formatter(), "{subpath[0][2]}/sv/{subdir[0][0]}/{subdir[0][0]}.junction.clustered.bedpe.gz")
def parse_sv(input_file, output_file):

    dir_name = os.path.dirname(output_file)
    sample_name = os.path.basename(dir_name)

    if not os.path.isdir(dir_name): os.mkdir(dir_name)
    yaml_flag = False
    sv_sampleConf = {"target": {}, "matched_control": {}, "non_matched_control_panel": {}}
    for complist in sample_conf.sv_detection:
        if sample_name == complist[0]:

            # tumor:exist, matched_normal:exist, non-matched-normal:exist
            if complist[0] != None and complist[1] != None and complist[2] != None:
                sv_sampleConf["target"]["label"] = sample_name
                sv_sampleConf["target"]["path_to_bam"] = input_file
                sv_sampleConf["target"]["path_to_output_dir"] = dir_name
                sv_sampleConf["matched_control"]["use"] = True
                sv_sampleConf["matched_control"]["path_to_bam"] = run_conf.project_root + '/bam/' + complist[1] + '/' + complist[1] + '.markdup.bam'
                sv_sampleConf["non_matched_control_panel"]["use"] = True
                sv_sampleConf["non_matched_control_panel"]["matched_control_label"] = complist[1]
                sv_sampleConf["non_matched_control_panel"]["data_path"] = run_conf.project_root +"/sv/non_matched_control_panel/"+ complist[2] +".merged.junction.control.bedpe.gz"
                yaml_flag = True
                break

            # tumor:exist, matched_normal:exist, non-matched-normal:none
            elif complist[0] != None and complist[1] != None and complist[2] == None:
                sv_sampleConf["target"]["label"] = sample_name
                sv_sampleConf["target"]["path_to_bam"] = input_file
                sv_sampleConf["target"]["path_to_output_dir"] = dir_name
                sv_sampleConf["matched_control"]["use"] = True
                sv_sampleConf["matched_control"]["path_to_bam"] = run_conf.project_root + '/bam/' + complist[1] + '/' + complist[1] + '.markdup.bam'
                sv_sampleConf["non_matched_control_panel"]["use"] = False
                yaml_flag = True
                break

            # tumor:exist, matched_normal:none, non-matched-normal:exist
            elif complist[0] != None and complist[1] == None and complist[2] != None:
                sv_sampleConf["target"]["label"] = sample_name
                sv_sampleConf["target"]["path_to_bam"] = input_file
                sv_sampleConf["target"]["path_to_output_dir"] = dir_name
                sv_sampleConf["matched_control"]["use"] = False
                sv_sampleConf["matched_control"]["path_to_bam"] = None
                sv_sampleConf["non_matched_control_panel"]["use"] = True
                sv_sampleConf["non_matched_control_panel"]["matched_control_label"] = None
                sv_sampleConf["non_matched_control_panel"]["data_path"] = run_conf.project_root +"/sv/non_matched_control_panel/"+ complist[2] +".merged.junction.control.bedpe.gz"
                yaml_flag = True
                break

    if not yaml_flag:
        sv_sampleConf["target"]["label"] = sample_name
        sv_sampleConf["target"]["path_to_bam"] = input_file
        sv_sampleConf["target"]["path_to_output_dir"] = dir_name
        sv_sampleConf["matched_control"]["use"] = False
        sv_sampleConf["matched_control"]["path_to_bam"] = None
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
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH")}
    sv_parse.task_exec(arguments, run_conf.project_root + '/log', run_conf.project_root + '/script')


# merge SV
@follows( parse_sv )
@transform(merge_bedpe_list, formatter(".+/(?P<NAME>.+).control.yaml"), "{subpath[0][2]}/sv/non_matched_control_panel/{NAME[0]}.merged.junction.control.bedpe.gz")
def merge_sv(input_files,  output_file):

    arguments = {"genomon_sv": genomon_conf.get("SOFTWARE", "genomon_sv"),
                 "control_conf": input_files[0],
                 "bedpe": output_file,
                 "param_conf": task_conf.get("genomon_sv", "param_file"),
                 "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),   
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH")}
    sv_merge.task_exec(arguments, run_conf.project_root + '/log', run_conf.project_root + '/script')


# filt SV
@follows( merge_sv )
@transform(filt_bedpe_list, formatter(), "{subpath[0][2]}/sv/{subdir[0][0]}/{subdir[0][0]}.genomonSV.result.txt")
def filt_sv(input_files,  output_file):

    dir_name = os.path.dirname(output_file)
    sample_name = os.path.basename(dir_name)
    sample_yaml = run_conf.project_root + "/sv/config/" + sample_name + ".yaml"

    arguments = {"genomon_sv": genomon_conf.get("SOFTWARE", "genomon_sv"),
                 "sample_conf": sample_yaml,
                 "param_conf": task_conf.get("genomon_sv", "param_file"),
                 "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),   
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH")}
    sv_filt.task_exec(arguments, run_conf.project_root + '/log', run_conf.project_root + '/script')


# summary
@follows( link_import_bam )
@follows( markdup )
@follows( filt_sv )
@follows( identify_mutations )
@transform(summary_bamstats_list, formatter(), "{subpath[0][2]}/summary/{subdir[0][0]}/{subdir[0][0]}.bamstats")
def bam_stats(input_file, output_file):
    dir_name = os.path.dirname(output_file)
    if not os.path.exists(dir_name): os.makedirs(dir_name)
      
    arguments = {"PCAP": genomon_conf.get("SOFTWARE", "PCAP"),
                 "PERL5LIB": genomon_conf.get("ENV", "PERL5LIB"),
                 "input": input_file,
                 "output": output_file}
    
    r_bamstats.task_exec(arguments, run_conf.project_root + '/log', run_conf.project_root + '/script')


@follows( link_import_bam )
@follows( markdup )
@follows( filt_sv )
@follows( identify_mutations )
@transform(summary_coverage_list, formatter(), "{subpath[0][2]}/summary/{subdir[0][0]}/{subdir[0][0]}.coverage")
def coverage(input_file, output_file):

    dir_name = os.path.dirname(output_file)
    if not os.path.exists(dir_name): os.makedirs(dir_name)
    sample_name = os.path.basename(dir_name)
    depth_output_file = dir_name+'/'+sample_name+'.depth'

    incl_bed_file = ""
    genome_file = ""
    data_type = ""
    if task_conf.get("coverage", "wgs_flag") == "True":
        genome_file = genomon_conf.get("REFERENCE", "hg19_genome")
        incl_bed_file = output_file + "genome.bed"
        incl_bed_w = task_conf.get("coverage", "wgs_incl_bed_width")
        r_coverage.create_incl_bed_wgs(genome_file, incl_bed_file, long(incl_bed_w), "")
        data_type = "wgs"

    arguments = {"data_type": data_type,
                 "i_bed_lines": task_conf.get("coverage", "wgs_i_bed_lines"),
                 "i_bed_size": task_conf.get("coverage", "wgs_i_bed_width"),
                 "incl_bed_file": incl_bed_file,
                 "genome_file": genome_file,
                 "gaptxt": genomon_conf.get("REFERENCE", "gaptxt"),
                 "bait_file": genomon_conf.get("REFERENCE", "bait_file"),
                 "BEDTOOLS": genomon_conf.get("SOFTWARE", "bedtools"),
                 "SAMTOOLS": genomon_conf.get("SOFTWARE", "samtools"),
                 "LD_LIBRARY_PATH": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "input": input_file,
                 "output": depth_output_file}

    r_coverage.task_exec(arguments, run_conf.project_root + '/log', run_conf.project_root + '/script')
    
    r_coverage.calc_coverage(depth_output_file, task_conf.get("coverage", "coverage"), output_file)
    
    os.unlink(dir_name+'/'+sample_name+'.depth')
    os.unlink(dir_name+'/'+sample_name+'.depth.input_bed')


@follows( bam_stats )
@follows( coverage )
@transform(summary_merge_list, formatter(), "{subpath[0][2]}/summary/{subdir[0][0]}/{subdir[0][0]}.tsv")
def merge(input_files, output_file):

    for f in input_files:
        if not os.path.exists(f):
            raise

    excel_file = os.path.splitext(output_file)[0] + ".xls"
    r_merge.mkxls(input_files, excel_file)
    r_merge.Excel2TSV(excel_file, output_file)


#####################
# post analysis stage
@active_if(task_conf.getboolean("post_analysis", "enable"))
@follows(identify_mutations)
@transform(pa_capt_list_mutation, formatter(), 
           [run_conf.project_root + "/post_analysis/mutation/bam/{subdir[0][0]}/{subdir[0][0]}.markdup.pickup.bam",
            run_conf.project_root + "/post_analysis/mutation/capture_script/{subdir[0][0]}.bat"])
def post_analysis_capt_mutation(input_file, output_file):

    ID = os.path.basename(input_file).replace("_genomon_mutations.result.txt", "")

    config = r_pa_capture.script_template_config.format(
                 mode = "mutation",
                 capture_max = task_conf.get("capture_igv", "capture_max"),
                 capture_width = task_conf.get("capture_igv", "capture_width"),
                 pickup_width = task_conf.get("capture_bam", "pickup_width"),
                 markdup_bam_suffix = ".markdup.bam",
                 pickup_bam_suffix = ".markdup.pickup.bam",
                 sept = "\\t",
                 header = "True",
                 col_pos_chr1 = "0",
                 col_pos_start = "1",
                 col_pos_chr2 = "0",
                 col_pos_end = "2",
                 samtools = genomon_conf.get("SOFTWARE", "samtools"),
                 bedtools = genomon_conf.get("SOFTWARE", "bedtools"),
                 )

    bam_script = ""
    if os.path.exists("%s/post_analysis/mutation/bam/%s/%s.markdup.pickup.bam" % (run_conf.project_root, ID, ID)) == False:
        bam_script = r_pa_capture.script_template_bam.format(
                 mode = "mutation",
                 input_file = input_file,
                 id = ID,
                 output_file = run_conf.project_root + '/post_analysis/mutation/bam_script/pickup.%s.sh' % (ID),
                 output_bam_dir = run_conf.project_root + '/post_analysis/mutation/bam/',
                 output_log_dir = run_conf.project_root + '/log/',
                 genomon_root = run_conf.project_root,
                 config = config
                 )
                 
    igv_script = ""
    if os.path.exists(run_conf.project_root + "/post_analysis/mutation/capture_script/capture.bat") == False:        
        igv_script = r_pa_capture.script_template_igv.format(
                 mode = "mutation",
                 input_file = input_file,
                 id = ID,
                 output_file = run_conf.project_root + '/post_analysis/mutation/capture_script/%s.bat' % (ID),
                 output_dir = run_conf.project_root + '/post_analysis/mutation/capture',
                 genomon_root = run_conf.project_root,
                 config = config
                 )
                 
    arguments = {"mode": "mutation",
                 "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),
                 "bam_script": bam_script,
                 "igv_script": igv_script
                 }
                 
    r_pa_capture.task_exec(arguments, run_conf.project_root + '/log', run_conf.project_root + '/script')

@active_if(task_conf.getboolean("post_analysis", "enable"))
@follows(post_analysis_capt_mutation)
@merge(pa_merge_list_mutation, 
       [run_conf.project_root + "/post_analysis/mutation/capture_script/capture.bat",
        run_conf.project_root + "/post_analysis/merge.mutation.csv"])
def post_analysis_merge_mutation(input_files, output_file):

    config = r_pa_merge.script_template_config.format(
                 mode = "mutation",
                 sept = "\\t",
                 header = "True",
                 suffix = "_genomon_mutations.result.txt",
                 filters = task_conf.get("merge_format_mutation", "filters"),
                 )

    merge_result = ""
    if os.path.exists(run_conf.project_root + "/post_analysis/merge.mutation.csv") == False:
        li = ""
        for item in input_files:
            if len(item) > 0:
                if len(li) > 0:
                    li += ","
                li += item
                
        merge_result = r_pa_merge.script_template_result.format(
                 mode = "mutation",
                 input_files = li,
                 output_file = run_conf.project_root + "/post_analysis/merge.mutation.csv",
                 config = config
                 )
                 
    merge_igv = ""
    if os.path.exists(run_conf.project_root + "/post_analysis/mutation/capture_script/capture.bat") == False:
        li = ""
        for item in input_files:
            if len(item) > 0:
                if len(li) > 0:
                    li += ","

                ID = os.path.basename(item).replace("_genomon_mutations.result.txt", "")
                li += run_conf.project_root + "/post_analysis/mutation/capture_script/" + ID + ".bat"
       
        merge_igv = r_pa_merge.script_template_igv.format(
                 input_files = li,
                 output_file = run_conf.project_root + "/post_analysis/mutation/capture_script/capture.bat",
                 config = config
                 )

    arguments = {"pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),
                 "merge_result": merge_result,
                 "merge_igv": merge_igv
                }
                 
    r_pa_merge.task_exec(arguments, run_conf.project_root + '/log', run_conf.project_root + '/script')

@active_if(task_conf.getboolean("post_analysis", "enable"))
@follows(filt_sv)
@transform(pa_capt_list_sv, formatter(), 
           [run_conf.project_root + "/post_analysis/sv/bam/{subdir[0][0]}/{subdir[0][0]}.markdup.pickup.bam",
            run_conf.project_root + "/post_analysis/sv/capture_script/{subdir[0][0]}.bat"])
def post_analysis_capt_sv(input_file,  output_file):

    ID = os.path.basename(input_file).replace(".genomonSV.result.txt", "")
    
    config = r_pa_capture.script_template_config.format(
                 mode = "sv",
                 capture_max = task_conf.get("capture_igv", "capture_max"),
                 capture_width = task_conf.get("capture_igv", "capture_width"),
                 pickup_width = task_conf.get("capture_bam", "pickup_width"),
                 markdup_bam_suffix = ".markdup.bam",
                 pickup_bam_suffix = ".markdup.pickup.bam",
                 sept = "\\t",
                 header = "False",
                 col_pos_chr1 = "0",
                 col_pos_start = "1",
                 col_pos_chr2 = "3",
                 col_pos_end = "4",
                 samtools = genomon_conf.get("SOFTWARE", "samtools"),
                 bedtools = genomon_conf.get("SOFTWARE", "bedtools"),
                 )

    bam_script = ""
    if os.path.exists("%s/post_analysis/sv/bam/%s/%s.markdup.pickup.bam" % (run_conf.project_root, ID, ID)) == False:
        bam_script = r_pa_capture.script_template_bam.format(
                 mode = "sv",
                 input_file = input_file,
                 id = ID,
                 output_file = run_conf.project_root + '/post_analysis/sv/bam_script/pickup.%s.sh' % (ID),
                 output_bam_dir = run_conf.project_root + '/post_analysis/sv/bam/',
                 output_log_dir = run_conf.project_root + '/log/',
                 genomon_root = run_conf.project_root,
                 config = config
                 )
                 
    igv_script = ""
    if os.path.exists(run_conf.project_root + "/post_analysis/sv/capture_script/capture.bat") == False:        
        igv_script = r_pa_capture.script_template_igv.format(
                 mode = "sv",
                 input_file = input_file,
                 id = ID,
                 output_file = run_conf.project_root + '/post_analysis/sv/capture_script/%s.bat' % (ID),
                 output_dir = run_conf.project_root + '/post_analysis/sv/capture/',
                 genomon_root = run_conf.project_root,
                 config = config
                 )
                 
    arguments = {"pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),
                 "bam_script": bam_script,
                 "igv_script": igv_script,
                 }
                 
    r_pa_capture.task_exec(arguments, run_conf.project_root + '/log', run_conf.project_root + '/script')

@active_if(task_conf.getboolean("post_analysis", "enable"))
@follows(post_analysis_capt_sv)
@merge(pa_merge_list_sv, 
       [run_conf.project_root + "/post_analysis/sv/capture_script/capture.bat", 
        run_conf.project_root + "/post_analysis/merge.sv.csv"])
def post_analysis_merge_sv(input_files, output_file):

    config = r_pa_merge.script_template_config.format(
                 mode = "sv",
                 sept = "\\t",
                 header = "False",
                 suffix = ".genomonSV.result.txt",
                 filters = task_conf.get("merge_format_sv", "filters"),
                 )

    merge_result = ""
    if os.path.exists(run_conf.project_root + "/post_analysis/merge.sv.csv") == False:
        li = ""
        for item in input_files:
            if len(item) > 0:
                if len(li) > 0:
                    li += ","
                li += item
                
        merge_result = r_pa_merge.script_template_result.format(
                 mode = "sv",
                 input_files = li,
                 output_file = run_conf.project_root + "/post_analysis/merge.sv.csv",
                 config = config
                 )
                 
    merge_igv = ""
    if os.path.exists(run_conf.project_root + "/post_analysis/sv/capture_script/capture.bat") == False:
        li = ""
        for item in input_files:
            if len(item) > 0:
                if len(li) > 0:
                    li += ","

                ID = os.path.basename(item).replace(".genomonSV.result.txt", "")
                li += run_conf.project_root + "/post_analysis/sv/capture_script/" + ID + ".bat"
                
        merge_igv = r_pa_merge.script_template_igv.format(
                 input_files = li,
                 output_file = run_conf.project_root + "/post_analysis/sv/capture_script/capture.bat",
                 config = config
                 )

    arguments = {"pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),
                 "merge_result": merge_result,
                 "merge_igv": merge_igv
                }
                 
    r_pa_merge.task_exec(arguments, run_conf.project_root + '/log', run_conf.project_root + '/script')

@active_if(task_conf.getboolean("post_analysis", "enable"))
@follows(write_summary)
@merge(pa_summary_list, run_conf.project_root + "/post_analysis/merge.summary.csv")
def post_analysis_merge_summary(input_files, output_file):

    li = ""
    for item in input_files:
        if len(item) > 0:
            if len(li) > 0:
                li += ","
            li += item

    arguments = {"mode":"summary",
                 "input_files":li,
                 "output_file":output_file,
                 "sept":"\\t",
                 "header":"True",
                 "suffix":".tsv",
                 "filters":task_conf.get("merge_format_sumary", "filters"),
                 "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH")}

    r_pa_merge.task_exec(arguments, run_conf.project_root + '/log', run_conf.project_root + '/script')
    
