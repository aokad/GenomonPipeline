import os
import shutil
import glob
from ruffus import *
from genomon_pipeline.config.run_conf import *
from genomon_pipeline.config.genomon_conf import *
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
from genomon_pipeline.dna_resource.qc_bamstats import *
from genomon_pipeline.dna_resource.qc_coverage import *
from genomon_pipeline.dna_resource.qc_merge import *
from genomon_pipeline.dna_resource.post_analysis import *
from genomon_pipeline.dna_resource.paplot import *

# set task classes
bamtofastq = Bam2Fastq(genomon_conf.get("bam2fastq", "qsub_option"), run_conf.drmaa)
fastq_splitter = Fastq_splitter(genomon_conf.get("split_fastq", "qsub_option"), run_conf.drmaa)
bwa_align = Bwa_align(genomon_conf.get("bwa_mem", "qsub_option"), run_conf.drmaa)
markduplicates = Markduplicates(genomon_conf.get("markduplicates", "qsub_option"), run_conf.drmaa)
mutation_call = Mutation_call(genomon_conf.get("mutation_call", "qsub_option"), run_conf.drmaa)
mutation_merge = Mutation_merge(genomon_conf.get("mutation_merge", "qsub_option"), run_conf.drmaa)
sv_parse = SV_parse(genomon_conf.get("sv_parse", "qsub_option"), run_conf.drmaa)
sv_merge = SV_merge(genomon_conf.get("sv_merge", "qsub_option"), run_conf.drmaa)
sv_filt = SV_filt(genomon_conf.get("sv_filt", "qsub_option"), run_conf.drmaa)
r_qc_bamstats = Res_QC_Bamstats(genomon_conf.get("qc_bamstats", "qsub_option"), run_conf.drmaa)
r_qc_coverage = Res_QC_Coverage(genomon_conf.get("qc_coverage", "qsub_option"), run_conf.drmaa)
r_qc_merge = Res_QC_Merge(genomon_conf.get("qc_merge", "qsub_option"), run_conf.drmaa)
r_pa_plot = Res_PA_Plot(genomon_conf.get("pa_plot", "qsub_option"), run_conf.drmaa)
r_post_analysis = Res_PostAnalysis(genomon_conf.get("post_analysis", "qsub_option"), run_conf.drmaa)

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
     if os.path.exists(run_conf.project_root + '/mutation/' + complist[0] + '/' + complist[0] + '.genomon_mutation.result.filt.txt'): continue
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
    tmp_list.append(run_conf.project_root + '/sv/control_panel/' + control_panel_name + ".control_info.txt")
    for sample in sample_conf.control_panel[control_panel_name]:
        tmp_list.append(run_conf.project_root+ "/sv/"+ sample +"/"+ sample +".junction.clustered.bedpe.gz")
    merge_bedpe_list.append(tmp_list)

# generate input list of 'SV filt'
filt_bedpe_list = []
for complist in sample_conf.sv_detection:
    if os.path.exists(run_conf.project_root + '/sv/' + complist[0] +'/'+ complist[0] +'.genomonSV.result.filt.txt'): continue
    filt_bedpe_list.append(run_conf.project_root+ "/sv/"+ complist[0] +"/"+ complist[0] +".junction.clustered.bedpe.gz")

# generate input list of 'qc'
qc_bamstats_list = []
qc_coverage_list = []
qc_merge_list = []
for sample in sample_conf.qc:
    if os.path.exists(run_conf.project_root + '/qc/' + sample + '/' + sample + '.genomonQC.result.txt'): continue
    qc_merge_list.append(
        [run_conf.project_root + '/qc/' + sample + '/' + sample + '.bamstats',
         run_conf.project_root + '/qc/' + sample + '/' + sample + '.coverage'])
    if not os.path.exists(run_conf.project_root + '/qc/' + sample + '/' + sample + '.bamstats'):
        qc_bamstats_list.append(run_conf.project_root + '/bam/' + sample +'/'+ sample +'.markdup.bam')
    if not os.path.exists(run_conf.project_root + '/qc/' + sample + '/' + sample + '.coverage'):
        qc_coverage_list.append(run_conf.project_root + '/bam/' + sample +'/'+ sample +'.markdup.bam')


sample_conf_name, ext = os.path.splitext(os.path.basename(run_conf.sample_conf_file))

# generate input list of 'post analysis for mutation'
pa_output_files_mutation = r_post_analysis.output_files("mutation", sample_conf.mutation_call, run_conf.project_root + '/post_analysis/' + sample_conf_name, genomon_conf)
run_pa = False
for key in pa_output_files_mutation:
    if not os.path.exists(pa_output_files_mutation[key]): run_pa = True
pa_files_mutation = []
pa_samples_mutation = {}
if run_pa == True:
    pa_samples_mutation = r_post_analysis.sample_split_case(sample_conf.mutation_call)
    for complist in sample_conf.sv_detection:
        pa_files_mutation.append(run_conf.project_root + '/mutation/' + complist[0] +'/'+ complist[0] +'.genomon_mutation.result.filt.txt')
        
# generate input list of 'post analysis for SV'
pa_output_files_sv = r_post_analysis.output_files("sv", sample_conf.sv_detection, run_conf.project_root + '/post_analysis/' + sample_conf_name, genomon_conf)
run_pa = False
for key in pa_output_files_sv:
    if not os.path.exists(pa_output_files_sv[key]): run_pa = True
pa_files_sv = []
pa_samples_sv = {}
if run_pa == True:
    pa_samples_sv = r_post_analysis.sample_split_case(sample_conf.sv_detection)
    for complist in sample_conf.sv_detection:
        pa_files_sv.append(run_conf.project_root + '/sv/' + complist[0] +'/'+ complist[0] +'.genomonSV.result.filt.txt')
        
# generate input list of 'post analysis for qc'
pa_output_files_qc = r_post_analysis.output_files("qc", sample_conf.qc, run_conf.project_root + '/post_analysis/' + sample_conf_name, genomon_conf)
run_pa = False
for key in pa_output_files_qc:
    if not os.path.exists(pa_output_files_qc[key]): run_pa = True
pa_files_qc = []
pa_samples_qc = []
if run_pa == True:
    for sample in sample_conf.qc:
        pa_samples_qc.append(sample)
        pa_files_qc.append(run_conf.project_root + '/qc/' + sample + '/' + sample + '.genomonQC.result.txt')

# generate input list of paplot
paplot_files_qc = []
paplot_files_sv = []
paplot_files_mutation = []
paplot_files_collate = []
if not os.path.exists(run_conf.project_root + '/paplot/' + sample_conf_name + '/index.html'):
    # mutation
    if pa_output_files_mutation.has_key("filt_case1"):
        paplot_files_mutation.append(pa_output_files_mutation["filt_case1"])
    if pa_output_files_mutation.has_key("filt_case2") and genomon_conf.getboolean("pa_plot", "include_unpanel"):
        paplot_files_mutation.append(pa_output_files_mutation["filt_case2"])
    if pa_output_files_mutation.has_key("filt_case3") and genomon_conf.getboolean("pa_plot", "include_unpair"):
        paplot_files_mutation.append(pa_output_files_mutation["filt_case3"])
    if pa_output_files_mutation.has_key("filt_case4") and genomon_conf.getboolean("pa_plot", "include_unpanel") and genomon_conf.getboolean("pa_plot", "include_unpair"):
        paplot_files_mutation.append(pa_output_files_mutation["filt_case4"])
            
    # sv
    if pa_output_files_sv.has_key("filt_case1"):
        paplot_files_sv.append(pa_output_files_sv["filt_case1"])
    if pa_output_files_sv.has_key("filt_case2") and genomon_conf.getboolean("pa_plot", "include_unpanel"):
        paplot_files_sv.append(pa_output_files_sv["filt_case2"])
    if pa_output_files_sv.has_key("filt_case3") and genomon_conf.getboolean("pa_plot", "include_unpair"):
        paplot_files_sv.append(pa_output_files_sv["filt_case3"])
    if pa_output_files_sv.has_key("filt_case4") and genomon_conf.getboolean("pa_plot", "include_unpanel") and genomon_conf.getboolean("pa_plot", "include_unpair"):
        paplot_files_sv.append(pa_output_files_sv["filt_case4"])

    # qc
    if pa_output_files_qc.has_key("all"):
        paplot_files_qc.append(pa_output_files_qc["all"])
        
    paplot_files_collate.extend(paplot_files_qc)
    paplot_files_collate.extend(paplot_files_sv)
    paplot_files_collate.extend(paplot_files_mutation)
    print paplot_files_collate
# prepare output directories
if not os.path.isdir(run_conf.project_root): os.mkdir(run_conf.project_root)
if not os.path.isdir(run_conf.project_root + '/script'): os.mkdir(run_conf.project_root + '/script')
if not os.path.isdir(run_conf.project_root + '/script/sv_merge'): os.mkdir(run_conf.project_root + '/script/sv_merge')
if not os.path.isdir(run_conf.project_root + '/script/post_analysis'): os.mkdir(run_conf.project_root + '/script/post_analysis')
if not os.path.isdir(run_conf.project_root + '/script/paplot'): os.mkdir(run_conf.project_root + '/script/paplot')
if not os.path.isdir(run_conf.project_root + '/log'): os.mkdir(run_conf.project_root + '/log')
if not os.path.isdir(run_conf.project_root + '/log/sv_merge'): os.mkdir(run_conf.project_root + '/log/sv_merge')
if not os.path.isdir(run_conf.project_root + '/log/post_analysis'): os.mkdir(run_conf.project_root + '/log/post_analysis')
if not os.path.isdir(run_conf.project_root + '/log/paplot'): os.mkdir(run_conf.project_root + '/log/paplot')
if not os.path.isdir(run_conf.project_root + '/fastq'): os.mkdir(run_conf.project_root + '/fastq')
if not os.path.isdir(run_conf.project_root + '/bam'): os.mkdir(run_conf.project_root + '/bam')
if not os.path.isdir(run_conf.project_root + '/mutation'): os.mkdir(run_conf.project_root + '/mutation')
if not os.path.isdir(run_conf.project_root + '/mutation/control_panel'): os.mkdir(run_conf.project_root + '/mutation/control_panel')
if not os.path.isdir(run_conf.project_root + '/sv'): os.mkdir(run_conf.project_root + '/sv')
if not os.path.isdir(run_conf.project_root + '/sv/non_matched_control_panel'): os.mkdir(run_conf.project_root + '/sv/non_matched_control_panel')
if not os.path.isdir(run_conf.project_root + '/sv/control_panel'): os.mkdir(run_conf.project_root + '/sv/control_panel')
if not os.path.isdir(run_conf.project_root + '/qc'): os.mkdir(run_conf.project_root + '/qc')
for sample in sample_conf.qc:
    if not os.path.isdir(run_conf.project_root + '/qc/' + sample): os.mkdir(run_conf.project_root + '/qc/' + sample)
if (genomon_conf.getboolean("post_analysis", "enable") == True):
    if not os.path.exists(run_conf.project_root + '/post_analysis'): os.mkdir(run_conf.project_root + '/post_analysis')
    if not os.path.exists(run_conf.project_root + '/post_analysis/' + sample_conf_name): os.mkdir(run_conf.project_root + '/post_analysis/' + sample_conf_name)
if not os.path.isdir(run_conf.project_root + '/config'): os.mkdir(run_conf.project_root + '/config')

genomon_conf_name, ext = os.path.splitext(os.path.basename(run_conf.genomon_conf_file))
shutil.copyfile(run_conf.genomon_conf_file, run_conf.project_root + '/config/' + genomon_conf_name +'_'+ run_conf.analysis_timestamp + ext)

for outputfiles in (bam2fastq_output_list, linked_fastq_list):
    for outputfile in outputfiles:
        sample = os.path.basename(os.path.dirname(outputfile[0][0]))
        fastq_dir = run_conf.project_root + '/fastq/' + sample
        bam_dir = run_conf.project_root + '/bam/' + sample
        if not os.path.isdir(fastq_dir): os.mkdir(fastq_dir)
        if not os.path.isdir(bam_dir): os.mkdir(bam_dir)

for target_sample_dict in (sample_conf.bam_import, sample_conf.fastq, sample_conf.bam_tofastq):
    for sample in target_sample_dict:
        script_dir = run_conf.project_root + '/script/' + sample
        log_dir = run_conf.project_root + '/log/' + sample
        if not os.path.isdir(script_dir): os.mkdir(script_dir)
        if not os.path.isdir(log_dir): os.mkdir(log_dir)

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
        control_conf = run_conf.project_root + '/sv/control_panel/' + control_panel_name + ".control_info.txt"
        with open(control_conf,  "w") as out_handle:
            for sample in sample_conf.control_panel[control_panel_name]:
                out_handle.write(sample+ "\t"+ run_conf.project_root+ "/sv/"+ sample +"/"+ sample+ "\n")

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
    bamtofastq.task_exec(arguments, run_conf.project_root + '/log/' + sample, run_conf.project_root + '/script/'+ sample)


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

    sample_name = os.path.basename(target_dir)

    for oo in output_files:
        os.unlink(oo)

    input_prefix, ext = os.path.splitext(input_files[0][0])
    arguments = {"lines": genomon_conf.get("split_fastq", "split_fastq_line_number"),
                 "fastq_filter": genomon_conf.get("split_fastq", "fastq_filter"),
                 "target_dir": target_dir,
                 "ext": ext}
    
    fastq_splitter.task_exec(arguments, run_conf.project_root + '/log/' + sample_name, run_conf.project_root + '/script/'+ sample_name, 2)
   
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
    split_lines = genomon_conf.get("split_fastq", "split_fastq_line_number")

    ans_quotient = all_line_num / int(split_lines)
    ans_remainder = all_line_num % int(split_lines)
    max_task_id = ans_quotient if ans_remainder == 0 else ans_quotient + 1
    
    arguments = {"input_dir": input_dir,
                 "output_dir": output_dir,
                 "sample_name": sample_name,
                 "bwa": genomon_conf.get("SOFTWARE", "bwa"),
                 "bwa_params": genomon_conf.get("bwa_mem", "bwa_params"),
                 "ref_fa":genomon_conf.get("REFERENCE", "ref_fasta"),
                 "biobambam": genomon_conf.get("SOFTWARE", "biobambam")}

    bwa_align.task_exec(arguments, run_conf.project_root + '/log/' + sample_name , run_conf.project_root + '/script/' + sample_name, max_task_id) 

    for task_id in range(max_task_id):
        num = str(task_id).zfill(4)
        os.unlink(input_dir +'/1_'+str(num)+'.fastq_split')
        os.unlink(input_dir +'/2_'+str(num)+'.fastq_split')
        os.unlink(output_dir+'/'+sample_name+'_'+str(num)+'.bwa.sam')


# merge sorted bams into one and mark duplicate reads with biobambam
@collate(map_dna_sequence, formatter(), "{subpath[0][2]}/bam/{subdir[0][0]}/{subdir[0][0]}.markdup.bam", "{subpath[0][2]}/bam/{subdir[0][0]}")
def markdup(input_files, output_file, output_dir):

    sample_name = os.path.basename(output_dir)

    output_prefix, ext = os.path.splitext(output_file)

    input_bam_files = ""
    for input_file in input_files:
        input_bam_files = input_bam_files + " I=" + input_file

    arguments = {"biobambam": genomon_conf.get("SOFTWARE", "biobambam"),
                 "out_prefix": output_prefix,
                 "input_bam_files": input_bam_files,
                 "out_bam": output_file}

    markduplicates.task_exec(arguments, run_conf.project_root + '/log/' + sample_name , run_conf.project_root + '/script/'+ sample_name)

    for input_file in input_files:
        os.unlink(input_file)
        os.unlink(input_file + ".bai")


# identify mutations
@follows( markdup )
@follows( link_import_bam )
@subdivide(markdup_bam_list, formatter(), "{subpath[0][2]}/mutation/{subdir[0][0]}/{subdir[0][0]}.genomon_mutation.result.filt.txt", "{subpath[0][2]}/mutation/{subdir[0][0]}")
def identify_mutations(input_file, output_file, output_dir):

    sample_name = os.path.basename(output_dir)

    active_inhouse_normal_flag = False
    if genomon_conf.has_option("annotation", "active_inhouse_normal_flag"):
        active_inhouse_normal_flag = genomon_conf.get("annotation", "active_inhouse_normal_flag")

    inhouse_normal_tabix_db = ""
    if genomon_conf.has_option("REFERENCE", "inhouse_normal_tabix_db"):
        inhouse_normal_tabix_db = genomon_conf.get("REFERENCE", "inhouse_normal_tabix_db")

    active_inhouse_tumor_flag = False
    if genomon_conf.has_option("annotation", "active_inhouse_tumor_flag"):
        active_inhouse_tumor_flag = genomon_conf.get("annotation", "active_inhouse_tumor_flag")

    inhouse_tumor_tabix_db = ""
    if genomon_conf.has_option("REFERENCE", "inhouse_tumor_tabix_db"):
        inhouse_tumor_tabix_db = genomon_conf.get("REFERENCE", "inhouse_tumor_tabix_db")

    active_HGMD_flag = False
    if genomon_conf.has_option("annotation", "active_HGMD_flag"):
        active_HGMD_flag = genomon_conf.get("annotation", "active_HGMD_flag")
        
    HGMD_tabix_db = ""
    if genomon_conf.has_option("REFERENCE", "HGMD_tabix_db"):
        HGMD_tabix_db = genomon_conf.get("REFERENCE", "HGMD_tabix_db")

    arguments = {
        # fisher mutation
        "fisher": genomon_conf.get("SOFTWARE", "fisher"),
        "map_quality": genomon_conf.get("fisher_mutation_call", "map_quality"),
        "base_quality": genomon_conf.get("fisher_mutation_call", "base_quality"),
        "min_allele_freq": genomon_conf.get("fisher_mutation_call", "disease_min_allele_frequency"),
        "max_allele_freq": genomon_conf.get("fisher_mutation_call", "control_max_allele_frequency"),
        "min_depth": genomon_conf.get("fisher_mutation_call", "min_depth"),
        "min_variant_read": genomon_conf.get("fisher_mutation_call", "min_variant_read"),
        "fisher_thres": genomon_conf.get("fisher_mutation_call", "fisher_thres_hold"),
        "post_10_q": genomon_conf.get("fisher_mutation_call", "post_10_q"),
        # realignment filter
        "mutfilter": genomon_conf.get("SOFTWARE", "mutfilter"),
        "realign_score_diff": genomon_conf.get("realignment_filter","score_diff"),
        "realign_window_size": genomon_conf.get("realignment_filter","window_size"),
        "realign_max_depth": genomon_conf.get("realignment_filter","max_depth"),
        # indel filter
        "indel_search_length": genomon_conf.get("indel_filter","search_length"),
        "indel_neighbor": genomon_conf.get("indel_filter","neighbor"),
        "indel_base_quality": genomon_conf.get("indel_filter","base_quality"),
        "indel_min_depth": genomon_conf.get("indel_filter","min_depth"),
        "indel_min_mismatch": genomon_conf.get("indel_filter","max_mismatch"),
        "indel_min_allele_freq": genomon_conf.get("indel_filter","max_allele_freq"),
        # breakpoint filter
        "bp_max_depth": genomon_conf.get("breakpoint_filter","max_depth"),
        "bp_min_clip_size": genomon_conf.get("breakpoint_filter","min_clip_size"),
        "bp_junc_num_thres": genomon_conf.get("breakpoint_filter","junc_num_thres"),
        "bp_map_quality": genomon_conf.get("breakpoint_filter","map_quality"),
        # simplerepeat filter
        "simple_repeat_db":genomon_conf.get("REFERENCE", "simple_repeat_tabix_db"),
        # EB filter
        "EBFilter": genomon_conf.get("SOFTWARE", "ebfilter"),
        "eb_map_quality": genomon_conf.get("eb_filter","map_quality"),
        "eb_base_quality": genomon_conf.get("eb_filter","base_quality"),
        "control_bam_list": input_file[2],
        # original_annotations
        "mutanno": genomon_conf.get("SOFTWARE", "mutanno"),
        "active_inhouse_normal_flag": active_inhouse_normal_flag,
        "inhouse_normal_database":inhouse_normal_tabix_db,
        "active_inhouse_tumor_flag": active_inhouse_tumor_flag,
        "inhouse_tumor_database":inhouse_tumor_tabix_db,
        "active_HGVD_flag": genomon_conf.get("annotation", "active_HGVD_flag"),
        "HGVD_database":genomon_conf.get("REFERENCE", "HGVD_tabix_db"),
        "active_HGMD_flag": active_HGMD_flag,
        "HGMD_database": HGMD_tabix_db,
        # annovar
        "active_annovar_flag": genomon_conf.get("annotation", "active_annovar_flag"),
        "annovar": genomon_conf.get("SOFTWARE", "annovar"),
        "table_annovar_params": genomon_conf.get("annotation", "table_annovar_params"),
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

    mutation_call.task_exec(arguments, run_conf.project_root + '/log/' + sample_name, run_conf.project_root + '/script/' + sample_name, max_task_id)
    
    arguments = {
        "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
        "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),   
        "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
        "control_bam": input_file[1],
        "control_bam_list": input_file[2],
        "active_annovar_flag": genomon_conf.get("annotation", "active_annovar_flag"),
        "active_HGVD_flag": genomon_conf.get("annotation", "active_HGVD_flag"),
        "active_HGMD_flag": active_HGMD_flag,
        "active_inhouse_normal_flag": active_inhouse_normal_flag,
        "active_inhouse_tumor_flag": active_inhouse_tumor_flag,
        "filecount": max_task_id,
        "mutil": genomon_conf.get("SOFTWARE", "mutil"),
        "eb_pval": genomon_conf.get("eb_filter","ebcall_pval-log10_thres"),
        "fish_pval": genomon_conf.get("fisher_mutation_call","fisher_pval-log10_thres"),
        "realign_pval": genomon_conf.get("realignment_filter","fisher_pval-log10_thres"),
        "tcount": genomon_conf.get("realignment_filter","disease_min_mismatch"),
        "ncount": genomon_conf.get("realignment_filter","control_max_mismatch"),
        "post10q": genomon_conf.get("fisher_mutation_call","post_10_q_thres"),
        "r_post10q": genomon_conf.get("realignment_filter","post_10_q_thres"),
        "meta_info_em": get_meta_info(["fisher", "mutfilter", "ebfilter", "mutil", "mutanno"]),
        "meta_info_e":  get_meta_info(["fisher", "mutfilter", "ebfilter", "mutil"]),
        "meta_info_m": get_meta_info(["fisher", "mutfilter", "mutil", "mutanno"]),
        "meta_info":   get_meta_info(["fisher", "mutfilter", "mutil"]),
        "out_prefix": output_dir + '/' + sample_name}

    mutation_merge.task_exec(arguments, run_conf.project_root + '/log/' + sample_name, run_conf.project_root + '/script/' + sample_name)

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
    if not os.path.isdir(dir_name): os.mkdir(dir_name)
    sample_name = os.path.basename(dir_name)

    arguments = {"genomon_sv": genomon_conf.get("SOFTWARE", "genomon_sv"),
                 "input_bam": input_file,
                 "output_prefix": output_file.replace(".junction.clustered.bedpe.gz", ""),
                 "param": genomon_conf.get("sv_parse", "params"),
                 "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),   
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "htslib": genomon_conf.get("SOFTWARE", "htslib")}

    sv_parse.task_exec(arguments, run_conf.project_root + '/log/' + sample_name , run_conf.project_root + '/script/' + sample_name)


# merge SV
@follows( parse_sv )
@transform(merge_bedpe_list, formatter(".+/(?P<NAME>.+).control_info.txt"), "{subpath[0][2]}/sv/non_matched_control_panel/{NAME[0]}.merged.junction.control.bedpe.gz")
def merge_sv(input_files,  output_file):

    arguments = {"genomon_sv": genomon_conf.get("SOFTWARE", "genomon_sv"),
                 "control_info": input_files[0],
                 "merge_output_file": output_file,
                 "param": genomon_conf.get("sv_merge", "params"),
                 "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),   
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "htslib": genomon_conf.get("SOFTWARE", "htslib")}

    sv_merge.task_exec(arguments, run_conf.project_root + '/log/sv_merge', run_conf.project_root + '/script/sv_merge')


# filt SV
@follows( merge_sv )
@transform(filt_bedpe_list, formatter(), "{subpath[0][2]}/sv/{subdir[0][0]}/{subdir[0][0]}.genomonSV.result.filt.txt")
def filt_sv(input_files,  output_file):

    dir_name = os.path.dirname(output_file)
    sample_name = os.path.basename(dir_name)
    #sample_yaml = run_conf.project_root + "/sv/config/" + sample_name + ".yaml"

    filt_param = ""

    for complist in sample_conf.sv_detection:
        if sample_name == complist[0]:

            if complist[1] != None:
                filt_param = filt_param + " --matched_control_bam " + run_conf.project_root + "/bam/" + complist[1] + '/' + complist[1] + ".markdup.bam"

            if complist[2] != None:
                filt_param = filt_param + " --non_matched_control_junction " + run_conf.project_root +"/sv/non_matched_control_panel/"+ complist[2] +".merged.junction.control.bedpe.gz"
                if complist[1] != None:
                    filt_param = filt_param + " --matched_control_label " + complist[1]

            break

    filt_param = filt_param.lstrip(' ') + ' ' + genomon_conf.get("sv_filt", "params")

    arguments = {"genomon_sv": genomon_conf.get("SOFTWARE", "genomon_sv"),
                 "input_bam": run_conf.project_root + "/bam/" + sample_name + '/' + sample_name + ".markdup.bam",
                 "output_prefix": run_conf.project_root + "/sv/" + sample_name + '/' + sample_name,
                 "reference_genome": genomon_conf.get("REFERENCE", "ref_fasta"),
                 "annotation_dir": genomon_conf.get("sv_filt", "annotation_dir"),
                 "param": filt_param,
                 "meta_info": get_meta_info(["genomon_sv", "sv_utils"]),
                 "sv_utils": genomon_conf.get("SOFTWARE", "sv_utils"),
                 "sv_utils_annotation_dir": genomon_conf.get("sv_filt", "sv_utils_annotation_dir"),
                 "sv_utils_param": genomon_conf.get("sv_filt", "sv_utils_params"),
                 "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),   
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "htslib": genomon_conf.get("SOFTWARE", "htslib")}

    sv_filt.task_exec(arguments, run_conf.project_root + '/log/' + sample_name, run_conf.project_root + '/script/' + sample_name)


# qc
@follows( link_import_bam )
@follows( markdup )
@follows( filt_sv )
@follows( identify_mutations )
@transform(qc_bamstats_list, formatter(), "{subpath[0][2]}/qc/{subdir[0][0]}/{subdir[0][0]}.bamstats")
def bam_stats(input_file, output_file):
    dir_name = os.path.dirname(output_file)
    sample_name = os.path.basename(dir_name)
    
    arguments = {"bamstats": genomon_conf.get("SOFTWARE", "bamstats"),
                 "PERL5LIB": genomon_conf.get("ENV", "PERL5LIB"),
                 "input": input_file,
                 "output": output_file}
    
    r_qc_bamstats.task_exec(arguments, run_conf.project_root + '/log/' + sample_name, run_conf.project_root + '/script/' + sample_name)


@follows( link_import_bam )
@follows( markdup )
@follows( filt_sv )
@follows( identify_mutations )
@transform(qc_coverage_list, formatter(), "{subpath[0][2]}/qc/{subdir[0][0]}/{subdir[0][0]}.coverage")
def coverage(input_file, output_file):
    dir_name = os.path.dirname(output_file)
    sample_name = os.path.basename(dir_name)
    depth_output_file = dir_name+'/'+sample_name+'.depth'

    incl_bed_file = ""
    genome_file = ""
    data_type = "exome"
    if genomon_conf.get("qc_coverage", "wgs_flag") == "True":
        genome_file = genomon_conf.get("REFERENCE", "hg19_genome")
        incl_bed_file = output_file + "genome.bed"
        incl_bed_w = genomon_conf.get("qc_coverage", "wgs_incl_bed_width")
        r_qc_coverage.create_incl_bed_wgs(genome_file, incl_bed_file, long(incl_bed_w), "")
        data_type = "wgs"

    arguments = {"data_type": data_type,
                 "i_bed_lines": genomon_conf.get("qc_coverage", "wgs_i_bed_lines"),
                 "i_bed_size": genomon_conf.get("qc_coverage", "wgs_i_bed_width"),
                 "incl_bed_file": incl_bed_file,
                 "genome_file": genome_file,
                 "gaptxt": genomon_conf.get("REFERENCE", "gaptxt"),
                 "bait_file": genomon_conf.get("REFERENCE", "bait_file"),
                 "BEDTOOLS": genomon_conf.get("SOFTWARE", "bedtools"),
                 "SAMTOOLS": genomon_conf.get("SOFTWARE", "samtools"),
                 "LD_LIBRARY_PATH": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "input": input_file,
                 "output": depth_output_file}

    r_qc_coverage.task_exec(arguments, run_conf.project_root + '/log/' + sample_name , run_conf.project_root + '/script/' + sample_name)
    
    r_qc_coverage.calc_coverage(depth_output_file, genomon_conf.get("qc_coverage", "coverage"), output_file)
    
    os.unlink(dir_name+'/'+sample_name+'.depth')
    os.unlink(dir_name+'/'+sample_name+'.depth.input_bed')


@follows( bam_stats )
@follows( coverage )
@collate(qc_merge_list, formatter(), "{subpath[0][2]}/qc/{subdir[0][0]}/{subdir[0][0]}.genomonQC.result.txt")
def merge_qc(input_files, output_file):

    r_qc_merge.write_qc(input_files[0], output_file, get_meta_info(["genomon_pipeline"]))

#####################
# post analysis stage
@active_if(genomon_conf.getboolean("post_analysis", "enable"))
@active_if(len(pa_output_files_mutation) > 0)
@follows(filt_sv)
@follows(identify_mutations)
@collate(pa_files_mutation, formatter(), pa_output_files_mutation.values())
def post_analysis_mutation(input_files, output_file):
        
    arguments = {"pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),
                 "genomon_pa":  genomon_conf.get("SOFTWARE", "genomon_pa"),
                 "mode": "mutation",
                 "genomon_root": run_conf.project_root,
                 "output_dir": run_conf.project_root + "/post_analysis/" + sample_conf_name,
                 "sample_sheet": os.path.abspath(run_conf.sample_conf_file),
                 "config_file": genomon_conf.get("post_analysis", "config_file"),
                 "samtools": genomon_conf.get("SOFTWARE", "samtools"),
                 "bedtools": genomon_conf.get("SOFTWARE", "bedtools"),
                 "input_file_case1": ",".join(pa_samples_mutation["case1"]),
                 "input_file_case2": ",".join(pa_samples_mutation["case2"]),
                 "input_file_case3": ",".join(pa_samples_mutation["case3"]),
                 "input_file_case4": ",".join(pa_samples_mutation["case4"])
                }
                 
    r_post_analysis.task_exec(arguments, run_conf.project_root + '/log/post_analysis', run_conf.project_root + '/script/post_analysis')
    
@active_if(genomon_conf.getboolean("post_analysis", "enable"))
@active_if(len(pa_output_files_sv) > 0)
@follows(filt_sv)
@follows(identify_mutations)
@collate(pa_files_sv, formatter(), pa_output_files_sv.values())
def post_analysis_sv(input_files, output_file):

    arguments = {"pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),
                 "genomon_pa":  genomon_conf.get("SOFTWARE", "genomon_pa"),
                 "mode": "sv",
                 "genomon_root": run_conf.project_root,
                 "output_dir": run_conf.project_root + "/post_analysis/" + sample_conf_name,
                 "sample_sheet": os.path.abspath(run_conf.sample_conf_file),
                 "config_file": genomon_conf.get("post_analysis", "config_file"),
                 "samtools": genomon_conf.get("SOFTWARE", "samtools"),
                 "bedtools": genomon_conf.get("SOFTWARE", "bedtools"),
                 "input_file_case1": ",".join(pa_samples_sv["case1"]),
                 "input_file_case2": ",".join(pa_samples_sv["case2"]),
                 "input_file_case3": ",".join(pa_samples_sv["case3"]),
                 "input_file_case4": ",".join(pa_samples_sv["case4"])
                }
                 
    r_post_analysis.task_exec(arguments, run_conf.project_root + '/log/post_analysis', run_conf.project_root + '/script/post_analysis')

@active_if(genomon_conf.getboolean("post_analysis", "enable"))
@active_if(len(pa_output_files_qc) > 0)
@follows(merge_qc)
@collate(pa_files_qc, formatter(), pa_output_files_qc.values())
def post_analysis_qc(input_files, output_file):

    arguments = {"pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),
                 "genomon_pa":  genomon_conf.get("SOFTWARE", "genomon_pa"),
                 "mode": "qc",
                 "genomon_root": run_conf.project_root,
                 "output_dir": run_conf.project_root + "/post_analysis/" + sample_conf_name,
                 "sample_sheet": os.path.abspath(run_conf.sample_conf_file),
                 "config_file": genomon_conf.get("post_analysis", "config_file"),
                 "samtools": genomon_conf.get("SOFTWARE", "samtools"),
                 "bedtools": genomon_conf.get("SOFTWARE", "bedtools"),
                 "input_file_case1": ",".join(pa_samples_qc),
                 "input_file_case2": "",
                 "input_file_case3": "",
                 "input_file_case4": ""
                }
                 
    r_post_analysis.task_exec(arguments, run_conf.project_root + '/log/post_analysis', run_conf.project_root + '/script/post_analysis')
    
@active_if(genomon_conf.getboolean("pa_plot", "enable"))
@active_if(len(paplot_files_collate) > 0)
@follows(post_analysis_mutation)
@follows(post_analysis_sv)
@follows(post_analysis_qc)
@collate(paplot_files_collate, formatter(), run_conf.project_root + '/paplot/' + sample_conf_name + '/index.html')
def post_analysis_plot(input_file, output_file):
    
    if not os.path.isdir(run_conf.project_root + '/paplot/'): os.mkdir(run_conf.project_root + '/paplot/')
    if not os.path.isdir(run_conf.project_root + '/paplot/' + sample_conf_name): os.mkdir(run_conf.project_root + '/paplot/' + sample_conf_name)

    remark = genomon_conf.get("pa_plot", "remarks")
    remark += "<ul>"
    
    for item in genomon_conf.get("pa_plot", "software").split(","):
        key = item.split(":")[0].strip(" ").rstrip(" ")
        name = item.split(":")[1].strip(" ").rstrip(" ")
        try:
            version = get_version(key).split("-")
        except Exception:
            print ("[WARNING] paplot: %s is not defined." % (key))
            continue
        
        remark += "<li>" + name + " " + version[-1] + "</li>"

    remark += "</ul>"
            
    arguments = {"pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),
                 "pa_plot":  genomon_conf.get("SOFTWARE", "pa_plot"),
                 "inputs_qc": ",".join(paplot_files_qc),
                 "inputs_sv": ",".join(paplot_files_sv),
                 "inputs_mutation": ",".join(paplot_files_mutation),
                 "output_dir": run_conf.project_root + "/paplot/" + sample_conf_name,
                 "title": genomon_conf.get("pa_plot", "title"),
                 "remarks": remark,
                 "config_file": genomon_conf.get("pa_plot", "config_file"),
                }
                 
    r_pa_plot.task_exec(arguments, run_conf.project_root + '/log/paplot', run_conf.project_root + '/script/paplot')


