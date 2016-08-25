import os
import shutil
from ruffus import *
from genomon_pipeline.config.run_conf import *
from genomon_pipeline.config.genomon_conf import *
from genomon_pipeline.config.sample_conf import *
from genomon_pipeline.rna_resource.star_align import *
from genomon_pipeline.rna_resource.fusionfusion import *
from genomon_pipeline.rna_resource.post_analysis import *
from genomon_pipeline.rna_resource.paplot import *
from genomon_pipeline.rna_resource.fusion_count import *
from genomon_pipeline.rna_resource.fusion_merge import *
from genomon_pipeline.rna_resource.genomon_expression import *
from genomon_pipeline.dna_resource.bamtofastq import *

# set task classes
bamtofastq = Bam2Fastq(genomon_conf.get("bam2fastq", "qsub_option"), run_conf.drmaa)
star_align = Star_align(genomon_conf.get("star_align", "qsub_option"), run_conf.drmaa)
fusionfusion = Fusionfusion(genomon_conf.get("fusionfusion", "qsub_option"), run_conf.drmaa)
fusion_count = Fusion_count(genomon_conf.get("fusion_count_control", "qsub_option"), run_conf.drmaa)
fusion_merge = Fusion_merge(genomon_conf.get("fusion_merge_control", "qsub_option"), run_conf.drmaa)
genomon_expression = Genomon_expression(genomon_conf.get("genomon_expression", "qsub_option"), run_conf.drmaa)
r_pa_plot = Res_PA_Plot(genomon_conf.get("pa_plot", "qsub_option"), run_conf.drmaa)
r_post_analysis = Res_PostAnalysis(genomon_conf.get("post_analysis", "qsub_option"), run_conf.drmaa)

# generate list of linked_fastq file path
linked_fastq_list = []
for sample in sample_conf.fastq:
    linked_fastq_list.append([run_conf.project_root + '/fastq/' + sample + '/1_1.fastq',
                              run_conf.project_root + '/fastq/' + sample + '/1_2.fastq'])

# generate output list of 'bam2fastq'
bam2fastq_output_list = []
for sample in sample_conf.bam_tofastq:
    bam2fastq_output_list.append([run_conf.project_root + '/fastq/' + sample + '/1_1.fastq',
                                  run_conf.project_root + '/fastq/' + sample + '/1_2.fastq'])

sample_list_fastq = sample_conf.fastq

sample_conf_name, ext = os.path.splitext(os.path.basename(run_conf.sample_conf_file))

# generate input list of 'post-analysis'
run_pa = False
if not os.path.exists(run_conf.project_root + '/post_analysis/' + sample_conf_name + '/merge_fusionfusion.txt'):
    run_pa = True
elif not os.path.exists(run_conf.project_root + '/post_analysis/' + sample_conf_name + '/merge_starqc.txt'):
    run_pa = True
else:
    for sample in sample_list_fastq:
        if not os.path.exists(run_conf.project_root + '/fusion/' + sample + '/fusion_fusion.result.txt'):
            run_pa = True
            break

pa_files = []
if run_pa == True:
    for sample in sample_list_fastq:
        pa_files.append(run_conf.project_root + '/fusion/' + sample + '/fusion_fusion.result.txt')

# generate input list of 'paplot'
run_paplot = False
if not os.path.exists(run_conf.project_root + '/paplot/' + sample_conf_name + '/index.html'):
    run_paplot = True
elif run_pa == True:
    run_paplot = True

paplot_files = []
if run_paplot == True:
    paplot_files.append(run_conf.project_root + '/post_analysis/' + sample_conf_name + '/merge_fusionfusion.txt')
    paplot_files.append(run_conf.project_root + '/post_analysis/' + sample_conf_name + '/merge_starqc.txt')

# prepare output directories
if not os.path.isdir(run_conf.project_root): os.makedirs(run_conf.project_root)
if not os.path.isdir(run_conf.project_root + '/script'): os.mkdir(run_conf.project_root + '/script')
if not os.path.isdir(run_conf.project_root + '/script/fusion_merge'): os.mkdir(run_conf.project_root + '/script/fusion_merge')
if not os.path.isdir(run_conf.project_root + '/script/post_analysis'): os.mkdir(run_conf.project_root + '/script/post_analysis')
if not os.path.isdir(run_conf.project_root + '/script/paplot'): os.mkdir(run_conf.project_root + '/script/paplot')
if not os.path.isdir(run_conf.project_root + '/log'): os.mkdir(run_conf.project_root + '/log')
if not os.path.isdir(run_conf.project_root + '/log/fusion_merge'): os.mkdir(run_conf.project_root + '/log/fusion_merge')
if not os.path.isdir(run_conf.project_root + '/log/post_analysis'): os.mkdir(run_conf.project_root + '/log/post_analysis')
if not os.path.isdir(run_conf.project_root + '/log/paplot'): os.mkdir(run_conf.project_root + '/log/paplot')
if not os.path.isdir(run_conf.project_root + '/fastq'): os.mkdir(run_conf.project_root + '/fastq')
if not os.path.isdir(run_conf.project_root + '/star'): os.mkdir(run_conf.project_root + '/star')
if not os.path.isdir(run_conf.project_root + '/fusion'): os.mkdir(run_conf.project_root + '/fusion')
if not os.path.isdir(run_conf.project_root + '/fusion/control_panel'): os.mkdir(run_conf.project_root + '/fusion/control_panel')
if not os.path.isdir(run_conf.project_root + '/config'): os.mkdir(run_conf.project_root + '/config')
if not os.path.isdir(run_conf.project_root + '/expression'): os.mkdir(run_conf.project_root + '/expression')
if (genomon_conf.getboolean("post_analysis", "enable") == True):
    if not os.path.exists(run_conf.project_root + '/post_analysis'): os.mkdir(run_conf.project_root + '/post_analysis')
    if not os.path.exists(run_conf.project_root + '/post_analysis/' + sample_conf_name): os.mkdir(run_conf.project_root + '/post_analysis/' + sample_conf_name)

for target_sample_dict in (sample_conf.bam_import, sample_conf.fastq, sample_conf.bam_tofastq):
    for sample in target_sample_dict:
        script_dir = run_conf.project_root + '/script/' + sample
        log_dir = run_conf.project_root + '/log/' + sample
        if not os.path.isdir(script_dir): os.mkdir(script_dir)
        if not os.path.isdir(log_dir): os.mkdir(log_dir)

genomon_conf_name, ext = os.path.splitext(os.path.basename(run_conf.genomon_conf_file))
sample_conf_name, ext = os.path.splitext(os.path.basename(run_conf.sample_conf_file))
shutil.copyfile(run_conf.genomon_conf_file, run_conf.project_root + '/config/' + genomon_conf_name +'_'+ run_conf.analysis_timestamp + ext)
shutil.copyfile(run_conf.sample_conf_file, run_conf.project_root + '/config/' + sample_conf_name +'_'+ run_conf.analysis_timestamp + ext)

expression_bams = []
# generate input list of genomon expression
for tumor_sample in sample_conf.expression:
    expression_bams.append(run_conf.project_root + '/star/' + tumor_sample + '/' + tumor_sample + '.Aligned.sortedByCoord.out.bam')

fusionfusion_bams = []
fusion_control_panel = []
chimeric_ctrl_sams = []
# generate input list of fusionfusion
for complist in sample_conf.fusionfusion:

    tumor_sample = complist[0]
    control_panel_name = complist[1]

    # generate input list of 'fusionfusion'
    tumor_bam = run_conf.project_root + '/star/' + tumor_sample + '/' + tumor_sample + '.Aligned.sortedByCoord.out.bam'
    merged_count = run_conf.project_root + '/fusion/control_panel/' + control_panel_name + ".merged.Chimeric.count" if control_panel_name != None else None
    fusionfusion_bams.append([tumor_bam, merged_count])
    
    if control_panel_name != None:
        # generate input list of 'fusion count'
        for panel_sample in sample_conf.control_panel[control_panel_name]:
            if panel_sample not in chimeric_ctrl_sams:
                chimeric_ctrl_sams.append(run_conf.project_root + '/star/' + panel_sample + '/' + panel_sample + '.Chimeric.out.sam')

        # generate input list of 'fusion merge control'
        control_panel_file = run_conf.project_root + '/fusion/control_panel/' + control_panel_name + ".Chimeric_count.list"
        fusion_control_panel.append(control_panel_file)
        with open(control_panel_file,  "w") as out_handle:
            for panel_sample in sample_conf.control_panel[control_panel_name]:
                out_handle.write(run_conf.project_root + '/fusion/' + panel_sample + '/' + panel_sample + '.Chimeric.count' + "\n")


# link the import bam to project directory
@originate(sample_conf.bam_import.keys())
def link_import_bam(sample):
    bam = sample_conf.bam_import[sample]
    link_dir = run_conf.project_root + '/star/' + sample
    bam_prefix, ext = os.path.splitext(bam)
    input_dir_name = os.path.dirname(bam)
    sample_name = os.path.basename(input_dir_name)
    input_chimeric_sam = input_dir_name + '/' + sample_name + ".Chimeric.out.sam"
    
    if not os.path.isdir(link_dir): os.mkdir(link_dir)
    if (not os.path.exists(link_dir +'/'+ sample +'.Aligned.sortedByCoord.out.bam')) and (not os.path.exists(link_dir +'/'+ sample +'.Aligned.sortedByCoord.out.bam.bai')): 
        os.symlink(bam, link_dir +'/'+ sample +'.Aligned.sortedByCoord.out.bam')
        if (os.path.exists(bam +'.bai')):
            os.symlink(bam +'.bai', link_dir +'/'+ sample +'.Aligned.sortedByCoord.out.bam.bai')
        elif (os.path.exists(bam_prefix +'.bai')):
            os.symlink(bam_prefix +'.bai', link_dir +'/'+ sample +'.Aligned.sortedByCoord.out.bam.bai')
    if (not os.path.exists(link_dir +'/'+ sample +'.Chimeric.out.sam')): 
        os.symlink(input_chimeric_sam, link_dir +'/'+ sample +'.Chimeric.out.sam')

# convert bam to fastq
@originate(bam2fastq_output_list)
def bam2fastq(outputfiles):
    sample = os.path.basename(os.path.dirname(outputfiles[0]))
    output_dir = run_conf.project_root + '/fastq/' + sample
            
    arguments = {"biobambam": genomon_conf.get("SOFTWARE", "biobambam"),
                 "input_bam": sample_conf.bam_tofastq[sample],
                 "f1_name": outputfiles[0],
                 "f2_name": outputfiles[1],
                 "o1_name": output_dir + '/unmatched_first_output.txt',
                 "o2_name": output_dir + '/unmatched_second_output.txt',
                 "t": output_dir + '/temp.txt',
                 "s": output_dir + '/single_end_output.txt'}

    if not os.path.isdir(output_dir): os.mkdir(output_dir)
    bamtofastq.task_exec(arguments, run_conf.project_root + '/log/' + sample, run_conf.project_root + '/script/'+ sample)

# link the input fastq files
@originate(linked_fastq_list, sample_list_fastq)
def link_input_fastq(output_file, sample_list_fastq):
    sample = os.path.basename(os.path.dirname(output_file[0]))
    link_dir = run_conf.project_root + '/fastq/' + sample

    # Todo
    # 1. should compare the timestamps between input and linked file
    # 2. check md5sum ?
    if not os.path.isdir(link_dir): os.mkdir(link_dir)
    if not os.path.exists(link_dir + '/1_1.fastq'): os.symlink(sample_list_fastq[sample][0][0], link_dir + '/1_1.fastq')
    if not os.path.exists(link_dir + '/1_2.fastq'): os.symlink(sample_list_fastq[sample][1][0], link_dir + '/1_2.fastq')


@transform([bam2fastq,link_input_fastq], formatter(), "{subpath[0][2]}/star/{subdir[0][0]}/{subdir[0][0]}.Aligned.sortedByCoord.out.bam")
def task_star_align(input_files, output_file):

    dir_name = os.path.dirname(output_file)
    sample_name = os.path.basename(dir_name)

    arguments = {"star": genomon_conf.get("SOFTWARE", "STAR"),
                 "star_genome": genomon_conf.get("REFERENCE", "star_genome"),
                 "additional_params": genomon_conf.get("star_align", "star_params"),
                 "samtools": genomon_conf.get("SOFTWARE", "samtools"),
                 "samtools_sort_params": genomon_conf.get("star_align", "samtools_sort_params"),
                 "fastq1": input_files[0],
                 "fastq2": input_files[1],
                 "out_prefix": dir_name + '/' + sample_name + '.'}

    if not os.path.isdir(dir_name): os.mkdir(dir_name)
    star_align.task_exec(arguments, run_conf.project_root + '/log/' + sample_name , run_conf.project_root + '/script/' + sample_name)


@follows( link_import_bam )
@follows( task_star_align )
@transform(chimeric_ctrl_sams, formatter(), "{subpath[0][2]}/fusion/{subdir[0][0]}/{subdir[0][0]}.Chimeric.count")
def task_fusion_count(input_file, output_file):

    input_dir_name = os.path.dirname(input_file)
    sample_name = os.path.basename(input_dir_name)
    output_dir_name = os.path.dirname(output_file) 

    arguments = {"chimera_utils": genomon_conf.get("SOFTWARE", "chimera_utils"),
                 "chimeric_sam": input_file,
                 "output": output_file,
                 "additional_params": genomon_conf.get("fusion_count_control", "params"),
                 "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),   
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH")}

    if not os.path.isdir(output_dir_name): os.mkdir(output_dir_name)
    fusion_count.task_exec(arguments, run_conf.project_root + '/log/' + sample_name, run_conf.project_root + '/script/' + sample_name)


@follows( task_fusion_count )
@transform(fusion_control_panel, formatter(".+/(?P<NAME>.+).Chimeric.count.list"), "{subpath[0][2]}/fusion/control_panel/{NAME[0]}.merged.Chimeric.count")
def task_fusion_merge(input_file, output_file):

    arguments = {"chimera_utils": genomon_conf.get("SOFTWARE", "chimera_utils"),
                 "htslib": genomon_conf.get("SOFTWARE", "htslib"),
                 "count_list": input_file,
                 "output": output_file,
                 "additional_params": genomon_conf.get("fusion_merge_control", "params"),
                 "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),   
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH")}

    fusion_merge.task_exec(arguments, run_conf.project_root + '/log/fusion_merge', run_conf.project_root + '/script/fusion_merge')


@follows( task_fusion_merge )
@transform(fusionfusion_bams, formatter(), "{subpath[0][2]}/fusion/{subdir[0][0]}/star.fusion.result.txt")
def task_fusionfusion(input_file, output_file):

    input_dir_name = os.path.dirname(input_file[0])
    sample_name = os.path.basename(input_dir_name)
    input_chimeric_sam = input_dir_name + '/' + sample_name + ".Chimeric.out.sam"
    output_dir_name = os.path.dirname(output_file) 

    params = ""
    if input_file[1] != None:
        params = "--pooled_control_file " + input_file[1] + " "

    arguments = {"fusionfusion": genomon_conf.get("SOFTWARE", "fusionfusion"),
                 "blat": genomon_conf.get("SOFTWARE", "blat"),
                 "ref_fa":genomon_conf.get("REFERENCE", "ref_fasta"),
                 "chimeric_sam": input_chimeric_sam,
                 "output_prefix": output_dir_name,
                 "annotation_dir": genomon_conf.get("fusionfusion", "annotation_dir"),
                 "additional_params": params + genomon_conf.get("fusionfusion", "params"),
                 "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),   
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH")}

    if not os.path.isdir(output_dir_name): os.mkdir(output_dir_name)
    fusionfusion.task_exec(arguments, run_conf.project_root + '/log/' + sample_name, run_conf.project_root + '/script/' + sample_name)


@follows( link_import_bam )
@follows( task_star_align )
@transform(expression_bams, formatter(), "{subpath[0][2]}/expression/{subdir[0][0]}/{subdir[0][0]}.sym2fkpm.txt")
def task_genomon_expression(input_file, output_file):

    input_dir_name = os.path.dirname(input_file)
    sample_name = os.path.basename(input_dir_name)
    output_dir_name = os.path.dirname(output_file)  

    arguments = {"genomon_expression": genomon_conf.get("SOFTWARE", "genomon_expression"),
                 "bedtools": genomon_conf.get("SOFTWARE", "bedtools"),
                 "input_bam": input_file,
                 "output_prefix": output_dir_name + '/' + sample_name,
                 "annotation_file": genomon_conf.get("genomon_expression", "annotation_file"),
                 "additional_params": genomon_conf.get("genomon_expression", "params"),
                 "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH")}  

    if not os.path.isdir(output_dir_name): os.mkdir(output_dir_name)
    genomon_expression.task_exec(arguments, run_conf.project_root + '/log/' + sample_name, run_conf.project_root + '/script/' + sample_name)



@active_if(genomon_conf.getboolean("post_analysis", "enable"))
@active_if(len(pa_files) > 0)
@follows(task_fusionfusion)
@collate(pa_files, formatter(), paplot_files)
def post_analysis(input_files, output_file):
    
    arguments = {"pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),
                 "genomon_pa":  genomon_conf.get("SOFTWARE", "genomon_pa"),
                 "genomon_root": run_conf.project_root,
                 "output_dir": run_conf.project_root + "/post_analysis/" + sample_conf_name,
                 "sample_sheet": os.path.abspath(run_conf.sample_conf_file),
                 "config_file": genomon_conf.get("post_analysis", "config_file"),
                 "input_file_case1": ",".join(sample_list_fastq.keys()),
                }
                 
    r_post_analysis.task_exec(arguments, run_conf.project_root + '/log/post_analysis', run_conf.project_root + '/script/post_analysis')
    
    
@active_if(genomon_conf.getboolean("pa_plot", "enable"))
@active_if(len(paplot_files) > 0)
@follows(post_analysis)
@collate(paplot_files, formatter(), run_conf.project_root + '/paplot/' + sample_conf_name + '/index.html')
def pa_plot(input_file, output_file):
    
    if not os.path.isdir(run_conf.project_root + '/paplot/'): os.mkdir(run_conf.project_root + '/paplot/')
    if not os.path.isdir(run_conf.project_root + '/paplot/' + sample_conf_name): os.mkdir(run_conf.project_root + '/paplot/' + sample_conf_name)

    # software version in index.html
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

    # input files
    qc_file = run_conf.project_root + '/post_analysis/' + sample_conf_name + '/merge_starqc.txt'
    sv_file = run_conf.project_root + '/post_analysis/' + sample_conf_name + '/merge_fusionfusion.txt'
    print paplot_files
    print qc_file
    print sv_file
    if (qc_file in paplot_files) == False: qc_file = ""
    if (sv_file in paplot_files) == False: sv_file = ""

    arguments = {"pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),
                 "pa_plot":  genomon_conf.get("SOFTWARE", "pa_plot"),
                 "inputs_qc": qc_file,
                 "inputs_sv": sv_file,
                 "output_dir": run_conf.project_root + "/paplot/" + sample_conf_name,
                 "title": genomon_conf.get("pa_plot", "title"),
                 "remarks": remark,
                 "config_file": genomon_conf.get("pa_plot", "config_file"),
                }
                 
    r_pa_plot.task_exec(arguments, run_conf.project_root + '/log/paplot', run_conf.project_root + '/script/paplot')


