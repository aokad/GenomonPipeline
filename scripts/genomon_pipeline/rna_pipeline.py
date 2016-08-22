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


# set task classes
star_align = Star_align(genomon_conf.get("star_align", "qsub_option"), run_conf.drmaa)
fusionfusion = Fusionfusion(genomon_conf.get("fusionfusion", "qsub_option"), run_conf.drmaa)
r_pa_plot = Res_PA_Plot(genomon_conf.get("pa_plot", "qsub_option"), run_conf.drmaa)
r_post_analysis = Res_PostAnalysis(genomon_conf.get("post_analysis", "qsub_option"), run_conf.drmaa)

# generate list of linked_fastq file path
linked_fastq_list = []
for sample in sample_conf.fastq:
    linked_fastq_list.append([run_conf.project_root + '/fastq/' + sample + '/1_1.fastq',
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
if not os.path.isdir(run_conf.project_root + '/script/post_analysis'): os.mkdir(run_conf.project_root + '/script/post_analysis')
if not os.path.isdir(run_conf.project_root + '/script/paplot'): os.mkdir(run_conf.project_root + '/script/paplot')
if not os.path.isdir(run_conf.project_root + '/log'): os.mkdir(run_conf.project_root + '/log')
if not os.path.isdir(run_conf.project_root + '/log/post_analysis'): os.mkdir(run_conf.project_root + '/log/post_analysis')
if not os.path.isdir(run_conf.project_root + '/log/paplot'): os.mkdir(run_conf.project_root + '/log/paplot')
if not os.path.isdir(run_conf.project_root + '/fastq'): os.mkdir(run_conf.project_root + '/fastq')
if not os.path.isdir(run_conf.project_root + '/star'): os.mkdir(run_conf.project_root + '/star')
if not os.path.isdir(run_conf.project_root + '/fusion'): os.mkdir(run_conf.project_root + '/fusion')
if not os.path.isdir(run_conf.project_root + '/config'): os.mkdir(run_conf.project_root + '/config')
if (genomon_conf.getboolean("post_analysis", "enable") == True):
    if not os.path.exists(run_conf.project_root + '/post_analysis'): os.mkdir(run_conf.project_root + '/post_analysis')
    if not os.path.exists(run_conf.project_root + '/post_analysis/' + sample_conf_name): os.mkdir(run_conf.project_root + '/post_analysis/' + sample_conf_name)

for sample in sample_conf.fastq:
        script_dir = run_conf.project_root + '/script/' + sample
        log_dir = run_conf.project_root + '/log/' + sample
        if not os.path.isdir(script_dir): os.mkdir(script_dir)
        if not os.path.isdir(log_dir): os.mkdir(log_dir)

genomon_conf_name, ext = os.path.splitext(os.path.basename(run_conf.genomon_conf_file))
shutil.copyfile(run_conf.genomon_conf_file, run_conf.project_root + '/config/' + genomon_conf_name +'_'+ run_conf.analysis_timestamp + ext)

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


@transform(link_input_fastq, formatter(), "{subpath[0][2]}/star/{subdir[0][0]}/{subdir[0][0]}.Aligned.sortedByCoord.out.bam")
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


@transform(task_star_align, formatter(), "{subpath[0][2]}/fusion/{subdir[0][0]}/star.fusion.result.txt")
def task_fusionfusion(input_file, output_file):

    input_dir_name = os.path.dirname(input_file)
    sample_name = os.path.basename(input_dir_name)
    input_chimeric_sam = input_dir_name + '/' + sample_name + ".Chimeric.out.sam"
    output_dir_name = os.path.dirname(output_file) 

    arguments = {"fusionfusion": genomon_conf.get("SOFTWARE", "fusionfusion"),
                 "chimeric_sam": input_chimeric_sam,
                 "output_prefix": output_dir_name,
                 "param_file": genomon_conf.get("fusionfusion", "param_file"),
                 "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),   
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH")}

    if not os.path.isdir(output_dir_name): os.mkdir(output_dir_name)
    fusionfusion.task_exec(arguments, run_conf.project_root + '/log/' + sample_name, run_conf.project_root + '/script/' + sample_name)


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


