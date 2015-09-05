import os
from ruffus import *
from genomon_pipeline.run_param import *
from genomon_pipeline.genomon_param import *
from genomon_pipeline.task_param import *
from genomon_pipeline.sample_list import *
from genomon_pipeline.rna_resource.star_align import *

def rna_pipeline_run():

    star_align = Star_align(task_param.get("star_align", "qsub_option"), run_param.project_root + '/script')
    
    # generate list of linked_fastq file path
    linked_fastq_list = []
    for sample in sample_list.fastq:
        linked_fastq_list.append([run_param.project_root + '/fastq/' + sample + '/1_1.fastq',
                                  run_param.project_root + '/fastq/' + sample + '/1_2.fastq'])

    sample_list_fastq = sample_list.fastq

    if not os.path.isdir(run_param.project_root): os.mkdir(run_param.project_root)
    if not os.path.isdir(run_param.project_root + '/script'): os.mkdir(run_param.project_root + '/script')
    if not os.path.isdir(run_param.project_root + '/log'): os.mkdir(run_param.project_root + '/log')
    if not os.path.isdir(run_param.project_root + '/fastq'): os.mkdir(run_param.project_root + '/fastq')
    if not os.path.isdir(run_param.project_root + '/star'): os.mkdir(run_param.project_root + '/star')
    if not os.path.isdir(run_param.project_root + '/fusion'): os.mkdir(run_param.project_root + '/fusion')

    @originate(linked_fastq_list, sample_list_fastq)
    def link_input_fastq(output_file, sample_list_fastq):
        sample = os.path.basename(os.path.dirname(output_file[0]))
        link_dir = run_param.project_root + '/fastq/' + sample
        if not os.path.isdir(link_dir): os.mkdir(link_dir)
        os.symlink(sample_list_fastq[sample][0][0], link_dir + '/1_1.fastq')
        os.symlink(sample_list_fastq[sample][1][0], link_dir + '/1_2.fastq')

    @transform(link_input_fastq, formatter(), "{subpath[0][2]}/star/{subdir[0][0]}/{subdir[0][0]}.bam")
    def task_star_align(input_files, output_file):
        dir_name = os.path.dirname(output_file)

        arguments = {"star": genomon_param.get("SOFTWARE", "STAR"),
                  "star_genome": genomon_param.get("REFERENCE", "star_genome"),
                  "additional_params": task_param.get("star_align", "star_params"),
                  "fastq1": input_files[0],
                  "fastq2": input_files[1],
                  "out_prefix": os.path.dirname(output_file),
                  "log": run_param.project_root + '/log'}

        if not os.path.isdir(dir_name): os.mkdir(dir_name)
        star_align.task_exec(arguments)
        # open(output_file, 'w')

    """
    @transform(task_star_align, formatter(), "{subpath[0][2]}/fusion/{subdir[0][0]}/{subdir[0][0]}.fusion.txt")
    def task_fusionfusion(input_files, output_file):
        dir_name = os.path.dirname(output_file)
        if not os.path.isdir(dir_name): os.mkdir(dir_name)
        open(output_file, 'w')
    """    

    pipeline_run(verbose = 3, multithread = 10)

