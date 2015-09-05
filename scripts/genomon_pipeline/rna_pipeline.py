import os
from ruffus import *
from genomon_pipeline.config.run_conf import *
from genomon_pipeline.config.genomon_conf import *
from genomon_pipeline.config.task_conf import *
from genomon_pipeline.config.sample_conf import *
from genomon_pipeline.rna_resource.star_align import *

def rna_pipeline_run():

    star_align = Star_align(task_conf.get("star_align", "qsub_option"), run_conf.project_root + '/script')
    
    # generate list of linked_fastq file path
    linked_fastq_list = []
    for sample in sample_conf.fastq:
        linked_fastq_list.append([run_conf.project_root + '/fastq/' + sample + '/1_1.fastq',
                                  run_conf.project_root + '/fastq/' + sample + '/1_2.fastq'])

    sample_list_fastq = sample_conf.fastq

    if not os.path.isdir(run_conf.project_root): os.mkdir(run_conf.project_root)
    if not os.path.isdir(run_conf.project_root + '/script'): os.mkdir(run_conf.project_root + '/script')
    if not os.path.isdir(run_conf.project_root + '/log'): os.mkdir(run_conf.project_root + '/log')
    if not os.path.isdir(run_conf.project_root + '/fastq'): os.mkdir(run_conf.project_root + '/fastq')
    if not os.path.isdir(run_conf.project_root + '/star'): os.mkdir(run_conf.project_root + '/star')
    if not os.path.isdir(run_conf.project_root + '/fusion'): os.mkdir(run_conf.project_root + '/fusion')

    @originate(linked_fastq_list, sample_list_fastq)
    def link_input_fastq(output_file, sample_list_fastq):
        sample = os.path.basename(os.path.dirname(output_file[0]))
        link_dir = run_conf.project_root + '/fastq/' + sample
        if not os.path.isdir(link_dir): os.mkdir(link_dir)
        os.symlink(sample_list_fastq[sample][0][0], link_dir + '/1_1.fastq')
        os.symlink(sample_list_fastq[sample][1][0], link_dir + '/1_2.fastq')

    @transform(link_input_fastq, formatter(), "{subpath[0][2]}/star/{subdir[0][0]}/{subdir[0][0]}.bam")
    def task_star_align(input_files, output_file):
        dir_name = os.path.dirname(output_file)

        arguments = {"star": genomon_conf.get("SOFTWARE", "STAR"),
                  "star_genome": genomon_conf.get("REFERENCE", "star_genome"),
                  "additional_params": task_conf.get("star_align", "star_params"),
                  "fastq1": input_files[0],
                  "fastq2": input_files[1],
                  "out_prefix": os.path.dirname(output_file),
                  "log": run_conf.project_root + '/log'}

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

