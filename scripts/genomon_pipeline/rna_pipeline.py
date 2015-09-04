import os
from ruffus import *
from genomon_pipeline.run_param import *
from genomon_pipeline.sample_list import *


def rna_pipeline_run():

    # generate list of linked_fastq file path
    linked_fastq_list = []
    for sample in sample_list.fastq:
        linked_fastq_list.append([run_param.project_root + '/fastq/' + sample + '/1_1.fastq',
                                  run_param.project_root + '/fastq/' + sample + '/1_2.fastq'])

    sample_list_fastq = sample_list.fastq

    if not os.path.isdir(run_param.project_root): os.mkdir(run_param.project_root)
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
    def map_rna_star(input_files, output_file):
        dir_name = os.path.dirname(output_file)
        if not os.path.isdir(dir_name): os.mkdir(dir_name)
        open(output_file, 'w')

    @transform(map_rna_star, formatter(), "{subpath[0][2]}/fusion/{subdir[0][0]}/{subdir[0][0]}.fusion.txt")
    def fusionfusion(input_files, output_file):
        dir_name = os.path.dirname(output_file)
        if not os.path.isdir(dir_name): os.mkdir(dir_name)
        open(output_file, 'w')

    pipeline_run(verbose = 6)

