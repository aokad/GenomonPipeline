import os
from ruffus import *
# import subprocess

project_root = "/home/yshira/project/rnatest"

sample_fastq = {'MCF-7': ['/home/yshira/Genomon/rna_pipeline/MCF-7/sequence1.txt', '/home/yshira/Genomon/rna_pipeline/MCF-7/sequence2.txt'],
                'MCF-8': ['/home/yshira/Genomon/rna_pipeline/MCF-8/sequence1.txt', '/home/yshira/Genomon/rna_pipeline/MCF-8/sequence2.txt']}


def linked_fastq(project_root, sample_fastq):
    linked_fastq_list = []
    for sample in sample_fastq:
        linked_fastq_list.append([project_root + '/fastq/' + sample + '/1_1.fastq',
                                  project_root + '/fastq/' + sample + '/1_2.fastq'])

    return linked_fastq_list


# print linked_fastq(project_root, sample_fastq)

if not os.path.isdir(project_root): os.mkdir(project_root)
if not os.path.isdir(project_root + '/fastq'): os.mkdir(project_root + '/fastq')
if not os.path.isdir(project_root + '/star'): os.mkdir(project_root + '/star')
if not os.path.isdir(project_root + '/fusion'): os.mkdir(project_root + '/fusion')

@originate (linked_fastq(project_root, sample_fastq), sample_fastq )
def link_input_fastq(output_file, sample_fastq):
    sample = os.path.basename(os.path.dirname(output_file[0]))
    link_dir = project_root + '/fastq/' + sample
    if not os.path.isdir(link_dir): os.mkdir(link_dir)
    os.symlink(sample_fastq[sample][0], link_dir + '/1_1.fastq')
    os.symlink(sample_fastq[sample][1], link_dir + '/1_2.fastq')

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
