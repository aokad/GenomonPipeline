import os
from ruffus import *
from genomon_pipeline.config.run_conf import *
from genomon_pipeline.config.genomon_conf import *
from genomon_pipeline.config.task_conf import *
from genomon_pipeline.config.sample_conf import *
from genomon_pipeline.rna_resource.star_align import *
from genomon_pipeline.rna_resource.mapsplice2_align import *
from genomon_pipeline.rna_resource.tophat2_align import *
from genomon_pipeline.rna_resource.fusionfusion import *
from genomon_pipeline.rna_resource.star_fusion import *
from genomon_pipeline.rna_resource.tophat_fusion import *



# set task classes
star_align = Star_align(task_conf.get("star_align", "qsub_option"), run_conf.project_root + '/script')
mapsplice2_align = Mapsplice2_align(task_conf.get("mapsplice2_align", "qsub_option"), run_conf.project_root + '/script')
tophat2_align = TopHat2_align(task_conf.get("tophat2_align", "qsub_option"), run_conf.project_root + '/script')
fusionfusion = Fusionfusion(task_conf.get("fusionfusion", "qsub_option"), run_conf.project_root + '/script')
star_fusion = Star_fusion(task_conf.get("star_fusion", "qsub_option"), run_conf.project_root + '/script')
tophat_fusion = TopHat_fusion(task_conf.get("tophat_fusion", "qsub_option"), run_conf.project_root + '/script')

# generate list of linked_fastq file path
linked_fastq_list = []
for sample in sample_conf.fastq:
    linked_fastq_list.append([run_conf.project_root + '/fastq/' + sample + '/1_1.fastq',
                              run_conf.project_root + '/fastq/' + sample + '/1_2.fastq'])

sample_list_fastq = sample_conf.fastq

# prepare output directories
if not os.path.isdir(run_conf.project_root): os.makedirs(run_conf.project_root)
if not os.path.isdir(run_conf.project_root + '/script'): os.mkdir(run_conf.project_root + '/script')
if not os.path.isdir(run_conf.project_root + '/log'): os.mkdir(run_conf.project_root + '/log')
if not os.path.isdir(run_conf.project_root + '/fastq'): os.mkdir(run_conf.project_root + '/fastq')
if not os.path.isdir(run_conf.project_root + '/star'): os.mkdir(run_conf.project_root + '/star')
if not os.path.isdir(run_conf.project_root + '/mapsplice2'): os.mkdir(run_conf.project_root + '/mapsplice2')
if not os.path.isdir(run_conf.project_root + '/tophat2'): os.mkdir(run_conf.project_root + '/tophat2')
if not os.path.isdir(run_conf.project_root + '/tophat_fusion'): os.mkdir(run_conf.project_root + '/tophat_fusion')
if not os.path.isdir(run_conf.project_root + '/fusionfusion'): os.mkdir(run_conf.project_root + '/fusionfusion')
if not os.path.isdir(run_conf.project_root + '/star_fusion'): os.mkdir(run_conf.project_root + '/star_fusion')


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
                 "ref_gtf": genomon_conf.get("REFERENCE", "ref_gtf"),
                 "additional_params": task_conf.get("star_align", "star_params"),
                 "samtools": genomon_conf.get("SOFTWARE", "samtools"),
                 "fastq1": input_files[0],
                 "fastq2": input_files[1],
                 "out_prefix": dir_name + '/' + sample_name + '.',
                 "log": run_conf.project_root + '/log'}

    if not os.path.isdir(dir_name): os.mkdir(dir_name)
    star_align.task_exec(arguments)


@transform(link_input_fastq, formatter(), "{subpath[0][2]}/mapsplice2/{subdir[0][0]}/alignments.bam")
def task_mapsplice2_align(input_files, output_file):
    
    dir_name = os.path.dirname(output_file)
    sample_name = os.path.basename(dir_name)
    
    arguments = {"mapsplice2": genomon_conf.get("SOFTWARE", "mapsplice2"),
                 "ref_gtf": genomon_conf.get("REFERENCE", "ref_gtf"),
                 "ref_fasta": genomon_conf.get("REFERENCE", "mapsplice2_ref"),
                 "bow_ind_mp2": genomon_conf.get("REFERENCE", "bowtie_ind_mp2"),
                 "additional_params": task_conf.get("mapsplice2_align", "mapsplice2_params"),
                 "python": genomon_conf.get("SOFTWARE", "python"), 
                 "fastq1": input_files[0],
                 "fastq2": input_files[1],
                 "out_dir": dir_name,
                 "log": run_conf.project_root + '/log'}

    if not os.path.isdir(dir_name): os.mkdir(dir_name)
    mapsplice2_align.task_exec(arguments)


@transform(link_input_fastq, formatter(), "{subpath[0][2]}/tophat2/{subdir[0][0]}/accepted_hits.bam")
def task_tophat2_align(input_files, output_file):
    
    dir_name = os.path.dirname(output_file)
    sample_name = os.path.basename(dir_name)
    
    arguments = {"tophat2_dir": genomon_conf.get("SOFTWARE", "tophat2_dir"),
                 "ref_gtf": genomon_conf.get("REFERENCE", "ref_gtf"),
                 "bowtie2_database": genomon_conf.get("REFERENCE", "bowtie2_db"),
                 "samtools_path": os.path.dirname(genomon_conf.get("SOFTWARE", "samtools")),
                 "bowtie_path": os.path.dirname(genomon_conf.get("SOFTWARE", "bowtie2")),
                 "additional_params": task_conf.get("tophat2_align", "tophat2_params"),
                 "fastq1": input_files[0],
                 "fastq2": input_files[1],
                 "output_dir": dir_name,
                 "log": run_conf.project_root + '/log'}

    if not os.path.isdir(dir_name): os.mkdir(dir_name)
    tophat2_align.task_exec(arguments)


@collate([task_mapsplice2_align, task_star_align, task_tophat2_align], formatter(), "{subpath[0][2]}/fusionfusion/{subdir[0][0]}/fusion_fusion.result.txt")
def task_fusionfusion(input_files, output_file):

    input_dir_name = os.path.dirname(input_files[1])
    sample_name = os.path.basename(input_dir_name)
    input_chimeric_sam = input_dir_name + '/' + sample_name + ".Chimeric.out.sam"
    output_dir_name = os.path.dirname(output_file) 

    arguments = {"fusionfusion": genomon_conf.get("SOFTWARE", "fusionfusion"),
                 "ms2_bam": input_files[0],
                 "star_chimeric_sam": input_chimeric_sam,
                 "th2_bam": input_files[2],
                 "output_prefix": output_dir_name,
                 "param_file": task_conf.get("fusionfusion", "param_file"),
                 "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),   
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "log": run_conf.project_root + '/log'}

    if not os.path.isdir(output_dir_name): os.mkdir(output_dir_name)
    fusionfusion.task_exec(arguments)


@transform(task_star_align, formatter(), "{subpath[0][2]}/star_fusion/{subdir[0][0]}/{subdir[0][0]}.fusion_candidates.txt")
def task_star_fusion(input_file, output_file):

    input_dir_name = os.path.dirname(input_file)
    sample_name = os.path.basename(input_dir_name)
    output_dir_name = os.path.dirname(output_file)

    arguments = {"star_fusion": genomon_conf.get("SOFTWARE", "STAR-Fusion"),
                 "chimeric_sam": input_dir_name + '/' + sample_name + ".Chimeric.out.sam",
                 "chimeric_junction": input_dir_name + '/' + sample_name + ".Chimeric.out.junction",
                 "gtf_file": genomon_conf.get("REFERENCE", "ref_gtf"),
                 "out_prefix": output_dir_name + '/' + sample_name, 
                 "additional_params": task_conf.get("star_fusion", "star_fusion_params"),
                 "environment_variables": genomon_conf.get("ENV", "PERL5LIB"),
                 "log": run_conf.project_root + '/log'}

    if not os.path.isdir(output_dir_name): os.mkdir(output_dir_name)
    star_fusion.task_exec(arguments)

"""
@transform(link_input_fastq, formatter(), "{subpath[0][2]}/tophat2/{subdir[0][0]}")
def task_tophat2_align(input_files, output_file):
    
    dir_name = output_file
    sample_name = os.path.basename(dir_name)
    
    arguments = {"tophat2_dir": genomon_conf.get("SOFTWARE", "tophat2_dir"),
                 "bowtie2_database": genomon_conf.get("REFERENCE", "bowtie2_db"),
                 "samtools_path": os.path.dirname(genomon_conf.get("SOFTWARE", "samtools")),
                 "bowtie_path": os.path.dirname(genomon_conf.get("SOFTWARE", "bowtie2")),
                 "additional_params": task_conf.get("tophat2_align", "tophat2_params"),
                 "output_dir": dir_name,
                 "log": run_conf.project_root + '/log'}
    
    if not os.path.isdir(dir_name): os.mkdir(dir_name)
    tophat2_align.task_exec(arguments)
"""


