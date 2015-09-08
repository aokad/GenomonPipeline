import os
from ruffus import *
from genomon_pipeline.config.run_conf import *
from genomon_pipeline.config.genomon_conf import *
from genomon_pipeline.config.task_conf import *
from genomon_pipeline.config.sample_conf import *
from genomon_pipeline.dna_resource.fastq_splitter import *
from genomon_pipeline.dna_resource.bwa_align import *
from genomon_pipeline.dna_resource.markduplicates import *


# set task classes
fastq_splitter = Fastq_splitter(task_conf.get("split_fast", "qsub_option"), run_conf.project_root + '/script')
bwa_align = Bwa_align(task_conf.get("bwa_mem", "qsub_option"), run_conf.project_root + '/script')
markduplicates = Markduplicates(task_conf.get("markduplicates", "qsub_option"), run_conf.project_root + '/script')

# generate list of linked_fastq file path
linked_fastq_list = []
for sample in sample_conf.fastq:
    linked_fastq_list.append([run_conf.project_root + '/fastq/' + sample + '/1_1.fastq',
                              run_conf.project_root + '/fastq/' + sample + '/1_2.fastq'])

sample_list_fastq = sample_conf.fastq

# prepare output directories
if not os.path.isdir(run_conf.project_root): os.mkdir(run_conf.project_root)
if not os.path.isdir(run_conf.project_root + '/script'): os.mkdir(run_conf.project_root + '/script')
if not os.path.isdir(run_conf.project_root + '/log'): os.mkdir(run_conf.project_root + '/log')
if not os.path.isdir(run_conf.project_root + '/fastq'): os.mkdir(run_conf.project_root + '/fastq')
if not os.path.isdir(run_conf.project_root + '/bam'): os.mkdir(run_conf.project_root + '/bam')

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


##################
#  split stage
@subdivide(link_input_fastq, formatter(), "{path[0]}/*_*.fastq_split", "{path[0]}")
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
@transform(split_files, formatter(".+/1_(?P<NAME>[0-9]+).fastq_split"), add_inputs("{path[0]}/2_{NAME[0]}.fastq_split"), "{subpath[0][2]}/bam/{subdir[0][0]}/{NAME[0]}.bamsorted.bam")
def map_dna_sequence(input_files, output_file):
    
    dir_name = os.path.dirname(output_file)
    output_bwa_sam = output_file.replace('.bamsorted.bam', '.bwa.sam')
         
    arguments = {"bwa": genomon_conf.get("SOFTWARE", "bwa"),
                 "bwa_params": task_conf.get("bwa_mem", "bwa_params"),
                 "ref_fa":genomon_conf.get("REFERENCE", "ref_fasta"),
                 "fastq1": input_files[0],
                 "fastq2": input_files[1],
                 "sam": output_bwa_sam,
                 "biobambam": genomon_conf.get("SOFTWARE", "biobambam"),
                 "bam": output_file,
                 "log": run_conf.project_root + '/log'}

    if not os.path.isdir(dir_name): os.mkdir(dir_name)
    bwa_align.task_exec(arguments) 
 

###################
# merge stage
@collate(map_dna_sequence, formatter(".+/([0-9]+).bamsorted.bam"), "{path[0]}/ga.bammarkdup.bam")
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




