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

markdup_bam_list = []
for complist in sample_conf.compare:
    markdup_bam_list.append([run_conf.project_root + '/bam/' + complist[0] + '/' + complist[0] + '.bammarkdup.bam',
                             run_conf.project_root + '/bam/' + complist[1] + '/' + complist[1] + '.bammarkdup.bam'])

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
@collate(map_dna_sequence, formatter(".+/(?P<SAMPLE>.+)/([0-9]+).bamsorted.bam"), "{path[0]}/{SAMPLE[0]}.bammarkdup.bam")
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
# mutation calling stage
@follows( markdup )
@transform(markdup_bam_list, suffix(".bammarkdup.bam"), ".candidate_mutations.tsv")
def identify_mutations(input_files, output_file):
    print "        fisher %s -> %s" % (input_files, output_file)
 
    interval_list = genomon_conf.get("REFERENCE", "interval_list")
    num_lines = sum(1 for line in open(interval_list))

    print num_lines

    # # For test run-----------------------------------------------------------------
    # inputT_prefix, ext = os.path.splitext(input_files[0])
    # inputN_prefix, ext = os.path.splitext(input_files[1])
    # cmd = "samtools view -h -b -o %s %s 17" % (inputT_prefix+".chr17.bam", input_files[0])
    # subprocess.call( cmd , shell=True)
    # cmd = "samtools view -h -b -o %s %s 17" % (inputN_prefix+".chr17.bam", input_files[1])
    # subprocess.call( cmd , shell=True)
    # # -----------------------------------------------------------------------------
 
    # output_prefix, ext = os.path.splitext(output_file)
 
    # map_quality = 30
    # base_quality = 15
    # min_allele = 0.08
    # max_allele = 0.1
    # min_depth = 10
    # min_variant = 4
    # input_disease_bam = inputT_prefix+".chr17.bam"  # test run
    # input_ctrl_bam = inputN_prefix+".chr17.bam"     # test run
    # # input_disease_bam = input_files[0]
    # # input_ctrl_bam = input_files[1]
    # cmd ="fisher comparison -o %s --ref_fa %s --mapping_quality %d --base_quality %d --min_allele_freq %f --max_allele_freq %f --min_depth %d -2 %s -1 %s --samtools_path %s" % (output_prefix+"fisher.tmp", ref_fasta, map_quality, base_quality, min_allele, max_allele, min_depth , input_ctrl_bam, input_disease_bam, samtools_path)
    # subprocess.call( cmd , shell=True)






