import os
from ruffus import *
from random import randint
import subprocess

bwa_path = "/home/w3varann/tools/bwa-0.7.8/bwa"
ref_fasta = "/home/w3varann/database/GRCh37/GRCh37.fa"
biobambam_path = "/home/kchiba/work_genomonprj/TEST/subdivide/biobambam-0.0.185-release-20150116094537-x86_64-etch-linux-gnu/bin"
samtools_path = "/home/w3varann/tools/samtools-1.2/samtools"
originate_fastq = [['Genomon_german21_test/data/analysis/021T/1_1.fastq', 'Genomon_german21_test/data/analysis/021T/1_2.fastq'], ['Genomon_german21_test/data/analysis/021N/1_1.fastq','Genomon_german21_test/data/analysis/021N/1_2.fastq']]
originate_bam = [['/home/kchiba/work_genomonprj/TEST/subdivide/Genomon_german21_test/data/analysis/021T/ga.bammarkdup.bam', '/home/kchiba/work_genomonprj/TEST/subdivide/Genomon_german21_test/data/analysis/021N/ga.bammarkdup.bam']]

@originate(originate_fastq)
def create_files(output_file): pass

###################
# split stage
@subdivide(create_files, formatter(),
    # Output parameter: Glob matches any number of output file names
    "{path[0]}/*_*.fastq_split",
    # Extra parameter:  Append to this for output file names
    "{path[0]}/")
def subdivide_files(input_files, output_files, output_file_name_stem):

    for oo in output_files:
        os.unlink(oo)

    pair_id = 0
    for input_file in input_files:
        pair_id += 1

        output_file_prefix = "%s/%d_" % (output_file_name_stem, pair_id)
        cmd = "split -a 3 -d -l 4000000 %s %s" % (input_file, output_file_prefix) 
        subprocess.call( cmd , shell=True )
        cmd = "ls -1 %s/%d_[0-9][0-9][0-9] | while read filename; do mv $filename $filename.fastq_split; done" % (output_file_name_stem, pair_id)
        subprocess.call( cmd , shell=True )


###################
# mapping stage
@transform(subdivide_files, formatter(".+/1_(?P<NAME>[0-9]+).fastq_split"), add_inputs("{path[0]}/2_{NAME[0]}.fastq_split"), "{path[0]}/{NAME[0]}.bamsorted.bam")
def map_dna_sequence(input_files, output_file):
    print "        alignment %s -> %s" % (input_files, output_file)

    output_bwa_sam = output_file.replace('.bamsorted.bam', '.bwa.sam')
    min_score = 0
    
    cmd = "%s mem -T %d %s %s %s > %s" % (bwa_path, min_score, ref_fasta, input_files[0], input_files[1], output_bwa_sam)
    #  -R '{read_group}' {additional_params}
    subprocess.call( cmd , shell=True )

    cmd = "%s/bamsort index=1 level=1 inputthreads=2 outputthreads=2 calmdnm=1 calmdnmrecompindentonly=1 calmdnmreference=%s tmpfile=%s inputformat=sam indexfilename=%s I=%s O=%s " % (biobambam_path, ref_fasta, output_file+".tmp", output_file+".bai", output_bwa_sam, output_file) 
    subprocess.call( cmd , shell=True )


###################
# merge stage
@collate(map_dna_sequence, formatter(".+/([0-9]+).bamsorted.bam"), "{path[0]}/ga.bammarkdup.bam")
def merge_sequences(input_files, output_file):
    print "        merge %s -> %s" % (input_files, output_file)

    output_prefix, ext = os.path.splitext(output_file)

    input_bam_files = ""
    for input_file in input_files:
        input_bam_files = input_bam_files + " I=" + input_file

    cmd = "%s/bammarkduplicates M=%s tmpfile=%s markthreads=2 rewritebam=1 rewritebamlevel=1 index=1 md5=1 %s O=%s" % (biobambam_path, output_prefix+".metrics", output_prefix+"tmp", input_bam_files, output_file)
    subprocess.call( cmd , shell=True )


###################
# set of starting bam files
@originate(originate_bam)
def pair_files(output_file): pass


###################
# fisher stage
@transform(pair_files, suffix(".bammarkdup.bam"), ".candidate_mutations.tsv")
def identify_mutations(input_files, output_file):
    print "        fisher %s -> %s" % (input_files, output_file)

    # For test run-----------------------------------------------------------------
    inputT_prefix, ext = os.path.splitext(input_files[0])
    inputN_prefix, ext = os.path.splitext(input_files[1])
    cmd = "samtools view -h -b -o %s %s 17" % (inputT_prefix+".chr17.bam", input_files[0])
    # subprocess.call( cmd , shell=True)
    cmd = "samtools view -h -b -o %s %s 17" % (inputN_prefix+".chr17.bam", input_files[1])
    # subprocess.call( cmd , shell=True)
    # -----------------------------------------------------------------------------

    output_prefix, ext = os.path.splitext(output_file)

    map_quality = 30
    base_quality = 15
    min_allele = 0.08
    max_allele = 0.1
    min_depth = 10
    min_variant = 4
    input_disease_bam = inputT_prefix+".chr17.bam"  # test run
    input_ctrl_bam = inputN_prefix+".chr17.bam"     # test run
    # input_disease_bam = input_files[0]
    # input_ctrl_bam = input_files[1]
    cmd ="fisher comparison -o %s --ref_fa %s --mapping_quality %d --base_quality %d --min_allele_freq %f --max_allele_freq %f --min_depth %d -2 %s -1 %s --samtools_path %s" % (output_prefix+"fisher.tmp", ref_fasta, map_quality, base_quality, min_allele, max_allele, min_depth , input_ctrl_bam, input_disease_bam, samtools_path)
    subprocess.call( cmd , shell=True)


pipeline_run()


