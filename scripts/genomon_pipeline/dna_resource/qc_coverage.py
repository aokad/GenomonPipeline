#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Res_QC_Coverage(Stage_task):

    task_name = "qc_coverage"

    script_template = """
#!/bin/bash
#
# Set SGE
#
#$ -S /bin/bash         # set shell in UGE
#$ -cwd                 # execute at the submitted dir
pwd                     # print current working directory
hostname                # print hostname
date                    # print date
set -xv
set -eu

export LD_LIBRARY_PATH={LD_LIBRARY_PATH}

input_bed={output}.input_bed

if [ {data_type} = "wgs" ]
then
    ########## WGS ##########
    
    # for hg19, create gap text (bedtools shuffle -excl option file)
    gaptxt_cut={output}.gap.txt
    cut -f 2,3,4 {gaptxt} | cut -c 4- > $gaptxt_cut
    
    # create temp bed (bedtools shuffle -i option file)
    temp_bed={output}_{i_bed_lines}_{i_bed_size}.bed
    
    echo '' > $temp_bed
    for J in `seq 1 {i_bed_lines}`
    do
        printf "1\t0\t%s\n" {i_bed_size} >> $temp_bed
    done
    
    # bedtools shuffle
    {BEDTOOLS} shuffle -i $temp_bed -g {genome_file} -incl {incl_bed_file} -excl $gaptxt_cut > $input_bed
    
    # depth
    if [ -e {output}.tmp ]; then
        rm {output}.tmp
    fi

    cat $input_bed | while read line; do (
        set -- $line
        let start=$2+1
        
        {SAMTOOLS} view -F 3072 -f 2 -b -h {input} $1:$start-$3 > {output}.tmp.bam
        {SAMTOOLS} index {output}.tmp.bam
        {SAMTOOLS} depth -r $1:$start-$3 {output}.tmp.bam >> {output}.tmp

    ) </dev/null; done
 
    rm $temp_bed {incl_bed_file} $gaptxt_cut {output}.tmp.bam {output}.tmp.bam.bai

else
    ########## exome ##########

    # merge bed (bedtools shuffle -incl option file)
    total_l=`cat {bait_file} | wc -l`
    header_l=`grep ^@ {bait_file} | wc -l`
    data_l=`expr $total_l - $header_l`
    tail -$data_l {bait_file} > {output}.noheader.bed
    {BEDTOOLS} sort -i {output}.noheader.bed > {output}.sort.bed
    {BEDTOOLS} merge -i {output}.sort.bed > {output}.merge.bed
    cut -c 4- {output}.merge.bed > $input_bed

    rm {output}.noheader.bed {output}.sort.bed {output}.merge.bed

    # depth
    {SAMTOOLS} view -F 3072 -f 2 -b -h {input} | {SAMTOOLS} depth -b $input_bed - > {output}.tmp
fi

mv {output}.tmp {output}

"""


    def __init__(self, qsub_option, script_dir):
        super(Res_QC_Coverage, self).__init__(qsub_option, script_dir)

    # 
    # create -incl BED, for bedtools shuffle
    #
    def create_incl_bed_wgs(self, input_genome, output_bed, width, suffix):
        import os
        import math
       
        if not os.path.exists(input_genome):
            print "Not exist file, " + input_genome
            return
    
        # read genome file to dict.
        f = open(input_genome, "r")
        genomes = {}    
        for line in f:
            cells = line.split("\t")
            if len(cells) < 2:
                break
            genomes.update({cells[0]:long(cells[1].rstrip('rn'))})
    
        f.close()
        
        # write BED file
        f = open(output_bed, "w") 
        for i in range(24):
            if i == 22:
                title = "chrX"
            elif i == 23:
                title = "chrY"        
            else:
                title = "chr" + str(i+1)
            
            size = genomes[title]
            start = 0
            end = width
            
            line_num = long(math.ceil(float(size)/float(width)))
            for j in range(line_num):
                if end > size:
                    end = size
                
                if i == 22:
                    head = suffix + "X"
                elif i == 23:
                    head = suffix + "Y"
                else:
                    head = suffix + str(i+1)
    
                f.write("{0}\t{1}\t{2}\n".format(head, start, end))
                start = end
                end += width
    
        f.close()
        
    #
    # Calculate coverage
    #
        
    def calc_coverage(self, depth_file, coverage_depth, output):
        import pandas
        import math

        if not os.path.exists(depth_file + ".input_bed"):
            print "Not exist file, " + depth_file + ".input_bed"
            return
            
        bed = pandas.read_csv(depth_file + ".input_bed", header = None, sep='\t', usecols = [1,2])
        all_count = sum(bed[2] - bed[1])
    
        # depth average & coverage
        depth_reader = pandas.read_csv(depth_file, header = None, sep='\t', usecols = [2], chunksize = 10000)
        
        all_sum = 0
        coverage = {}
        for depth in depth_reader:
            all_sum += sum(depth[2])
            
            for num in coverage_depth.split(','):
                filt = depth[(depth[2] >= int(num))]
                if (num in coverage) == True:
                    count = coverage[num][0] + len(filt)
                else:
                    count = len(filt)
    
                coverage[num] = [count, 0]
        
        for num in coverage_depth.split(','):
            ratio = float(coverage[num][0])/float(all_count)
            coverage[num] = [coverage[num][0], ratio]
                
        ave = float(all_sum)/float(all_count)
        
        # depth std
        depth_reader = pandas.read_csv(depth_file, header = None, sep='\t', usecols = [2], dtype = 'float', chunksize = 10000)
        dist2 = 0
        for depth in depth_reader:
            dist2 += sum((depth[2] - ave)**2)
        
        std = math.sqrt(dist2/float(all_count))
    
        #
        # Output result
        #
        header_string =  "total_depth\tbait_size\taverage_depth\tdepth_stdev"
        data_string = "{0}\t{1}\t{2}\t{3}".format(all_sum, all_count, ave, std)
        
        num_header_string = ""
        ratio_header_string = ""
        num_string = ""
        ratio_string = ""
        
        for num in coverage_depth.split(','):
        
            num_header_string += "\t{0}x".format(num)
            ratio_header_string += "\t{0}x_ratio".format(num)
            num_string += "\t{0}".format(coverage[num][0])
            ratio_string += "\t{0}".format(coverage[num][1])
            
        header_string += ratio_header_string + num_header_string
        data_string += ratio_string + num_string
        
        f = open(output, "w")
        f.write(header_string)
        f.write('\n')
        f.write(data_string)
        f.close()
