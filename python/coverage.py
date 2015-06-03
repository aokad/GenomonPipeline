#! /usr/bin/python
"""

Calculate coverage


"""

import sys
import os
from datetime import datetime
import argparse
import random
import subprocess
from bedparse import BedExtract

#
# Subroutines
#

def PrintHeader( myself, arg ):
    now = datetime.now()

    #print '#' * 84
    #print '# Summary'
    #print '# Generated by %(my)s' % { 'my': myself }
    #print '# %(y)d.%(m)d.%(d)d.%(h)d:%(m)d' % { 'y': now.year, 'm': now.month, 'd': now.day, 'h': now.hour, 'm': now.minute }
    #print '#' * 84 + ''
    print 'file:\t{input_bam} '.format( input_bam = arg.input_bam )
    print 'bin:\t{bin_size} '.format( bin_size = arg.bin_size )
    print 'sample_num:\t{sample_num} '.format( sample_num = arg.sample_num )
    print 'coverage:\t{coverage} '.format( coverage = arg.coverage_depth )


#
# Main
#
def main():

#
# Argument parse
#
    argvs = sys.argv
    myself = argvs[ 0 ]
    argc = len(argvs)

    parser = argparse.ArgumentParser( description = "Calculate coverage" )
    parser.add_argument( '-i', '--input_bam', help = "Input BAM file", type = str )
    parser.add_argument( '-t', '--coverage_tmp', help = "Temporary output file", type = str, default = './tmp.txt' )
    parser.add_argument( '-f', '--ref_fasta', help = "Genome FASTA file", type = str, default = '/home/w3varann/database/hg19/hg19.fa' )
    parser.add_argument( '-g', '--genome_size', help = "Genome size file", type = str, default = '/home/w3varann/database/hg19/hg19.chrom.sizes' )
    parser.add_argument( '-e', '--genome_bed', help = "Genome bed file", type = str, default = '/home/w3varann/database/hg19.fa/hg19_minus_gap.bed' )
    parser.add_argument( '-b', '--bin_size', help = "Input BAM file", type = int, default = 1000 )
    parser.add_argument( '-n', '--sample_num', help = "Number of samples to pick", type = int, default = 1000 )
    parser.add_argument( '-s', '--samtools', help = "Path to samtools", type = str, default = '/home/w3varann/tools/samtools-1.2/samtools' )
    parser.add_argument( '-c', '--coverage_depth', help = "List of coverage depth", type = str, default = '2,10,20,30,40,50,100' )
    parser.add_argument( '-r', '--chr_str', help = "Add 'chr' in bed", dest='chr_str', action = 'store_true' )
    parser.set_defaults( chr_str=False )

    arg = parser.parse_args()
    if not arg.input_bam:
        print parser.print_help();
        sys.exit( 1 )

    #
    # Print header
    #
    PrintHeader( myself, arg )

    try:
        #
        # Parse genomesize file
        #
        gen_bed = BedExtract( arg.genome_bed, arg.genome_size, bin_size = arg.bin_size, chr_str = arg.chr_str )

        #
        # Make index file if not exist
        #
        if not os.path.exists( arg.input_bam + '.bai' ):
            samtools_index_cmd = '{samtools} index {bam}'.format( 
                                    samtools = arg.samtools,
                                    bam = arg.input_bam,
                                    )
            process = subprocess.Popen( samtools_index_cmd,
                              shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)
            std_out, std_err = process.communicate()
            p_return_code = process.returncode

        #
        # Calculate coverage
        #
        cov = 0
        sum = 0
        coverage = {}
        reads = 0
        for i in range( 0, arg.sample_num ):
            #
            # Run samtools mpileup
            #
            mpileup_file = arg.coverage_tmp + str( i )
            f_first = True
#            f_first = False
            while f_first:
#            while f_first or 0 == os.path.getsize( mpileup_file ):
                f_first = False
                position = gen_bed.get_random_interval()
                #print "{0}:{1}".format( position[ 0 ], position[ 1 ] )

                samtools_mpileup_cmd = '{samtools} mpileup -f {fa} -r {chr}:{start}-{end} {bam} > {mpileup}'.format(
                                            samtools = arg.samtools,
                                            fa = arg.ref_fasta,
                                            bam = arg.input_bam,
                                            chr = position[ 0 ], start = position[ 1 ], end = position[ 1 ] + arg.bin_size - 1,
                                            mpileup = mpileup_file
                                        )
                process = subprocess.Popen( samtools_mpileup_cmd,
                                  shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
                std_out, std_err = process.communicate()
                p_return_code = process.returncode
                #print os.path.getsize( mpileup_file )

            #
            # Calculate coverage
            #
            with open( mpileup_file ) as f:
                for line in f:
                    line_list = line.split( "\t" )
                    reads += 1
                    if line_list[ 2 ] != 'N':
                        cov += 1
                        for num in arg.coverage_depth.split( ',' ):
                            sum += int( line_list[ 3 ] )
                            if int( line_list[ 3 ] ) >= int( num ):
                                if num in coverage:
                                    coverage[ num ] += 1
                                else:
                                    coverage[ num ] = 1

            os.remove( mpileup_file )

        #
        # Output result
        #
        data_string =  "non-N_total_depth\ttotal_bases\tnon-N_bases"
        for num in arg.coverage_depth.split( ',' ):
            data_string += "\t{0}x\t{0}x_ratio".format( num )
        print data_string

        data_string = "\n{0}\t{1}\t{2}".format( sum, reads, cov )
        for num in arg.coverage_depth.split( ',' ):
            if num in coverage:
                data_string += "\t{num}\t{ratio}".format(
                                cov = num,
                                num = coverage[ num ],
                                ratio = float( coverage[ num ] )/float( cov ) ),
            else:
                data_string += "\t0\t0"

        print data_string

    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print( "Unexpected error: {error}".format( error = sys.exc_info()[0] ) )
        print("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )



if __name__ == "__main__":
    main()

