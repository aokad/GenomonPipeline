"""
wgs_pipeline.py

"""

import sys
import os
from datetime import datetime
from ruffus import *
from runtask import RunTask
from glob import glob
import difflib
from sample import Sample

#####################################################################
#
# Private modules
#
from __main__ import *
import genomon_rc as res
import wgs_resource as wgs_res
from utils import *

#####################################################################
#
# Subroutines
#

#####################################################################
#
# File update check
#
def check_file_exists_for_split_fastq(
    input_file1,
    input_file2,
    output_file1,
    output_file2,
    ):

    """
    Checks if output file exists for split_fastq
        glob output files with prefix and suffix

    """

    ( output_prefix, output_suffix ) = os.path.splitext( output_file1 )
    outfile_list = glob( "{prefix}*{suffix}".format( prefix = output_prefix, suffix = output_suffix ) )
    if not outfile_list:
        return True, "Missing file {outprefix}*{outsuffix} for {input}.".format(
                            outprefix = output_prefix,
                            outsuffix = output_suffix,
                            input = input_file1 )
    else:
        in_time = os.path.getmtime( input_file1 )
        out_time = os.path.getmtime( outfile_list[ 0 ] )
        if in_time > out_time:
            return True, "{outprefix}*{outsuffix} is older than {input}.".format(
                                outprefix = output_prefix,
                                outsuffix = output_suffix,
                                input = input_file1 )
        else:
            return False, "File {outprefix}*{outsuffix} exits for {input}.".format(
                                outprefix = output_prefix,
                                outsuffix = output_suffix,
                                input = input_file1 )

def check_file_exists_for_merge_bam(
    input_file1,
    input_file2,
    output_file1,
    output_file2,
    ):

    """
    Checks if output file exists for merge_bam

    """

    ( input_prefix, input_suffix ) = os.path.splitext( input_file1 )
    if not os.path.exists( output_file1 ):
        return True, "Missing file {outputfile} for {input_prefix}*{input_suffix}.".format(
                            outputfile = output_file1,
                            input_prefix = input_prefix,
                            input_suffix = input_suffix )
    else:
        return False, "File {outputfile} exists for {input_prefix}*{input_suffix}.".format(
                            outputfile = output_file1,
                            input_prefix = input_prefix,
                            input_suffix = input_suffix )

def check_file_exists_for_bwa_mem(
        input_file1,
        input_file2,
        output_file1,
        output_file2
    ):

    """
    Checks if output file exists for bwa_mem

    """

    ( output_prefix, output_suffix ) = os.path.splitext( output_file1 )
    if not glob( "{prefix}*{suffix}".format( prefix = output_prefix, suffix = output_suffix ) ):
        return True, "Missing file {outputfile} for {inputfile}.".format(
                            outputfile = output_file1,
                            inputfile = input_file1 )
    else:
        return False, "File {output} exits for {input}.".format( output = output_file1, input = input_file1 )

def check_file_exists_for_markduplicates(
        input_file1,
        input_file2,
        output_file1,
        output_file2
    ):

    """
    Checks if output file exists for bwa_mem

    """

    if Geno.job.get( 'use_biobambam' ):
        if not os.path.exists( output_file1 ):
            return True, "Missing file {outputfile} for {inputfile}.".format(
                                outputfile = output_file1,
                                inputfile = input_file1 )
        else:
            ( input_prefix, input_suffix ) = os.path.splitext( input_file1 )
            input_file_list = glob( "{prefix}*_bamsorted{suffix}".format( prefix = input_prefix, suffix = input_suffix ) )
            if len( input_file_list ) > 0:
                in_time = os.path.getmtime( input_file_list[ 0 ] )
                out_time = os.path.getmtime( output_file1 )
                if in_time > out_time:
                    return True, "{output} is older than {input}.".format( output = output_file1, input = input_file1 )
                else:
                    return False, "File {output} exits for {input}.".format( output = output_file1, input = input_file1 )
            else:
                return False, "Input files do not exist."
    else:
        return check_file_exists_for_input_output( input_file1, input_file2, output_file1, output_file2 )


def check_file_exists_for_input_output(
        input_file1,
        input_file2,
        output_file1,
        output_file2
    ):

    """
    Checks if output file exists

    """

    if not os.path.exists( output_file1 ):
        return True, "Missing file {outputfile} for {inputfile}.".format(
                            outputfile = output_file1,
                            inputfile = input_file1 )

    else:
        in_time = os.path.getmtime( input_file1 )
        out_time = os.path.getmtime( output_file1 )
        if in_time > out_time:
            return True, "{output} is older than {input}.".format( output = output_file1, input = input_file1 )
        else:
            return False, "File {output} exits for {input}.".format( output = output_file1, input = input_file1 )


#####################################################################
#
# Data preparation subrotines for parallelization
#

#
# For STAGE 1 bam2fastq
#
def generate_params_for_bam2fastq():
    """
    Generate parameter list for bam2fastq

    """

    Sample.make_param( 'bam2fastq', '.fastq', 'fastq', 1, 2 )
    for param in Sample.param( 'bam2fastq' ):
        yield param


#
# For STAGE 2 split_fastq
#
def generate_params_for_split_fastq():
    """
    Generate parameter list for spilt_fastq_files

    """
    
    Sample.make_param( 'split_fastq', '.fastq', 'fastq', 2, 2 )
    for param in Sample.param( 'split_fastq' ):
        yield param

#
# For STAGE 3 cutadapt
#
def generate_params_for_cutadapt():
    """
    Generate parameter list for cutadapt

    """
    
    Sample.make_param( 'cutadapt', '_cutadapt.fastq', 'fastq', 2, 2 )
    for param in Sample.param( 'cutadapt' ):
        yield param

#
# For STAGE 4 bwa_mem
#

    
def generate_params_for_bwa_mem():
    """
    Generate parameter list for bwa_mem

    """

    Sample.make_param( 'bwa_mem', '.bam', 'bam', 2, 1 ) 
    for param in Sample.param( 'bwa_mem' ):
        yield param

#
# For STAGE 5 merge_bam        
#        
def generate_params_for_merge_bam():
    """
    Generate parameter list for merge_bam

    """
    
    Sample.make_param( 'merge_bam', '.bam', 'bam', 1, 1 )
    for param in Sample.param( 'merge_bam' ):
        yield param

#
# For STAGE 6 markduplicates
#
def generate_params_for_markduplicates ():
    """
    Generate parameter list for markduplicates

    """
    
    Sample.make_param( 'markduplicates', '_markdup.bam', 'bam', 1, 1 )
    for param in Sample.param( 'markduplicates' ):
        yield param

#
# For STAGE 7 fisher_mutation_call
#
def generate_params_for_fisher_mutation_call():
    """
    Generate parameter list for fisher_mutation_call

    """
    
    Sample.make_param( 'fisher_mutation_call', '.txt', 'mutation', 2, 1 )
    for param in Sample.param( 'fisher_mutation_call' ):
        yield param


#
# Extract compressed file
#
def extract_fastq( input_file_list, file_ext ):
    """
        Extract compressed files

        file_type:  gz
                    bz2
    """
    #
    # Make shell script for array job
    #
    array_in  = "IN_FILE=(\n"
    array_out = "OUT_FILE=(\n"
    id = 0
    for infile in input_file_list:
        id += 1
        output = infile.replace( file_ext, '' )
        array_in  += "[{id}]=\"{file}\"\n".format( id = id, file = infile )
        array_out += "[{id}]=\"{file}\"\n".format( id = id, file = output )
    array_in  += ")\n"
    array_out += ")\n"

    input_file = '${IN_FILE[$SGE_TASK_ID]}'
    output_file = '${OUT_FILE[$SGE_TASK_ID]}'

    shell_script_full_path = make_script_file_name( 'extract_fastq', Geno )
    shell_script_file = open( shell_script_full_path, 'w' )
    shell_script_file.write( res.extract_gz.format(
                                    log = Geno.dir[ 'log' ],
                                    array_data = array_in + array_out,
                                    input_file = input_file,
                                    output_file = output_file ) )
    shell_script_file.close()

    #
    # Run
    #
    return_code = Geno.RT.run_arrayjob(
                            Geno.job.get( 'job_queue' )[ function_name ],
                            Geno.job.get( 'memory' )[ function_name ],
                            shell_script_full_path,
                            '1-{id}:1'.format( id = id ) )
    if return_code != 0:
        log.error( "{function}: runtask failed" )
        raise


#####################################################################
#
#   Main functions
#

#
# Stage 1: bam2fastq
#
def bamtofastq(
    input_file1,
    input_file2,
    output_file1,
    output_file2
    ):
    """
    Split fastq file into small pieces for parallele processing

    """
    try:
        function_name = whoami()
        log.info( "#{function}".format( function = function_name ) )

        #
        # Make sure files exist.
        #
        if not os.path.isfile( input_file ):
            log.error( "file: {file} does not exist.".format( file=input_file ) )
            return 0

        #
        # Make shell script
        #
        shell_script_full_path = make_script_file_name( function_name, Geno )
        shell_script_file = open( shell_script_full_path, 'w' )

        file_type = Geno.job.get( 'input_file_type' )
        if 'paired_bam' == file_type:
            shell_script_file.write( wgs_res.bamtofastq_p.format(
                                            log = Geno.dir[ 'log' ],
                                            bamfile = input_file1,
                                            outfastq1 = output_file1,
                                            outfastq2 = output_file2,
                                            tmpfastq = output_file1 + '.tmp',
                                            bamtofastq = Geno.conf.get( 'SOFTWARE', 'bamtofastq' )
                                            ) )
        elif 'single_bam' == file_type:
            shell_script_file.write( wgs_res.bamtofastq_s.format(
                                            log = Geno.dir[ 'log' ],
                                            bamfile = input_file1,
                                            outfastq1 = output_file1,
                                            tmpfastq = output_file1 + '.tmp',
                                            bamtofastq = Geno.conf.get( 'SOFTWARE', 'bamtofastq' )
                                            ) )
        shell_script_file.close()

        #
        # Run
        #
        return_code = Geno.RT.runtask(
                        Geno.job.get( 'job_queue' )[ function_name ],
                        Geno.job.get( 'memory' )[ function_name ],
                        shell_script_full_path )
        if return_code != 0:
            log.error( "{function}: runtask failed" )
            raise

    except IOError as (errno, strerror):
        log.error( "{function}: I/O error({num}): {error}".format(num = errno, error = strerror, function = whoami() ) )
        return_code = False

    except ValueError:
        log.error( "{function}: ValueError".format( function = whoami() ) )
        return_code = False

    except:
        log.error( "{function}: Unexpected error: {error}".format( function = whoami(), error = sys.exc_info()[0] ) )


    return True

#
# Stage 2: split_fastq
#
def split_fastq(
    input_file1,
    input_file2,
    output_file1,
    output_file2
    ):
    """
    Split fastq file into small pieces for parallele processing

    input file types:
        *.fastq
        *.fastq.gz
        *.fastq.bz2

    """
    try:
        function_name = whoami()
        log.info( "#{function}".format( function = function_name ) )

        #
        # Make sure files exist.
        #
        if not os.path.isfile( input_file1 ):
            log.error( "file: {file} does not exist.".format( file=input_file1 ) )
            return 0

        #
        # Make data for array job 
        #
        array_in  = "IN_FILE=(\n"
        array_out = "OUT_FILE=(\n"
        id = 0
        for infile, outfile in ( [ input_file1, output_file1 ],
                                 [ input_file2, output_file2 ] ):
            id += 1
            out_prefix = make_sample_file_name( outfile, "{dir}/{base}_" )
            array_in  += "[{id}]=\"{file}\"\n".format( id = id, file = infile )
            array_out += "[{id}]=\"{file}\"\n".format( id = id, file = out_prefix )
        array_in  += ")\n"
        array_out += ")\n"

        input_file = '${IN_FILE[$SGE_TASK_ID]}'
        output_prefix = '${OUT_FILE[$SGE_TASK_ID]}'

        output_suffix = '.fastq'
        suffix_len = len( output_suffix )

        #
        # Make shell script for array job
        #
        shell_script_full_path = make_script_file_name( function_name, Geno )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( wgs_res.splitfile.format(
                                        log = Geno.dir[ 'log' ],
                                        fastq_filter = Geno.job.get( 'fastq_filter' ),
                                        array_data = array_in + array_out,
                                        lines_per_file = Geno.job.get( 'split_fastq_line_number' ),
                                        input_file = input_file,
                                        suffix_len = suffix_len,
                                        output_suffix = output_suffix,
                                        output_prefix = output_prefix ) )
        shell_script_file.close()

        #
        # Run
        #
        return_code = Geno.RT.run_arrayjob(
                            Geno.job.get( 'job_queue' )[ function_name ],
                            Geno.job.get( 'memory' )[ function_name ],
                            shell_script_full_path,
                            '1-{id}:1'.format( id = id ) )
        if return_code != 0:
            log.error( "{function}: runtask failed" )
            raise


    except IOError as (errno, strerror):
        log.error( "{function}: I/O error({num}): {error}".format(num = errno, error = strerror, function = whoami() ) )
        return_code = False

    except ValueError:
        log.error( "{function}: ValueError".format( function = whoami() ) )
        return_code = False

    except:
        log.error( "{function}: Unexpected error: {error}".format( function = whoami(), error = sys.exc_info()[0] ) )


    return True


#
# Stage 3: cutadapt
#
def cutadapt(
    input_file1,
    input_file2,
    output_file1,
    output_file2
    ):
    """
       Apply cutadapt to remove adaptor sequences from reads

       Input file type:
        *.fastq
        *.fastq.gz
        *.fastq.bz2

       cutadapt accepts *.fastq and *.fastq.gz.
       *.fastq.bz2 needs to be decompressed

    """

    try:
        function_name = whoami()
        log.info( "#{function}".format( function = function_name ) )

        file_ext = os.path.splitext( input_file1 )[ 1 ]
        if file_ext == '.bz2':
            extract_fastq( ( input_file1, input_file2 ), file_ext )

        #
        # Make data for array job
        #
        if not os.path.isfile( input_file ):
            log.error( "file: {file} does not exist.".format( file=input_file ) )
            return 0


        output_file1 = make_sample_file_name( output_file1, "{dir}/{base}_" ) + '_cutadapt.fastq'
        output_file2 = make_sample_file_name( output_file2, "{dir}/{base}_" ) + '_cutadapt.fastq'

        
        array_data1 = "IN_FILE=(\n[1]=\"{file1}\"\n[2]=\"{file2}\"\n)\n".format(
                        file1 = input_file1,
                        file2 = input_file2 )
        array_data2 = "OUT_FILE=(\n[1]=\"{file1}\"\n[2]=\"{file2}\"\n)\n".format(
                        file1 = output_file1,
                        file2 = output_file2 )

        input_file = '${IN_FILE[$SGE_TASK_ID]}'
        output_file = '${OUT_FILE[$SGE_TASK_ID]}'

        suffix_len = len( output_suffix )

        #
        # Make shell script
        #
        shell_script_full_path = make_script_file_name( function_name, Geno )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( wgs_res.cutadapt.format(
                                        log = Geno.dir[ 'log' ],
                                        infastq = input_file,
                                        outfastq = output_file,
                                        array_data = array_data1 + array_data2,
                                        tmpoutfastq = output_file + '.tmp',
                                        optadapters = '-a ' + ' -a '.join( Geno.job.get( 'adaptor' ) ),
                                        casavacode = 2,
                                        cutadapt = Geno.conf.get( 'SOFTWARE', 'cutadapt' ),
                                        scriptdir = Geno.dir[ 'script' ],
                                        ) )
        shell_script_file.close()

        #
        # Run
        #
        return_code = Geno.RT.run_arrayjob(
                            Geno.job.get( 'job_queue' )[ function_name ],
                            Geno.job.get( 'memory' )[ function_name ],
                            shell_script_full_path,
                            '1-2:1' )
        if return_code != 0:
            log.error( "{function}: runtask failed" )
            raise

    except IOError as (errno, strerror):
        log.error( "{function}: I/O error({num}): {error}".format(function = whoami(), num = errno, error = strerror) )
        return_code = False

    except ValueError:
        log.error( "{function}: ValueError".format( function = whoami() ) )
        return_code = False

    except:
        log.error( "{function}: Unexpected error: {error}".format( function = whoami(), error = sys.exc_info()[0] ) )


    return True

#
# Stage 4: bwa_mem
#
def bwa_mem(
    input_file1,
    input_file2,
    output_file1,
    output_file2,
    use_biobambam
    ):
    """
       Align sequence reads in FASTQ file.
       Convert the result SAM file to BAM file.
    """

    try:
        function_name = whoami()
        log.info( "#{function}".format( function = function_name ) )

        file_ext = os.path.splitext( input_file1 )[ 1 ]

        #
        # Make data for array job 
        #
        fastq_list = []
        for file_name in ( input_file1, input_file2 ):
            ( input_prefix, input_suffix ) = os.path.splitext( file_name )
            fastq_list.append( sorted( glob( "{input_prefix}*{input_suffix}".format(
                                        input_prefix = input_prefix,
                                        input_suffix = input_suffix ) ) ) )
        str1 = "FILE1=(\n"
        str2 = "FILE2=(\n"
        str3 = "FILE3=(\n"
        id = 0
        for fastq1, fastq2 in zip( fastq_list[ 0 ], fastq_list[ 1 ] ):
            output_file = "{dir}/{base}.bam".format(
                                dir = os.path.dirname( output_file1),
                                base = os.path.splitext( os.path.basename( fastq1 ) )[ 0 ] )
            id += 1
            if file_ext == '.bz2':
                str1 += " [{id}]=\"'<bzip2 -dc {file}'\"\n".format( id = id, file = fastq1 )
                str2 += " [{id}]=\"'<bzip2 -dc {file}'\"\n".format( id = id, file = fastq2 )
            else:
                str1 += " [{id}]=\"{file}\"\n".format( id = id, file = fastq1 )
                str2 += " [{id}]=\"{file}\"\n".format( id = id, file = fastq2 )
            str3 += " [{id}]=\"{file}\"\n".format( id = id, file = output_file )
        str1 += " )\n"
        str2 += " )\n"
        str3 += " )\n"

        
        #
        # Make shell script for array job
        #
        shell_script_full_path = make_script_file_name( function_name, Geno )
        shell_script_file = open( shell_script_full_path, 'w' )
        if use_biobambam:
            bwa_mem_resource = wgs_res.bwa_mem_biobambam
        else:
            bwa_mem_resource = wgs_res.bwa_mem
        shell_script_file.write( bwa_mem_resource.format(
                                        log = Geno.dir[ 'log' ],
                                        array_data = str1 + str2 + str3,
                                        fastq1 = "${FILE1[$SGE_TASK_ID]}",
                                        fastq2 = "${FILE2[$SGE_TASK_ID]}",
                                        bam = "${FILE3[$SGE_TASK_ID]}",
                                        rg_id = Geno.job.get( 'rg_id' ),
                                        sample_desc = Geno.job.get( 'sample_desc' ),
                                        library = Geno.job.get( 'library' ),
                                        platform = Geno.job.get( 'platform' ),
                                        platform_unit = Geno.job.get( 'platform_unit' ),
                                        seq_center = Geno.job.get( 'seq_center' ),
                                        pred_med_insert = Geno.job.get( 'pred_med_insert' ),
                                        min_score = Geno.job.get( 'min_score' ),
                                        ref_fa = Geno.conf.get( 'REFERENCE', 'ref_fasta' ),
                                        bwa = Geno.conf.get( 'SOFTWARE', 'bwa' ),
                                        samtools = Geno.conf.get( 'SOFTWARE', 'samtools' ),
                                        biobambam = Geno.conf.get( 'SOFTWARE', 'biobambam' )
                    )
                )
        shell_script_file.close()

        #
        # Run
        #
        return_code = Geno.RT.run_arrayjob(
                            Geno.job.get( 'job_queue' )[ function_name ],
                            Geno.job.get( 'memory' )[ function_name ],
                            shell_script_full_path,
                            "1-{id}:1".format( id = id ) )
        if return_code != 0:
            log.error( "{function}: runtask failed" )
            raise
            


    except IOError as (errno, strerror):
        log.error( "{function}: I/O error({num}): {error}".format(function = whoami(), num = errno, error = strerror) )
        return_code = False

    except ValueError:
        log.error( "{function}: ValueError".format( function = whoami() ) )
        return_code = False

    except:
        log.error( "{function}: Unexpected error: {error}".format( function = whoami(), error = sys.exc_info()[0] ) )


    return True

#
# Stage 5: merge_bam
#
def merge_bam(
    input_file1,
    input_file2,
    output_file1,
    output_file2,
    use_biobambam
    ):
    """
       Merge split bam files

    """

    try:
        function_name = whoami()
        log.info( "#{function}".format( function = function_name ) )

        #
        # Make data for array job 
        #
        ( input_prefix, input_suffix ) = os.path.splitext( input_file1 )

        #
        # Make shell script
        #
        if use_biobambam:
            input_file_name = "{input_prefix}*_bamsorted{input_suffix}".format(
                                        input_prefix = input_prefix,
                                        input_suffix = input_suffix )
            bam_merge_resource = wgs_res.biobambam_merge_bam
            input_files = ''
            for file_name in glob( input_file_name ):
                input_files += "I={file_name} ".format( file_name = file_name )

        else:
            input_files = "{input_prefix}*_sorted{input_suffix}".format(
                                        input_prefix = input_prefix,
                                        input_suffix = input_suffix )
            bam_merge_resource = wgs_res.samtools_merge_bam

        shell_script_full_path = make_script_file_name( function_name, Geno )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( bam_merge_resource.format(
                                        log = Geno.dir[ 'log' ],
                                        input_bam_files = input_files,
                                        output_bam_file = output_file1,
                                        samtools = Geno.conf.get( 'SOFTWARE', 'samtools' ),
                                        biobambam = Geno.conf.get( 'SOFTWARE', 'biobambam' )
                                        ) )
        shell_script_file.close()

        #
        # Run
        #
        return_code = Geno.RT.runtask(
                            Geno.job.get( 'job_queue' )[ 'merge_bam' ],
                            Geno.job.get( 'memory' )[ 'merge_bam' ],
                            shell_script_full_path )
        if return_code != 0:
            log.error( "{function}: runtask failed" )
            raise


    except IOError as (errno, strerror):
        log.error( "{function}: I/O error({num}): {error}".format( function = whoami(), num = errno, error = strerror) )
        return_code = False

    except ValueError:
        log.error( "{function}: ValueError".format( function = whoami() ) )
        return_code = False

    except:
        log.error( "{function}: Unexpected error: {error}".format( function = whoami(), error = sys.exc_info()[0] ) )


    return True

#
# Stage 6: markduplicates
#
def markduplicates(
    input_file1,
    input_file2,
    output_file1,
    output_file2,
    use_biobambam
    ):
    """
       Mutaion calling

    """

    try:
        function_name = whoami()
        log.info( "#{function}".format( function = function_name ) )

        #
        # Make shell script
        #
        if use_biobambam:
            #
            # Use biobambam bammarkduplicates for bamsort results
            #
            #
            # Make shell script for array job
            #
            ( input_prefix, input_suffix ) = os.path.splitext( input_file1 )
            input_file_list =  glob( "{input_prefix}*_bamsorted{input_suffix}".format(
                                        input_prefix = input_prefix,
                                        input_suffix = input_suffix ) )
            input_files =''
            for infile in input_file_list:
                input_files += "I={file} ".format( file = infile )

            shell_script_full_path = make_script_file_name( function_name, Geno )
            shell_script_file = open( shell_script_full_path, 'w' )
            shell_script_file.write( wgs_res.biobambam_markduplicates.format(
                                            log = Geno.dir[ 'log' ],
                                            input_bam_files = input_files,
                                            output_bam = output_file1,
                                            biobambam = Geno.conf.get( 'SOFTWARE', 'biobambam' )
                                        )
                                    )
            shell_script_file.close()

        else:
            #
            # Use Picard MarkDuplicates for samtools merge result
            #
            shell_script_full_path = make_script_file_name( function_name, Geno )
            shell_script_file = open( shell_script_full_path, 'w' )
            shell_script_file.write( wgs_res.markduplicates.format(
                                            log = Geno.dir[ 'log' ],
                                            input_bam = input_file1,
                                            output_bam = output_file1,
                                            picard = Geno.conf.get( 'SOFTWARE', 'picard' )
                                        )
                                    )
            shell_script_file.close()

        #
        # Run
        #
        return_code = Geno.RT.runtask(
                            Geno.job.get( 'job_queue' )[ function_name ],
                            Geno.job.get( 'memory' )[ function_name ],
                            shell_script_full_path )

        if return_code != 0:
            log.error( "{function}: runtask failed" )
            raise


    except IOError as (errno, strerror):
        log.error( "{function}: I/O error({num}): {error}".format( function = whoami(), num = errno, error = strerror) )
        return_code = False

    except ValueError:
        log.error( "{function}: ValueError".format( function = whoami() ) )
        return_code = False

    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        log.error( "{function}: Unexpected error: {error}".format( function = whoami(), error = sys.exc_info()[0] ) )
        log.error("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )


    return True
#
# Stage 7: fisher_mutation_call
#
def fisher_mutation_call(
    input_file1,
    input_file2,
    output_file1,
    output_file2
    ):
    """
       Mutaion calling

    """

    try:
        function_name = whoami()
        log.info( "#{function}".format( function = function_name ) )

        #
        # Make shell script
        #
        shell_script_full_path = make_script_file_name( function_name, Geno )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( wgs_res.fisher_mutation_call.format(
                                        log = Geno.dir[ 'log' ],
                                        ref_fa = Geno.conf.get( 'REFERENCE', 'ref_fasta' ),
                                        control_input_bam = input_file1,
                                        disease_input_bam = '',
                                        output_txt = output_file1,
                                        max_indel = Geno.job.get( 'max_indel' ),
                                        max_distance = Geno.job.get( 'max_distance' ),
                                        base_quality = Geno.job.get( 'base_quality' ),
                                        map_quality = Geno.job.get( 'map_quality' ),
                                        mismatch_rate = Geno.job.get( 'mismatch_rate' ),
                                        min_depth = Geno.job.get( 'min_depth' ),
                                        script_dir = Geno.dir[ 'script' ]
                                    )
                                )
        shell_script_file.close()

        #
        # Run
        #
        return_code = Geno.RT.runtask(
                            Geno.job.get( 'job_queue' )[ function_name ],
                            Geno.job.get( 'memory' )[ function_name ],
                            shell_script_full_path )
        if return_code != 0:
            log.error( "{function}: runtask failed" )
            raise


    except IOError as (errno, strerror):
        log.error( "{function}: I/O error({num}): {error}".format( function = whoami(), num = errno, error = strerror) )
        return_code = False

    except ValueError:
        log.error( "{function}: ValueError".format( function = whoami() ) )
        return_code = False

    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        log.error( "{function}: Unexpected error: {error}".format( function = whoami(), error = sys.exc_info()[0] ) )
        log.error("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )


    return True

#####################################################################
#
#   STAGE 0 data preparation
#
Sample = Sample()

#####################################################################
#
#   STAGE 1 bam2fastq
#   in:     bam
#   out:    fastq
#
@active_if( 'bam2fastq' in Geno.job.get( 'tasks' )[ 'WGS' ] )
@parallel( generate_params_for_bam2fastq )
@check_if_uptodate( check_file_exists_for_input_output )
def stage_1( input_file, output_file ):
    bam2fastq( input_file, output_file )

#####################################################################
#
#   STAGE 2 split_fastq
#   in:     fastq
#   out:    fastq * X
#
@follows( stage_1 )
@active_if( 'split_fastq' in Geno.job.get( 'tasks' )[ 'WGS' ] )
@files( generate_params_for_split_fastq )
@check_if_uptodate( check_file_exists_for_split_fastq )
def stage_2( input_file1, input_file2, output_file1, output_file2 ):
    split_fastq( input_file1, input_file2, output_file1, output_file2 )

#####################################################################
#
#   STAGE 3 cutadapt
#   in:     fastq
#   out:    fastq
#
@follows( stage_2 )
@active_if( 'cutadapt' in Geno.job.get( 'tasks' )[ 'WGS' ] )
@files( generate_params_for_cutadapt )
@check_if_uptodate( check_file_exists_for_input_output )
def stage_3( input_file1, input_file2, output_file1, output_file2 ):
    cutadapt( input_file1, input_file2, output_file1, output_file2 )

#####################################################################
#
#   STAGE 4 bwa_mem
#
#   in:     fastq1, fastq2
#   out:    bam
#
@follows( stage_3 )
@active_if( 'bwa_mem' in Geno.job.get( 'tasks' )[ 'WGS' ] )
@files( generate_params_for_bwa_mem )
@check_if_uptodate( check_file_exists_for_bwa_mem )
def stage_4(  input_file1, input_file2, output_file1, output_file2 ):
    if Geno.job.get( 'use_biobambam' ):
        bwa_mem( input_file1, input_file2, output_file1, output_file2, True )
    else:
        bwa_mem( input_file1, input_file2, output_file1, output_file2, False )


#####################################################################
#
#   STAGE 5 merge
#
#   in:     bam x {number}
#   out:    bam
#
@follows( stage_4 )
@active_if( 'merge_bam' in Geno.job.get( 'tasks' )[ 'WGS' ] )
@files( generate_params_for_merge_bam )
@check_if_uptodate( check_file_exists_for_merge_bam )
def stage_5( input_file1, input_file2, output_file1, output_file2 ):
    if Geno.job.get( 'use_biobambam' ):
        if 'markdulicates' in Geno.job.get( 'tasks' )[ 'WGS']:
            return True
        else:
            merge_bam( input_file1, input_file2, output_file1, output_file2, True )
    else:
        merge_bam( input_file1, input_file2, output_file1, output_file2, False )


#####################################################################
#
#   STAGE 6 markduplicates
#
#   in:     bam x {number}
#   out:    bam
#
@follows( stage_5 )
@active_if( 'markduplicates' in Geno.job.get( 'tasks' )[ 'WGS' ] )
@files( generate_params_for_markduplicates )
@check_if_uptodate( check_file_exists_for_markduplicates )
def stage_6( input_file1, input_file2, output_file1, output_file2 ):
    if Geno.job.get( 'use_biobambam' ):
        markduplicates( input_file1, input_file2, output_file1, output_file2, True )
    else:
        markduplicates( input_file1, input_file2, output_file1, output_file2, False )

#####################################################################
#
#   STAGE 7 fisher_mutation call
#
#   in:     bam
#   out:    vcf
#
@follows( stage_6 )
@active_if ( 'fisher_mutation_call' in Geno.job.get( 'tasks' )[ 'WGS' ] )
@files( generate_params_for_fisher_mutation_call )
@check_if_uptodate( check_file_exists_for_input_output )
def stage_7(  input_file1, input_file2, output_file1, output_file2 ):
    fisher_mutation_call(  input_file1, input_file2, output_file1, output_file2 )


#####################################################################
#
#   LAST STAGE 
#
@follows( stage_7 )
def last_function():
    log.info( "Genomon pipline has finished successflly!" )
    return True

