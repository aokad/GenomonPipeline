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
    log.debug( "# {function}\n".format( function = whoami() ) )
    log.debug( "in1:   {input}\n".format( input=input_file1) )
    log.debug( "in2:   {input}\n".format( input=input_file2) )
    log.debug( "out1:  {output}\n".format( output=output_file1) )
    log.debug( "out2:  {output}\n".format( output=output_file2) )


    ( output_prefix, output_suffix ) = os.path.splitext( output_file1 )
    outfile_list = glob( "{prefix}*{suffix}".format( prefix = output_prefix, suffix = output_suffix ) )
    if not outfile_list:
        log.debug( "Missing file {outprefix}*{outsuffix} for {input}.".format(
                            outprefix = output_prefix,
                            outsuffix = output_suffix,
                            input = input_file1 ) )
        return True, "Missing file {outprefix}*{outsuffix} for {input}.".format(
                            outprefix = output_prefix,
                            outsuffix = output_suffix,
                            input = input_file1 )
    else:
        in_time = os.path.getmtime( input_file1 )
        out_time = os.path.getmtime( outfile_list[ 0 ] )
        if in_time > out_time:
            log.debug( "{outprefix}*{outsuffix} is older than {input}.".format(
                                outprefix = output_prefix,
                                outsuffix = output_suffix,
                                input = input_file1 ) )
            return True, "{outprefix}*{outsuffix} is older than {input}.".format(
                                outprefix = output_prefix,
                                outsuffix = output_suffix,
                                input = input_file1 )
        else:
            log.debug( "File {outprefix}*{outsuffix} exits for {input}.".format(
                                outprefix = output_prefix,
                                outsuffix = output_suffix,
                                input = input_file1 ) )
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
    Checks if output file exists for mutation_call

    """
    log.debug( "# {function}\n".format( function = whoami() ) )
    log.debug( "in1:   {input}\n".format( input=input_file1) )
    log.debug( "in2:   {input}\n".format( input=input_file2) )
    log.debug( "out1:  {output}\n".format( output=output_file1) )
    log.debug( "out2:  {output}\n".format( output=output_file2) )

    ( output_prefix, output_suffix ) = os.path.splitext( output_file1 )
    if not glob( "{prefix}*{suffix}".format( prefix = output_prefix, suffix = output_suffix ) ):
        log.debug("Missing file {outputfile} for {inputfile}.".format(
                            outputfile = output_file1,
                            inputfile = input_file1 ) )
        return True, "Missing file {outputfile} for {inputfile}.".format(
                            outputfile = output_file1,
                            inputfile = input_file1 )
    else:
        log.debug( "File {output} exits for {input}.".format( output = output_file1, input = input_file1 ) )
        return False, "File {output} exits for {input}.".format( output = output_file1, input = input_file1 )

def check_file_exists_for_input_output(
        input_file1,
        input_file2,
        output_file1,
        output_file2
    ):

    """
    Checks if output file exists for mutation_call

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

    Sample.make_param( 'bam2fastq', '.fastq', 1, 2 )
    for param in Sample.current():
        yield param


#
# For STAGE 2 split_fastq
#
def generate_params_for_split_fastq():
    """
    Generate parameter list for spilt_fastq_files

    """

    
    Sample.make_param( 'split_fastq', '.fastq', 2, 2 )
    for param in Sample.current():
        yield param

#
# For STAGE 3 cutadapt
#
def generate_params_for_cutadapt():
    """
    Generate parameter list for cutadapt

    """
    
    Sample.make_param( 'cutadapt', '_cutadapt.fastq', 2, 2 )
    for param in Sample.current():
        yield param

#
# For STAGE 4 bwa_mem
#

    
def generate_params_for_bwa_mem():
    """
    Generate parameter list for bwa_mem

    """

    Sample.make_param( 'bwa_mem', '.bam', 2, 1 ) 
    for parameter in Sample.current():
        yield parameter

#
# For STAGE 5 merge_bam        
#        
def generate_params_for_merge_bam():
    """
    Generate parameter list for merge_bam

    """
    
    Sample.make_param( 'merge_bam', '.bam', 1, 1 )
    for parameter in Sample.current():
        yield parameter
#
# For STAGE 6 mutation_call
#
def generate_params_for_mutation_call():

    
    Sample.make_param( 'mutation_call', '.txt', 1, 1 )
    for parameter in Sample.current():
        yield parameter

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
        log.info( "# {function}\n".format( function = whoami() ) )
        log.info( "in1:   {input}\n".format( input=input_file1) )
        log.info( "in2:   {input}\n".format( input=input_file2) )
        log.info( "out1:  {output}\n".format( output=output_file1) )
        log.info( "out2:  {output}\n".format( output=output_file2) )

        #
        # Make sure files exist.
        #
        if not os.path.isfile( input_file ):
            log.error( "file: {file} does not exist.".format( file=input_file ) )
            return 0

        #
        # Make shell script
        #
        function_name = whoami()
        shell_script_full_path = make_script_file_name( function_name )
        shell_script_file = open( shell_script_full_path, 'w' )

        file_type = Geno.job.get( 'input_file_type' )
        if 'paired_bam' == file_type:
            shell_script_file.write( res.bamtofastq_p.format(
                                            log = Geno.dir[ 'log' ],
                                            bamfile = input_file1,
                                            outfastq1 = output_file1,
                                            outfastq2 = output_file2,
                                            tmpfastq = output_file1 + '.tmp',
                                            bamtofastq = Geno.conf.get( 'SOFTWARE', 'bamtofastq' )
                                            ) )
        elif 'single_bam' == file_type:
            shell_script_file.write( res.bamtofastq_s.format(
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
        Geno.RT.runtask( Geno.job.get( 'job_type' )[ function_name ],
                         Geno.job.get( 'memory' )[ function_name ],
                         shell_script_full_path )

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

    """
    try:
        log.info( "# {function}\n".format( function = whoami() ) )
        log.info( "in1:   {input}\n".format( input=input_file1) )
        log.info( "in2:   {input}\n".format( input=input_file2) )
        log.info( "out1:  {output}\n".format( output=output_file1) )
        log.info( "out2:  {output}\n".format( output=output_file2) )

        #
        # Make sure files exist.
        #
        if not os.path.isfile( input_file1 ):
            log.error( "file: {file} does not exist.".format( file=input_file1 ) )
            return 0

        #
        # Make data for array job 
        #
        output_prefix1 = make_sample_file_name( output_file1, "{dir}/{base}_" )
        output_prefix2 = make_sample_file_name( output_file2, "{dir}/{base}_" )
        array_data1 = "IN_FILE=(\n[1]=\"{file1}\"\n[2]=\"{file2}\"\n)\n".format(
                        file1 = input_file1,
                        file2 = input_file2 )
        array_data2 = "OUT_FILE=(\n[1]=\"{file1}\"\n[2]=\"{file2}\"\n)\n".format(
                        file1 = output_prefix1,
                        file2 = output_prefix2 )
        input_file = '${IN_FILE[$SGE_TASK_ID]}'
        output_prefix = '${OUT_FILE[$SGE_TASK_ID]}'

        output_suffix = '.fastq'
        suffix_len = len( output_suffix )
        function_name = whoami()

        #
        # Make shell script for array job
        #
        shell_script_full_path = make_script_file_name( function_name )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( res.splitfile.format(
                                        log = Geno.dir[ 'log' ],
                                        array_data = array_data1 + array_data2,
                                        lines_per_file = Geno.job.get( 'split_fastq_line_number' ),
                                        input_file = input_file,
                                        suffix_len = suffix_len,
                                        output_suffix = output_suffix,
                                        output_prefix = output_prefix ) )
        shell_script_file.close()

        #
        # Run
        #
        Geno.RT.run_arrayjob( Geno.job.get( 'job_type' )[ function_name ],
                              Geno.job.get( 'memory' )[ function_name ],
                              shell_script_full_path,
                              '1-2:1' )


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
# Stage 2: cutadapt
#
def cutadapt(
    input_file1,
    input_file2,
    output_file1,
    output_file2
    ):
    """
       Apply cutadapt to remove adaptor sequences from reads
    """

    try:
        log.info( "# {function}\n".format( function = whoami() ) )
        log.info( "in1:   {input}\n".format( input=input_file1) )
        log.info( "in2:   {input}\n".format( input=input_file2) )
        log.info( "out1:  {output}\n".format( output=output_file1) )
        log.info( "out2:  {output}\n".format( output=output_file2) )

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
        function_name = whoami()

        #
        # Make shell script
        #
        shell_script_full_path = make_script_file_name( function_name )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( res.cutadapt.format(
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
        Geno.RT.run_arrayjob( Geno.job.get( 'job_type' )[ function_name ],
                              Geno.job.get( 'memory' )[ function_name ],
                              shell_script_full_path,
                              '1-2:1' )


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
# Stage 2: bwa_mem
#
def bwa_mem(
    input_file1,
    input_file2,
    output_file1,
    output_file2
    ):
    """
       Align sequence reads in FASTQ file.
       Convert the result SAM file to BAM file.
    """

    try:
        log.info( "# {function}\n".format( function = whoami() ) )
        log.info( "in1:   {input}\n".format( input=input_file1) )
        log.info( "in2:   {input}\n".format( input=input_file2) )
        log.info( "out1:  {output}\n".format( output=output_file1) )
        log.info( "out2:  {output}\n".format( output=output_file2) )

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
            bam = "{dir}/{file}.bam".format( dir = Geno.dir[ 'bam' ],
                                             file = os.path.splitext( os.path.basename( fastq1 ) )[ 0 ] )
            id += 1
            str1 += " [{id}]=\"{file}\"\n".format( id = id, file = fastq1 )
            str2 += " [{id}]=\"{file}\"\n".format( id = id, file = fastq2 )
            str3 += " [{id}]=\"{file}\"\n".format( id = id, file = bam )
        str1 += " )\n"
        str2 += " )\n"
        str3 += " )\n"

        
        #
        # Make shell script for array job
        #
        function_name = whoami()
        shell_script_full_path = make_script_file_name( function_name )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( res.bwa_mem.format(
                                        log = Geno.dir[ 'log' ],
                                        bwa = Geno.conf.get( 'SOFTWARE', 'bwa' ),
                                        hg19_fa = Geno.conf.get( 'REFERENCE', 'hg19_fasta' ),
                                        array_data = str1 + str2 + str3,
                                        fastq1 = "${FILE1[$SGE_TASK_ID]}",
                                        fastq2 = "${FILE2[$SGE_TASK_ID]}",
                                        samtools = Geno.conf.get( 'SOFTWARE', 'samtools' ),
                                        bam = "${FILE3[$SGE_TASK_ID]}" ) )
        shell_script_file.close()

        #
        # Run
        #
        Geno.RT.run_arrayjob( Geno.job.get( 'job_type' )[ function_name ],
                              Geno.job.get( 'memory' )[ function_name ],
                              shell_script_full_path,
                              "1-{id}:1".format( id = id ) )


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
# Stage 3: merge_bam
#
def merge_bam(
    input_file1,
    input_file2,
    output_file1,
    output_file2
    ):
    """
       Merge split bam files

    """

    try:
        log.info( "# {function}\n".format( function = whoami() ) )
        log.info( "in1:   {input}\n".format( input=input_file1) )
        log.info( "in2:   {input}\n".format( input=input_file2) )
        log.info( "out1:  {output}\n".format( output=output_file1) )
        log.info( "out2:  {output}\n".format( output=output_file2) )

        #
        # Make data for array job 
        #
        ( input_prefix, input_suffix ) = os.path.splitext( input_file1 )
        input_files = "{input_prefix}*{input_suffix}".format(
                                    input_prefix = input_prefix,
                                    input_suffix = input_suffix )

        #
        # Make shell script
        #
        function_name = whoami()

        shell_script_full_path = make_script_file_name( function_name )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( res.merge_bam.format(
                                        log = Geno.dir[ 'log' ],
                                        input_bam_files = input_files,
                                        output_bam_file = output_file1,
                                        samtools = Geno.conf.get( 'SOFTWARE', 'samtools' ) ) )
        shell_script_file.close()

        #
        # Run
        #
        Geno.RT.runtask( Geno.job.get( 'job_type' )[ 'merge_bam' ],
                         Geno.job.get( 'memory' )[ 'merge_bam' ],
                         shell_script_full_path )


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
# Stage 4: mutation_call
#
def mutation_call(
    input_file1,
    input_file2,
    output_file1,
    output_file2
    ):
    """
       Mutaion calling

    """

    try:
        log.info( "# {function}\n".format( function = whoami() ) )
        log.info( "in1:   {input}\n".format( input=input_file1) )
        log.info( "in2:   {input}\n".format( input=input_file2) )
        log.info( "out1:  {output}\n".format( output=output_file1) )
        log.info( "out2:  {output}\n".format( output=output_file2) )

        #
        # Make shell script
        #
        function_name = whoami()
        shell_script_full_path = make_script_file_name( function_name )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( res.mutation_call.format(
                                        log = Geno.dir[ 'log' ],
                                        input_bam = input_file1,
                                        outputpath = Geno.dir[ 'fisher' ],
                                        num = 1,
                                        region = 1,
                                        type = '',
                                        map_qual_thres = Geno.job.get( 'map_cutoff' ),
                                        base_qual_thres = Geno.job.get( 'base_quality_cutoff' ),
                                        max_indel = Geno.job.get( 'max_indel' ),
                                        genref = Geno.conf.get( 'REFERENCE', 'hg19_fastq' ),
                                        samtools  =  Geno.conf.get( 'SOFTWARE', 'samtools' ),
                                        script_dir = Geno.dir[ 'script' ],
                                        R = Geno.conf.get( 'SOFTWARE', 'R' )
                                    )
                                )
        shell_script_file.close()

        #
        # Run
        #
        Geno.RT.runtask( Geno.job.get( 'job_type' )[ function_name ],
                         Geno.job.get( 'memory' )[ function_name ],
                         shell_script_full_path )


    except IOError as (errno, strerror):
        log.error( "{function}: I/O error({num}): {error}".format( function = whoami(), num = errno, error = strerror) )
        return_code = False

    except ValueError:
        log.error( "{function}: ValueError".format( function = whoami() ) )
        return_code = False

    except:
        log.error( "{function}: Unexpected error: {error}".format( function = whoami(), error = sys.exc_info()[0] ) )


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
    bwa_mem( input_file1, input_file2, output_file1, output_file2 )

#####################################################################
#
#   STAGE 5 merge
#
#   in:     bam * X
#   out:    bam
#
@follows( stage_4 )
@active_if( 'merge_bam' in Geno.job.get( 'tasks' )[ 'WGS' ] )
@files( generate_params_for_merge_bam )
@check_if_uptodate( check_file_exists_for_merge_bam )
def stage_5( input_file1, input_file2, output_file1, output_file2 ):
    merge_bam( input_file1, input_file2, output_file1, output_file2 )

#####################################################################
#
#   STAGE 6 mutation call
#
#   in:     bam
#   out:    vcf
#
@follows( stage_5 )
@active_if ( 'mutation_call' in Geno.job.get( 'tasks' )[ 'WGS' ] )
@files( generate_params_for_mutation_call )
@check_if_uptodate( check_file_exists_for_input_output )
def stage_6(  input_file1, input_file2, output_file1, output_file2 ):
    mutation_call(  input_file1, input_file2, output_file1, output_file2 )


#####################################################################
#
#   LAST STAGE 
#
@follows( stage_6 )
def last_function():
    log.info( "Genomon pipline has finished successflly!" )
    return True
