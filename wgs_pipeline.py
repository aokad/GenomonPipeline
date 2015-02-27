"""
wgs_pipeline.py

"""

import sys
import os
from datetime import datetime
from ruffus import *
from runtask import RunTask
from glob import glob

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
def check_file_exists_for_bwa_mem(
        input_file1,
        input_file2,
        output_file
        ):

    """
    Checks if output file exists for bwa_mem

    """

    if not os.path.exists( output_file ):
        return True, "Missing file {output} for {input}.".format( output = output_file, input = input_file1 )
    else:
        in_time = os.path.getmtime( input_file1 )
        out_time = os.path.getmtime( output_file )
        if in_time > out_time:
            return True, "{output} is older than {input}.".format( output = output_file, input = input_file1 )
        else:
            return False, "File {output} exits for {input}.".format( output = output_file, input = input_file1 )

def check_file_exists_for_split_fastq_files(
    input_file,
    output_prefix,
    output_suffix,
    ):

    """
    Checks if output file exists for split_fastq_files
        glob output files with prefix and suffix

    """

    if not glob( "{prefix}*{suffix}".format( prefix = output_prefix, suffix = output_suffix ) ):
        return True, "Missing file {outprefix}*{outsuffix} for {input}.".format(
                            outprefix = output_prefix,
                            outsuffix = output_suffix,
                            input = input_file )
    else:
        return False, "File {outprefix}*{outsuffix} exists for {input}.".format(
                            outprefix = output_prefix,
                            outsuffix = output_suffix,
                            input = input_file )

def check_file_exists_for_merge_bam(
    input_prefix,
    input_suffix,
    output_file
    ):

    """
    Checks if output file exists for split_fastq_files

    """

    if not os.path.exists( output_file ):
        return True, "Missing file {outputfile} for {input_prefix}*{input_suffix}.".format(
                            outputfile = output_file,
                            input_prefix = input_prefix,
                            input_suffix = input_suffix )
    else:
        return False, "File {outputfile} exists for {input_prefix}*{input_suffix}.".format(
                            outputfile = output_file,
                            input_prefix = input_prefix,
                            input_suffix = input_suffix )
def check_file_exists_for_mutation_call(
    input_file,
    output_file
    ):

    """
    Checks if output file exists for split_fastq_files

    """

    if not os.path.exists( output_file ):
        return True, "Missing file {outputfile} for {input_file}.".format(
                            outputfile = output_file,
                            inputfile = input_file )
    else:
        return False, "File {outputfile} exists for {inputfile}.".format(
                            outputfile = output_file,
                            inputfile = input_file )


#####################################################################
#
# Data preparation subrotines for parallelization
#

#
# Get starting files
#
def get_starting_files():
    """
    Get the list of starting files from
    specified job configuration file.

    1) Get file name from job configuration files
    2) Globl files
    3) Create list for each pair of fastq files

    """

    file_type   = Geno.job.get( 'input_file_type' )
    bed_file    = Geno.job.get( 'bed_file' )

    file_list = []

    try:
        #
        # A) Paired-end fastq files
        #
        if file_type == 'paired':
            fastq_file = []
            for pair_id in Geno.job.get( 'pair_id' ):
                fastq_file.append( Geno.job.get( 'fastq_file' ).format( pair_id = pair_id ) )

            if ( fastq_file[ 0 ] != None and
                 fastq_file[ 1 ] != None  ):

                glob_file_list = []
                for pair_id in [ 0, 1 ]:
                    glob_file_list.append( sorted( glob( Geno.data + '/' + fastq_file[ pair_id ] ) ) )

                for i in range( len( glob_file_list[ 0 ] ) ):
                    file_list.append( [
                                    glob_file_list[ 0 ][ i ],
                                    glob_file_list[ 1 ][ i ]
                                    ] )

        #
        # B) Single-end fastq files
        #
        elif file_type == 'single':
            fastq_file = Geno.job.get( 'fastq_file' )
            fastq_dir  = Geno.job.get( 'fastq_dir' )

            if fastq_file != None:
                file_list = sorted( glob( Geno.data + fastq_file ) )

        #
        # C) bam files
        #
        else:
            bam_files = Geno.job.get( 'bam_file' )
            file_list = sorted( glob( Geno.data + bam_files ) )

    except:
        log.error( "{function} failed.".format( function = whoami() ) )
        log.error( "Unexpected error: {1}".format( sys.exc_info()[0] ) )

    return file_list

#
# For STAGE 1 split fastq
#
def generate_parameters_for_split_fastq_files( starting_files ):
    """
    Generate parameter list for STAGE 1 spilt_fastq_files

    Create a list of the following
        input_file
        output_prefix
        output_suffix

    """
    input_fastq_files = []

    if 'paired' == Geno.job.get( 'input_file_type' ):
        for file1, file2 in starting_files:
                for input_file in ( file1, file2 ):
                    tmp_name = os.path.splitext( os.path.basename( input_file ) )
                    basename = tmp_name[ 0 ]
                    output_prefix = "{dir}/out/fastq/{file}".format(
                                            dir = Geno.results,
                                            file = basename 
                                            )
                    output_suffix = tmp_name[ 1 ]
                    input_fastq_files.append( [ input_file, output_prefix, output_suffix ] )
    else:
        for file_name in starting_files:
            tmp_name = os.path.splitext( os.path.basename( file_name ) )
            basename = tmp_name[ 0 ]
            output_prefix = "{dir}/out/fastq/{file}".format(
                                    dir = Geno.results,
                                    file = basename
                                    )
            output_suffix = tmp_name[ 1 ]
            input_fastq_files.append( [ file_name, output_prefix, output_suffix ] )


    return input_fastq_files


#
# For STAGE 2 fastq to bam
#

def generate_parameters_for_bwa_mem_single():
    """
    Generate parameter list for STAGE 2 bwa_mem

    input_file1
    input_file2
    output_file

    """
    global bwa_mem_parameters
    bwa_mem_parameters = []

    for parameter in split_fastq_files_parameters:
        split_files = sorted( glob( "{prefix}*{suffix}".format(
                                        prefix = parameter[ 1 ],
                                        suffix = parameter[ 2 ] 
                        ) ) )
        for split_file in split_files:
            bam_file = "{dir}/out/bam/{file}.bam".format(
                                dir = Geno.results,
                                file = os.path.splitext( os.path.basename( split_files ))[ 0 ]
                )
            out_param = [ split_file,
                          '',
                          bam_file ]
            bwa_mem_parameters.append( out_param )

    return bwa_mem_parameters

def generate_parameters_for_bwa_mem_pair():
    """
    Generate parameter list for STAGE 2 bwa_mem

    input_file1
    input_file2
    output_file
    job_type

    """
    global bwa_mem_parameters
    bwa_mem_parameters = []

    for id1 in range( 0, len( split_fastq_files_parameters ), 2 ):
        split_files1 = sorted( glob( "{filename}*".format( filename = split_fastq_files_parameters[ id1 ][ 1 ] ) ) )
        split_files2 = sorted( glob( "{filename}*".format( filename = split_fastq_files_parameters[ id1 + 1 ][ 1 ] ) ) )

        for file1, file2 in zip( split_files1, split_files2 ):
            bam_file = "{dir}/out/bam/{file}.bam".format(
                                dir = Geno.results,
                                file = os.path.splitext( os.path.basename( file1 ))[ 0 ]
                )

            out_param = [ file1,
                          file2,
                          bam_file ]
            bwa_mem_parameters.append( out_param )
    return bwa_mem_parameters

def generate_parameters_for_bwa_mem():
    """
    Generate parameter list for STAGE 2 bwa_mem

    """
    seq_type = Geno.job.get( 'input_file_type' )

    if 'paired' == seq_type:
        parameter_list = generate_parameters_for_bwa_mem_pair()
    elif 'single' == seq_type:
        parameter_list = generate_parameters_for_bwa_mem_single()


    for parameter in parameter_list:
        yield parameter

#
# For STAGE 3 merge_bam        
#        
def generate_parameters_for_merge_bam():
    """
    Generate parameter list for STAGE 2 bwa_mem

    """
    global starting_files
    global merge_bam_parameters
    merge_bam_parameters = []
    for parameter in starting_files:
        base, ext = os.path.splitext( os.path.basename( parameter[ 0 ] ) )
        basename = "{dir}/out/bam/{basename}".format(
                            dir = Geno.results,
                            basename = base )
        output_bam = basename + '.bam'
        out_param = [ basename,
                      '.bam',
                      output_bam ]
        merge_bam_parameters.append( out_param )
        yield out_param

#
# For STAGE 4 mutation_call
#
def generate_parameters_for_mutation_call():
    pass

#####################################################################
#
#   Main functions
#
def split_fastq_files(
    input_file,
    output_prefix,
    output_suffix
    ):
    """
    Split fastq file into small pieces for parallele processing

    """
    try:
        log.info( "# {function}\n".format( function = whoami() ) + 
                     "in:  {input}\n".format( input=input_file) )

        #
        # Make sure files exist.
        #
        if not os.path.isfile( input_file ):
            log.error( "file: {file} does not exist.".format( file=input_file ) )
            return 0

        #
        # Make shell script
        #
        now = datetime.now()
        suffix_len = len( output_suffix )
        shell_script_name = res.shell_script_format.format(
                                        name='split_fastq',
                                        year=now.year,
                                        month=now.month,
                                        day=now.day,
                                        hour=now.hour,
                                        min=now.minute,
                                        msecond=now.microsecond )
        shell_script_full_path = "{results}/script/{file}.sh".format(
                                        results = Geno.results,
                                        file = shell_script_name )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( res.splitfile.format(
                                        results = Geno.results,
                                        lines_per_file = Geno.job.get( 'split_fastq_line_number' ),
                                        input_file = input_file,
                                        suffix_len = suffix_len,
                                        output_suffix = output_suffix,
                                        output_prefix = output_prefix ) )
        shell_script_file.close()

        #
        # Run
        #
        Geno.RT.runtask( Geno.job.get( 'job_type' )[ 'split_fastq' ],
                         Geno.job.get( 'memory' )[ 'split_fastq' ],
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

def bwa_mem(
    input_file1,            # 1st parameter is Input
    input_file2,            # 2nd parameter is Input
    output_file             # 3rd parameter is Output
    ):
    """
       Align sequence reads in FASTQ file.
       Convert the result SAM file to BAM file.
    """

    try:
        log.info( "# {function}\n".format( function = whoami() ) +
                     "in1: {input}\n".format( input=input_file1) +
                     "in2: {input}\n".format( input=input_file2) +
                     "out: {output}\n".format( output=output_file) )

        #
        # Make sure files exist.
        #
        if ( not os.path.isfile( input_file1 ) or
             not os.path.isfile( input_file2 ) ) :
            log.error( "file: {file} does not exist.".format( file=input_file1 ) )
            log.error( "file: {file} does not exist.".format( file=input_file2 ) )
            return 0

        #
        # Make shell script
        #
        now = datetime.now()
        shell_script_name = res.shell_script_format.format(
                                        name='fastq2bam',
                                        year=now.year,
                                        month=now.month,
                                        day=now.day,
                                        hour=now.hour,
                                        min=now.minute,
                                        msecond=now.microsecond )
        shell_script_full_path = "{results}/script/{file}.sh".format(
                                        results = Geno.results,
                                        file = shell_script_name )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( res.fastq2bam.format(
                                        results = Geno.results,
                                        bwa = Geno.conf.get( 'SOFTWARE', 'bwa' ),
                                        hg19_fa = Geno.conf.get( 'REFERENCE', 'hg19_fasta' ),
                                        fastq1 = input_file1,
                                        fastq2 = input_file2,
                                        samtools = Geno.conf.get( 'SOFTWARE', 'samtools' ),
                                        bam = output_file ) )
        shell_script_file.close()

        #
        # Run
        #
        Geno.RT.runtask( Geno.job.get( 'job_type' )[ 'split_fastq' ],
                         Geno.job.get( 'memory' )[ 'split_fastq' ],
                         shell_script_full_path )


    except IOError as (errno, strerror):
        log.error( "{function}: I/O error({num}): {error}".format(function = whoami(), num = errno, error = strerror) )
        return_code = False

    except ValueError:
        log.error( "{function}: ValueError".format( function = whoami() ) )
        return_code = False

    except:
        log.error( "{function}: Unexpected error: {error}".format( function = whoami(), error = sys.exc_info()[0] ) )


    return True

def merge_bam(
    input_prefix,
    input_suffix,
    output_file
    ):
    """
       Merge split bam files

    """

    try:
        log.info( "# {function}\n".format( function = whoami() ) +
                     "in: {prefix}*{suffix}\n".format( prefix=input_prefix, suffix = input_suffix ) +
                     "out: {output}\n".format( output=output_file) )

        #
        # Make shell script
        #
        now = datetime.now()
        shell_script_name = res.shell_script_format.format(
                                        name='merge_bam',
                                        year=now.year,
                                        month=now.month,
                                        day=now.day,
                                        hour=now.hour,
                                        min=now.minute,
                                        msecond=now.microsecond )
        shell_script_full_path = "{results}/script/{file}.sh".format(
                                        results = Geno.results,
                                        file = shell_script_name )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( res.merge_bam.format(
                                        results = Geno.results,
                                        input_bam_files = "{prefix}*{suffix}".format(
                                                        prefix = input_prefix,
                                                        suffix = input_suffix
                                                        ),
                                        output_bam_file = output_file,
                                        samtools = Geno.conf.get( 'SOFTWARE', 'samtools' )))
        shell_script_file.close()

        #
        # Run
        #
        Geno.RT.runtask( Geno.job.get( 'job_type' )[ 'split_fastq' ],
                         Geno.job.get( 'memory' )[ 'split_fastq' ],
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

def mutation_call(
    input_file,
    output_file
    ):
    """
       Mutaion calling

    """

    try:
        log.info( "# {function}\n".format( function = whoami() ) +
                     "in: {prefix}*{suffix}\n".format( prefix=input_prefix, suffix = input_suffix ) +
                     "out: {output}\n".format( output=output_file) )

        #
        # Make shell script
        #
        now = datetime.now()
        shell_script_name = res.shell_script_format.format(
                                        name='mutation_call',
                                        year=now.year,
                                        month=now.month,
                                        day=now.day,
                                        hour=now.hour,
                                        min=now.minute,
                                        msecond=now.microsecond )
        shell_script_full_path = "{results}/script/{file}.sh".format(
                                        results = Geno.results,
                                        file = shell_script_name )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( res.mutation_call.format(
                                        results = Geno.results,
                                        input_bam_files = input_file,
                                        samtools = Geno.conf.get( 'SOFTWARE', 'samtools' )))
        shell_script_file.close()

        #
        # Run
        #
        Geno.RT.runtask( Geno.job.get( 'job_type' )[ 'split_fastq' ],
                         Geno.job.get( 'memory' )[ 'split_fastq' ],
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
#   Global variables
#
global starting_files
global split_fastq_files_parameters
global bwa_mem_parameters
global merge_bam_parameters

#####################################################################
#
#   STAGE 0 data preparation
#
starting_files = get_starting_files()
split_fastq_files_parameters = generate_parameters_for_split_fastq_files( starting_files )

#####################################################################
#
#   STAGE 1 split fastq
#
@parallel( split_fastq_files_parameters )
@check_if_uptodate( check_file_exists_for_split_fastq_files )
def stage_1( input_file, output_prefix, output_suffix ):
    split_fastq_files( input_file, output_prefix, output_suffix )

#####################################################################
#
#   STAGE 2 fastq to bam
#
@files( generate_parameters_for_bwa_mem )
@check_if_uptodate( check_file_exists_for_bwa_mem )
@follows( stage_1 )
def stage_2( input_file1, input_file2, output_file ):
    bwa_mem( input_file1, input_file2, output_file )

#####################################################################
#
#   STAGE 3 merge
#
@files( generate_parameters_for_merge_bam )
@check_if_uptodate( check_file_exists_for_merge_bam )
@follows( stage_2 )
def stage_3( input_prefix, input_suffix, output_file ):
    merge_bam( input_prefix, input_suffix, output_file )

#####################################################################
#
#   STAGE 4 mutation call
#
# @files( generate_parameters_for_mutation_call )
# @check_if_uptodate( check_file_exists_for_mutation_call )
# @follows( stage_3 )
# def stage_4( input_file, output_file ):
#     mutation_call( input_file, output_file )


#####################################################################
#
#   LAST STAGE 
#
@follows( stage_3 )
def last_function():
    log.info( "Genomon pipline has finished successflly!" )
    return True
