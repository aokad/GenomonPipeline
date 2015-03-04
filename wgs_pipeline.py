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
def check_file_exists_for_split_fastq(
    input_file,
    output_prefix,
    output_suffix,
    ):

    """
    Checks if output file exists for split_fastq
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

def check_file_exists_for_cutadapt(
        input_file,
        output_file
        ):

    """
    Checks if output file exists for cutadapt

    """

    if not os.path.exists( output_file ):
        return True, "Missing file {output} for {input}.".format( output = output_file, input = input_file )
    else:
        in_time = os.path.getmtime( input_file )
        out_time = os.path.getmtime( output_file )
        if in_time > out_time:
            return True, "{output} is older than {input}.".format( output = output_file, input = input_file )
        else:
            return False, "File {output} exits for {input}.".format( output = output_file, input = input_file )

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

def check_file_exists_for_merge_bam(
    input_prefix,
    input_suffix,
    output_file
    ):

    """
    Checks if output file exists for merge_bam

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
    Checks if output file exists for mutation_call

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
# For STAGE 1 split_fastq
#
def generate_parameters_for_split_fastq( starting_files ):
    """
    Generate parameter list for STAGE 1 spilt_fastq_files

    Create a list of the following
        input_file
        output_prefix
        output_suffix

    """
    input_fastq_files = []
    input_file_type = Geno.job.get( 'input_file_type' )
    if 'paired_fastq' == input_file_type:
        for file1, file2 in starting_files:
                for input_file in ( file1, file2 ):
                    tmp_name = os.path.splitext( os.path.basename( input_file ) )
                    basename = tmp_name[ 0 ]
                    output_prefix = "{dir}/{file}".format(
                                            dir = Geno.dir[ 'fastq' ],
                                            file = basename 
                                            )
                    output_suffix = tmp_name[ 1 ]
                    input_fastq_files.append( [ input_file, output_prefix, output_suffix ] )
    elif 'single_fastq' == input_file_type:
        for file_name in starting_files:
            tmp_name = os.path.splitext( os.path.basename( file_name ) )
            basename = tmp_name[ 0 ]
            output_prefix = "{dir}/{file}".format(
                                    dir = Geno.dir[ 'fastq' ],
                                    file = basename
                                    )
            output_suffix = tmp_name[ 1 ]
            input_fastq_files.append( [ file_name, output_prefix, output_suffix ] )
    elif 'bam' == input_file_type:
        pass

    return input_fastq_files


#
# For STAGE 2 cutadapt
#
def generate_parameters_for_cutadapt():
    """
    Generate parameter list for STAGE 1 spilt_fastq_files

    Create a list of the following
        input_file
        output_file

    """
    input_file_type = Geno.job.get( 'input_file_type' )
    if 'paired_fastq' == input_file_type:
        for file1, file2 in starting_files:
                for input_file in ( file1, file2 ):
                    tmp_name = os.path.splitext( os.path.basename( input_file ) )
                    basename = tmp_name[ 0 ]
                    glob_input_file_name = "{dir}/{file}*{ext}".format(
                                                dir = Geno.dir[ 'fastq' ],
                                                file = basename, 
                                                ext = tmp_name[ 1 ] )

                    for tmp_file in glob( glob_input_file_name ):
                        tmp_basename = os.path.basename( tmp_file )
                        output_file = "{dir}/cutadapt_{file}".format(
                                                dir = Geno.dir[ 'fastq' ],
                                                file = tmp_basename )
                        yield [ tmp_file, output_file ]
    elif 'sinle_fastq' == input_file_type:
        for file_name in starting_files:
            tmp_name = os.path.splitext( os.path.basename( file_name ) )
            basename = tmp_name[ 0 ]
            glob_input_file_name = "{dir}/{file}*{ext}".format(
                                        dir = Geno.dir[ 'fastq' ],
                                        file = basename,
                                        ext = tmp_name[ 1 ] )
            for tmp_file in glob( glob_input_file_name ):
                tmp_basename = os.path.basename( tmp_file )
                output_file = "{dir}/cutadapt_{file}".format(
                                        dir = Geno.dir[ 'fastq' ],
                                        file = tmp_basename )
                yield [ tmp_file, output_file ]



#
# For STAGE 2 bwa_mem
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

    for parameter in split_fastq_parameters:
        split_files = sorted( glob( "{prefix}*{suffix}".format(
                                        prefix = parameter[ 1 ],
                                        suffix = parameter[ 2 ] 
                        ) ) )
        for split_file in split_files:
            bam_file = "{dir}/{file}.bam".format(
                                dir = Geno.dir[ 'bam' ],
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

    for id1 in range( 0, len( split_fastq_parameters ), 2 ):
        split_files1 = sorted( glob( "{filename}*".format( filename = split_fastq_parameters[ id1 ][ 1 ] ) ) )
        split_files2 = sorted( glob( "{filename}*".format( filename = split_fastq_parameters[ id1 + 1 ][ 1 ] ) ) )

        for file1, file2 in zip( split_files1, split_files2 ):
            bam_file = "{dir}/{file}.bam".format(
                                dir = Geno.dir[ 'bam' ],
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
    input_file_type = Geno.job.get( 'input_file_type' )

    if 'paired_fastq' == input_file_type:
        parameter_list = generate_parameters_for_bwa_mem_pair()
    elif 'single_fastq' == input_file_type:
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
        basename = "{dir}/{basename}".format(
                            dir = Geno.dir[ 'bam' ],
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

#
# Stage 1: split_fastq
#
def split_fastq(
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
        function_name = 'split_fastq'
        suffix_len = len( output_suffix )
        shell_script_name = res.file_timestamp_format.format(
                                        name=function_name,
                                        year=now.year,
                                        month=now.month,
                                        day=now.day,
                                        hour=now.hour,
                                        min=now.minute,
                                        msecond=now.microsecond )
        shell_script_full_path = "{script}/{file}.sh".format(
                                        script = Geno.dir[ 'script' ],
                                        file = shell_script_name )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( res.splitfile.format(
                                        log = Geno.dir[ 'log' ],
                                        lines_per_file = Geno.job.get( 'split_fastq_line_number' ),
                                        input_file = input_file,
                                        suffix_len = suffix_len,
                                        output_suffix = output_suffix,
                                        output_prefix = output_prefix ) )
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
# Stage 2: cutadapt
#
def cutadapt(
    input_file,
    output_file
    ):
    """
       Apply cutadapt to remove adaptor sequences from reads
    """

    try:
        log.info( "# {function}\n".format( function = whoami() ) +
                     "in: {input}\n".format( input=input_file) +
                     "out: {output}\n".format( output=output_file) )

        #
        # Make sure files exist.
        #
        if not os.path.isfile( input_file ):
            log.error( "file: {file} does not exist.".format( file=input_file ) )
            return 0

        #
        # Make shell script
        #
        function_name = 'cutadapt'
        now = datetime.now()
        shell_script_name = res.file_timestamp_format.format(
                                        name=function_name,
                                        year=now.year,
                                        month=now.month,
                                        day=now.day,
                                        hour=now.hour,
                                        min=now.minute,
                                        msecond=now.microsecond )
        shell_script_full_path = "{script}/{file}.sh".format(
                                        script = Geno.dir[ 'script' ],
                                        file = shell_script_name )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( res.cutadapt.format(
                                        log = Geno.dir[ 'log' ],
                                        infastq = input_file,
                                        outfastq = output_file,
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
        Geno.RT.runtask( Geno.job.get( 'job_type' )[ function_name ],
                         Geno.job.get( 'memory' )[ function_name ],
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
#
# Stage 2: bwa_mem
#
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
        function_name = 'bwa_mem'
        now = datetime.now()
        shell_script_name = res.file_timestamp_format.format(
                                        name=function_name,
                                        year=now.year,
                                        month=now.month,
                                        day=now.day,
                                        hour=now.hour,
                                        min=now.minute,
                                        msecond=now.microsecond )
        shell_script_full_path = "{script}/{file}.sh".format(
                                        script = Geno.dir[ 'script' ],
                                        file = shell_script_name )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( res.bwa_mem.format(
                                        log = Geno.dir[ 'log' ],
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
        Geno.RT.runtask( Geno.job.get( 'job_type' )[ function_name ],
                         Geno.job.get( 'memory' )[ function_name ],
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

#
# Stage 3: merge_bam
#
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
        shell_script_name = res.file_timestamp_format.format(
                                        name='merge_bam',
                                        year=now.year,
                                        month=now.month,
                                        day=now.day,
                                        hour=now.hour,
                                        min=now.minute,
                                        msecond=now.microsecond )
        shell_script_full_path = "{script}/{file}.sh".format(
                                        script = Geno.dir[ 'script' ],
                                        file = shell_script_name )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( res.merge_bam.format(
                                        log = Geno.dir[ 'log' ],
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
def fisher_mutation_call(
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
        function_name = 'fisher_mutation_call'
        shell_script_name = res.file_timestamp_format.format(
                                        name=function_name,
                                        year=now.year,
                                        month=now.month,
                                        day=now.day,
                                        hour=now.hour,
                                        min=now.minute,
                                        msecond=now.microsecond )
        shell_script_full_path = "{script}/{file}.sh".format(
                                        script = Geno.dir[ 'script' ],
                                        file = shell_script_name )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( res.mutation_call.format(
                                        log = Geno.dir[ 'log' ],
                                        input_bam = input_file,
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
#   Global variables
#
global starting_files
global split_fastq_parameters
global bwa_mem_parameters
global merge_bam_parameters

#####################################################################
#
#   STAGE 0 data preparation
#
starting_files = get_starting_files()
split_fastq_parameters = generate_parameters_for_split_fastq( starting_files )

#####################################################################
#
#   STAGE 1 split_fastq
#
@active_if ( 'split_fastq' in Geno.job.get( 'tasks' )[ 'WGS' ] )
@parallel( split_fastq_parameters )
@check_if_uptodate( check_file_exists_for_split_fastq )
def stage_1( input_file, output_prefix, output_suffix ):
    split_fastq( input_file, output_prefix, output_suffix )

#####################################################################
#
#   STAGE 2 cutadapt
#
@active_if ( 'cutadapt' in Geno.job.get( 'tasks' )[ 'WGS' ] )
@files( generate_parameters_for_cutadapt )
@check_if_uptodate( check_file_exists_for_cutadapt )
@follows( stage_1 )
def stage_2( input_file, output_file ):
    cutadapt( input_file, output_file )

#####################################################################
#
#   STAGE 3 fastq to bam
#
@active_if ( 'bwa_mem' in Geno.job.get( 'tasks' )[ 'WGS' ] )
@files( generate_parameters_for_bwa_mem )
@check_if_uptodate( check_file_exists_for_bwa_mem )
@follows( stage_2 )
def stage_3( input_file1, input_file2, output_file ):
    bwa_mem( input_file1, input_file2, output_file )

#####################################################################
#
#   STAGE 4 merge
#
@active_if ( 'merge_bam' in Geno.job.get( 'tasks' )[ 'WGS' ] )
@files( generate_parameters_for_merge_bam )
@check_if_uptodate( check_file_exists_for_merge_bam )
@follows( stage_3 )
def stage_4( input_prefix, input_suffix, output_file ):
    merge_bam( input_prefix, input_suffix, output_file )

#####################################################################
#
#   STAGE 5 mutation call
#
#@active_if ( 'fisher_mutation_call' in Geno.job.get( 'tasks' )[ 'WGS' ] )
# @files( generate_parameters_for_mutation_call )
# @check_if_uptodate( check_file_exists_for_mutation_call )
# @follows( stage_4 )
# def stage_5( input_file, output_file ):
#     fisher_mutation_call( input_file, output_file )


#####################################################################
#
#   LAST STAGE 
#
@follows( stage_4 )
def last_function():
    log.info( "Genomon pipline has finished successflly!" )
    return True
