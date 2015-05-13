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

#####################################################################
#
# Private modules
#
from __main__ import *
from resource import genomon_rc as res
from resource import wgs_resource as wgs_res
from utils import *
from sample import Sample

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

def check_file_exists_for_merge(
    input_file_list,
    output_file
    ):

    """
    Checks if output file exists for merge_bam

    """

    if not os.path.exists( output_file ):
        return True, "Missing file {outputfile} for input file.".format(
                            outputfile = output_file )
    else:
        return False, "File {outputfile} exists for input file.".format(
                            outputfile = output_file )

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
        input_file_list,
        output_file
    ):

    """
    Checks if output file exists for bwa_mem

    """

    if Geno.job.get_job( 'use_biobambam' ):
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

def check_file_exists_for_input2_output(
        input_file1,
        input_file2,
        output_file
    ):

    """
    Checks if output file exists

    """

    if not os.path.exists( output_file ):
        return True, "Missing file {outputfile} for {inputfile}.".format(
                            outputfile = output_file,
                            inputfile = input_file1 )

    else:
        in_time = os.path.getmtime( input_file1 )
        out_time = os.path.getmtime( output_file )
        if in_time > out_time:
            return True, "{output} is older than {input}.".format( output = output_file, input = input_file1 )
        else:
            return False, "File {output} exits for {input}.".format( output = output_file, input = input_file1 )


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
    input_file_list = {}
    for param in Sample.param( 'merge_bam' ):
        dir_name = os.path.dirname( param[ 0 ] )
        if not ( dir_name in input_file_list ):
            input_file_list[ dir_name ] = []
        input_file_list[ dir_name ].append( param[ 0 ] )

    for dir_name in input_file_list.keys():
        return_list = [ input_file_list[ dir_name ], dir_name + '/' + Geno.job.get_job( 'sample_name' ) + '_merged.bam' ]
        yield return_list
        

#
# For STAGE 6 markduplicates
#
def generate_params_for_markduplicates ():
    """
    Generate parameter list for markduplicates

    """
    
    Sample.make_param( 'markduplicates', '.bam', 'bam', 1, 1 )
    input_file_list = {}
    for param in Sample.param( 'markduplicates' ):
        dir_name = os.path.dirname( param[ 0 ] )
        if not ( dir_name in input_file_list ):
            input_file_list[ dir_name ] = []
        input_file_list[ dir_name ].append( param[ 0 ] )

    for dir_name in input_file_list.keys():
        return_list = [ input_file_list[ dir_name ], dir_name + '/' + Geno.job.get_job( 'sample_name' ) + '_markdup.bam' ]
        yield return_list

#
# For STAGE 7 fisher_mutation_call
#
def generate_params_for_fisher_mutation_call():
    """
    Generate parameter list for fisher_mutation_call

    """
    
    Sample.make_param( 'fisher_mutation_call', '.txt', 'mutation', 2, 1 )
    input_file_list = {}
    for param in Sample.param( 'fisher_mutation_call' ):
        dir_name = os.path.dirname( param[ 0 ] )
        if not ( dir_name in input_file_list ):
            if 'markduplicates' in Geno.job.get_job( 'tasks' )[ 'WGS' ]:
                input_bam =  dir_name + '/' + Geno.job.get_job( 'sample_name' ) + '_markdup.bam'
            elif 'merge_bam' in Geno.job.get_job( 'tasks' )[ 'WGS' ]:
                input_bam =  dir_name + '/' + Geno.job.get_job( 'sample_name' ) + '_merged.bam'
            else:
                input_bam =  param[ 0 ]

            input_file_list[ dir_name ] = input_bam
            mutation_dir_name = os.path.dirname( param[ 2 ] )
            return_list = [ input_bam,
                            'None',
                            mutation_dir_name + '/' + Geno.job.get_job( 'sample_name' ) + '.txt' ]
            yield return_list


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
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ],
                            id_start = 1,
                            id_end = id )
    Geno.status.save_status( 'extract_fastq', input_file1, return_code )
    if return_code != 0:
        with log_mutex:
            log.error( "{function}: runtask failed",format( function = 'extract_fastq' ) )
        raise


#####################################################################
#
#   Main functions
#

#
# Stage 1: bam2fastq
#
def bam2fastq(
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
        with log_mutex:
            log.info( "#{function}".format( function = function_name ) )

        #
        # Make sure files exist.
        #
        if not os.path.isfile( input_file ):
            with log_mutex:
                log.error( "file: {file} does not exist.".format( file=input_file ) )
            raise

        #
        # Make shell script
        #
        shell_script_full_path = make_script_file_name( function_name, Geno )
        shell_script_file = open( shell_script_full_path, 'w' )

        file_type = Geno.job.get_job( 'input_file_type' )
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
                        shell_script_full_path,
                        Geno.job.get_job( 'cmd_options' )[ function_name ],
                        shell_script_full_path )

        Geno.status.save_status( function_name, input_file1, return_code )

        if return_code != 0:
            with log_mutex:
                log.error( "{function}: runtask failed".format( function = function_name ) )
            raise

    except IOError as (errno, strerror):
        with log_mutex:
            log.error( "{function}: I/O error({num}): {error}".format(num = errno, error = strerror, function = whoami() ) )
        return_value = False

    except ValueError:
        with log_mutex:
            log.error( "{function}: ValueError".format( function = whoami() ) )
        return_value = False

    except:
        with log_mutex:
            log.error( "{function}: Unexpected error: {error}".format( function = whoami(), error = sys.exc_info()[0] ) )
        return_value = False

    else:
        return_value = True


    return return_value

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
        with log_mutex:
            log.info( "#{function}".format( function = function_name ) )

        #
        # Make sure files exist.
        #
        if not os.path.isfile( input_file1 ):
            with log_mutex:
                log.error( "file: {file} does not exist.".format( file=input_file1 ) )
            raise

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
                                        fastq_filter = Geno.job.get_param( 'split_fastq', 'fastq_filter' ),
                                        array_data = array_in + array_out,
                                        lines_per_file = Geno.job.get_param( 'split_fastq', 'split_fastq_line_number' ),
                                        input_file = input_file,
                                        suffix_len = suffix_len,
                                        output_suffix = output_suffix,
                                        output_prefix = output_prefix ) )
        shell_script_file.close()

        #
        # Run
        #
        return_code = Geno.RT.run_arrayjob(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ],
                            id_start = 1,
                            id_end = id )
        Geno.status.save_status( function_name, input_file1, return_code )
        if return_code != 0:
            with log_mutex:
                log.error( "{function}: runtask failed".format( function = function_name ) )
            raise


    except IOError as (errno, strerror):
        with log_mutex:
            log.error( "{function}: I/O error({num}): {error}".format(num = errno, error = strerror, function = whoami() ) )
        return_value = False

    except ValueError:
        with log_mutex:
            log.error( "{function}: ValueError".format( function = whoami() ) )
        return_value = False

    except:
        with log_mutex:
            log.error( "{function}: Unexpected error: {error}".format( function = whoami(), error = sys.exc_info()[0] ) )
        return_value = False

    else:
        return_value = True

    return return_value


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
        with log_mutex:
            log.info( "#{function}".format( function = function_name ) )

        file_ext = os.path.splitext( input_file1 )[ 1 ]
        if file_ext == '.bz2':
            extract_fastq( ( input_file1, input_file2 ), file_ext )

        #
        # Make data for array job
        #
        if not os.path.isfile( input_file ):
            with log_mutex:
                log.error( "file: {file} does not exist.".format( file=input_file ) )
            raise

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
                                        optadapters = '-a ' + ' -a '.join( Geno.job.get_param( 'cutadapt', 'adaptor' ) ),
                                        casavacode = 2,
                                        cutadapt = Geno.conf.get( 'SOFTWARE', 'cutadapt' ),
                                        scriptdir = Geno.dir[ 'script' ],
                                        ) )
        shell_script_file.close()

        #
        # Run
        #
        return_code = Geno.RT.run_arrayjob(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ],
                            id_start = 1,
                            id_end = 2 )
        Geno.status.save_status( function_name, input_file1, return_code )
        if return_code != 0:
            with log_mutex:
                log.error( "{function}: runtask failed".format( function = function_name ) )
            raise


    except IOError as (errno, strerror):
        with log_mutex:
            log.error( "{function}: I/O error({num}): {error}".format(function = whoami(), num = errno, error = strerror) )
        return_value = False

    except ValueError:
        with log_mutex:
            log.error( "{function}: ValueError".format( function = whoami() ) )
        return_value = False

    except:
        with log_mutex:
            log.error( "{function}: Unexpected error: {error}".format( function = whoami(), error = sys.exc_info()[0] ) )
        return_value = False

    else:
        return_value = True

    return return_value

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
        with log_mutex:
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
                                        read_group = Geno.job.get_param( 'bwa_mem', 'bwa_read_group' ),
                                        min_score = Geno.job.get_param( 'bwa_mem', 'min_score' ),
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
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ],
                            id_start = 1,
                            id_end = id )
        Geno.status.save_status( function_name, input_file1, return_code )
        if return_code != 0:
            with log_mutex:
                log.error( "{function}: runtask failed".format( function = function_name ) )
            raise
            

    except IOError as (errno, strerror):
        with log_mutex:
            log.error( "{function}: I/O error({num}): {error}".format(function = whoami(), num = errno, error = strerror) )
        return_value = False

    except ValueError:
        with log_mutex:
            log.error( "{function}: ValueError".format( function = whoami() ) )
        return_value = False

    except:
        with log_mutex:
            log.error( "{function}: Unexpected error: {error}".format( function = whoami(), error = sys.exc_info()[0] ) )
        return_value = False

    else:
        return_value = True

    return return_value

#
# Stage 5: merge_bam
#
def merge_bam(
    input_file_list,
    output_file,
    use_biobambam
    ):
    """
       Merge split bam files

    """

    try:
        function_name = whoami()
        with log_mutex:
            log.info( "#{function}".format( function = function_name ) )

        #
        # Make data for array job 
        #
        input_files = ''
        for input_file in input_file_list:
            ( input_prefix, input_suffix ) = os.path.splitext( input_file )

            #
            # Make shell script
            #
            if use_biobambam:
                input_file_name = "{input_prefix}*_bamsorted{input_suffix}".format(
                                            input_prefix = input_prefix,
                                            input_suffix = input_suffix )
                bam_merge_resource = wgs_res.biobambam_merge_bam
                for file_name in glob( input_file_name ):
                    input_files += "I={file_name} ".format( file_name = file_name )

            else:
                input_files += "{input_prefix}*_sorted{input_suffix} ".format(
                                            input_prefix = input_prefix,
                                            input_suffix = input_suffix )

        if use_biobambam:
            bam_merge_resource = wgs_res.biobambam_merge_bam
        else:
            bam_merge_resource = wgs_res.samtools_merge_bam

        shell_script_full_path = make_script_file_name( function_name, Geno )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( bam_merge_resource.format(
                                        log = Geno.dir[ 'log' ],
                                        input_bam_files = input_files,
                                        output_bam_file = output_file,
                                        samtools = Geno.conf.get( 'SOFTWARE', 'samtools' ),
                                        biobambam = Geno.conf.get( 'SOFTWARE', 'biobambam' )
                                        ) )
        shell_script_file.close()

        #
        # Run
        #
        return_code = Geno.RT.runtask(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ] )
        Geno.status.save_status( function_name, input_file, return_code )
        if return_code != 0:
            with log_mutex:
                log.error( "{function}: runtask failed".format( function = function_name ) )
            raise


    except IOError as (errno, strerror):
        with log_mutex:
            log.error( "{function}: I/O error({num}): {error}".format( function = whoami(), num = errno, error = strerror) )
        return_value = False

    except ValueError:
        with log_mutex:
            log.error( "{function}: ValueError".format( function = whoami() ) )
        return_value = False

    except:
        with log_mutex:
            log.error( "{function}: Unexpected error: {error}".format( function = whoami(), error = sys.exc_info()[0] ) )
        return_value = False

    else:
        return_value = True

    return return_value

#
# Stage 6: markduplicates
#
def markduplicates(
    input_file_list,
    output_file,
    use_biobambam
    ):
    """
       Mutaion calling

    """

    try:
        function_name = whoami()
        with log_mutex:
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
            for input_file in input_file_list:
                ( input_prefix, input_suffix ) = os.path.splitext( input_file )
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
                                            output_bam = output_file,
                                            biobambam = Geno.conf.get( 'SOFTWARE', 'biobambam' )
                                        )
                                    )
            shell_script_file.close()

        else:
            #
            # Use Picard MarkDuplicates for samtools merge result
            #
            tmp_options = Geno.job.get_job( 'cmd_options' )[ 'markduplicates' ]
            tmp_memory = int( tmp_options[ tmp_options.find( 's_vmem=' ) + len('s_vmem=') : tmp_options.find('G') ] )
            input_file = output_file.replace( 'markdup', 'merged' )

            if tmp_memory > 3:
                java_memory = str( tmp_memory - 2 ) + 'G'
            else:
                java_memory = '1G'

            shell_script_full_path = make_script_file_name( function_name, Geno )
            shell_script_file = open( shell_script_full_path, 'w' )
            shell_script_file.write( wgs_res.markduplicates.format(
                                            log = Geno.dir[ 'log' ],
                                            input_bam = input_file,
                                            output_bam = output_file,
                                            memory = java_memory,
                                            samtools = Geno.conf.get( 'SOFTWARE', 'samtools' ),
                                            picard = Geno.conf.get( 'SOFTWARE', 'picard' )
                                        )
                                    )
            shell_script_file.close()

        #
        # Run
        #
        return_code = Geno.RT.runtask(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ] )
        Geno.status.save_status( function_name, input_file_list[ 0 ], return_code )
        if return_code != 0:
            with log_mutex:
                log.error( "{function}: runtask failed".format( function = function_name ) )
            raise


    except IOError as (errno, strerror):
        with log_mutex:
            log.error( "{function}: I/O error({num}): {error}".format( function = whoami(), num = errno, error = strerror) )
        return_value = False

    except ValueError:
        with log_mutex:
            log.error( "{function}: ValueError".format( function = whoami() ) )
        return_value = False

    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        log.error( "{function}: Unexpected error: {error}".format( function = whoami(), error = sys.exc_info()[0] ) )
        log.error("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )
        return_value = False

    else:
        return_value = True

    return return_value

#
# Stage 7: fisher_mutation_call
#
def fisher_mutation_call(
    control_input_file,
    disease_input_file,
    output_file,
    ):
    """
       Mutaion calling

    """

    try:
        function_name = whoami()
        with log_mutex:
            log.info( "#{function}".format( function = function_name ) )

        #
        # Make shell script
        #
        shell_script_full_path = make_script_file_name( function_name, Geno )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( wgs_res.fisher_mutation_call.format(
                                        log = Geno.dir[ 'log' ],
                                        ref_fa = Geno.conf.get( 'REFERENCE', 'ref_fasta' ),
                                        control_input_bam = control_input_file,
                                        disease_input_bam = disease_input_file,
                                        output_txt = output_file,
                                        max_indel = Geno.job.get_param( 'fisher_mutation_call', 'max_indel' ),
                                        max_distance = Geno.job.get_param( 'fisher_mutation_call', 'max_distance' ),
                                        base_quality = Geno.job.get_param( 'fisher_mutation_call', 'base_quality' ),
                                        map_quality = Geno.job.get_param( 'fisher_mutation_call', 'map_quality' ),
                                        mismatch_rate = Geno.job.get_param( 'fisher_mutation_call', 'mismatch_rate' ),
                                        min_depth = Geno.job.get_param( 'fisher_mutation_call', 'min_depth' ),
                                        samtools = Geno.conf.get( 'SOFTWARE', 'samtools' ),
                                        python = Geno.conf.get( 'SOFTWARE', 'python' ),
                                        script_dir = Geno.dir[ 'script' ]
                                    )
                                )
        shell_script_file.close()

        #
        # Run
        #
        return_code = Geno.RT.run_arrayjob(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ],
                            id_start = 1,
                            id_end = wgs_res.interval_num )
        Geno.status.save_status( function_name, control_input_file, return_code )
        if return_code != 0:
            with log_mutex:
                log.error( "{function}: runtask failed".format( function = function_name ) )
            raise

        #
        # Make shell script
        #
        shell_script_full_path = make_script_file_name( 'merge_fisher_result', Geno )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( wgs_res.merge_fisher_result.format(
                                        log = Geno.dir[ 'log' ],
                                        script_dir = Geno.dir[ 'script' ],
                                        output_txt = output_file ) )
        shell_script_file.close()

        #
        # Run
        #
        return_code = Geno.RT.runtask(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ] )
        Geno.status.save_status( 'merge_fisher_result', control_input_file, return_code )
        if return_code != 0:
            with log_mutex:
                log.error( "{function}: runtask failed".format( function = 'merge_fisher_result' ) )
            raise

    except IOError as (errno, strerror):
        with log_mutex:
            log.error( "{function}: I/O error({num}): {error}".format( function = whoami(), num = errno, error = strerror) )
        return_value = False

    except ValueError:
        with log_mutex:
            log.error( "{function}: ValueError".format( function = whoami() ) )
        return_value = False

    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        with log_mutex:
            log.error( "{function}: Unexpected error: {error}".format( function = whoami(), error = sys.exc_info()[0] ) )
            log.error("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )
        return_value = False

    else:
        return_value = True

    return return_value

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
@active_if( 'bam2fastq' in Geno.job.get_job( 'tasks' )[ 'WGS' ] )
@parallel( generate_params_for_bam2fastq )
@check_if_uptodate( check_file_exists_for_input_output )
def stage_1( input_file, output_file ):
    return_value =  bam2fastq( input_file, output_file )
    if not return_value:
        raise

#####################################################################
#
#   STAGE 2 split_fastq
#   in:     fastq
#   out:    fastq * X
#
@follows( stage_1 )
@active_if( 'split_fastq' in Geno.job.get_job( 'tasks' )[ 'WGS' ] )
@files( generate_params_for_split_fastq )
@check_if_uptodate( check_file_exists_for_split_fastq )
def stage_2( input_file1, input_file2, output_file1, output_file2 ):
    return_value = split_fastq( input_file1, input_file2, output_file1, output_file2 )
    if not return_value:
        raise

#####################################################################
#
#   STAGE 3 cutadapt
#   in:     fastq
#   out:    fastq
#
@follows( stage_2 )
@active_if( 'cutadapt' in Geno.job.get_job( 'tasks' )[ 'WGS' ] )
@files( generate_params_for_cutadapt )
@check_if_uptodate( check_file_exists_for_input_output )
def stage_3( input_file1, input_file2, output_file1, output_file2 ):
    return_value = cutadapt( input_file1, input_file2, output_file1, output_file2 )
    if not return_value:
        raise

#####################################################################
#
#   STAGE 4 bwa_mem
#
#   in:     fastq1, fastq2
#   out:    bam
#
@follows( stage_3 )
@active_if( 'bwa_mem' in Geno.job.get_job( 'tasks' )[ 'WGS' ] )
@files( generate_params_for_bwa_mem )
@check_if_uptodate( check_file_exists_for_bwa_mem )
def stage_4(  input_file1, input_file2, output_file1, output_file2 ):
    return_value = bwa_mem( input_file1,
                            input_file2,
                            output_file1,
                            output_file2, 
                            Geno.job.get_job( 'use_biobambam' ) )
    if not return_value:
        raise

#####################################################################
#
#   STAGE 5 merge
#
#   in:     bam x {number}
#   out:    bam
#
@follows( stage_4 )
@active_if( 'merge_bam' in Geno.job.get_job( 'tasks' )[ 'WGS' ] )
@files( generate_params_for_merge_bam )
@check_if_uptodate( check_file_exists_for_merge)
def stage_5( input_file_list, output_file ):
    if ( Geno.job.get_job( 'use_biobambam' ) and
         'markduplicates' in Geno.job.get_job( 'tasks' )[ 'WGS'] ):
        return_value = True
    else:
        return_value = merge_bam( input_file_list,
                                  output_file,
                                  Geno.job.get_job( 'use_biobambam' ) )

    if not return_value:
        raise

#####################################################################
#
#   STAGE 6 markduplicates
#
#   in:     bam x {number}
#   out:    bam
#
@follows( stage_5 )
@active_if( 'markduplicates' in Geno.job.get_job( 'tasks' )[ 'WGS' ] )
@files( generate_params_for_markduplicates )
@check_if_uptodate( check_file_exists_for_merge )
def stage_6( input_file_list, output_file ):
    return_value =  markduplicates( input_file_list,
                                    output_file,
                                    Geno.job.get_job( 'use_biobambam' ) )

    if not return_value:
        raise

#####################################################################
#
#   STAGE 7 fisher_mutation_call
#
#   in:     bam
#   out:    txt
#
@follows( stage_6 )
@active_if ( 'fisher_mutation_call' in Geno.job.get_job( 'tasks' )[ 'WGS' ] )
@files( generate_params_for_fisher_mutation_call )
@check_if_uptodate( check_file_exists_for_input2_output )
def stage_7(  input_file1, input_file2, output_file ):
    return_value = fisher_mutation_call(  input_file1, input_file2, output_file )
    if not return_value:
        raise


#####################################################################
#
#   LAST STAGE 
#
@follows( stage_7 )
def last_function():
    with log_mutex:
        log.info( "Genomon pipline has finished successflly!" )
    return True

