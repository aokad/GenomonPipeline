"""
star_pipeline.py

"""

import sys
import os
import shutil
from datetime import datetime
from ruffus import *
from runtask import RunTask

#####################################################################
#
# Private modules
#
from __main__ import *
from resource import genomon_rc as res
from resource import star_resource as star_res
from utils import *
from sample import Sample

#####################################################################
#
# Subroutines
#
def save_status_of_this_process( process_name, output_file, return_code ):

    use_subdir = ( Geno.job.get_job( 'sample_subdir' ) != None )
    Geno.status.save_status( process_name, output_file, return_code, use_subdir = use_subdir )

def get_status_of_this_process( process_name, output_file ):

    use_subdir = ( Geno.job.get_job( 'sample_subdir' ) != None )
    exit_status = Geno.status.check_exit_status(
                    process_name,
                    output_file,
                    use_subdir = use_subdir ) 

    return exit_status
#####################################################################
#
def check_file_exists(input_file, output_file):
    if not os.path.exists(output_file):
        return True, "Missing file %s" % output_file
    else:
        in_time = os.path.getmtime( input_file )
        out_time = os.path.getmtime( output_file )
        if in_time > out_time:
            return True, "Missing file %s" % output_file
        else:
            return False, "File %s exists" % output_file

def check_file_exists_for_star_genome( fasta, star_genome_dir ):

    exit_status = get_status_of_this_process( 'star_genome', star_genome_dir )
    if ( not os.path.exists( star_genome_dir + '/Genome' ) or
         not os.path.exists( star_genome_dir + '/SA' ) or
         not os.path.exists( star_genome_dir + '/SAindex' ) or
         exit_status != 0 ):
        return True, "Missing file %s" % output_file
    else:
        return False, "File %s exists" % output_file

def check_file_exists_for_star( input_file1, input_file2, output_prefix ):
    exit_status = get_status_of_this_process( 'star', output_prefix )
    if exit_status != 0:
        return True, "Missing file %s" % output_prefix
    else:
        output_tmp = output_prefix + '_Chimeric.out.junction'
        return check_file_exists( input_file1, output_tmp )

def check_file_exists_for_star_fusion( input_file1, input_file2, output_prefix ):
    exit_status = get_status_of_this_process( 'star_fusion', output_prefix )
    if exit_status != 0:
        return True, "Missing file %s" % output_prefix
    else:
        output_tmp = output_prefix + '.fusion_candidates.txt'
        return check_file_exists( input_file1, output_tmp )


#####################################################################
#
def generate_params_for_star_genome( ):
    yield Geno.conf.get( 'REFERENCE', 'ref_fasta' ), \
          Geno.conf.get( 'REFERENCE', 'ref_gtf' ), \
          Geno.conf.get( 'REFERENCE', 'star_genome' )
    
def generate_params_for_star( ):
    global Sample
    Sample.make_param( 'star', None, '', 'star', 2, 1 )
    for infile1, infile2, outfile1, outfil2 in Sample.param( 'star' ):
        yield infile1, infile2, outfile1

def generate_params_for_star_fusion( ):
    global Sample
    Sample.make_param( 'star_fusion', None, '', 'star_fusion', 1, 1 )
    for infile1, infile2, outdir1, outdir2 in Sample.param( 'star_fusion' ):
        yield ( infile1 + '_Chimeric.out.sam',
                infile1 + '_Chimeric.out.junction',
                outdir1 )

#####################################################################
#
#   STAGE 0 data preparation
#
Sample = Sample()

#####################################################################
#
#   STAGE 1 fastq to bam by star
#
@active_if( 'star_genome' in Geno.job.get_job( 'tasks' )[ 'STAR' ] )
@files( generate_params_for_star_genome )
@check_if_uptodate( check_file_exists_for_star_genome )
def star_genome(
        fasta,
        gtf,
        output_dir
        ):
    """
        Stage 1: star_genome: make genome data for star

    """
    return_code = True

    try:
        function_name = whoami()
        log.info( "#{function}".format( function = function_name ) )

        #
        # Make shell script
        #
        shell_script_full_path = make_script_file_name( function_name, Geno )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( star_res.star.format(
                                        log = Geno.dir[ 'log' ],
                                        ref_fasta = fasta,
                                        ref_gtf = gtf,
                                        out_dir = output_dir,
                                        star = Geno.conf.get( 'SOFTWARE', 'star' ),
                                        additional_params = Geno.job.get_param( 'star_genome', 'additional_params' )
                                    )
                                )
        shell_script_file.close()

        #
        # Run
        #
        return_code = Geno.RT.runtask(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ] )

        if return_code != 0:
            log.error( "{function}: runtask failed".format( function = function_name ) )
            raise

        save_status_of_this_process( function_name, output_dir, return_code )

    except IOError as (errno, strerror):
        log.error( "{function}: I/O error({num}): {error}".format(
                        function = whoami(),
                        num = errno,
                        error = strerror)
                )
        return_code = False

    except ValueError:
        log.error( "{function}: ValueError".format(
                        function = whoami()
                    )
                )
        return_code = False

    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        log.error( "{function}: Unexpected error: {error}.".format(
                    function = whoami()) , error = sys.exc_info()[0] )
        log.error("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )
        return_code = False


    return return_code

#####################################################################
#
#   STAGE 2 fastq to bam by star
#
@follows( star_genome )
@active_if( 'star' in Geno.job.get_job( 'tasks' )[ 'STAR' ] )
@files( generate_params_for_star )
@check_if_uptodate( check_file_exists_for_star )
def star(
        input_file1,
        input_file2,
        output_prefix,
        ):
    """
        Stage 2: star

    """
    return_code = True

    try:
        function_name = whoami()
        log.info( "#{function}".format( function = function_name ) )

        if input_file1[ -3: ] == '.gz' or input_file1[ -4: ] == '.gz2':


            #
            # Extract file
            #
            if ( len( input_file2 ) >= 3 and ( input_file2[ -3: ] == '.gz' or input_file2[ -4: ] == '.gz2' ) ):
                id = 2
                fastq_file2 = input_file2[ :-3 ]
            else:
                id = 1
                fastq_file2 = ''

            fastq_file1 = input_file1[ :-3 ]
            array_data  = "IN_FILE=(\n[1]=\"{file1}\"\n[2]=\"{file2}\"\n)\n".format(
                                        file1 = input_file1, file2 = input_file2 ) + \
                          "OUT_FILE=(\n[1]=\"{file1}\"\n[2]=\"{file2}\"\n)\n".format(
                                        file1 = fastq_file1, file2 = fastq_file2 )

            if ( get_status_of_this_process( function_name + '_decomp', output_prefix ) != 0 or
                 not os.path.exists( fastq_file1 ) or
                 ( fastq_file2 != '' and not os.path.exists( fastq_file2 ) ) ):
                #
                # Make shell script
                #
                shell_script_full_path = make_script_file_name( function_name + '_decomp', Geno )
                shell_script_file = open( shell_script_full_path, 'w' )
                shell_script_file.write( star_res.extract_fastq.format(
                                                log = Geno.dir[ 'log' ],
                                                array_data = array_data,
                                                input_file = "${IN_FILE[$SGE_TASK_ID]}",
                                                output_file = "${OUT_FILE[$SGE_TASK_ID]}"
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

                if return_code != 0:
                    log.error( "{function}: runtask failed".format( function = function_name ) )
                    raise

                save_status_of_this_process( function_name + '_decomp' , output_prefix, return_code )

        else:
            fastq_file1 = input_file1
            fastq_file2 = input_file2

        #
        # Run Star
        # Make shell script
        #
        shell_script_full_path = make_script_file_name( function_name, Geno )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( star_res.star_map.format(
                                        log = Geno.dir[ 'log' ],
                                        star_genome = Geno.conf.get( 'REFERENCE', 'star_genome' ),
                                        fastq1 = fastq_file1,
                                        fastq2 = fastq_file2,
                                        out_prefix = output_prefix + '_',
                                        star = Geno.conf.get( 'SOFTWARE', 'STAR' ),
                                        additional_params = Geno.job.get_param( 'star', 'additional_params' )
                                    )
                                )
        shell_script_file.close()

        #
        # Run
        #
        return_code = Geno.RT.runtask(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ] )

        if return_code != 0:
            log.error( "{function}: runtask failed".format( function = function_name ) )
            raise

        save_status_of_this_process( function_name, output_prefix, return_code )

    except IOError as (errno, strerror):
        log.error( "{function}: I/O error({num}): {error}".format(
                        function = whoami(),
                        num = errno,
                        error = strerror)
                )
        return_code = False

    except ValueError:
        log.error( "{function}: ValueError".format(
                        function = whoami()
                    )
                )
        return_code = False

    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        log.error( "{function}: Unexpected error: {error}.".format(
                    function = whoami()) , error = sys.exc_info()[0] )
        log.error("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )
        return_code = False


    return return_code

#####################################################################
#
#   STAGE 3
#
@follows( star )
@active_if( 'star_fusion' in Geno.job.get_job( 'tasks' )[ 'STAR' ] )
@files( generate_params_for_star_fusion )
@check_if_uptodate( check_file_exists_for_star_fusion )
def star_fusion(
        sam,
        junction,
        output_prefix
        ):
    """
        Stage 3: star_fusion

    """
    return_code = True

    try:
        function_name = whoami()
        log.info( "#{function}".format( function = function_name ) )

        #
        # Make shell script
        #
        shell_script_full_path = make_script_file_name( function_name, Geno )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( star_res.star_fusion.format(
                                        log = Geno.dir[ 'log' ],
                                        chimeric_sam = sam,
                                        chimeric_junction = junction,
                                        out_prefix = output_prefix,
                                        gtf_file = Geno.conf.get( 'REFERENCE', 'ref_gtf' ),
                                        star_fusion = Geno.conf.get( 'SOFTWARE', 'STAR-Fusion' ),
                                        additional_params = Geno.job.get_param( 'star_fusion', 'additional_params' )
                                    )
                                )
        shell_script_file.close()

        #
        # Run
        #
        return_code = Geno.RT.runtask(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ] )

        if return_code != 0:
            log.error( "{function}: runtask failed".format( function = function_name ) )
            raise

        save_status_of_this_process( function_name, output_prefix, return_code )

    except IOError as (errno, strerror):
        log.error( "{function}: I/O error({num}): {error}".format(
                        function = whoami(),
                        num = errno,
                        error = strerror)
                )
        return_code = False

    except ValueError:
        log.error( "{function}: ValueError".format(
                        function = whoami()
                    )
                )
        return_code = False

    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        log.error( "{function}: Unexpected error: {error}.".format(
                    function = whoami()) , error = sys.exc_info()[0] )
        return_code = False


    return return_code

#####################################################################
#
#   LAST STAGE 
#

@follows( star_fusion )
def last_function():
    log.info( "Genomon pipline has finished successflly!" )
    return True
