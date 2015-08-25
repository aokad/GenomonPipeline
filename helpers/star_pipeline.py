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

use_subdir = ( Geno.job.get_job( 'sample_name' ) != None )

#####################################################################
#
# Subroutines
#

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

def check_file_exists_for_star_genome( fasta, gtf, star_genome_dir ):

    exit_status = get_status_of_this_process( 'star_genome', star_genome_dir, Geno, use_subdir )
    if ( not os.path.exists( star_genome_dir + '/Genome' ) or
         not os.path.exists( star_genome_dir + '/SA' ) or
         not os.path.exists( star_genome_dir + '/SAindex' ) or
         exit_status != 0 ):
        return True, "Missing file %s" % star_genome_dir
    else:
        return False, "File %s exists" % star_genome_dir

def check_file_exists_for_star( input_file1, input_file2, output_prefix ):
    exit_status = get_status_of_this_process( 'star', output_prefix, Geno, use_subdir )
    if exit_status != 0:
        return True, "Missing file %s" % output_prefix
    else:
        output_tmp = output_prefix + '_Log.final.out'
        return check_file_exists( input_file1, output_tmp )

def check_file_exists_for_star_fusion( input_file1, input_file2, output_prefix ):
    exit_status = get_status_of_this_process( 'star_fusion', output_prefix, Geno, use_subdir )
    if exit_status != 0:
        return True, "Missing file %s" % output_prefix
    else:
        output_tmp = output_prefix + '.junction_breakpts_to_genes.txt'
        print Geno.job.get_job( 'tasks' )
        return check_file_exists( input_file1, output_tmp )

def check_file_exists_for_fusionfusion( input_file1, input_file2, output_prefix ):
    exit_status = get_status_of_this_process( 'fusionfusion', output_prefix, Geno, use_subdir )
    if exit_status != 0:
        return True, "Missing file %s" % output_prefix
    else:
        output_tmp = output_prefix + '/star.fusion.result.txt'
        print output_tmp
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
    for infile1, infile2, outfile1, outfile2 in Sample.param( 'star' ):
        yield infile1, infile2, outfile1

def generate_params_for_star_fusion( ):
    global Sample
    Sample.make_param( 'star_fusion', None, '', 'star_fusion', 1, 1 )
    for infile1, infile2, outdir1, outdir2 in Sample.param( 'star_fusion' ):
        yield ( infile1 + '_Chimeric.out.sam',
                infile1 + '_Chimeric.out.junction',
                outdir1 )

def generate_params_for_fusionfusion( ):
    global Sample
    Sample.make_param( 'fusionfusion', None, '', 'fusionfusion', 1, 1 )
    for infile1, infile2, outdir1, outdir2 in Sample.param( 'fusionfusion' ):
        yield ( infile1 + '_Chimeric.out.sam',
                infile1 + '_Chimeric.out.sam',
                outdir1 )

#####################################################################
#
# Functions
#
def star_genome(
        fasta,
        gtf,
        output_dir
        ):
    """
        Stage 1: star_genome: make genome data for star

    """
    return_value = True

    try:
        function_name = whoami()
        log.info( "#{function}".format( function = function_name ) )

        #
        # Make shell script
        #
        args = {
            'log': Geno.dir[ 'log' ],
            'ref_fasta': fasta,
            'ref_gtf': gtf,
            'out_dir': output_dir,
            'star': Geno.conf.get( 'SOFTWARE', 'star' ),
            'scriptdir': Geno.dir[ 'script' ],
            'additional_params': Geno.job.get_param( 'star_genome', 'additional_params' )
        }
        shell_script_full_path = make_script_file( function_name, star_res.star_genome, Geno, **args )
        
        #
        # Run
        #
        return_code = Geno.RT.runtask(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ] )

        if return_code != 0:
            log.error( "{function}: runtask failed".format( function = function_name ) )
            raise

        save_status_of_this_process( function_name, output_dir, return_code, Geno, use_subdir )

    except IOError as (errno, strerror):
        log.error( "{function}: I/O error({num}): {error}".format(
                        function = whoami(),
                        num = errno,
                        error = strerror)
                )
        return_value = False

    except ValueError:
        log.error( "{function}: ValueError".format(
                        function = whoami()
                    )
                )
        return_value = False

    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        log.error( "{function}: Unexpected error: {error}.".format(
                    function = whoami() , error = sys.exc_info()[0] ) )
        log.error("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )
        return_value = False


    return return_value

def star(
        input_file1,
        input_file2,
        output_prefix,
        ):
    """
        Stage 2: star

    """
    return_value = True

    try:
        function_name = whoami()
        log.info( "#{function}".format( function = function_name ) )

        if input_file1[ -3: ] == '.gz' or input_file1[ -4: ] == '.bz2':
            #
            # Extract file
            #
            if input_file1[ -3: ] == '.gz':
                fastq_file1 = input_file1[ :-3 ]
            elif input_file1[ -4: ] == '.bz2':
                fastq_file1 = input_file1[ :-4 ]

            if ( len( input_file2 ) >= 3 ): # if input_file1 has a filename in it.
                id = 2
                if ( input_file2[ -3: ] == '.gz' ):
                    fastq_file2 = input_file2[ :-3 ]

                elif ( input_file2[ -4: ] == '.bz2' ):
                    fastq_file2 = input_file2[ :-4 ]

            else:
                id = 1
                fastq_file2 = ''

            array_data  = "IN_FILE=(\n[1]=\"{file1}\"\n[2]=\"{file2}\"\n)\n".format(
                                        file1 = input_file1, file2 = input_file2 ) + \
                          "OUT_FILE=(\n[1]=\"{file1}\"\n[2]=\"{file2}\"\n)\n".format(
                                        file1 = fastq_file1, file2 = fastq_file2 )

            if ( get_status_of_this_process( function_name + '_decomp', output_prefix, Geno, use_subdir ) != 0 or
                 not os.path.exists( fastq_file1 ) or
                 ( fastq_file2 != '' and not os.path.exists( fastq_file2 ) ) ):
                #
                # Make shell script
                #
                args = {
                    'log': Geno.dir[ 'log' ],
                    'array_data': array_data,
                    'input_file': "${IN_FILE[$SGE_TASK_ID]}",
                    'output_file': "${OUT_FILE[$SGE_TASK_ID]}",
                    'scriptdir': Geno.dir[ 'script' ]
                    }
                make_script_file( function_name + '_decomp', star_res.extract_fastq, Geno, **args)

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

                save_status_of_this_process( function_name + '_decomp' , output_prefix, return_code, Geno, use_subdir )

        else:
            fastq_file1 = input_file1
            fastq_file2 = input_file2

        #
        # Run Star
        # Make shell script
        #
        args = {
            'log': Geno.dir[ 'log' ],
            'star_genome': Geno.conf.get( 'REFERENCE', 'star_genome' ),
            'fastq1': fastq_file1,
            'fastq2': fastq_file2,
            'out_prefix': output_prefix + '_',
            'star': Geno.conf.get( 'SOFTWARE', 'STAR' ),
            'additional_params': Geno.job.get_param( 'star', 'additional_params' ),
            'scriptdir': Geno.dir[ 'script' ]
        }
        shell_script_full_path = make_script_file( function_name, star_res.star_map, Geno, **args )

        #
        # Run
        #
        return_code = Geno.RT.runtask(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ] )

        if return_code != 0:
            log.error( "{function}: runtask failed".format( function = function_name ) )
            raise

        save_status_of_this_process( function_name, output_prefix, return_code, Geno, use_subdir )

    except IOError as (errno, strerror):
        log.error( "{function}: I/O error({num}): {error}".format(
                        function = whoami(),
                        num = errno,
                        error = strerror)
                )
        return_value = False

    except ValueError:
        log.error( "{function}: ValueError".format(
                        function = whoami()
                    )
                )
        return_value = False

    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        log.error( "{function}: Unexpected error: {error}.".format(
                    function = whoami() , error = sys.exc_info()[0] ) )
        log.error("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )
        return_value = False


    return return_value

def star_fusion(
        sam,
        junction,
        output_prefix
        ):
    """
        Stage 3: star_fusion

    """
    return_value = True

    try:
        function_name = whoami()
        log.info( "#{function}".format( function = function_name ) )

        #
        # Make shell script
        #
        arg = { log: Geno.dir[ 'log' ],
                chimeric_sam: sam,
                chimeric_junction: junction,
                out_prefix: output_prefix,
                gtf_file: Geno.conf.get( 'REFERENCE', 'ref_gtf' ),
                environment_variables: Geno.conf.get( 'ENV', 'PERL5LIB' ),
                star_fusion: Geno.conf.get( 'SOFTWARE', 'STAR-Fusion' ),
                additional_params: Geno.job.get_param( 'star_fusion', 'additional_params' ),
                scriptdir: Geno.dir[ 'script' ] }
        shell_script_full_path = make_script_file( function_name, star_res.star_fusion, Geno, **arg)

        #
        # Run
        #
        return_code = Geno.RT.runtask(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ] )

        if return_code != 0:
            log.error( "{function}: runtask failed".format( function = function_name ) )
            raise

        save_status_of_this_process( function_name, output_prefix, return_code, Geno, use_subdir )

    except IOError as (errno, strerror):
        log.error( "{function}: I/O error({num}): {error}".format(
                        function = whoami(),
                        num = errno,
                        error = strerror)
                )
        return_value = False

    except ValueError:
        log.error( "{function}: ValueError".format(
                        function = whoami()
                    )
                )
        return_value = False

    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        log.error( "{function}: Unexpected error: {error}.".format(
                    function = whoami() , error = sys.exc_info()[0] ) )
        return_value = False


    return return_value

def fusionfusion(
        sam,
        junction,
        output_prefix
        ):
    """
        Stage 4: fusion_fusion

    """
    return_value = True

    try:
        function_name = whoami()
        log.info( "#{function}".format( function = function_name ) )

        #
        # Make shell script
        #
        args =  {
            'log': Geno.dir[ 'log' ],
            'pythonhome': Geno.conf.get( 'ENV', 'PYTHONHOME' ),
            'ld_library_path': Geno.conf.get( 'ENV', 'LD_LIBRARY_PATH'),
            'pythonpath': Geno.conf.get( 'ENV', 'PYTHONPATH' ),
            'chimeric_sam': sam,
            'output_prefix': output_prefix,
            # 'out_dir': output_dir,
            'fusion_fusion': Geno.conf.get( 'SOFTWARE', 'fusionfusion' ),
            'param_file': Geno.job.get_param( 'fusionfusion', 'param_file' ),
            'scriptdir': Geno.dir[ 'script' ]
        }

        shell_script_full_path = make_script_file( function_name, star_res.fusionfusion, Geno, **args )

        #
        # Run
        #
        return_code = Geno.RT.runtask(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ] )

        if return_code != 0:
            log.error( "{function}: runtask failed".format( function = function_name ) )
            raise

        save_status_of_this_process( function_name, output_prefix, return_code, Geno, use_subdir )

    except IOError as (errno, strerror):
        log.error( "{function}: I/O error({num}): {error}".format(
                        function = whoami(),
                        num = errno,
                        error = strerror)
                )
        return_value = False

    except ValueError:
        log.error( "{function}: ValueError".format(
                        function = whoami()
                    )
                )
        return_value = False

    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        log.error( "{function}: Unexpected error: {error}.".format(
                    function = whoami() , error = sys.exc_info()[0] ) )
        return_value = False


    return return_value

#####################################################################
#
#   STAGE 0 data preparation
#
Sample = Sample()
if Geno.input_file_list:
    Sample.set_sample_list( Geno.input_file_list )

#####################################################################
#
#   STAGE 1 fastq to bam by star
#
@active_if( 'star_genome' in Geno.job.get_job( 'tasks' )[ 'RNA' ] )
@files( generate_params_for_star_genome )
@check_if_uptodate( check_file_exists_for_star_genome )
def stage_1( fasta, gtf, output_dir ):
    if not star_genome( fasta, gtf, output_dir ) and Geno.options.stop_pipeline:
        raise Exception( 'star_genome failed.' )

#####################################################################
#
#   STAGE 2 fastq to bam by star
#
@follows( stage_1 )
@active_if( 'star' in Geno.job.get_job( 'tasks' )[ 'RNA' ] )
@files( generate_params_for_star )
@check_if_uptodate( check_file_exists_for_star )
def stage_2( input_file1, input_file2, output_prefix ):
    if not star( input_file1, input_file2, output_prefix ) and Geno.options.stop_pipeline:
        raise Exception( 'star failed.' )

#####################################################################
#
#   STAGE 3
#
@follows( stage_2 )
@active_if( 'star_fusion' in Geno.job.get_job( 'tasks' )[ 'RNA' ] )
@files( generate_params_for_star_fusion )
@check_if_uptodate( check_file_exists_for_star_fusion )
def stage_3( sam, junction, output_prefix ):
    if not star_fusion( sam, junction, output_prefix ) and Geno.options.stop_pipeline:
        raise Exception( 'star_fusion failed.' )

#####################################################################
#
#   STAGE 4
#
@follows( stage_2 )
@active_if( 'fusionfusion' in Geno.job.get_job( 'tasks' )[ 'RNA' ] )
@files( generate_params_for_fusionfusion )
@check_if_uptodate( check_file_exists_for_fusionfusion )
def stage_4( sam, junction, output_prefix ):
    if not fusionfusion( sam, junction, output_prefix ) and Geno.options.stop_pipeline:
        raise Exception( 'fusionfusion failed.' )

#####################################################################
#
#   LAST STAGE 
#

@follows( stage_4 )
def last_function():
    log.info( "Genomon pipline has finished successflly!" )
    return True
