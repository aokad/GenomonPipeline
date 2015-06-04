"""
wes_pipeline.py

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
from resource import rna_resource as rna_res
from utils import *
from sample import Sample

def check_file_exists(input_file, output_file):
    if not os.path.exists(output_file):
        return True, "Missing file %s" % output_file
    else:
        return False, "File %s exists" % output_file

#####################################################################
#
#   STAGE 0 data preparation
#
Sample = Sample()
Sample.make_param( 'tophat2', '.bam', 'bam', 1, 1 )
starting_file_list = []
for infile1, infile2, outfile1, outfil2 in Sample.param( 'tophat2' ):
    starting_file_list.append( ( infile1, outfile1 ) )


#####################################################################
#
#   STAGE 1 fastq to bam by tophat2
#
@active_if ( 'tophat2' in Geno.job.get_job( 'tasks' )[ 'RNA' ] )
@parallel( starting_file_list )
@check_if_uptodate( check_file_exists )
def tophat2(
        input_file,
        output_file,
        ):
    """
        Stage 1: tophat2

    """
    return_code = True

    try:
        function_name = whoami()
        log.info( "#{function}".format( function = function_name ) )

        bowtie_bin = Geno.conf.get( 'SOFTWARE', 'bowtie2' )
        bowtie_path = os.path.split( bowtie_bin )[0]

        #
        # Make shell script
        #
        shell_script_full_path = make_script_file_name( function_name, Geno )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( rna_res.tophat2.format(
                                        log = Geno.dir[ 'log' ],
                                        ref_fa = Geno.conf.get( 'REFERENCE', 'ref_fasta' ),
                                        input_fastq = input_file,
                                        output_file = output_file,
                                        ref_gtf = Geno.conf.get( 'REFERENCE', 'ref_gtf' ),
                                        bowtie2_database = Geno.conf.get( 'REFERENCE', 'bowtie2_db' ),
                                        bowtie_path = bowtie_path,
                                        tophat2 = Geno.conf.get( 'SOFTWARE', 'tophat2' )
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

        Geno.status.save_status( function_name, input_file, return_code )

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
        log.error( "{function}: Unexpected error: {error} ".format(
                    function = whoami()))
        log.error("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )
        return_code = False


    return return_code

#####################################################################
#
#   STAGE 2
#
@follows( tophat2 )
@active_if ( 'cufflinks' in Geno.job.get_job( 'tasks' )[ 'RNA' ] )
@transform( tophat2, suffix( "2.txt" ), "3.txt" )
@check_if_uptodate( check_file_exists )
def cufflinks(
        input_file,
        output_file
        ):
    """
        Stage 2

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
        shell_script_file.write( rna_res.cufflinks.format(
                                        log = Geno.dir[ 'log' ],
                                        bam_file = input_file,
                                        output_dir = output_dir,
                                        ref_gtf = Geno.conf.get( 'REFERENCE', 'gtf' ),
                                        cufflinks = Geno.conf.get( 'SOFTWARE', 'cufflinks' )
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
                            id_end = rna_res.interval_num )

        if return_code != 0:
            log.error( "{function}: runtask failed".format( function = function_name ) )
            raise

        Geno.status.save_status( function_name, input_file, return_code )

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
        log.error( "{function}: Unexpected error: {error} ".format(
                    function = whoami()))
        return_code = False


    return return_code

#####################################################################
#
#   STAGE 3
#
@follows( cufflinks )
@active_if ( 'cummeRbund' in Geno.job.get_job( 'tasks' )[ 'RNA' ] )
@transform( cufflinks, suffix( "3.txt" ), "4.txt" )
@check_if_uptodate( check_file_exists )
def cummeRbund(
        input_file,
        output_file
        ):
    """
        Stage 3

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
        shell_script_file.write( rna_res.fisher_mutation_call.format(
                                        log = Geno.dir[ 'log' ],
                                        ref_fa = Geno.conf.get( 'REFERENCE', 'ref_fasta' ),
                                        input_fastq = input_file,
                                        output_bam = output_file,
                                        ref_gtf = Geno.conf.get( 'REFERENCE', 'gtf' ),
                                        bowtie2_db = Geno.conf.get( 'REFERENCE', 'bowtie2' ),
                                        tophat2 = Geno.conf.get( 'SOFTWARE', 'tophat2' ),
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
                            id_end = rna_res.interval_num )
        if return_code != 0:
            log.error( "{function}: runtask failed".format( function = function_name ) )
            raise

        Geno.status.save_status( function_name, input_file, return_code )

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
        log.error( "{function}: Unexpected error: {error} ".format(
                    function = whoami()))
        return_code = False


    return return_code


#####################################################################
#
#   LAST STAGE 
#

@follows( cummeRbund )
def last_function():
    log.info( "Genomon pipline has finished successflly!" )
    return True
