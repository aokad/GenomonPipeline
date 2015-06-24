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
# Check file exits.
#
def check_file_exists( function_name, input_file, output_dir, output_file):
    if not os.path.exists(output_file):
        return True, "Missing file %s" % output_file
    else:
        exit_status = get_status_of_this_process( function_name, output_dir )
        in_time = os.path.getmtime( input_file )
        out_time = os.path.getmtime( output_file )
        if exit_status != 0 or in_time > out_time:
            return True, "Missing file %s" % output_file
        else:
            return False, "File %s exists" % output_file

def check_file_exists_for_tophat2( input_file, output_dir ):
    return check_file_exists( 'tophat2',
                              input_file,
                              output_dir,
                              output_dir + '/accepted_hits.bam' )

def check_file_exists_for_cufflinks( input_file, output_dir ):
    return check_file_exists( 'cufflinks',
                              input_file,
                              output_dir,
                              output_dir + '/genes.fpkm_tracking' )

def check_file_exists_for_cuffdiff( input_file_list1, input_file_list2, output_dir1, output_dir2 ):
    for input_file in input_file_list1 + input_file_list2:
        
        output_file = output_dir2 + '/gene_exp.diff'
        if check_file_exists( 'cuffdiff', input_file, output_dir2, output_file ):
            return True, "Missing file %s" % output_file

    return False, "File %s exists" % output_file

def check_file_exists_for_cummeRbund( input_dir, output_dir ):
    input_file = input_dir + '/gene_exp.diff'
    output_file = output_dir + '/' + os.path.split( output_dir )[ 1 ] + '_disp.png'
    return check_file_exists( 'cummeRbund',
                              input_file,
                              output_dir,
                              output_file )

#####################################################################
#
# generate parameters
#
# control_disease_pairs:
#   #Normal:     Disease
#   s_B1N:              s_B1T
#   "s_B2Na, s_B2Nb":   s_B2T       --> merge s_B2Na + s_B2Nb
#   s_B2Na:             s_B2T
#   s_B2Nb:             s_B2T       --> compare 3 files separately.
#   s_B3N: " s_B3Ta,s_B3Tb "        --> merge s_B3Ta + s_B3Tb
#   s_B3N:
#                       - s_B3Ta
#                       - s_B3Tb    --> compare 3 files separately.
#   Normal:  # Normal only
#                       - s_B4N
#                       - s_B5N
#   Disease: # Disease only
#                       - s_B6T
#                       - s_B7T
#
def generate_params_for_tophat2( ):
    global Sample
    Sample.make_param( 'tophat2', '', 'bam', 1, 1 )
    for infile1, infile2, outfile1, outfil2 in Sample.param( 'tophat2' ):
        yield ( infile1, outfile1 )

def generate_params_for_cufflinks( ):
    global Sample
    Sample.make_param( 'cufflinks', '', 'cufflinks', 1, 1 )
    for infile1, infile2, outdir1, outdir2 in Sample.param( 'cufflinks' ):
        yield ( infile1, outdir1 )

def generate_params_for_cuffdiff( ):
    global Sample
    Sample.make_param( 'cuffdiff', '', 'cuffdiff', 1, 1 )

    #
    # Make sure that 'control_disease_pairs' is defined,
    #            and 'sample_subdir' is defined.
    #
    ctrl_dis_pairs = Geno.job.get_job( 'control_disease_pairs' )
    if not Sample.subdir_exists() or ctrl_dis_pairs == None:
        raise

    #
    # Get sample_subdir list
    #
    list_id = 0
    data_dict = {}
    param_list = Sample.param( 'cuffdiff' )
    for infile1, infile2, outdir1, outdir2 in param_list:
        tmp_dir = os.path.basename( os.path.split( infile2 )[0] )
        data_dict[ tmp_dir ] = list_id
        list_id += 1

    #
    # Create parameters for comparison
    #
    for data_type in ctrl_dis_pairs.keys():
        if data_type != 'Normal' and data_type != 'Disease':
            normal_list = data_type.replace( ' ', '' ).split( ',' )
            disease_list = ctrl_dis_pairs[ data_type ].replace( ' ', '' ).split( ',' )
            normal_dir_list = []
            disease_dir_list = []
            for normal in normal_list:
                if normal in data_dict.keys():
                    normal_dir_list.append( param_list[ data_dict[ normal ] ][ 0 ] )
            for disease in disease_list:
                if disease in data_dict.keys():
                    disease_dir_list.append( param_list[ data_dict[ disease ] ][ 0 ] )
            normal_outdir = param_list[ data_dict[ normal_list[ 0 ] ] ][ 3 ]
            disease_outdir = param_list[ data_dict[ disease_list[ 0 ] ] ][ 3 ]
            yield( normal_dir_list,
                   disease_dir_list,
                   os.path.split( normal_outdir )[ 0 ],
                   os.path.split( disease_outdir )[ 0 ] )


def generate_params_for_cummeRbund( ):
    global Sample
    Sample.make_param( 'cummeRbund', '', 'cummeRbund', 1, 1 )

    #
    # Make sure that 'control_disease_pairs' is defined,
    #            and 'sample_subdir' is defined.
    #
    ctrl_dis_pairs = Geno.job.get_job( 'control_disease_pairs' )
    if not Sample.subdir_exists() or ctrl_dis_pairs == None:
        raise

    #
    # Get sample_subdir list
    #
    list_id = 0
    data_dict = {}
    param_list = Sample.param( 'cummeRbund' )
    for infile1, infile2, outdir1, outdir2 in param_list:
        tmp_dir = os.path.basename( os.path.split( infile2 )[0])
        data_dict[ tmp_dir ] = list_id
        list_id += 1

    #
    # Create parameters for comparison
    #
    for data_type in ctrl_dis_pairs.keys():
        if data_type != 'Normal' and data_type != 'Disease':
            normal_list = data_type.replace( ' ', '' ).split( ',' )
            disease_list = ctrl_dis_pairs[ data_type ].replace( ' ', '' ).split( ',' )

            for normal in normal_list:
                if normal in data_dict.keys():
                    normal_dir = os.path.split( param_list[ data_dict[ normal ] ][ 0 ] )[ 0 ]
                    break
            for disease in disease_list:
                if disease in data_dict.keys():
                    disease_dir = os.path.split( param_list[ data_dict[ disease ] ][ 0 ] )[ 0 ]
                    break

            normal_outdir = os.path.split( param_list[ data_dict[ normal_list[ 0 ] ] ][ 3 ] )[ 0 ]
            disease_outdir = os.path.split( param_list[ data_dict[ disease_list[ 0 ] ] ][ 3 ] )[ 0 ]
            yield( disease_dir, disease_outdir )

#####################################################################
#
#   STAGE 0 data preparation
#
Sample = Sample()

#####################################################################
#
#   STAGE 1 fastq to bam by tophat2
#
@active_if( 'tophat2' in Geno.job.get_job( 'tasks' )[ 'RNA' ] )
@files( generate_params_for_tophat2 )
@check_if_uptodate( check_file_exists_for_tophat2 )
def tophat2(
        input_file,
        output_dir,
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
                                        output_dir = output_dir,
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
        runtask_return_code = Geno.RT.runtask(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ] )

        if runtask_return_code != 0:
            log.error( "{function}: runtask failed".format( function = function_name ) )
            raise

        save_status_of_this_process( function_name, output_dir, runtask_return_code )

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
#   STAGE 2 estimate RNA expression for each sample by cufflinks
#
@follows( tophat2 )
@active_if( 'cufflinks' in Geno.job.get_job( 'tasks' )[ 'RNA' ] )
@files( generate_params_for_cufflinks )
@check_if_uptodate( check_file_exists_for_cufflinks )
def cufflinks(
        input_file,
        output_dir
        ):
    """
        Stage 2: cufflinks

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
                                        bam_file = input_file + '/accepted_hits.bam',
                                        output_dir = output_dir,
                                        ref_gtf = Geno.conf.get( 'REFERENCE', 'ref_gtf' ),
                                        cufflinks = Geno.conf.get( 'SOFTWARE', 'cufflinks' )
                                    )
                                )
        shell_script_file.close()

        #
        # Run
        #
        runtask_return_code = Geno.RT.runtask(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ] )

        if runtask_return_code != 0:
            log.error( "{function}: runtask failed".format( function = function_name ) )
            raise

        save_status_of_this_process( function_name, output_dir, runtask_return_code )

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
#   STAGE 3 estimate RNA expression and compare normal vs disease
#
@follows( cufflinks )
@active_if( 'cuffdiff' in Geno.job.get_job( 'tasks' )[ 'RNA' ] )
@files( generate_params_for_cuffdiff )
@check_if_uptodate( check_file_exists_for_cuffdiff )
def cuffdiff(
        control_input_dir_list,
        disease_input_dir_list,
        control_output_dir,
        disease_output_dir
        ):
    """
        Stage 3: cuffdiff

    """
    return_code = True

    try:
        function_name = whoami()
        log.info( "#{function}".format( function = function_name ) )

        #
        # Make shell script
        #
        env_variable_str = "LD_LIBRARY_PATH=$LD_LIBRARY_PATH:{libmaus_PATH}".format(
                                    libmaus_PATH = Geno.conf.get( 'ENV', 'libmaus_PATH' ) )

        if Geno.job.get_job( 'use_biobambam' ):
            biobambam_format = "I="
        else:
            biobambam_format = ""

        control_merge_bam_flag = ( len( control_input_dir_list ) > 1 )
        disease_merge_bam_flag = ( len( disease_input_dir_list ) > 1 )

        control_bam_file_list = ''
        for input_dir in control_input_dir_list:
            control_bam_file_list +=  "{biobambam_format}{input_dir}/accepted_hits.bam ".format(
                                        input_dir = input_dir,
                                        biobambam_format = biobambam_format )

        disease_bam_file_list = ''
        for input_dir in disease_input_dir_list:
            disease_bam_file_list +=  "{biobambam_format}{input_dir}/accepted_hits.bam  ".format(
                                        input_dir = input_dir,
                                        biobambam_format = biobambam_format )


        bam_file_list_array = "INPUT_FILE=(\n"
        output_dir_array    = "OUT_DIR=(\n"
        merge_bam_flag      = "FLAG=(\n"

        bam_file_list_array += " [1]=\"{file_list}\"\n".format( file_list = control_bam_file_list )
        output_dir_array    += " [1]=\"{output_dir}\"\n".format( output_dir = control_output_dir )
        merge_bam_flag      += " [1]=\"{flag}\"\n".format( flag = control_merge_bam_flag )

        bam_file_list_array += " [2]=\"{file_list}\"\n".format( file_list = disease_bam_file_list )
        output_dir_array    += " [2]=\"{output_dir}\"\n".format(  output_dir = disease_output_dir )
        merge_bam_flag      += " [2]=\"{flag}\"\n".format(  flag = disease_merge_bam_flag )

        bam_file_list_array += ")\n"
        output_dir_array    += ")\n"
        merge_bam_flag      += ")\n"

        #
        # Make shell script
        #
        shell_script_full_path = make_script_file_name( function_name + '_merge', Geno )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( rna_res.cuffdiff_merge.format(
                                        log = Geno.dir[ 'log' ],
                                        env_variables = env_variable_str,
                                        use_biobambam = Geno.job.get_job( 'use_biobambam' ),
                                        array_data = bam_file_list_array + output_dir_array + merge_bam_flag,
                                        input_bam_files = "${INPUT_FILE[$SGE_TASK_ID]}",
                                        merge_bam_flag = "${FLAG[$SGE_TASK_ID]}",
                                        merged_bam_file = '${OUT_DIR[$SGE_TASK_ID]}/merged_accepted_hits.bam',
                                        out_dir = '${OUT_DIR[$SGE_TASK_ID]}',
                                        biobambam = Geno.conf.get( 'SOFTWARE', 'biobambam' ),
                                        samtools = Geno.conf.get( 'SOFTWARE', 'samtools' )
                                    )
                                )

        shell_script_file.close()

        runtask_return_code = Geno.RT.run_arrayjob(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ],
                            id_start = 1,
                            id_end = 2 )

        if runtask_return_code != 0:
            log.error( "{function}: runtask failed".format( function = function_name ) )
            raise

        #
        # Make shell script
        #
        data_labels = control_output_dir.split( '/' )[ -2 ]
        data_labels += ','
        data_labels += disease_output_dir.split( '/' )[ -2 ]
        shell_script_full_path = make_script_file_name( function_name, Geno )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( rna_res.cuffdiff.format(
                                        log = Geno.dir[ 'log' ],
                                        env_variables = env_variable_str,
                                        merged_control_bam_file = control_output_dir + '/merged_accepted_hits.bam',
                                        merged_disease_bam_file = disease_output_dir + '/merged_accepted_hits.bam',
                                        output_dir = disease_output_dir,
                                        data_labels = data_labels,
                                        ref_gtf = Geno.conf.get( 'REFERENCE', 'ref_gtf' ),
                                        cuffdiff = Geno.conf.get( 'SOFTWARE', 'cuffdiff' )
                                    )
                                )
        shell_script_file.close()

        runtask_return_code = Geno.RT.runtask(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ] )

        if runtask_return_code != 0:
            log.error( "{function}: runtask failed".format( function = function_name ) )
            raise

        save_status_of_this_process( function_name, disease_output_dir, runtask_return_code )

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
#   STAGE 4 run cummeRbund and plot data of control vs disease
#
@follows( cuffdiff )
@active_if ( 'cummeRbund' in Geno.job.get_job( 'tasks' )[ 'RNA' ] )
@files ( generate_params_for_cummeRbund )
@check_if_uptodate( check_file_exists_for_cummeRbund )
def cummeRbund(
        input_dir,
        output_dir
        ):
    """
        Stage 4: cummeRbund

    """
    return_code = True

    try:
        function_name = whoami()
        log.info( "#{function}".format( function = function_name ) )

        env_variable_str = "R_HOME={r_path}\nR_LIBS={r_path}\nR_LIBS_USER={r_path}".format(
                                    r_path = Geno.conf.get( 'ENV', 'R_LIBS' ) )
        #
        # Make shell script
        #
        sample_name = os.path.basename( input_dir )
        shell_script_full_path = make_script_file_name( function_name, Geno )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( rna_res.cummeRbund.format(
                                        log = Geno.dir[ 'log' ],
                                        env_variables = env_variable_str,
                                        sample_name = sample_name,
                                        input_dir = input_dir,
                                        output_dir = output_dir,
                                        R = Geno.conf.get( 'SOFTWARE', 'R' ),
                                        script_dir = Geno.dir[ 'script' ]
                                    )
                                )
        shell_script_file.close()

        #
        # Run
        #
        runtask_return_code = Geno.RT.runtask(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ] )
        if runtask_return_code != 0:
            log.error( "{function}: runtask failed".format( function = function_name ) )
            raise

        save_status_of_this_process( function_name, output_dir, runtask_return_code )

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
