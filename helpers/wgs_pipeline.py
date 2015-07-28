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
# File update check
#

def check_file_exists_for_bam2fastq(
        input_file1,
        input_file2,
        output_file1,
        output_file2
        ):
    return check_file_exists_for_input_output(
                'bam2fastq',
                input_file1,
                input_file2,
                output_file1,
                output_file2)

def check_file_exists_for_cutadapt(
        input_file1,
        input_file2,
        output_file1,
        output_file2
        ):
    return check_file_exists_for_input_output(
                'cutadapt',
                input_file1,
                input_file2,
                output_file1,
                output_file2)

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

        exit_status = get_status_of_this_process( 'split_fastq', output_file1 )
        if exit_status != 0 and in_time > out_time:
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
    input_file_list,
    output_file
    ):
    if ( Geno.job.get_job( 'use_biobambam' ) and
         'markduplicates' in Geno.job.get_job( 'tasks' )[ 'WGS'] ):
        return_code = False, "BioBambam does not need merge."
    else:
        return_code = check_file_exists_for_merge(
            'merge_bam',
            input_file_list,
            output_file )

    return return_code

def check_file_exists_for_markduplicates(
    input_file_list,
    output_file
    ):
    return check_file_exists_for_merge(
        'markduplicates',
        input_file_list,
        output_file )

def check_file_exists_for_merge(
    process_name,
    input_file_list,
    output_file
    ):

    """
    Checks if output file exists for merge_bam

    """

    exit_status = get_status_of_this_process( process_name, output_file )

    if not os.path.exists( output_file ) or exit_status != 0:
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

    exit_status = get_status_of_this_process( 'bwa_mem', output_file1 )
    ( output_prefix, output_suffix ) = os.path.splitext( output_file1 )
    if Geno.job.get_job( 'use_biobambam' ):
        glob_filename = "{prefix}*_bamsorted{suffix}".format( prefix = output_prefix, suffix = output_suffix ) 
    else:
        glob_filename = "{prefix}*_sorted{suffix}".format( prefix = output_prefix, suffix = output_suffix ) 

    if exit_status != 0 or not glob( glob_filename ):
        return True, "Missing file {outputfile} for {inputfile}.".format(
                            outputfile = output_file1,
                            inputfile = input_file1 )
    else:
        return False, "File {output} exits for {input}.".format( output = output_file1, input = input_file1 )

def check_file_exists_for_input_output(
        process_name,
        input_file1,
        input_file2,
        output_file1,
        output_file2
    ):

    """
    Checks if output file exists

    """

    exit_status = get_status_of_this_process( process_name, output_file1 )
    if exit_status != 0 or not os.path.exists( output_file1 ):
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

def check_file_exists_for_fisher_mutation_call(
        control_file_list,
        disease_file_list,
        control_output_dir,
        output_file
    ):
    """
    Check if output file exists for fisher_mutation_call

    """
    exit_status = get_status_of_this_process( 'fisher_mutation_call', output_file )

    if exit_status != 0 or not os.path.exists( output_file ):
        return True, "Missing file {outputfile} for {inputfile}.".format(
                            outputfile = output_file,
                            inputfile = control_file_list[ 0 ] )

    else:
        in_time = os.path.getmtime( disease_file_list[ 0 ] )
        out_time = os.path.getmtime( output_file )
        if in_time > out_time:
            return True, "{output} is older than {input}.".format( output = output_file, input = control_file_list[ 0 ] )
        else:
            return False, "File {output} exits for {input}.".format( output = output_file, input = control_file_list[ 0 ] )


def check_file_exists_for_bam_stats(
        input_file1,
        input_file2,
        output_file
    ):

    """
    Checks if output file exists for a general input x 2 and output x 1 case

    """
    exit_status = get_status_of_this_process( 'bam_stats_calc', output_file ) + \
                  get_status_of_this_process( 'bam_stats_merge', output_file )

    summary_dir_name = os.path.dirname( output_file )
    out_sum_file =  summary_dir_name + '/' + Geno.job.get_job( 'sample_name' ) + '.txt'
    if exit_status != 0 or not os.path.exists( summary_dir_name + '/' + Geno.job.get_job( 'sample_name' ) + '.txt' ):
        return True, "Missing file {outputfile} for {inputfile}.".format(
                            outputfile = out_sum_file,
                            inputfile = input_file1 )

    else:
        in_time = os.path.getmtime( input_file1 )
        out_time = os.path.getmtime( out_sum_file )
        if in_time > out_time:
            return True, "{output} is older than {input}.".format( output = out_sum_file, input = input_file1 )
        else:
            return False, "File {output} exits for {input}.".format( output = out_sum_file, input = input_file1 )

def check_file_exists_for_itd_detection(
        control_file_list,
        tumor_file_list,
        ctrl_output_dir_list,
        tumor_output_dir_list
    ):
    """
    """

    in_file_list = control_file_list if control_file_list != None else [] + \
                tumor_file_list if tumor_file_list != None else []
    out_dir_list = ctrl_output_dir_list if ctrl_output_dir_list != None else [] + \
                tumor_output_dir_list if tumor_output_dir_list != None else []
    exit_status = get_status_of_this_process( 'itd_detection', out_dir_list[ 0 ] )

    for in_file, out_dir in zip( in_file_list,  out_dir_list ) :
        if exit_status != 0 or not os.path.exists( out_dir ):
            return True, "Missing file {outputfile} for inputfile.".format(
                                outputfile = out_dir,
                                inputfile = in_file )

        else:
            in_time = os.path.getmtime( in_file )
            out_time = os.path.getmtime( out_dir )
            if in_time > out_time:
                return True, "{output} is older than {input}.".format( output = out_dir, input = in_file )

    return False, "Output files  exits."

def check_file_exists_for_annotation(
        input_file,
        output_file
    ):
    """
    """
    exit_status = get_status_of_this_process( 'annotation', output_file )
    if Geno.job.get_param( 'annotation', 'use_table_annovar' ):
        output_tmp = glob( output_file + '.*_multianno.txt' )
        file_exists = output_tmp != []
        if file_exists:
            output_tmp = output_tmp[ 0 ]
    else:
        output_tmp = output_file + '.genome_summary.csv' 
        file_exists = os.path.exists( output_tmp )

    if exit_status != 0 or not file_exists:
        return True, "Missing file {outputfile} for {inputfile}.".format(
                            outputfile = output_tmp,
                            inputfile = input_file )

    else:
        in_time = os.path.getmtime( input_file )
        out_time = os.path.getmtime( output_tmp )
        if in_time > out_time:
            return True, "{output} is older than {input}.".format( output = output_tmp, input = input_file )
        else:
            return False, "File {output} exits for {input}.".format( output = output_tmp, input = input_file )

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

    Sample.make_param( 'bam2fastq', None, '.fastq', 'fastq', 1, 2 )
    for param in Sample.param( 'bam2fastq' ):
        yield param


#
# For STAGE 2 split_fastq
#
def generate_params_for_split_fastq():
    """
    Generate parameter list for spilt_fastq_files

    """
    
    Sample.make_param( 'split_fastq', None, '.fastq', 'fastq', 2, 2 )
    for param in Sample.param( 'split_fastq' ):
        yield param

#
# For STAGE 3 cutadapt
#
def generate_params_for_cutadapt():
    """
    Generate parameter list for cutadapt

    """
    
    Sample.make_param( 'cutadapt', None, '_cutadapt.fastq', 'fastq', 2, 2 )
    for param in Sample.param( 'cutadapt' ):
        yield param

#
# For STAGE 4 bwa_mem
#

    
def generate_params_for_bwa_mem():
    """
    Generate parameter list for bwa_mem

    """

    Sample.make_param( 'bwa_mem', None, '.bam', 'bam', 2, 1 ) 
    for param in Sample.param( 'bwa_mem' ):
        yield param

#
# For STAGE 5 merge_bam        
#        
def generate_params_for_merge_bam():
    """
    Generate parameter list for merge_bam

    """
    
    Sample.make_param( 'merge_bam', None, '.bam', 'bam', 1, 1 )
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
    
    Sample.make_param( 'markduplicates', None, '.bam', 'bam', 1, 1 )
    input_file_list = {}
    for param in Sample.param( 'markduplicates' ):
        dir_name = os.path.dirname( param[ 0 ] )
        if not ( dir_name in input_file_list ):
            input_file_list[ dir_name ] = []
        input_file_list[ dir_name ].append( param[ 0 ] )

    if Geno.job.get_job( 'use_biobambam' ):
        for dir_name in input_file_list.keys():
            return_list = [ input_file_list[ dir_name ], dir_name + '/' + Geno.job.get_job( 'sample_name' ) + '_markdup.bam' ]
            yield return_list
    else:
        for dir_name in input_file_list.keys():
            return_list = [ dir_name +'/' + Geno.job.get_job( 'sample_name' ) + '_merged.bam',
                            dir_name + '/' + Geno.job.get_job( 'sample_name' ) + '_markdup.bam' ]
            yield return_list

#
# For STAGE 7 bam_stats
#
def generate_params_for_bam_stats ():
    """
    Generate parameter list for bam_stats

    """
    
    Sample.make_param( 'bam_stats', None, '.txt', 'summary', 1, 1 )
    input_file_list = {}
    for param in Sample.param( 'bam_stats' ):
        dir_name = os.path.dirname( param[ 0 ] )
        if not ( dir_name in input_file_list ):
            if 'markduplicates' in Geno.job.get_job( 'tasks' )[ 'WGS' ]:
                input_bam =  dir_name + '/' + Geno.job.get_job( 'sample_name' ) + '_markdup.bam'
            elif 'merge_bam' in Geno.job.get_job( 'tasks' )[ 'WGS' ]:
                input_bam =  dir_name + '/' + Geno.job.get_job( 'sample_name' ) + '_merged.bam'
            else:
                input_bam =  param[ 0 ]

            input_file_list[ dir_name ] = input_bam
            summary_dir_name = os.path.dirname( param[ 2 ] )
            return_list = [ input_bam,
                            'None',
                            summary_dir_name + '/' + Geno.job.get_job( 'sample_name' ) + '.txt' ]
            yield return_list


#
# For STAGE 8 fisher_mutation_call
#

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
def make_bam_filename_for_markdup_result( filename ):

    dir_name = os.path.dirname( filename )
    if 'markduplicates' in Geno.job.get_job( 'tasks' )[ 'WGS' ]:
        input_bam =  dir_name + '/' + Geno.job.get_job( 'sample_name' ) + '_markdup.bam'
    elif 'merge_bam' in Geno.job.get_job( 'tasks' )[ 'WGS' ]:
        input_bam =  dir_name + '/' + Geno.job.get_job( 'sample_name' ) + '_merged.bam'
    else:
        input_bam =  filename

    return input_bam

def generate_params_for_fisher_mutation_call():
    """
    Generate parameter list for fisher_mutation_call

    """
    
    Sample.make_param( 'fisher_mutation_call', None, '.txt', 'mutation', 2, 1 )

    ctrl_dis_pairs = Geno.job.get_job( 'control_disease_pairs' )

    #
    # No comparison. Run fisher exact test on each sample.
    #
    if ctrl_dis_pairs == None:
        input_file_list = {}
        for param in Sample.param( 'fisher_mutation_call' ):
            dir_name = os.path.dirname( param[ 0 ] )
            if not ( dir_name in input_file_list ):
                #
                # Set input_bam file name
                #
                input_bam = make_bam_filename_for_markdup_result( param[ 0 ] )

                #
                # Make parameters
                #
                input_file_list[ dir_name ] = input_bam
                mutation_dir_name = os.path.dirname( param[ 2 ] )
                return_list = [
                                [ None ],
                                [ input_bam ],
                                None ,
                                mutation_dir_name + '/' + Geno.job.get_job( 'sample_name' ) + '.txt'
                              ]
                yield return_list

    #
    # Comparison. Run fisher exact test on tumor, nomral pair.
    #
    else:
        #
        # Get the list of subdirs
        #
        list_id = 0
        data_dict = {}
        param_list = Sample.param( 'fisher_mutation_call' )
        for infile1, infile2, outdir1, outdir2 in param_list:
            tmp_dir = os.path.basename( os.path.split( infile2 )[ 0 ] )
            data_dict[ tmp_dir ] = list_id
            list_id += 1

        #   
        # Create parameters for comparison
        #
        for data_type in ctrl_dis_pairs.keys():
            #
            # Normal only and Disease only cases are not going to be processed.
            #
            if data_type != 'Normal' and data_type != 'Disease':
                #
                # Make list
                #

                #
                # Normal
                # Always merge in the following case.
                # "s_B2Na, s_B2Nb":     s_B2T
                #
                # Does not have to merge.
                # s_B2Na:   s_B2T
                # s_B2Nb:   s_B2T
                #
                normal_list = data_type.replace( ' ', '' ).split( ',' )
                normal_dir_list = []
                normal_outfile = None
                for normal_sample in normal_list:
                    if normal_sample in data_dict.keys():
                        normal_input = make_bam_filename_for_markdup_result(
                                param_list[ data_dict[ normal_sample ] ][ 0 ] )
                        normal_dir_list.append( normal_input )
                        if normal_outfile == None:
                            mutation_dir_name = os.path.dirname( param_list[ data_dict[ normal_sample ] ][ 2 ] )
                            normal_outfile = mutation_dir_name + '/' + Geno.job.get_job( 'sample_name' ) + '.txt'

                # Disease
                # s_B2N:    - "s_B2Ta, s_B2Tb"
                #           - s_B2Tc
                if isinstance( ctrl_dis_pairs[ data_type ], list ):
                    disease_list = ctrl_dis_pairs[ data_type ]
                else:
                    disease_list = [ ctrl_dis_pairs[ data_type ] ]

                for disease_data in disease_list:
                    disease_merge_list = []
                    disease_outfile = None
                    if -1 != disease_data.find( ',' ):
                        for disease_sample in disease_data.replace( ' ', '' ).split( ',' ):
                            disease_input = make_bam_filename_for_markdup_result( 
                                    param_list[ data_dict[ disease_sample ] ][ 0 ] )
                            disease_merge_list.append( disease_input )
                            if disease_outfile == None:
                                mutation_dir_name = os.path.dirname( param_list[ data_dict[ disease_sample ] ][ 3 ] )
                                disease_outfile = mutation_dir_name + '/' + Geno.job.get_job( 'sample_name' ) + '.txt'
                    else:
                        disease_input = make_bam_filename_for_markdup_result( 
                                param_list[ data_dict[ disease_data ] ][ 0 ] )
                        disease_merge_list.append( disease_input )
                        mutation_dir_name = os.path.dirname( param_list[ data_dict[ disease_data ] ][ 3 ] )
                        disease_outfile =  mutation_dir_name + '/' + Geno.job.get_job( 'sample_name' ) + '.txt'
                    #
                    # Return parameters
                    #
                    yield( normal_dir_list,
                           disease_merge_list,
                           os.path.split( normal_outfile )[ 0 ],
                           disease_outfile )

def generate_params_for_itd_detection( ):
    """
        Generate parameters for itd_detection

    """
    
    param_return = False
    task_name_list = ( 'markduplicates', 'merge_bam', 'bwa_mem' )
    id = 0
    while( not param_return ):
        param_return = Sample.make_param( 'itd_detection', task_name_list[ id ], '.txt', 'itd', 1, 1 )
        id += 1

    ctrl_dis_pairs = Geno.job.get_job( 'control_disease_pairs' )

    #
    # No comparison. Run Genomon-ITDetector on each sample.
    #
    if ctrl_dis_pairs == None:
        #
        # Make parameters
        #
        input_file_list = []
        output_dir_list = []
        for param in Sample.param( 'itd_detection' ):
            input_file_list.append( make_bam_filename_for_markdup_result( param[ 0 ] ) )
            output_dir_list.append( os.path.split( param[ 2 ] )[ 0 ] + '/'+ Geno.job.get_job( 'sample_name' ) )

        return_list = [
                        None,
                        input_file_list,
                        None,
                        output_dir_list,
                      ]
        yield return_list

    #
    # Comparison. Run Genomon-ITDetector exact test on tumor, nomral pair.
    #
    else:
        #
        # Get the list of subdirs
        #
        list_id = 0
        data_dict = {}
        param_list = Sample.param( 'itd_detection' )
        for infile1, infile2, outdir1, outdir2 in param_list:
            tmp_dir = os.path.basename( os.path.split( infile2 )[0] )
            data_dict[ tmp_dir ] = list_id
            list_id += 1

        #   
        # Create parameters for comparison
        #
        normal_file_list = []
        disease_file_list = []
        normal_outfile_list = []
        disease_outfile_list = []
        for data_type in ctrl_dis_pairs.keys():
            #
            # Normal only and Disease only cases are not going to be processed.
            #
            if data_type == 'Normal':
                for normal_dir in ctrl_dis_pairs[ 'Normal' ]:
                    normal_input = make_bam_filename_for_markdup_result( 
                            param_list[ data_dict[ normal_dir ] ][ 0 ] )
                    normal_file_list.append( normal_input )

                    itd_dir_name = os.path.dirname( param_list[ data_dict[ normal_dir ] ][ 2 ] )
                    normal_outfile_list.append( itd_dir_name )

            elif data_type == 'Disease':
                for tumor_dir in ctrl_dis_pairs[ 'Disease' ]:
                    disease_input = make_bam_filename_for_markdup_result(
                            param_list[ data_dict[ tumor_dir ] ][ 0 ] )
                    disease_file_list.append( disease_input )

                    itd_dir_name = os.path.dirname( param_list[ data_dict[ tumor_dir ] ][ 2 ] )
                    disease_outfile_list.append( itd_dir_name )

            else:
                #
                # Make list
                #

                #
                # Normal
                # case 1) "s_B2Na, s_B2Nb":     s_B2T
                # case 2) s_B2Na:   s_B2T
                #
                normal_list = data_type.replace( ' ', '' ).split( ',' )
                for normal_sample in normal_list:
                    if normal_sample in data_dict.keys():
                        normal_input = make_bam_filename_for_markdup_result(
                                param_list[ data_dict[ normal_sample ] ][ 0 ] )
                        normal_file_list.append( normal_input )
                        itd_dir_name = os.path.dirname( param_list[ data_dict[ normal_sample ] ][ 2 ] )
                        normal_outfile_list.append( itd_dir_name )

                # Disease
                # s_B2N:    - "s_B2Ta, s_B2Tb"
                #           - s_B2Tc
                if isinstance( ctrl_dis_pairs[ data_type ], list ):
                    disease_list = ctrl_dis_pairs[ data_type ]
                else:
                    disease_list = [ ctrl_dis_pairs[ data_type ] ]

                for disease_data in disease_list:
                    if -1 != disease_data.find( ',' ):
                        for disease_sample in disease_data.replace( ' ', '' ).split( ',' ):
                            disease_input = make_bam_filename_for_markdup_result( 
                                    param_list[ data_dict[ disease_sample ] ][ 0 ] )
                            disease_file_list.append( disease_input )
                            itd_dir_name = os.path.dirname( param_list[ data_dict[ disease_sample ] ][ 3 ] )
                            disease_outfile_list.append( itd_dir_name )
                    else:
                        disease_input = make_bam_filename_for_markdup_result( 
                                param_list[ data_dict[ disease_data ] ][ 0 ] )
                        disease_file_list.append( disease_input )
                        itd_dir_name = os.path.dirname( param_list[ data_dict[ disease_data ] ][ 3 ] )
                        disease_outfile_list.append( itd_dir_name ) 
        #
        # Return parameters
        #
        yield( normal_file_list,
               disease_file_list,
               normal_outfile_list,
               disease_outfile_list )

def generate_params_for_annotation():
    """
        Generate parameters for ANNOVAR annotation

    """

    param_return = False
    Sample.make_param( 'annotation', 'fisher_mutation_call', '', 'annotation', 1, 1 )
    ctrl_dis_pairs = Geno.job.get_job( 'control_disease_pairs' )

    if ctrl_dis_pairs == None:
        for param in Sample.param( 'annotation' ):
            mutation_dir_name = os.path.dirname( param[ 0 ] )
            annotation_dir_name = os.path.dirname( param[ 2 ] )
            yield ( mutation_dir_name + '/' + Geno.job.get_job( 'sample_name' ) + '.txt',
                    annotation_dir_name + '/' + Geno.job.get_job( 'sample_name' ) ) 

    else:
        disease_subdir_list = []
        for subdir in ctrl_dis_pairs.keys():
            if subdir == 'Disease':
                disease_subdir_list += ctrl_dis_pairs[ 'Disease' ]
            elif subdir  != 'Normal':
                if isinstance( ctrl_dis_pairs[ subdir ], list ):
                    disease_subdir_list += ctrl_dis_pairs[ subdir ]
                elif -1 == ctrl_dis_pairs[ subdir ].find( ',' ):
                    disease_subdir_list.append( ctrl_dis_pairs[ subdir ] )
                else:
                    disease_subdir_list += ctrl_dis_pairs[ subdir ].replace( ' ', '' ).split( ',' )

        for param in Sample.param( 'annotation' ):
            mutation_dir_name = os.path.dirname( param[ 0 ] )
            annotation_dir_name = os.path.dirname( param[ 3 ] )
            if os.path.basename( mutation_dir_name ) in disease_subdir_list:
                yield ( mutation_dir_name + '/' + Geno.job.get_job( 'sample_name' ) + '.txt',
                        annotation_dir_name + '/' + Geno.job.get_job( 'sample_name' ) ) 



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
        runtask_return_code = Geno.RT.runtask(
                        shell_script_full_path,
                        Geno.job.get_job( 'cmd_options' )[ function_name ] )

        if runtask_return_code != 0:
            with log_mutex:
                log.error( "{function}: runtask failed".format( function = function_name ) )
            raise

        save_status_of_this_process( function_name, output_file1, runtask_return_code )

    except IOError as (errno, strerror):
        with log_mutex:
            log.error( "{function}: I/O error({num}): {error}".format(num = errno, error = strerror, function = whoami() ) )
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
        runtask_return_code = Geno.RT.run_arrayjob(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ],
                            id_start = 1,
                            id_end = id )
        if runtask_return_code != 0:
            with log_mutex:
                log.error( "{function}: runtask failed".format( function = function_name ) )
            raise

        save_status_of_this_process( function_name, output_file1, runtask_return_code )

    except IOError as (errno, strerror):
        with log_mutex:
            log.error( "{function}: I/O error({num}): {error}".format(num = errno, error = strerror, function = whoami() ) )
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
        runtask_return_code = Geno.RT.run_arrayjob(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ],
                            id_start = 1,
                            id_end = 2 )
        if runtask_return_code != 0:
            with log_mutex:
                log.error( "{function}: runtask failed".format( function = function_name ) )
            raise

        save_status_of_this_process( function_name, output_file1, runtask_return_code )

    except IOError as (errno, strerror):
        with log_mutex:
            log.error( "{function}: I/O error({num}): {error}".format(function = whoami(), num = errno, error = strerror) )
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

        env_variable_str = "LD_LIBRARY_PATH=$LD_LIBRARY_PATH:{libmaus_PATH}".format(
                                        libmaus_PATH = Geno.conf.get( 'ENV', 'libmaus_PATH' ) )
        shell_script_file.write( bwa_mem_resource.format(
                                        log = Geno.dir[ 'log' ],
                                        array_data = str1 + str2 + str3,
                                        fastq1 = "${FILE1[$SGE_TASK_ID]}",
                                        fastq2 = "${FILE2[$SGE_TASK_ID]}",
                                        bam = "${FILE3[$SGE_TASK_ID]}",
                                        read_group = Geno.job.get_job(  'bam_read_group' ),
                                        min_score = Geno.job.get_param( 'bwa_mem', 'min_score' ),
                                        env_variables = env_variable_str,
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
        runtask_return_code = Geno.RT.run_arrayjob(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ],
                            id_start = 1,
                            id_end = id )
        if runtask_return_code != 0:
            with log_mutex:
                log.error( "{function}: runtask failed".format( function = function_name ) )
            raise
            
        save_status_of_this_process( function_name, output_file1, runtask_return_code )

    except IOError as (errno, strerror):
        with log_mutex:
            log.error( "{function}: I/O error({num}): {error}".format(function = whoami(), num = errno, error = strerror) )
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

        env_variable_str = "LD_LIBRARY_PATH=$LD_LIBRARY_PATH:{libmaus_PATH}".format(
                                libmaus_PATH = Geno.conf.get( 'ENV', 'libmaus_PATH' ) )

        shell_script_full_path = make_script_file_name( function_name, Geno )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( bam_merge_resource.format(
                                        log = Geno.dir[ 'log' ],
                                        input_bam_files = input_files,
                                        output_bam_file = output_file,
                                        env_variables = env_variable_str,
                                        samtools = Geno.conf.get( 'SOFTWARE', 'samtools' ),
                                        biobambam = Geno.conf.get( 'SOFTWARE', 'biobambam' )
                                        ) )
        shell_script_file.close()

        #
        # Run
        #
        runtask_return_code = Geno.RT.runtask(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ] )
        if runtask_return_code != 0:
            with log_mutex:
                log.error( "{function}: runtask failed".format( function = function_name ) )
            raise

        save_status_of_this_process( function_name, output_file, runtask_return_code )

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

            env_variable_str = "LD_LIBRARY_PATH=$LD_LIBRARY_PATH:{libmaus_PATH}".format(
                                    libmaus_PATH = Geno.conf.get( 'ENV', 'libmaus_PATH' ) )

            shell_script_full_path = make_script_file_name( function_name, Geno )
            shell_script_file = open( shell_script_full_path, 'w' )
            shell_script_file.write( wgs_res.biobambam_markduplicates.format(
                                            log = Geno.dir[ 'log' ],
                                            input_bam_files = input_files,
                                            output_bam = output_file,
                                            env_variables = env_variable_str,
                                            biobambam = Geno.conf.get( 'SOFTWARE', 'biobambam' )
                                        )
                                    )
            shell_script_file.close()

        else:
            #
            # Use Picard MarkDuplicates for samtools merge result
            #
            java_memory = Geno.job.get_param( 'markduplicates', 'java_memory' )

            shell_script_full_path = make_script_file_name( function_name, Geno )
            shell_script_file = open( shell_script_full_path, 'w' )
            shell_script_file.write( wgs_res.markduplicates.format(
                                            log = Geno.dir[ 'log' ],
                                            input_bam = input_file_list,
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
        runtask_return_code = Geno.RT.runtask(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ] )
        if runtask_return_code != 0:
            with log_mutex:
                log.error( "{function}: runtask failed".format( function = function_name ) )
            raise

        save_status_of_this_process( function_name, output_file, runtask_return_code )

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

#
# Stage 7: bam_stats
#
def bam_stats(
    control_input_file,
    disease_input_file,
    output_file,
    ):
    """
       BAM file statistics calculation

    """

    try:
        function_name = whoami()
        with log_mutex:
            log.info( "#{function}".format( function = function_name ) )

        #
        # Make shell script
        #
        bam_file_list = []
        if control_input_file != 'None':
            bam_file_list.append( control_input_file )
        if disease_input_file != 'None':
            bam_file_list.append( disease_input_file )

        if Geno.conf.get( 'REFERENCE', 'chr_str_in_fa' ) == 'True':
            chr_str_flag = '--chr_str'
        else:
            chr_str_flag = ''

        for bam_file in bam_file_list:
            #
            # Calculate bam file statistics
            #
            shell_script_name = function_name + '_calc'
            shell_script_full_path = make_script_file_name( shell_script_name , Geno )
            shell_script_file = open( shell_script_full_path, 'w' )
            shell_script_file.write( wgs_res.bam_stats_calc.format(
                                            log = Geno.dir[ 'log' ],
                                            bam_file = bam_file,
                                            output_txt = output_file,
                                            coverage = Geno.job.get_param( 'bam_stats', 'coverage' ),
                                            bed_file = Geno.job.get_param( 'bam_stats', 'bed_file' ),
                                            ref_fa = Geno.conf.get( 'REFERENCE', 'ref_fasta' ),
                                            genome_size = Geno.conf.get( 'REFERENCE', 'genome_size' ),
                                            pcap = Geno.conf.get( 'SOFTWARE', 'PCAP' ),
                                            samtools = Geno.conf.get( 'SOFTWARE', 'samtools' ),
                                            python = Geno.conf.get( 'SOFTWARE', 'python' ),
                                            chr_str_in_fa = chr_str_flag,
                                            script_dir = Geno.dir[ 'script' ]
                                        )
                                    )
            shell_script_file.close()

            calc_return_code = Geno.RT.run_arrayjob(
                                shell_script_full_path,
                                Geno.job.get_job( 'cmd_options' )[ function_name ],
                                id_start = 1,
                                id_end = 14 )

            if calc_return_code == 0:
                save_status_of_this_process( shell_script_name, output_file, calc_return_code )

            #
            # Merge results
            #
            shell_script_name = function_name + '_merge'
            shell_script_full_path = make_script_file_name( shell_script_name, Geno )
            shell_script_file = open( shell_script_full_path, 'w' )
            shell_script_file.write( wgs_res.bam_stats_merge.format(
                                            log = Geno.dir[ 'log' ],
                                            bam_file = bam_file,
                                            output_txt = output_file,
                                            python = Geno.conf.get( 'SOFTWARE', 'python' ),
                                            script_dir = Geno.dir[ 'script' ]
                                        )
                                    )
            shell_script_file.close()

            merge_return_code = Geno.RT.run_arrayjob(
                                shell_script_full_path,
                                Geno.job.get_job( 'cmd_options' )[ function_name ])

            #
            # Check return code
            #
            #  This is a calculation of statistics. Do not raise error even if it fails.
            if calc_return_code != 0:
                with log_mutex:
                    log.error( "{function}: runtask failed".format( function = function_name + '_calc' ) )
                raise
            if merge_return_code != 0:
                with log_mutex:
                    log.error( "{function}: runtask failed".format( function = function_name + '_merge' ) )
                raise

            if merge_return_code == 0:
                save_status_of_this_process( shell_script_name, output_file, merge_return_code )

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

#
# Stage 8: fisher_mutation_call
#
def fisher_mutation_call(
    control_input_file_list,
    disease_input_file_list,
    control_output_dir,
    disease_output_file,
    ):
    """
       Mutaion calling

    """

    try:
        function_name = whoami()
        with log_mutex:
            log.info( "#{function}".format( function = function_name ) )

        #
        # Merge bam files of multiple control sets or disease sets
        #

        #
        # There are control and tumor data to compare in fisher test.
        #
        if control_output_dir != None:
            #
            # Make shell script for merge bam files
            #
            env_variable_str = "LD_LIBRARY_PATH=$LD_LIBRARY_PATH:{libmaus_PATH}".format(
                                        libmaus_PATH = Geno.conf.get( 'ENV', 'libmaus_PATH' ) )

            control_merge_bam_flag = ( len( control_input_file_list ) > 1 )
            disease_merge_bam_flag = ( len( disease_input_file_list ) > 1 )

            if control_merge_bam_flag:
                samtools_control_bam_file_list = ' '.join( control_input_file_list )
                bambam_control_bam_file_list = 'I=' + ' I='.join( control_input_file_list )
            else:
                samtools_control_bam_file_list = control_input_file_list[ 0 ]
                bambam_control_bam_file_list = 'I=' + control_input_file_list[ 0 ]

            if disease_merge_bam_flag:
                samtools_disease_bam_file_list = ' '.join( disease_input_file_list )
                bambam_disease_bam_file_list = 'I=' + ' I='.join( disease_input_file_list )
            else:
                samtools_disease_bam_file_list = disease_input_file_list[ 0 ]
                bambam_disease_bam_file_list = 'I=' + disease_input_file_list[ 0 ]

            disease_output_dir = os.path.split( disease_output_file )[ 0 ]

            samtools_bam_file_list_array = "SAMTOOLS_INPUT_FILE=(\n"
            bambam_bam_file_list_array = "BAMBAM_INPUT_FILE=(\n"
            output_dir_array    = "OUT_FILE=(\n"
            merge_bam_flag      = "FLAG=(\n"

            samtools_bam_file_list_array += " [1]=\"{0}\"\n".format( samtools_control_bam_file_list )
            bambam_bam_file_list_array += " [1]=\"{0}\"\n".format( bambam_control_bam_file_list )
            output_dir_array    += " [1]=\"{0}\"\n".format( control_output_dir )
            merge_bam_flag      += " [1]=\"{0}\"\n".format( control_merge_bam_flag )

            samtools_bam_file_list_array += " [2]=\"{0}\"\n".format( samtools_disease_bam_file_list )
            bambam_bam_file_list_array += " [2]=\"{0}\"\n".format( bambam_disease_bam_file_list )
            output_dir_array    += " [2]=\"{0}\"\n".format(  disease_output_dir )
            merge_bam_flag      += " [2]=\"{0}\"\n".format(  disease_merge_bam_flag )

            samtools_bam_file_list_array += ")\n"
            bambam_bam_file_list_array += ")\n"
            output_dir_array    += ")\n"
            merge_bam_flag      += ")\n"

            #
            # Make shell script
            #
            shell_script_full_path = make_script_file_name( function_name + '_merge_bam', Geno )
            shell_script_file = open( shell_script_full_path, 'w' )
            shell_script_file.write( wgs_res.fisher_merge_bams.format(
                                            log = Geno.dir[ 'log' ],
                                            env_variables = env_variable_str,
                                            use_biobambam = Geno.job.get_job( 'use_biobambam' ),
                                            array_data = samtools_bam_file_list_array + bambam_bam_file_list_array +\
                                                         output_dir_array + merge_bam_flag,
                                            bambam_input_bam_files = "${BAMBAM_INPUT_FILE[$SGE_TASK_ID]}",
                                            samtools_input_bam_files = "${SAMTOOLS_INPUT_FILE[$SGE_TASK_ID]}",
                                            input_bam_file = "${SAMTOOLS_INPUT_FILE[$SGE_TASK_ID]}",
                                            merge_bam_flag = "${FLAG[$SGE_TASK_ID]}",
                                            merged_bam_file = '${{OUT_FILE[$SGE_TASK_ID]}}/{sample_name}.bam'.format(
                                                                sample_name = Geno.job.get_job( 'sample_name' ) ),
                                            out_dir = '${OUT_FILE[$SGE_TASK_ID]}',
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

            control_input_bam = control_output_dir + '/{sample_name}.bam'.format(
                                    sample_name = Geno.job.get_job( 'sample_name' ) )
            disease_input_bam = disease_output_dir + '/{sample_name}.bam'.format(
                                    sample_name = Geno.job.get_job( 'sample_name' ) )
        else:
        #
        # Not need to compare control and tumor in fisher test.
        # Run beta-distrubtion on allele frequencies in disease data.
        #
            control_input_bam = 'None'
            disease_input_bam = disease_input_file_list[ 0 ]


        #
        # Mutation call
        #
        #
        # Make shell script
        #
        shell_script_full_path = make_script_file_name( function_name, Geno )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( wgs_res.fisher_mutation_call.format(
                                        log = Geno.dir[ 'log' ],
                                        ref_fa = Geno.conf.get( 'REFERENCE', 'ref_fasta' ),
                                        control_input_bam = control_input_bam,
                                        disease_input_bam = disease_input_bam,
                                        output_txt = disease_output_file,
                                        remove_intermediate = True,
                                        max_indel = Geno.job.get_param( 'fisher_mutation_call', 'max_indel' ),
                                        max_distance = Geno.job.get_param( 'fisher_mutation_call', 'max_distance' ),
                                        base_quality = Geno.job.get_param( 'fisher_mutation_call', 'base_quality' ),
                                        map_quality = Geno.job.get_param( 'fisher_mutation_call', 'map_quality' ),
                                        mismatch_rate = Geno.job.get_param( 'fisher_mutation_call', 'mismatch_rate' ),
                                        min_depth = Geno.job.get_param( 'fisher_mutation_call', 'min_depth' ),
                                        samtools = Geno.conf.get( 'SOFTWARE', 'samtools' ),
                                        python = Geno.conf.get( 'SOFTWARE', 'python' ),
                                        script_dir = Geno.dir[ 'script' ]
                                 ) )
        shell_script_file.close()

        #
        # Run
        #
        runtask_return_code = Geno.RT.run_arrayjob(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ],
                            id_start = 1,
                            id_end = wgs_res.interval_num )

        if runtask_return_code != 0:
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
                                        output_txt = disease_output_file ) )
        shell_script_file.close()

        #
        # Run
        #
        runtask_return_code = Geno.RT.runtask(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ] )

        if runtask_return_code != 0:
            with log_mutex:
                log.error( "{function}: runtask failed".format( function = 'merge_fisher_result' ) )
            raise

        save_status_of_this_process( function_name, disease_output_file, runtask_return_code )

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

#
# Stage 9: itd_detection
#
def itd_detection(
    control_file_list,
    tumor_file_list,
    ctrl_output_file_list,
    tumor_output_file_list,
    ):
    """
        Genomon-ITDetector

    """

    try:
        function_name = whoami()
        with log_mutex:
            log.info( "#{function}".format( function = function_name ) )

        #
        # Preparation
        #
        # 1) create_ctrl_panel
        #   If create_ctrl_panel is set True,
        #       ctrl_panel is created from normal samples.
        #   If ctrl_panel_normal is set True and ctrl_panel_normal has a list of samples,
        #       ctrl_panel is created from the list of samples.
        #
        # 2) itd_ctrl_panel_files
        #   1. Get the file names.
        #   2. Cat the all files with extension 'inhouse_itd_list'
        #       into config/genomon_20150713_0927_049401_inhouse_itd.list
        #   4. soft-link the file to normal_inhouse_itd.list
        #   5. Cat the all files with extension 'inhouse_breakpoint_list'
        #       into config/genomon_20150713_0927_049401_inhouse_itd.list
        #   6. soft-link the file to normal_inhouse_breakpoint.list
        #   7. Pass the file name to detectITD.sh
        #       > detectITD.sh [input bam file] [output dir] [sample name] [inhouse file dir]
        #
        #   inhouse file dir:   normal_inhouse_itd.list
        #                       normal_inhouse_breakpoint.list
        #
        # 3) ctrl_panel_normal
        #   Get the list of normal samples.
        #   Run itd detection on the samples first.
        #   Get the inhouse data and add the list to normal_inhouse_itd.list and normal_inhouse_breakpint.list
        #
        normal_itd_file = Geno.dir[ 'config' ] + '/normal_inhouse_itd.list'
        normal_bp_file = Geno.dir[ 'config' ] + '/normal_inhouse_breakpoint.list'

        #
        # Make data for array job for create_ctrl_panel is True.
        #
        # First, look for 'ctrl_panel_normal'.
        # If there is no 'ctrl_panel_normal', use normal as normal samples to create ctrl_panel.
        in_file_list = []
        out_dir_list = []
        if Geno.job.get_param( 'itd_detection', 'create_ctrl_panel' ):
            ctrl_panel_normal_list = Geno.job.get_param( 'itd_detection', 'ctrl_panel_normal' )
            if ctrl_panel_normal_list:
                for ctrl_panel_normal in ctrl_panel_normal_list:
                    for i,x in [ (i,x) for i,x in enumerate( control_file_list ) if x.find( ctrl_panel_normal ) != -1 ]:
                        in_file_list.append( x )
                        out_dir_list.append( ctrl_output_file_list[ i ] )

                    for i, x in [ (i,x) for i,x in enumerate( tumor_file_list ) if x.find( ctrl_panel_normal ) != -1 ]:
                        in_file_list.append( x )
                        out_dir_list.append( tumor_output_file_list[ i ] )

            else:
                in_file_list = control_file_list if control_file_list != None else []
                out_dir_list = ctrl_output_file_list if ctrl_output_file_list != None else []

            #
            # Run itd_detection to create ctrl_panel
            #
            input_files = "FILE1=(\n"
            output_files = "FILE2=(\n"
            name_list = "NAME=(\n"
            id = 0
            for input_file, output_file in zip( in_file_list, out_dir_list ):
                id += 1
                name = os.path.basename( input_file )[ :-4 ]
                input_files += "[{id}]=\"{file}\"\n".format( id = id, file = input_file )
                output_files += "[{id}]=\"{file}\"\n".format( id = id, file = output_file )
                name_list += "[{id}]=\"{name}\"\n".format( id = id, name = name )
            input_files += ")\n"
            output_files += ")\n"
            name_list += ")\n"
            
            #
            # Make shell script for array job
            #
            shell_script_full_path = make_script_file_name( function_name, Geno )
            shell_script_file = open( shell_script_full_path, 'w' )
            shell_script_file.write( wgs_res.itd_detection.format(
                                            log = Geno.dir[ 'log' ],
                                            array = input_files + output_files + name_list,
                                            itd_inhouse_files = ' ',
                                            bam_file = "${FILE1[$SGE_TASK_ID]}",
                                            output_file = "${FILE2[$SGE_TASK_ID]}",
                                            name = "${NAME[$SGE_TASK_ID]}",
                                            itd_inhouse_dir = Geno.dir[ 'config' ],
                                            itd_detector = Geno.conf.get( 'SOFTWARE', 'itd_detector' ) ) )
            shell_script_file.close()

            #
            # Run
            #
            runtask_return_code = Geno.RT.run_arrayjob(
                                shell_script_full_path,
                                Geno.job.get_job( 'cmd_options' )[ function_name ],
                                id_start = 1,
                                id_end = id )

            if runtask_return_code != 0:
                with log_mutex:
                    log.error( "{function}: runtask failed".format( function = function_name ) )
                raise

        #
        # Make control panel file
        #

        #
        # Gather the inhouse files created from the firs itd_detection above.
        #
        inhouse_bp_list = []
        inhouse_itd_list = []
        for out_dir in out_dir_list:
            inhouse_bp_list.append( out_dir + '/inhouse_breakpoint.tsv' )
            inhouse_itd_list.append( out_dir + '/inhouse_itd.tsv' )

        #
        # Create inhouse file for detectITD.sh from 'itd_ctrl_panel_files' + on-fly-created file above.
        #
        timestamp = res.timestamp_format.format(
                                            year=Geno.now.year,
                                            month=Geno.now.month,
                                            day=Geno.now.day,
                                            hour=Geno.now.hour,
                                            min=Geno.now.minute,
                                            msecond=Geno.now.microsecond )
        itd_dest_file = "genomon_{timestamp}_itd.list".format(
                                    timestamp = timestamp )
        bp_dest_file = "genomon_{timestamp}_breakpoint.list".format(
                                    timestamp = timestamp )
        itd_dest = open( Geno.dir[ 'config' ] + '/' + itd_dest_file, 'w' )
        bp_dest = open( Geno.dir[ 'config' ] + '/' + bp_dest_file, 'w' )

        ctrl_panel_files = Geno.job.get_param( 'itd_detection', 'itd_ctrl_panel_files' )
        if len( ctrl_panel_files ) > 0 :
            for filename in ctrl_panel_files:
                if os.path.exists( filename ):
                    if filename[ -9: ] == '_itd.list':
                        shutil.copyfileobj( open( filename, 'r' ), itd_dest )
                    elif filename[ -16: ] == '_breakpoint.list':
                        shutil.copyfileobj( open( filename, 'r' ), bp_dest )

        for filename in inhouse_bp_list:
            if os.path.exists( filename ):
                shutil.copyfileobj( open( filename, 'r' ), bp_dest )
        for filename in inhouse_itd_list:
            if os.path.exists( filename ):
                shutil.copyfileobj( open( filename, 'r' ), itd_dest )

        itd_dest.close()
        bp_dest.close()

        if os.path.exists( normal_itd_file ):
            os.remove( normal_itd_file )
        if os.path.exists( normal_bp_file ):
            os.remove( normal_bp_file )
        os.symlink( itd_dest_file, normal_itd_file )
        os.symlink( bp_dest_file, normal_bp_file )

        #
        # Run the main itd detection
        #

        #
        # Remove ctrl_panel_normal from the list in sample_subdir
        #
        data_list = []
        if Geno.job.get_param( 'itd_detection', 'create_ctrl_panel' ) and ctrl_panel_normal_list:
            for input_file, output_dir in zip( control_file_list + tumor_file_list ,
                                                   ctrl_output_file_list + tumor_output_file_list ):
                for ctrl_panel_normal in ctrl_panel_normal_list:
                    if input_file.find( ctrl_panel_normal ) != -1:
                        data_list.append( ( input_file, output_dir ) )
        elif Geno.job.get_param( 'itd_detection', 'create_ctrl_panel' ):
            for input_file, output_dir in zip( tumor_file_list, tumor_output_file_list ):
                data_list.append( (input_file, output_dir ) )
        else:
            input_file_list = control_file_list if control_file_list != None else [] +\
                              tumor_file_list if tumor_file_list != None else [] 
            output_file_list = ctrl_output_file_list if ctrl_output_file_list != None else [] +\
                               tumor_output_file_list if tumor_output_file_list != None else [] 
            for input_file, output_dir in zip( input_file_list, output_file_list ):
                data_list.append( (input_file, output_dir ) )


        input_files = "FILE1=(\n"
        output_files = "FILE2=(\n"
        name_list = "NAME=(\n"
        id = 0
        for input_file, output_dir in data_list:
            id += 1
            name = os.path.basename( input_file )[ :-4 ]
            input_files += "[{id}]=\"{file}\"\n".format( id = id, file = input_file )
            output_files += "[{id}]=\"{file}\"\n".format( id = id, file = output_dir )
            name_list += "[{id}]=\"{name}\"\n".format( id = id, name = name )
        input_files += ")\n"
        output_files += ")\n"
        name_list += ")\n"
        
        #
        # Make shell script for array job
        #
        shell_script_full_path = make_script_file_name( function_name, Geno )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( wgs_res.itd_detection.format(
                                        log = Geno.dir[ 'log' ],
                                        array = input_files + output_files + name_list,
                                        itd_inhouse_files = itd_dest_file + ',' + bp_dest_file,
                                        bam_file = "${FILE1[$SGE_TASK_ID]}",
                                        output_file = "${FILE2[$SGE_TASK_ID]}",
                                        name = "${NAME[$SGE_TASK_ID]}",
                                        itd_inhouse_dir = Geno.dir[ 'config' ],
                                        itd_detector = Geno.conf.get( 'SOFTWARE', 'itd_detector' ) ) )
        shell_script_file.close()

        #
        # Run
        #
        runtask_return_code = Geno.RT.run_arrayjob(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ],
                            id_start = 1,
                            id_end = id )

        if runtask_return_code != 0:
            with log_mutex:
                log.error( "{function}: runtask failed".format( function = function_name ) )
            raise

        save_status_of_this_process( function_name, data_list[ 0 ][ 1 ], runtask_return_code )

    except IOError as (errno, strerror):
        with log_mutex:
            log.error( "{function}: I/O error({num}): {error}".format(function = whoami(), num = errno, error = strerror) )
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

#
# Stage 10: annotation
#
def annotation(
    input_file,
    output_file,
    ):
    """

    """

    try:
        function_name = whoami()
        with log_mutex:
            log.info( "#{function}".format( function = function_name ) )

        output_vcf = input_file[ :input_file.find( '.txt' ) ] + '.vcf'

        #
        # Make data for array job 
        #
        shell_script_full_path = make_script_file_name( function_name, Geno )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( wgs_res.annotation.format(
                                        log = Geno.dir[ 'log' ],
                                        input_file = input_file,
                                        output_prefix = output_file,
                                        output_vcf = output_vcf,
                                        output_in_vcf = Geno.job.get_param( 'annotation', 'output_in_vcf' ),
                                        use_table_annovar = Geno.job.get_param( 'annotation', 'use_table_annovar' ),
                                        summarize_annovar_params = Geno.job.get_param( 'annotation', 'summarize_annovar_params' ),
                                        table_annovar_params = Geno.job.get_param( 'annotation', 'table_annovar_params' ),
                                        ref_fa = Geno.conf.get( 'REFERENCE', 'ref_fasta' ),
                                        annovar = Geno.conf.get( 'SOFTWARE', 'annovar' ),
                                        scriptdir = Geno.dir[ 'script' ],
                                        python = Geno.conf.get( 'SOFTWARE', 'python' ),
                                        ) )
        shell_script_file.close()

        #
        # Run
        #
        runtask_return_code = Geno.RT.runtask(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ] )

        if runtask_return_code != 0:
            with log_mutex:
                log.error( "{function}: runtask failed".format( function = function_name ) )
            raise

        save_status_of_this_process( function_name, output_file, runtask_return_code )

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
@check_if_uptodate( check_file_exists_for_bam2fastq )
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
@check_if_uptodate( check_file_exists_for_cutadapt )
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
@check_if_uptodate( check_file_exists_for_merge_bam)
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
@check_if_uptodate( check_file_exists_for_markduplicates )
def stage_6( input_file_list, output_file ):
    return_value =  markduplicates( input_file_list,
                                    output_file,
                                    Geno.job.get_job( 'use_biobambam' ) )

    if not return_value:
        raise

#####################################################################
#
#   STAGE 7 bam statistics
#
#   in:     bam
#   out:    txt
#
@follows( stage_6 )
@active_if ( 'bam_stats' in Geno.job.get_job( 'tasks' )[ 'WGS' ] )
@files( generate_params_for_bam_stats )
@check_if_uptodate( check_file_exists_for_bam_stats )
def stage_7( input_file1, input_file2, output_file ):
    return_value = bam_stats(  input_file1, input_file2, output_file )
    if not return_value:
        raise

#####################################################################
#
#   STAGE 8 fisher_mutation_call
#
#   in:     bam
#   out:    txt
#
@follows( stage_6, stage_7 )
@active_if ( 'fisher_mutation_call' in Geno.job.get_job( 'tasks' )[ 'WGS' ] )
@files( generate_params_for_fisher_mutation_call )
@check_if_uptodate( check_file_exists_for_fisher_mutation_call )
def stage_8( control_file_list, disease_file_list, control_outdir, output_file ):
    return_value = fisher_mutation_call(  control_file_list, disease_file_list, control_outdir, output_file )
    if not return_value:
        raise

#####################################################################
#
#   STAGE 9 itd_detection
#
#   in:     bam
#   out:    txt
#
@follows( stage_8 )
@active_if ( 'itd_detection' in Geno.job.get_job( 'tasks' )[ 'WGS' ] )
@files( generate_params_for_itd_detection )
@check_if_uptodate( check_file_exists_for_itd_detection )
def stage_9( control_file_list, tumor_file_list, control_output_dir_list, tumor_output_dir_list ):
    return_value = itd_detection( control_file_list, tumor_file_list, control_output_dir_list, tumor_output_dir_list )
    if not return_value:
        raise

#####################################################################
#
#   STAGE 10 annotation
#
#   in:     txt
#   out:    xls or vcf
#
@follows( stage_9 )
@active_if ( 'annotation' in Geno.job.get_job( 'tasks' )[ 'WGS' ] )
@files( generate_params_for_annotation )
@check_if_uptodate( check_file_exists_for_annotation )
def stage_10(  input_file, output_file ):
    return_value = annotation(  input_file, output_file )
    if not return_value:
        raise

#####################################################################
#
#   LAST STAGE 
#
@follows( stage_10 )
def last_function():
    with log_mutex:
        log.info( "Genomon pipline has finished successflly!" )
    return True

