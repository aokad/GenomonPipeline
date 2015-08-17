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

use_subdir = ( Geno.job.get_job( 'sample_name' ) != None  or
               Geno.sample_list )

#####################################################################
#
# Subroutines
#
def make_sample_name( dir_name = None, filename = None, add_sample_group = False ):
    if use_subdir:
        if filename:
            subdir_name = os.path.basename( os.path.split( filename )[ 0 ] )
        elif dir_name:
            subdir_name = os.path.basename( dir_name )
        else:
            subdir_name = None
    else:
        subdir_name = None

    if add_sample_group and subdir_name:
        return_str = Geno.job.get_job( 'sample_group' ) + '_' + subdir_name
    elif subdir_name:
        return_str = subdir_name
    else:
        return_str = ''

    return return_str

def delete_intermediate_files( output_file ):
    #
    # Delete intermediate files 
    #
    # Replace large intermediate files by 0-byte-size files.
    #
    # starting from ( sample_name/R1.fastq, sample_name/R2.fastq )
    #
    # Files to delete
    # split fastq files: out/fastq/sample_name/R1_000000.fastq
    #                    out/fastq/sample_name/R2_000001.fastq
    # bwa_mem:           out/bam/sample_name/R1_000000_sorted.bam
    # bam_merge:         out/bam/sample_name/sample_name_merged.bam
    #
    pass

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

    ( output_prefix1, output_suffix1 ) = os.path.splitext( output_file1 )
    ( output_prefix2, output_suffix2 ) = os.path.splitext( output_file2 )
    outfile_list1 =  glob( "{prefix}*{suffix}".format( prefix = output_prefix1, suffix = output_suffix1 ) )
    outfile_list2 = glob( "{prefix}*{suffix}".format( prefix = output_prefix2, suffix = output_suffix2 ) )
    if not outfile_list1:
        return True, "Missing file {outprefix}*{outsuffix} for {input}.".format(
                            outprefix = output_prefix1,
                            outsuffix = output_suffix1,
                            input = input_file1 )
    elif not outfile_list2:
        return True, "Missing file {outprefix}*{outsuffix} for {input}.".format(
                            outprefix = output_prefix2,
                            outsuffix = output_suffix2,
                            input = input_file2 )
    else:
        exit_status = get_status_of_this_process( 'split_fastq', output_file1, Geno, use_subdir )

        in_time1 = os.path.getmtime( input_file1 )
        out_time1 = os.path.getmtime( outfile_list1[ 0 ] )

        in_time2 = os.path.getmtime( input_file2 )
        out_time2 = os.path.getmtime( outfile_list2[ 0 ] )

        if exit_status != 0 or in_time1 > out_time1 or in_time2 > out_time2:
            return True, "{outprefix}*{outsuffix} is older than {input}.".format(
                                outprefix = output_prefix1,
                                outsuffix = output_suffix1,
                                input = input_file1 )
        else:
            return False, "File {outprefix}*{outsuffix} exits for {input}.".format(
                                outprefix = output_prefix1,
                                outsuffix = output_suffix1,
                                input = input_file1 )


def check_file_exists_for_merge_bam(
    input_file_list,
    output_file
    ):
    if ( Geno.job.get_param( 'others', 'use_biobambam' ) and
         'markduplicates' in Geno.job.get_job( 'tasks' )[ 'DNA'] ):
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

    exit_status = get_status_of_this_process( process_name, output_file, Geno, use_subdir )

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

    exit_status = get_status_of_this_process( 'bwa_mem', output_file1, Geno, use_subdir )
    ( output_prefix, output_suffix ) = os.path.splitext( output_file1 )
    if Geno.job.get_param( 'others', 'use_biobambam' ):
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

    exit_status = get_status_of_this_process( process_name, output_file1, Geno, use_subdir )
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
    exit_status = get_status_of_this_process( 'fisher_mutation_call', output_file, Geno, use_subdir )

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
    exit_status = get_status_of_this_process( 'bam_stats_calc', output_file, Geno, use_subdir ) + \
                  get_status_of_this_process( 'bam_stats_merge', output_file, Geno, use_subdir )

    summary_dir_name = os.path.dirname( output_file )
    out_sum_file =  summary_dir_name + '/' + make_sample_name( filename = output_file ) + '.txt'
    if exit_status != 0 or not os.path.exists( out_sum_file ):
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
    exit_status = get_status_of_this_process( 'itd_detection', out_dir_list[ 0 ], Geno, use_subdir )

    for in_file, out_dir in zip( in_file_list,  out_dir_list ) :
        if exit_status != 0 or not os.path.exists( out_dir + '/itd_list.tsv' ):
            return True, "Missing file {outputfile} for inputfile.".format(
                                outputfile = out_dir,
                                inputfile = in_file )

        else:
            in_time = os.path.getmtime( in_file )
            out_time = os.path.getmtime( out_dir )
            if in_time > out_time:
                return True, "{output} is older than {input}.".format( output = out_dir, input = in_file )

    return False, "Output files  exits."


def check_file_exists_for_sv_parse(target_label, target_bam, target_outdir, match_use, match_bam):

    exit_status = get_status_of_this_process('sv_parse', target_label, Geno, use_subdir)
    input = target_bam
    output = target_outdir + "/" + target_label + ".junction.clustered.bedpe.gz"

    if exit_status != 0 or not os.path.exists(output):
        return True, "Missing file {outputfile} for {inputfile}.".format(
                            outputfile = output,
                            inputfile = input )
 
    else:
        in_time = os.path.getmtime(input)
        out_time = os.path.getmtime(output)
        if in_time > out_time:
            return True, "{outputfile} is older than {inputfile}.".format(outputfile = output, inputfile = input)
        else:
            return False, "File {outputfile} exits for {inputfile}.".format(outputfile = output, inputfile = input)


def check_file_exists_for_sv_filt(target_label, target_outdir):

    exit_status = get_status_of_this_process('sv_filt', target_label, Geno, use_subdir)
    input = target_outdir + "/" + target_label + ".junction.clustered.bedpe.gz"
    output = target_outdir + "/" + target_label + ".genomonSV.result.txt"

    if exit_status != 0 or not os.path.exists(output):
        return True, "Missing file {outputfile} for {inputfile}.".format(
                            outputfile = output,
                            inputfile = input )

    else:
        in_time = os.path.getmtime(input)
        out_time = os.path.getmtime(output)
        if in_time > out_time:
            return True, "{outputfile} is older than {inputfile}.".format(outputfile = output, inputfile = input)
        else:
            return False, "File {outputfile} exits for {inputfile}.".format(outputfile = output, inputfile = input)


def check_file_exists_for_mutation_filter(
        target_list,
        target_normal_bam,
        taraget_tumor_bam,
        output_list,
        output_dir
    ):

    exit_status = get_status_of_this_process( 'mutation_filter', output_list, Geno, use_subdir )

    if exit_status != 0 or not os.path.exists( output_list ):
        return True, "Missing file {output} for {input}.".format( output = output_list, input = target_list)
    else:
        in_time = os.path.getmtime( target_list )
        out_time = os.path.getmtime( output_list )
        if in_time > out_time:
            return True, "{output} is older than {input}.".format( output = output_list, input = target_list)
        else:
            return False, "File {output} exits for {input}.".format( output = output_list, input = target_list)


def check_file_exists_for_annotation(
        input_file,
        output_file
    ):
    """
    """
    exit_status = get_status_of_this_process( 'annotation', output_file, Geno, use_subdir )
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
        return_list = [ input_file_list[ dir_name ],
                        dir_name + '/' + make_sample_name( filename = input_file_list[ dir_name ][ 0 ] ) + '_merged.bam' ]
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

    if Geno.job.get_param( 'others', 'use_biobambam' ):
        for dir_name in input_file_list.keys():
            return_list = [ input_file_list[ dir_name ],
                            dir_name + '/' + make_sample_name( dir_name = dir_name ) + '_markdup.bam' ]
            yield return_list
    else:
        for dir_name in input_file_list.keys():
            return_list = [ dir_name +'/' + make_sample_name( dir_name = dir_name ) + '_merged.bam',
                            dir_name + '/' + make_sample_name( dir_name = dir_name ) + '_markdup.bam' ]
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
            if 'markduplicates' in Geno.job.get_job( 'tasks' )[ 'DNA' ]:
                input_bam =  dir_name + '/' + make_sample_name( filename = param[ 0 ] ) + '_markdup.bam'
            elif 'merge_bam' in Geno.job.get_job( 'tasks' )[ 'DNA' ]:
                input_bam =  dir_name + '/' + param[ 0 ] + '_merged.bam'
            else:
                input_bam =  param[ 0 ]

            input_file_list[ dir_name ] = input_bam
            summary_dir_name = os.path.dirname( param[ 2 ] )
            return_list = [ input_bam,
                            'None',
                            summary_dir_name + '/' + make_sample_name( filename = param[ 2 ] ) + '.txt' ]
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
    if 'markduplicates' in Geno.job.get_job( 'tasks' )[ 'DNA' ]:
        input_bam =  dir_name + '/' + make_sample_name( dir_name = dir_name ) + '_markdup.bam'
    elif 'merge_bam' in Geno.job.get_job( 'tasks' )[ 'DNA' ]:
        input_bam =  dir_name + '/' + make_sample_name( dir_name = dir_name ) + '_merged.bam'
    else:
        input_bam =  filename

    return input_bam

def generate_params_for_fisher_mutation_call():
    """
    Generate parameter list for fisher_mutation_call

    """
    
    param_return = False
    task_name_list = ( 'markduplicates', 'merge_bam', 'bwa_mem' )
    id = 0
    while( not param_return ):
        param_return = Sample.make_param( 'fisher_mutation_call', task_name_list[ id ], '.txt', 'mutation', 2, 1 )
        id += 1

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
                                mutation_dir_name + '/' + make_sample_name( dir_name = mutation_dir_name ) + '.txt'
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
            data_in_data_dict = False
            normal_list = data_type.replace( ' ', '' ).split( ',' )
            for data_tmp in normal_list:
                if data_tmp in data_dict.keys():
                    data_in_data_dict = True
                    break

            if data_type != 'Normal' and data_type != 'Disease' and data_in_data_dict:
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
                normal_dir_list = []
                normal_outfile = None
                for normal_sample in normal_list:
                    if normal_sample in data_dict.keys():
                        normal_input = make_bam_filename_for_markdup_result(
                                param_list[ data_dict[ normal_sample ] ][ 0 ] )
                        normal_dir_list.append( normal_input )
                        if normal_outfile == None:
                            mutation_dir_name = os.path.dirname( param_list[ data_dict[ normal_sample ] ][ 2 ] )
                            normal_outfile = mutation_dir_name + '/' + make_sample_name( dir_name = mutation_dir_name) + '.txt'

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
                                disease_outfile = mutation_dir_name + '/' + make_sample_name( dir_name = mutation_dir_name) + '.txt'
                    else:
                        disease_input = make_bam_filename_for_markdup_result( 
                                param_list[ data_dict[ disease_data ] ][ 0 ] )
                        disease_merge_list.append( disease_input )
                        mutation_dir_name = os.path.dirname( param_list[ data_dict[ disease_data ] ][ 3 ] )
                        disease_outfile =  mutation_dir_name + '/' + make_sample_name( dir_name = mutation_dir_name ) + '.txt'
                    #
                    # Return parameters
                    #
                    yield( normal_dir_list,
                           disease_merge_list,
                           os.path.split( normal_outfile )[ 0 ],
                           disease_outfile )

            elif data_type == 'Disease':
                for tumor_dir in ctrl_dis_pairs[ 'Disease' ]:
                    input_bam = make_bam_filename_for_markdup_result( param_list[ data_dict[ tumor_dir ] ][ 0 ] )
                    mutation_dir_name = os.path.dirname( param_list[ data_dict[ tumor_dir ] ][ 3 ] )

                    return_list = [
                                    [ None ],
                                    [ input_bam ],
                                    None ,
                                    mutation_dir_name + '/' + make_sample_name( dir_name = mutation_dir_name ) + '.txt'
                                  ]
                    yield return_list

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
            dir_name = os.path.split( param[ 2 ] )[ 0 ]
            output_dir_list.append( dir_name + '/'+ make_sample_name( dir_name = dir_name ) )

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

            elif data_type in data_dict.keys():
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
                    if -1 == disease_data.find( ',' ):
                        disease_sample = disease_data
                    else:
                        disease_sample = disease_data.replace( ' ', '' ).split( ',' )[ 0 ]

                    disease_input = make_bam_filename_for_markdup_result( 
                            param_list[ data_dict[ disease_sample ] ][ 0 ] )
                    disease_file_list.append( disease_input )
                    itd_dir_name = os.path.dirname( param_list[ data_dict[ disease_sample ] ][ 3 ] )
                    disease_outfile_list.append( itd_dir_name )

        #
        # Return parameters
        #
        yield( normal_file_list,
               disease_file_list,
               normal_outfile_list,
               disease_outfile_list )


def generate_params_for_sv_parse():


    Sample.make_param( 'sv_detection', 'markduplicates', '.txt', 'sv', 1, 1 )
    param_list = Sample.param( 'sv_detection' )

    # get the bam-path and output dir (for sv) for each sample
    sample_name2bampath = {}
    sample_name2outputdir = {}
    for input1, input2, output1, output2 in param_list:
        bam_path = make_bam_filename_for_markdup_result(input1)
        bam_dirname = os.path.dirname(bam_path)
        sample_name = os.path.basename(bam_dirname)
        out_path = os.path.dirname(output1)

        sample_name2bampath[sample_name] = bam_path
        sample_name2outputdir[sample_name] = out_path


    # on tumor-control relationships
    control2tumor = Geno.job.get_job( 'control_disease_pairs' )
    # reverse the key-value for ease of treatment
    tumor2control = dict([(v, k) for k, v in control2tumor.items()])


    for sample_name in sample_name2bampath:

        if sample_name in tumor2control:
            control_bam = sample_name2bampath[tumor2control[sample_name]]
            yield [sample_name, sample_name2bampath[sample_name], sample_name2outputdir[sample_name], True, control_bam]
        else:
            yield [sample_name, sample_name2bampath[sample_name], sample_name2outputdir[sample_name], False, None]
            

def generate_params_for_sv_filt():
 
    Sample.make_param( 'sv_detection', 'markduplicates', '.txt', 'sv', 1, 1 )
    param_list = Sample.param( 'sv_detection' )
    
    # get the bam-path and output dir (for sv) for each sample
    sample_name2outputdir = {}
    for input1, input2, output1, output2 in param_list:
        bam_path = make_bam_filename_for_markdup_result(input1)
        bam_dirname = os.path.dirname(bam_path)
        sample_name = os.path.basename(bam_dirname)
        out_path = os.path.dirname(output1)

        sample_name2outputdir[sample_name] = out_path



    # on tumor-control relationships
    control2tumor = Geno.job.get_job( 'control_disease_pairs' )
    # reverse the key-value for ease of treatment
    tumor2control = dict([(v, k) for k, v in control2tumor.items()])

    for sample_name in sample_name2outputdir:
        print sample_name
        if sample_name in tumor2control:
            yield [sample_name, sample_name2outputdir[sample_name]]


def generate_params_for_mutation_filter( ):

    ctrl_dis_pairs = Geno.job.get_job( 'control_disease_pairs' )
    for subdir in ctrl_dis_pairs.keys():
        in_mutation_txt  =  Geno.dir[ 'mutation' ]  + '/' + ctrl_dis_pairs[subdir]      + '/' + Geno.job.get_job( 'sample_name' ) + '.txt'
        in_normal_bam    =  Geno.dir[ 'bam' ]       + '/' + str(subdir)                 + '/' + Geno.job.get_job( 'sample_name' ) + '_markdup.bam'
        in_tumor_bam     =  Geno.dir[ 'bam' ]       + '/' + str(ctrl_dis_pairs[subdir]) + '/' + Geno.job.get_job( 'sample_name' ) + '_markdup.bam'
        out_mutation_txt =  Geno.dir[ 'mutfilter' ] + '/' + ctrl_dis_pairs[subdir]      + '/' + Geno.job.get_job( 'sample_name' ) + '_mutfilter.txt'
        out_mutation_dir =  Geno.dir[ 'mutfilter' ] + '/' + ctrl_dis_pairs[subdir] 

        yield (in_mutation_txt, in_tumor_bam, in_normal_bam, out_mutation_txt, out_mutation_dir) 


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
            yield ( mutation_dir_name + '/' + make_sample_name( dir_name = mutation_dir_name ) + '.txt',
                    annotation_dir_name + '/' + make_sample_name( dir_name = annotation_dir_name ) ) 

    else:
        disease_subdir_list = []
        for subdir in ctrl_dis_pairs.keys():
            if subdir == 'Disease':
                disease_subdir_list += ctrl_dis_pairs[ 'Disease' ]
            elif subdir  != 'Normal':
                if isinstance( ctrl_dis_pairs[ subdir ], list ):
                    #
                    # Case:
                    #   'Normal1':
                    #           - 'Tumor1'
                    #           - 'Tumor2'
                    disease_subdir_list += ctrl_dis_pairs[ subdir ]
                elif -1 == ctrl_dis_pairs[ subdir ].find( ',' ):
                    #
                    # Case:
                    #   'Normal1': 'Tumor1'
                    disease_subdir_list.append( ctrl_dis_pairs[ subdir ] )
                else:
                    #
                    # Case:
                    #   'Normal1': 'Tumor1,Tumor2'
                    # Output only the first sample  in 'Tumor1,Tumor2'
                    disease_subdir_list += [ ctrl_dis_pairs[ subdir ].replace( ' ', '' ).split( ',' )[ 0 ] ]

        for param in Sample.param( 'annotation' ):
            mutation_dir_name = os.path.dirname( param[ 0 ] )
            annotation_dir_name = os.path.dirname( param[ 3 ] )
            if os.path.basename( mutation_dir_name ) in disease_subdir_list:
                yield ( mutation_dir_name + '/' + make_sample_name( dir_name = mutation_dir_name ) + '.txt',
                        annotation_dir_name + '/' + make_sample_name( dir_name = annotation_dir_name ) ) 



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


    args = { 'log': Geno.dir[ 'log' ],
             'array_data': array_in + array_out,
             'input_file': input_file,
             'output_file': output_file,
             'scriptdir': Geno.dir[ 'script' ]  }
    shell_script_full_path = make_script_file( 'extract_fastq', res.extract_gz, Geno, **args )

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
        raise Exception( '{0} failed.'.format( function_name ) )


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
        if not os.path.isfile( input_file1 ):
            with log_mutex:
                log.error( "file: {file} does not exist.".format( file=input_file1 ) )
            raise Exception( '{0} failed.'.format( function_name ) )

        #
        # Make shell script
        #
        shell_script_full_path = make_script_file_name( function_name, Geno )
        shell_script_file = open( shell_script_full_path, 'w' )

        file_type = Geno.job.get_job( 'input_file_type' )

        if 'bam' == file_type:
            shell_script_file.write( wgs_res.bamtofastq_p.format(
                                            log = Geno.dir[ 'log' ],
                                            bamfile = input_file1,
                                            outfastq1 = output_file1,
                                            outfastq2 = output_file2,
                                            tmpfastq = output_file1 + '.tmp',
                                            biobambam = Geno.conf.get( 'SOFTWARE', 'biobambam' ),
                                            scriptdir = Geno.dir[ 'script' ]
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
            raise Exception( '{0} failed.'.format( function_name ) )

        save_status_of_this_process( function_name, output_file1, runtask_return_code, Geno, use_subdir )

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
            raise Exception( '{0} failed.'.format( function_name ) )

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
        args = { 'log': Geno.dir[ 'log' ],
                 'fastq_filter': Geno.job.get_param( 'split_fastq', 'fastq_filter' ),
                 'array_data': array_in + array_out,
                 'lines_per_file': Geno.job.get_param( 'split_fastq', 'split_fastq_line_number' ),
                 'input_file': input_file,
                 'suffix_len': suffix_len,
                 'output_suffix': output_suffix,
                 'output_prefix': output_prefix,
                 'scriptdir': Geno.dir[ 'script' ] }
        shell_script_full_path = make_script_file( function_name, wgs_res.splitfile, Geno, **args )

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
            raise Exception( '{0} failed.'.format( function_name ) )

        save_status_of_this_process( function_name, output_file1, runtask_return_code, Geno, use_subdir )

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
            raise Exception( '{0} failed.'.format( function_name ) )

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
        args = {
            'log': Geno.dir[ 'log' ],
            'infastq': input_file,
            'outfastq': output_file,
            'array_data': array_data1 + array_data2,
            'tmpoutfastq': output_file + '.tmp',
            'optadapters': '-a ' + ' -a '.join( Geno.job.get_param( 'cutadapt', 'adaptor' ) ),
            'casavacode': 2,
            'cutadapt': Geno.conf.get( 'SOFTWARE', 'cutadapt' ),
            'scriptdir': Geno.dir[ 'script' ]
        }
        shell_script_full_path = make_script_file( function_name, wgs_res.cutadapt, Geno, **args )

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
            raise Exception( '{0} failed.'.format( function_name ) )

        save_status_of_this_process( function_name, output_file1, runtask_return_code, Geno, use_subdir )

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

        if id > 0: # Error check: make sure that there are fastq files to align
            #
            # Make shell script for array job
            #
            if use_biobambam:
                bwa_mem_resource = wgs_res.bwa_mem_biobambam
            else:
                bwa_mem_resource = wgs_res.bwa_mem

            env_variable_str = "LD_LIBRARY_PATH=$LD_LIBRARY_PATH:{libmaus_PATH}".format(
                                            libmaus_PATH = Geno.conf.get( 'ENV', 'libmaus_PATH' ) )
            args = {
                'log': Geno.dir[ 'log' ],
                'array_data': str1 + str2 + str3,
                'fastq1': "${FILE1[$SGE_TASK_ID]}",
                'fastq2': "${FILE2[$SGE_TASK_ID]}",
                'bam': "${FILE3[$SGE_TASK_ID]}",
                'read_group': Geno.job.get_job(  'bam_read_group' ),
                'min_score': Geno.job.get_param( 'bwa_mem', 'min_score' ),
                'additional_params': Geno.job.get_param( 'bwa_mem', 'additional_params' ),
                'env_variables': env_variable_str,
                'ref_fa': Geno.conf.get( 'REFERENCE', 'ref_fasta' ),
                'bwa': Geno.conf.get( 'SOFTWARE', 'bwa' ),
                'samtools': Geno.conf.get( 'SOFTWARE', 'samtools' ),
                'biobambam': Geno.conf.get( 'SOFTWARE', 'biobambam' ),
                'scriptdir': Geno.dir[ 'script' ] 
            }
            shell_script_full_path = make_script_file( function_name, bwa_mem_resource, Geno, **args )

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
                raise Exception( '{0} failed.'.format( function_name ) )
                
            save_status_of_this_process( function_name, output_file1, runtask_return_code, Geno, use_subdir )

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

        args = {
            'log': Geno.dir[ 'log' ],
            'input_bam_files': input_files,
            'output_bam_file': output_file,
            'env_variables': env_variable_str,
            'samtools': Geno.conf.get( 'SOFTWARE', 'samtools' ),
            'biobambam': Geno.conf.get( 'SOFTWARE', 'biobambam' ),
            'scriptdir': Geno.dir[ 'script' ] 
        }
        shell_script_full_path = make_script_file( function_name, bam_merge_resource, Geno, **args )

        #
        # Run
        #
        runtask_return_code = Geno.RT.runtask(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ function_name ] )
        if runtask_return_code != 0:
            with log_mutex:
                log.error( "{function}: runtask failed".format( function = function_name ) )
            raise Exception( '{0} failed.'.format( function_name ) )

        save_status_of_this_process( function_name, output_file, runtask_return_code, Geno, use_subdir )

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
                                            biobambam = Geno.conf.get( 'SOFTWARE', 'biobambam' ),
                                            scriptdir = Geno.dir[ 'script' ] 
                                        ) )
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
                                            picard = Geno.conf.get( 'SOFTWARE', 'picard' ),
                                            scriptdir = Geno.dir[ 'script' ] 
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
            raise Exception( '{0} failed.'.format( function_name ) )

        save_status_of_this_process( function_name, output_file, runtask_return_code, Geno, use_subdir )

        #
        # Delete intermediate files 
        #
        # Replace large intermediate files by 0-byte-size files.
        #
        # starting from ( sample_name/R1.fastq, sample_name/R2.fastq )
        #
        # Files to delete
        # split fastq files: out/fastq/sample_name/R1_000000.fastq
        #                    out/fastq/sample_name/R2_000001.fastq
        # bwa_mem:           out/bam/sample_name/R1_000000_sorted.bam
        # bam_merge:         out/bam/sample_name/sample_name_merged.bam
        #
        if True:
            delete_intermediate_files( output_file )

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
            env_variable_str = "PERL5LIB=$PERL5LIB:{0}".format(
                                    Geno.conf.get( 'ENV', 'PERL5LIB' ) )

            shell_script_name = function_name + '_calc'
            shell_script_full_path = make_script_file_name( shell_script_name , Geno )
            shell_script_file = open( shell_script_full_path, 'w' )
            shell_script_file.write( wgs_res.bam_stats_calc.format(
                                            log = Geno.dir[ 'log' ],
                                            env_variables = env_variable_str,
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
                                            scriptdir = Geno.dir[ 'script' ]
                                        )
                                    )
            shell_script_file.close()

            calc_return_code = Geno.RT.run_arrayjob(
                                shell_script_full_path,
                                Geno.job.get_job( 'cmd_options' )[ function_name ],
                                id_start = 1,
                                id_end = 14 )

            if calc_return_code == 0:
                save_status_of_this_process( shell_script_name, output_file, calc_return_code, Geno, use_subdir )

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
                                            scriptdir = Geno.dir[ 'script' ]
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
                raise Exception( '{0} failed.'.format( function_name ) )
            if merge_return_code != 0:
                with log_mutex:
                    log.error( "{function}: runtask failed".format( function = function_name + '_merge' ) )
                raise Exception( '{0} failed.'.format( function_name ) )

            if merge_return_code == 0:
                save_status_of_this_process( shell_script_name, output_file, merge_return_code, Geno, use_subdir )

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
            bambam_bam_file_list_array   = "BAMBAM_INPUT_FILE=(\n"
            output_dir_array  = "OUT_FILE=(\n"
            sample_name_array = "SAMPLE_NAME=(\n"
            merge_bam_flag    = "FLAG=(\n"

            samtools_bam_file_list_array += " [1]=\"{0}\"\n".format( samtools_control_bam_file_list )
            bambam_bam_file_list_array   += " [1]=\"{0}\"\n".format( bambam_control_bam_file_list )
            output_dir_array  += " [1]=\"{0}\"\n".format( control_output_dir )
            sample_name_array += " [1]=\"{0}\"\n".format( make_sample_name( dir_name = control_output_dir ) )
            merge_bam_flag    += " [1]=\"{0}\"\n".format( control_merge_bam_flag )

            samtools_bam_file_list_array += " [2]=\"{0}\"\n".format( samtools_disease_bam_file_list )
            bambam_bam_file_list_array   += " [2]=\"{0}\"\n".format( bambam_disease_bam_file_list )
            output_dir_array  += " [2]=\"{0}\"\n".format(  disease_output_dir )
            sample_name_array += " [2]=\"{0}\"\n".format(  make_sample_name( dir_name = disease_output_dir ) )
            merge_bam_flag    += " [2]=\"{0}\"\n".format(  disease_merge_bam_flag )

            samtools_bam_file_list_array += ")\n"
            bambam_bam_file_list_array   += ")\n"
            output_dir_array  += ")\n"
            sample_name_array += ")\n"
            merge_bam_flag    += ")\n"

            #
            # Make shell script
            #
            shell_script_full_path = make_script_file_name( function_name + '_merge_bam', Geno )
            shell_script_file = open( shell_script_full_path, 'w' )
            shell_script_file.write( wgs_res.fisher_merge_bams.format(
                                            log = Geno.dir[ 'log' ],
                                            env_variables = env_variable_str,
                                            use_biobambam = Geno.job.get_param( 'others', 'use_biobambam' ),
                                            array_data = samtools_bam_file_list_array + bambam_bam_file_list_array +\
                                                         output_dir_array + merge_bam_flag + sample_name_array,
                                            bambam_input_bam_files = "${BAMBAM_INPUT_FILE[$SGE_TASK_ID]}",
                                            samtools_input_bam_files = "${SAMTOOLS_INPUT_FILE[$SGE_TASK_ID]}",
                                            input_bam_file = "${SAMTOOLS_INPUT_FILE[$SGE_TASK_ID]}",
                                            merge_bam_flag = "${FLAG[$SGE_TASK_ID]}",
                                            merged_bam_file = '${OUT_FILE[$SGE_TASK_ID]}/${SAMPLE_NAME[$SGE_TASK_ID]}.bam',
                                            out_dir = '${OUT_FILE[$SGE_TASK_ID]}',
                                            biobambam = Geno.conf.get( 'SOFTWARE', 'biobambam' ),
                                            samtools = Geno.conf.get( 'SOFTWARE', 'samtools' ),
                                            scriptdir = Geno.dir[ 'script' ] 
                                        ) )

            shell_script_file.close()

            runtask_return_code = Geno.RT.run_arrayjob(
                                shell_script_full_path,
                                Geno.job.get_job( 'cmd_options' )[ function_name ],
                                id_start = 1,
                                id_end = 2 )

            if runtask_return_code != 0:
                log.error( "{function}: runtask failed".format( function = function_name ) )
                raise Exception( '{0} failed.'.format( function_name ) )

            control_input_bam = control_output_dir + '/{sample_name}.bam'.format(
                                    sample_name = make_sample_name( dir_name = control_output_dir ) )
            disease_input_bam = disease_output_dir + '/{sample_name}.bam'.format(
                                    sample_name = make_sample_name( dir_name = disease_output_dir ) )
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
                                        chr_str_in_fa = Geno.conf.get( 'REFERENCE', 'chr_str_in_fa' ),
                                        remove_intermediate = True,
                                        max_indel = Geno.job.get_param( 'fisher_mutation_call', 'max_indel' ),
                                        max_distance = Geno.job.get_param( 'fisher_mutation_call', 'max_distance' ),
                                        base_quality = Geno.job.get_param( 'fisher_mutation_call', 'base_quality' ),
                                        map_quality = Geno.job.get_param( 'fisher_mutation_call', 'map_quality' ),
                                        mismatch_rate = Geno.job.get_param( 'fisher_mutation_call', 'mismatch_rate' ),
                                        min_depth = Geno.job.get_param( 'fisher_mutation_call', 'min_depth' ),
                                        samtools = Geno.conf.get( 'SOFTWARE', 'samtools' ),
                                        python = Geno.conf.get( 'SOFTWARE', 'python' ),
                                        scriptdir = Geno.dir[ 'script' ]
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
            raise Exception( '{0} failed.'.format( function_name ) )

        #
        # Make shell script
        #
        shell_script_full_path = make_script_file_name( 'merge_fisher_result', Geno )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( wgs_res.merge_fisher_result.format(
                                        log = Geno.dir[ 'log' ],
                                        output_txt = disease_output_file,
                                        scriptdir = Geno.dir[ 'script' ] 
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
                log.error( "{function}: runtask failed".format( function = 'merge_fisher_result' ) )
            raise Exception( '{0} failed.'.format( function_name ) )

        save_status_of_this_process( function_name, disease_output_file, runtask_return_code, Geno, use_subdir )

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
            normal_function_name = function_name + '_normal'
            shell_script_full_path = make_script_file_name( normal_function_name, Geno )
            shell_script_file = open( shell_script_full_path, 'w' )
            shell_script_file.write( wgs_res.itd_detection.format(
                                            log = Geno.dir[ 'log' ],
                                            array = input_files + output_files + name_list,
                                            itd_inhouse_files = ' ',
                                            bam_file = "${FILE1[$SGE_TASK_ID]}",
                                            output_file = "${FILE2[$SGE_TASK_ID]}",
                                            name = "${NAME[$SGE_TASK_ID]}",
                                            itd_inhouse_dir = Geno.dir[ 'config' ],
                                            itd_detector = Geno.conf.get( 'SOFTWARE', 'itd_detector' ),
                                            scriptdir = Geno.dir[ 'script' ] 
                                     ) )
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
                    log.error( "{function}: runtask failed. Return code: {code}.".format(
                                            function = normal_function_name,
                                            code = runtask_return_code ) )
                raise Exception( '{0} failed.'.format( function_name ) )

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
        # Remove ctrl_panel_normal from the list in sample_name
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
            input_file_list = control_file_list if control_file_list != None else []
            input_file_list += tumor_file_list if tumor_file_list != None else [] 
            output_file_list = ctrl_output_file_list if ctrl_output_file_list != None else []
            output_file_list += tumor_output_file_list if tumor_output_file_list != None else [] 
            for input_file, output_dir in zip( input_file_list, output_file_list ):
                data_list.append( ( input_file, output_dir ) )


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
        tumor_function_name = function_name + '_tumor'
        shell_script_full_path = make_script_file_name( tumor_function_name, Geno )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( wgs_res.itd_detection.format(
                                        log = Geno.dir[ 'log' ],
                                        array = input_files + output_files + name_list,
                                        itd_inhouse_files = itd_dest_file + ',' + bp_dest_file,
                                        bam_file = "${FILE1[$SGE_TASK_ID]}",
                                        output_file = "${FILE2[$SGE_TASK_ID]}",
                                        name = "${NAME[$SGE_TASK_ID]}",
                                        itd_inhouse_dir = Geno.dir[ 'config' ],
                                        itd_detector = Geno.conf.get( 'SOFTWARE', 'itd_detector' ),
                                        scriptdir = Geno.dir[ 'script' ] 
                                ) )
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
                log.error( "{function}: runtask failed".format( function = tumor_function_name ) )
            raise Exception( '{0} failed.'.format( function_name ) )

        save_status_of_this_process( function_name, data_list[ 0 ][ 1 ], runtask_return_code, Geno, use_subdir )

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
#  mutation filter
#
def mutation_filter(
        target_list,
        target_tumor_bam,
        target_normal_bam,
        output_list,
        output_dir
        ):
    """
        mutation filter: 

    """
    return_value = True

    try:
        function_name = whoami()
        log.info( "#{function}".format( function = function_name ) )

        #
        # Make shell script
        #
        shell_script_full_path = make_script_file_name( function_name, Geno )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( wgs_res.mutation_filter.format(
                                     log = Geno.dir[ 'log' ],
                                     scriptdir = Geno.dir[ 'script' ],
                                     mutfilter = Geno.conf.get( 'SOFTWARE', 'mutfilter' ),
                                     realignment_params = Geno.job.get_param( 'mutfilter', 'realignment_params' ),
                                     indel_params = Geno.job.get_param( 'mutfilter', 'indel_params' ),
                                     breakpoint_params = Geno.job.get_param( 'mutfilter', 'breakpoint_params' ),
                                     simplerepeat_params = Geno.job.get_param( 'mutfilter', 'simplerepeat_params' ),
                                     target_list = target_list,
                                     target_tumor_bam = target_tumor_bam,
                                     target_normal_bam = target_normal_bam,
                                     output_list = output_list,
                                     ref_fasta = Geno.conf.get( 'REFERENCE', 'ref_fasta' ),
                                     simple_repeat_db = Geno.conf.get( 'REFERENCE', 'simple_repeat_tabix_db' ),
                                     blat = Geno.conf.get( 'SOFTWARE', 'blat' ),
                                     tmp_out_realignment = output_dir + '/' + "tmp_realignment.txt",
                                     tmp_out_indel =  output_dir + '/' + "tmp_indel.txt",
                                     tmp_out_breakpoint =  output_dir + '/' + "tmp_breakpoint.txt"
                                     )  )
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

        save_status_of_this_process( function_name, output_list, return_code, Geno, use_subdir )

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
        log.error( "{function}: Unexpected error: {error}.".format(
                    function = whoami() , error = sys.exc_info()[0] ) )
        log.error("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )
        return_value = False

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
                                        python = Geno.conf.get( 'SOFTWARE', 'python' ),
                                        scriptdir = Geno.dir[ 'script' ] 
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
            raise Exception( '{0} failed.'.format( function_name ) )

        save_status_of_this_process( function_name, output_file, runtask_return_code, Geno, use_subdir )

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
if Geno.input_file_list:
    Sample.set_sample_list( Geno.input_file_list )

#####################################################################
#
#   STAGE 1 bam2fastq
#   in:     bam
#   out:    fastq
#
@active_if( 'bam2fastq' in Geno.job.get_job( 'tasks' )[ 'DNA' ] )
@check_if_uptodate( check_file_exists_for_bam2fastq )
@files( generate_params_for_bam2fastq )
def stage_1( input_file1, input_file2, output_file1, output_file2 ):
    return_value =  bam2fastq( input_file1, input_file2, output_file1, output_file2 )
    if not return_value:
        raise Exception( 'stage_1 failed.' )

#####################################################################
#
#   STAGE 2 split_fastq
#   in:     fastq
#   out:    fastq * X
#
@follows( stage_1 )
@active_if( 'split_fastq' in Geno.job.get_job( 'tasks' )[ 'DNA' ] )
@check_if_uptodate( check_file_exists_for_split_fastq )
@files( generate_params_for_split_fastq )
def stage_2( input_file1, input_file2, output_file1, output_file2 ):
    return_value = split_fastq( input_file1, input_file2, output_file1, output_file2 )
    if not return_value:
        raise Exception( 'stage_2 failed.' )

#####################################################################
#
#   STAGE 3 cutadapt
#   in:     fastq
#   out:    fastq
#
@follows( stage_2 )
@active_if( 'cutadapt' in Geno.job.get_job( 'tasks' )[ 'DNA' ] )
@check_if_uptodate( check_file_exists_for_cutadapt )
@files( generate_params_for_cutadapt )
def stage_3( input_file1, input_file2, output_file1, output_file2 ):
    return_value = cutadapt( input_file1, input_file2, output_file1, output_file2 )
    if not return_value:
        raise Exception( 'stage_3 failed' )

#####################################################################
#
#   STAGE 4 bwa_mem
#
#   in:     fastq1, fastq2
#   out:    bam
#
@follows( stage_3 )
@active_if( 'bwa_mem' in Geno.job.get_job( 'tasks' )[ 'DNA' ] )
@check_if_uptodate( check_file_exists_for_bwa_mem )
@files( generate_params_for_bwa_mem )
def stage_4(  input_file1, input_file2, output_file1, output_file2 ):
    return_value = bwa_mem( input_file1,
                            input_file2,
                            output_file1,
                            output_file2, 
                            Geno.job.get_param( 'others', 'use_biobambam' ) )
    if not return_value:
        raise Exception( 'stage_4 failed.' )

#####################################################################
#
#   STAGE 5 merge
#
#   in:     bam x {number}
#   out:    bam
#
@follows( stage_4 )
@active_if( 'merge_bam' in Geno.job.get_job( 'tasks' )[ 'DNA' ] )
@check_if_uptodate( check_file_exists_for_merge_bam)
@files( generate_params_for_merge_bam )
def stage_5( input_file_list, output_file ):
    if ( Geno.job.get_param( 'others', 'use_biobambam' ) and
         'markduplicates' in Geno.job.get_job( 'tasks' )[ 'DNA'] ):
        return_value = True
    else:
        return_value = merge_bam( input_file_list,
                                  output_file,
                                  Geno.job.get_param( 'others', 'use_biobambam' ) )

    if not return_value:
        raise Exception( 'stage_5 failed.' )

#####################################################################
#
#   STAGE 6 markduplicates
#
#   in:     bam x {number}
#   out:    bam
#
@follows( stage_5 )
@active_if( 'markduplicates' in Geno.job.get_job( 'tasks' )[ 'DNA' ] )
@check_if_uptodate( check_file_exists_for_markduplicates )
@files( generate_params_for_markduplicates )
def stage_6( input_file_list, output_file ):
    return_value =  markduplicates( input_file_list,
                                    output_file,
                                    Geno.job.get_param( 'others', 'use_biobambam' ) )

    if not return_value:
        raise Exception( 'stage_6 failed.' )

#####################################################################
#
#   STAGE 7 bam statistics
#
#   in:     bam
#   out:    txt
#
@follows( stage_6 )
@active_if ( 'bam_stats' in Geno.job.get_job( 'tasks' )[ 'DNA' ] )
@check_if_uptodate( check_file_exists_for_bam_stats )
@files( generate_params_for_bam_stats )
def stage_7( input_file1, input_file2, output_file ):
    return_value = bam_stats(  input_file1, input_file2, output_file )
    if not return_value:
        raise Exception( 'stage_7 failed.' )

#####################################################################
#
#   STAGE 8 fisher_mutation_call
#
#   in:     bam
#   out:    txt
#
@follows( stage_6, stage_7 )
@active_if ( 'fisher_mutation_call' in Geno.job.get_job( 'tasks' )[ 'DNA' ] )
@check_if_uptodate( check_file_exists_for_fisher_mutation_call )
@files( generate_params_for_fisher_mutation_call )
def stage_8( control_file_list, disease_file_list, control_outdir, output_file ):
    return_value = fisher_mutation_call(  control_file_list, disease_file_list, control_outdir, output_file )
    if not return_value:
        raise Exception( 'stage_8 failed.' )

#####################################################################
#
#   STAGE 9 itd_detection
#
#   in:     bam
#   out:    txt
#
@follows( stage_8 )
@active_if ( 'itd_detection' in Geno.job.get_job( 'tasks' )[ 'DNA' ] )
@check_if_uptodate( check_file_exists_for_itd_detection )
@files( generate_params_for_itd_detection )
def stage_9( control_file_list, tumor_file_list, control_output_dir_list, tumor_output_dir_list ):
    return_value = itd_detection( control_file_list, tumor_file_list, control_output_dir_list, tumor_output_dir_list )
    if not return_value:
        raise Exception( 'stage_9 failed.' )

#####################################################################
#
# STAGE 10 sv_detection (parse)
#
#  in:  bam
#  out: .gz 
#
@follows ( stage_6 )
@active_if( 'sv_detection' in Geno.job.get_job( 'tasks' )[ 'DNA' ] )
@check_if_uptodate( check_file_exists_for_sv_parse )
@parallel( generate_params_for_sv_parse )
def sv_detection_parse( target_label, target_bam, target_outdir, match_use, match_bam ):


    ##########
    # generate sample config yaml file
    sv_sampleConf = {"target": {}, "matched_control": {}, "non_matched_control_panel": {}}
    sv_sampleConf["target"]["label"] = target_label 
    sv_sampleConf["target"]["path_to_bam"] = target_bam
    sv_sampleConf["target"]["path_to_output_dir"] = target_outdir
    sv_sampleConf["matched_control"]["use"] = match_use 
    sv_sampleConf["matched_control"]["path_to_bam"] = match_bam
    sv_sampleConf["non_matched_control_panel"]["use"] = False
    hOUT = open(target_outdir + "/" + target_label + ".yaml", "w")
    print >> hOUT, yaml.dump(sv_sampleConf, default_flow_style = False)
    hOUT.close()
    ##########

    return_code = True

    try:
        function_name = whoami()
        log.info( "#{function}".format( function = function_name ) )

        # Make shell script
        shell_script_full_path = make_script_file_name( function_name, Geno )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( wgs_res.sv_parse_filt.format(
                                     log = Geno.dir[ 'log' ],
                                     pythonhome = Geno.conf.get( 'ENV', 'PYTHONHOME' ),
                                     ld_library_path = Geno.conf.get( 'ENV', 'LD_LIBRARY_PATH'),
                                     pythonpath = Geno.conf.get( 'ENV', 'PYTHONPATH' ),
                                     scriptdir = Geno.dir[ 'script' ],
                                     genomon_sv = Geno.conf.get( 'SOFTWARE', 'genomon_sv' ),
                                     method = "parse",
                                     sample_conf = target_outdir + "/" + target_label + ".yaml",
                                     param_conf = Geno.job.get_param( 'genomon_sv', 'param_file' )
                                     ) )
        shell_script_file.close()

        # Run
        return_code = Geno.RT.runtask(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ "sv_detection" ] )

        if return_code != 0:
            log.error( "{function}: runtask failed".format( function = function_name ) )
            raise

        save_status_of_this_process( "sv_parse", target_label, return_code, Geno, use_subdir )


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


@follows ( sv_detection_parse )
@active_if( 'sv_detection' in Geno.job.get_job( 'tasks' )[ 'DNA' ] )
@check_if_uptodate( check_file_exists_for_sv_filt )
@parallel( generate_params_for_sv_filt )
def sv_detection_filt( target_label, target_outdir ):

    return_code = True

    try:
        function_name = whoami()
        log.info( "#{function}".format( function = function_name ) )

        # Make shell script
        shell_script_full_path = make_script_file_name( function_name, Geno )
        shell_script_file = open( shell_script_full_path, 'w' )
        shell_script_file.write( wgs_res.sv_parse_filt.format(
                                     log = Geno.dir[ 'log' ],
                                     pythonhome = Geno.conf.get( 'ENV', 'PYTHONHOME' ),
                                     ld_library_path = Geno.conf.get( 'ENV', 'LD_LIBRARY_PATH'),
                                     pythonpath = Geno.conf.get( 'ENV', 'PYTHONPATH' ),
                                     scriptdir = Geno.dir[ 'script' ],
                                     genomon_sv = Geno.conf.get( 'SOFTWARE', 'genomon_sv' ),
                                     method = "filt",
                                     sample_conf = target_outdir + "/" + target_label + ".yaml",
                                     param_conf = Geno.job.get_param( 'genomon_sv', 'param_file' )
                                     ) )
        shell_script_file.close()
        
        # Run
        return_code = Geno.RT.runtask(
                            shell_script_full_path,
                            Geno.job.get_job( 'cmd_options' )[ "sv_detection" ] )
                            
        if return_code != 0:
            log.error( "{function}: runtask failed".format( function = function_name ) )
            raise
            
        save_status_of_this_process( "sv_filt", target_label, return_code, Geno, use_subdir )
        
        
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
#   mutation_filter
#
#   in:     txt
#   out:    xls or vcf
#
@follows( stage_9 )
@active_if ( 'mutation_filter' in Geno.job.get_job( 'tasks' )[ 'DNA' ] )
@files( generate_params_for_mutation_filter )
@check_if_uptodate( check_file_exists_for_mutation_filter )
def stage_mutation_filter(target_list, target_normal_bam, taraget_numor_bam, output_list, output_dir):
    return_value = mutation_filter(target_list, target_normal_bam, taraget_numor_bam, output_list, output_dir)

    if not return_value:
        raise


#####################################################################
#
#   STAGE 10 annotation
#
#   in:     txt
#   out:    xls or vcf
#
@follows( stage_mutation_filter )
@active_if ( 'annotation' in Geno.job.get_job( 'tasks' )[ 'DNA' ] )
@check_if_uptodate( check_file_exists_for_annotation )
@files( generate_params_for_annotation )
def stage_10(  input_file, output_file ):
    return_value = annotation(  input_file, output_file )
    if not return_value:
        raise Exception( 'stage_10 failed.' )

#####################################################################
#
#   LAST STAGE 
#
@follows( stage_10, sv_detection_filt )
def last_function():
    with log_mutex:
        log.info( "Genomon pipline has finished successflly!" )
    return True

