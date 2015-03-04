"""
Utilty functions

"""
import os
import inspect

#####################################################################
#
# Private modules
#
from __main__ import *

def whoami():
    return inspect.stack()[1][3]
def whosdaddy():
    return inspect.stack()[2][3]

########################################
def split_file( file_name, output_file_name, split_len ):
    """
    Split files by line number

    """

    try:
        input = open( file_name, 'r' )
        count = 0
        at = 0
        dest = None

        for line in input:
            if count % split_len == 0:
                if dest: dest.close()
                dest = open( output_file_name.format( num=at ) )
                at += 1
            dest.write( line )
            count += 1

        dest.close()
        input.close()

    except IOError as ( errno, strerror ):
        log.error( "split_file failed." )
        log.error( "IOError {0}]{1}".format( errno, strerror ) )
    except:
        log.error( "split_file failed." )
        log.error( "Unexpected error: {1}".format( sys.exc_info()[0] ) )

########################################
def get_starting_files():
    """
    Get the list of starting files from
    specified job configuration file.

    1) Get file name from job configuration files
    2) Glob files

    Reserved words
    ==============
    input_file_type:    paired_fastq|single_fastq|bam
    file_name:
    control_file_name:
    disease_filename
    pair_id
    bed_file
    
    """

    file_type = Geno.job.get( 'input_file_type' )
    bed_file  = Geno.job.get( 'bed_file' )

    file_list = []

    try:
        #
        # A) Paired-end fastq files
        #
        if file_type == 'paired_fastq':
            fastq_file = [ None, None, None ]
            pair_id_list = Geno.job.get( 'pair_id' )

            for pair_id in pair_id_list:
                fastq_file[ pair_id ] = ( Geno.job.get( 'file_name' ).format( pair_id = pair_id ) )
                if not fastq_file[ pair_id ]:
                    log.error( "{function} failed.".format( function = whoami() ) )
                    log.error( "job configuration file does not set file_name or pair_id properly." )
                    raise

            glob_file_list = []
            for pair_id in pair_id_list:
                glob_file_list.append( sorted( glob( Geno.dir[ 'data' ] + '/' + fastq_file[ pair_id ] ) ) )

            if 0 == len( glob_file_list ):
                log.error( "{function} failed.".format( function = whoami() ) )
                log.error( "Paired fastq files not found." )
                raise

            for i in range( len( glob_file_list[ 0 ] ) ):
                file_list.append( [
                                glob_file_list[ 0 ][ i ],
                                glob_file_list[ 1 ][ i ]
                                ] )

        #
        # B) Single-end fastq files
        #
        elif file_type == 'single_fastq':
            fastq_file = Geno.job.get( 'file_name' )

            if fastq_file != None:
                file_list = sorted( glob( Geno.dir[ 'data' ] + fastq_file ) )
            else:
                log.error( "{function} failed.".format( function = whoami() ) )
                log.error( "Single fastq files not found." )
                raise

        #
        # C) bam files
        #
        elif file_type == 'bam':
            bam_files = Geno.job.get( 'file_name' )
            file_list = sorted( glob( Geno.dir[ 'data' ] + bam_files ) )

            if 0 == len( file_list ):
                log.error( "{function} failed.".format( function = whoami() ) )
                log.error( "Bam files not found." )
                raise

    except:
        log.error( "{function} failed.".format( function = whoami() ) )
        log.error( "Unexpected error: {1}".format( sys.exc_info()[0] ) )

    return file_list

