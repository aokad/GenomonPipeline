#  Copyright Human Genome Center, Institute of Medical Science, the University of Tokyo
#  @since 2015

from __main__ import *
from utils import *

"""
    Genomon Sample data object

"""

class Sample():
    def __init__( self ):
        self.__sample_id_list = [ 'start' ]
        self.__current_sample_id = 0
        self.__param = []

    def __del__( self ):
        pass


    def current( self ):
        if len( self.__param ) > 0:
            return self.__param[ self.__current_sample_id ]
        else:
            return None

    def previous( self ):
        if self.__current_sample_id > 0:
            return self.__param[ self.__current_sample_id - 1 ]
        else:
            return None

    #
    # File name generator
    # -------------------
    # Pattern
    #   bam => fastq                            : out_ext = '.fastq', num_in = 1, num_out = 1
    #   bam => ( fastq1, fastq2 )               : out_ext = '.fastq', num_in = 1, num_out = 2
    #   fastq => bam                            : out_ext = '.bam',   num_in = 1, num_out = 1
    #   ( fastq1, fastq2 ) => bam               : out_ext = '.bam',   num_in = 2, num_out = 1
    #   ( fastq1, fastq2 ) => (fastq1, fastq2 ) : out_ext = '.fastq', num_in = 2, num_out = 2
    #
    def make_param(
            self,
            task_name,
            out_ext,
            num_in,
            num_out
            ):
        """
        Generate parameter list with output file extension

        Create a list of the following
            input_file
            output_file

        """
        #
        # First get starting files
        #
        if 0 == self.__current_sample_id:
            self.__param.append( [] )
            self.get_starting_files()
            self.__sample_id_list.append( 'start' )

        if not ( task_name in self.__sample_id_list ):
            self.__current_sample_id += 1
            self.__sample_id_list.append( task_name )
            self.__param.append( [] )
            param_list = self.previous()
            for infile1, infile2, outfile1, outfile2 in param_list:
                if num_in == 1 and num_out == 1:
                    real_outfile = make_sample_file_name( outfile1,
                                                         "{dir}/{base}{ext}",
                                                         dir = Geno.dir[ out_ext[1: ] ],
                                                         ext = out_ext )
                    return_list = [ outfile1, None, real_outfile, None ]
                    self.__param[ self.__current_sample_id ].append( return_list )

                elif num_in == 1 and num_out == 2:
                    real_outfile1 = make_sample_file_name( outfile1,
                                                         "{dir}/{base}_1{ext}",
                                                         dir = Geno.dir[ out_ext[1: ] ],
                                                         ext = out_ext )
                    real_outfile2 = make_sample_file_name( outfile1,
                                                         "{dir}/{base}_2{ext}",
                                                         dir = Geno.dir[ out_ext[1: ] ],
                                                         ext = out_ext )
                    return_list = [ outfile1, None, real_outfile1, real_outfile2 ]
                    self.__param[ self.__current_sample_id ].append( return_list )

                elif num_in == 2 and num_out == 1:
                    real_outfile = make_sample_file_name( outfile1,
                                                         "{dir}/{base}{ext}",
                                                         dir = Geno.dir[ out_ext[1: ] ],
                                                         ext = out_ext )
                    return_list = [ outfile1, outfile2, real_outfile, None ]
                    self.__param[ self.__current_sample_id ].append( return_list )

                elif num_in == 2 and num_out == 2:
                    real_outfile1 = make_sample_file_name( outfile1,
                                                         "{dir}/{base}{ext}",
                                                         dir = Geno.dir[ out_ext[1: ] ],
                                                         ext = out_ext )
                    real_outfile2 = make_sample_file_name( outfile2,
                                                         "{dir}/{base}{ext}",
                                                         dir = Geno.dir[ out_ext[1: ] ],
                                                         ext = out_ext )
                    return_list = [ outfile1, outfile2, real_outfile1, real_outfile2 ]
                    self.__param[ self.__current_sample_id ].append( return_list )

    def get_starting_files( self ):
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
                    self.current().append( [
                                    '',
                                    '',
                                    glob_file_list[ 0 ][ i ],
                                    glob_file_list[ 1 ][ i ]
                                    ] )

            #
            # B) Single-end fastq files
            #
            elif file_type == 'single_fastq':
                fastq_file = Geno.job.get( 'file_name' )

                if fastq_file != None:
                    for file_name in sorted( glob( Geno.dir[ 'data' ] + fastq_file ) ):
                        self.__param[ self.__current_sample_id ].append( [ '', '', file_name, '' ] )
                else:
                    log.error( "{function} failed.".format( function = whoami() ) )
                    log.error( "Single fastq files not found." )
                    raise

            #
            # C) bam files
            #
            elif file_type == 'bam':
                bam_files = Geno.job.get( 'file_name' )
                for file_name in sorted( glob( Geno.dir[ 'data' ] + bam_files ) ):
                    self.__param[ self.__current_sample_id ].append( [ '', '', file_name, '' ] )

                if 0 == len( self.current() ):
                    log.error( "{function} failed.".format( function = whoami() ) )
                    log.error( "Bam files not found." )
                    raise

        except:
            log.error( "{function} failed.".format( function = whoami() ) )
            log.error( "Unexpected error: {1}".format( sys.exc_info()[0] ) )

