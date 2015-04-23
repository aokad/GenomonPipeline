#  Copyright Human Genome Center, Institute of Medical Science, the University of Tokyo
#  @since 2015

from __main__ import *
from utils import *

"""
    Genomon Sample data object

"""

class Sample():
    def __init__( self ):
        self.__sample_id_list = [ ]
        self.__current_sample_id = 0
        self.__param = []
        self.__f_subdir = False

    def __del__( self ):
        pass

    def param ( self, task_id ):
        if len( self.__param ) > 0:
            for id in range( len( self.__sample_id_list ) ):
                if self.__sample_id_list[ id ] == task_id:
                    break
            if self.__sample_id_list[ id ] == task_id:
                return self.__param[ id ]
            else:
                return self.__param[ self.__current_sample_id ]
        else:
            return None

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
    #   bam  => txt                             : out_ext = '.txt',   num_in = 1, num_out = 1
    #   (control_bam, diseaes_bam )  => txt     : out_ext = '.txt',   num_in = 2, num_out = 1
    #
    def make_param(
            self,
            task_name,
            out_ext,
            out_type,
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
            self.__sample_id_list.append( 'start' )
            self.get_starting_files()

        if not ( task_name in self.__sample_id_list ):
            self.__current_sample_id += 1
            self.__sample_id_list.append( task_name )
            self.__param.append( [] )

            for infile1, infile2, outfile1, outfile2 in self.previous():
                file_ext = Geno.job.get( 'file_ext' )
                if self.__current_sample_id == 1 and file_ext:
                    file_type = Geno.job.get( 'input_file_type' )
                    if file_type in ( 'paired_fastq', 'single_fastq'):
                        replace_ext = '.fastq'
                    if file_type in ( 'bam'):
                        replace_ext = '.bam'
                    tmp_out1 = outfile1.replace( file_ext, replace_ext )
                    tmp_out2 = outfile2.replace( file_ext, replace_ext )
                else:
                    tmp_out1 = outfile1
                    tmp_out2 = outfile2

                if self.__f_subdir:
                    subdir = '/' + os.path.basename( os.path.dirname( tmp_out1 ) ) + '/'
                else:
                    subdir = '/'

                if num_in == 1 and num_out == 2:
                    real_outfile1 = make_sample_file_name( tmp_out1,
                                                         "{dir}{subdir}{base}_1{ext}",
                                                         dir = Geno.dir[ out_type ],
                                                         subdir = subdir,
                                                         ext = out_ext )
                    real_outfile2 = make_sample_file_name( tmp_out1,
                                                         "{dir}{subdir}{base}_2{ext}",
                                                         dir = Geno.dir[ out_type ],
                                                         subdir = subdir,
                                                         ext = out_ext )
                else:
                    real_outfile1 = make_sample_file_name( tmp_out1,
                                                         "{dir}{subdir}{base}{ext}",
                                                         dir = Geno.dir[ out_type ],
                                                         subdir = subdir,
                                                         ext = out_ext )
                    real_outfile2 = make_sample_file_name( tmp_out2,
                                                         "{dir}{subdir}{base}{ext}",
                                                         dir = Geno.dir[ out_type ],
                                                         subdir = subdir,
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
        pair_id

        sample_subdir:
        control_subdir:
        disease_subdir

        bed_file
        
        """

        file_type = Geno.job.get( 'input_file_type' )

        try:
            #
            # A) Paired-end fastq files
            #
            if file_type == 'paired_fastq':
                fastq_file = [ None, None, None ]
                pair_id_list = Geno.job.get( 'pair_id' )

                glob_file_list = []
                i = 0
                for pair_id in pair_id_list:
                    glob_file_list.append( self.filename_list( Geno.job.get( 'file_name' ).format( pair_id = pair_id ) ) )

                    if 0 == len( glob_file_list[ i ] ):
                        log.error( "{function} failed.".format( function = whoami() ) )
                        log.error( "Paired fastq files not found." )
                        raise

                    i += 1

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
                filename_list = self.filename_list( Geno.dir( 'file_name' ) )

                if 0 == len( filename_list ):
                    log.error( "{function} failed.".format( function = whoami() ) )
                    log.error( "Single fastq files not found." )
                    raise
                else:
                    for file_name in filename_list:
                        self.current().append( [ '', '', file_name, '' ] )

            #
            # C) bam files
            #
            elif file_type == 'bam':
                filename_list = self.filename_list( Geno.dir( 'file_name' ) )

                if 0 == len( filename_list ):
                    log.error( "{function} failed.".format( function = whoami() ) )
                    log.error( "Bam files not found." )
                    raise
                else:
                    for file_name in filename_list:
                        self.current().append( [ '', '', file_name, '' ] )


        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            log.error( "{function} failed.".format( function = whoami() ) )
            log.error("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )


    def filename_list( self, filename_format ):
        fastq_type_list = []
        fastq_file_list = []

        if not filename_format:
            return fastq_file_list

        sample_subdir  = Geno.job.get( 'sample_subdir' )
        control_subdir = Geno.job.get( 'control_subdir' )
        disease_subdir = Geno.job.get( 'disease_subdir' )

        if sample_subdir:
            self.__f_subdir = True
            fastq_type_list.append( "{subdir}/{filename}".format(
                                        subdir = sample_subdir,
                                        filename = filename_format ) )
        elif control_subdir and disease_subdir:
            self.__f_subdir = True
            fastq_type_list.append( "{subdir}/{filename}".format(
                                        subdir = control_subdir,
                                        filename = filename_format ) )
            fastq_type_list.append( "{subdir}/{filename}".format(
                                        subdir = disease_subdir,
                                        filename = filename_format ) )

        else:
            self.__f_subdir = False
            fastq_type_list.append( filename_format )

        for fastq_file in fastq_type_list:
            fastq_file_list = fastq_file_list + sorted( glob( Geno.dir[ 'data' ] + '/' + fastq_file ) )

        return fastq_file_list
