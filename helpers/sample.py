#  Copyright Human Genome Center, Institute of Medical Science, the University of Tokyo
#  @since 2015

from __main__ import *
from utils import *

"""
    Genomon Sample data object

"""

class Sample():
    def __init__( self ):
        self.__sample_list = None
        self.__sample_id_list = [ ]
        self.__current_sample_id = 0
        self.__param = []
        self.__f_subdir = False
        self.__f_filelist = False

    def __del__( self ):
        pass

    def set_sample_list( self, sample_list ):
        self.__sample_list = sample_list

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
    #   (control_bam, diseaes_bam )  => txt     : out_ext = '.txt',   num_in = 2, num_out = 1
    #
    def make_param(
            self,
            task_name,
            from_task_name,
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
            if not self.get_starting_files():
                return False

        if not ( task_name in self.__sample_id_list ):
            if from_task_name != None:
                param_list = self.param( from_task_name )
            else:
                param_list = self.current()

            if param_list == None:
                return False

            self.__current_sample_id += 1
            self.__sample_id_list.append( task_name )
            self.__param.append( [] )

            for infile1, infile2, outfile1, outfile2 in param_list:
                file_ext = Geno.job.get_job( 'file_ext' )
                if self.__current_sample_id == 1 and file_ext:
                    file_type = Geno.job.get_job( 'input_file_type' )
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
                elif self.__f_filelist and self.__current_sample_id == 1:
                    subdir = '/' + os.path.split( os.path.split( tmp_out1 )[ 0 ] )[ 1 ] + '_'
                elif self.__sample_list:
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


        return True

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

        sample_name:
        control_disease_pairs:

        bed_file
        
        """

        file_type = Geno.job.get_job( 'input_file_type' )

        try:

            if self.__sample_list:
                for fastq1_list, fastq2_list in self.__sample_list:
                    for fastq1, fastq2 in zip ( fastq1_list.split( ',' ), fastq2_list.split( ',' ) ):
                        self.current().append( ('', '', fastq1, fastq2 ) )

            #
            # A) Paired-end fastq files
            #
            elif file_type == 'paired_fastq':
                pair_id_list = Geno.job.get_job( 'pair_id' )

                glob_file_list = []
                for pair_id in pair_id_list:
                    glob_file_list.append( self.filename_list( pair_id = pair_id ) )

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
            #
            # C) bam files
            #
            elif file_type == 'single_fastq' or file_type == 'bam':
                filename_list = self.filename_list()

                if 0 == len( filename_list ):
                    log.error( "{function} failed.".format( function = whoami() ) )
                    log.error( "Files not found." )
                    raise
                else:
                    for file_name in filename_list:
                        self.current().append( [ '', '', file_name, '' ] )


        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            log.error( "{function} failed.".format( function = whoami() ) )
            log.error("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )
            return False

        return True

    def filename_list( self, pair_id = None ):
        input_type_list = []
        input_file_list = []

        #
        # Get resoruce from job yaml file
        #
        filename_format_tmp = Geno.job.get_job( 'file_name' )
        sample_name  = Geno.job.get_job( 'sample_name' )

        if not filename_format_tmp:
            return input_file_list

        #
        # Support single file format or list of file format
        #
        # example 1)
        #   file_name:  'R{pair_id}.fastq'
        #
        # example 2)
        #   file_name:  
        #               -'test1/R{pair_id}.fastq'
        #               -'test2/R{pair_id}.fastq'
        #
        if isinstance( filename_format_tmp, list ):
            self.__f_filelist = True
            file_format_list = filename_format_tmp
        else:
            file_format_list = [ filename_format_tmp ]

        for filename_format in file_format_list:
            if pair_id:
                filename_format = filename_format.format( pair_id = pair_id )

            if sample_name:
                self.__f_subdir = True
                input_type_list.append( "{subdir}/{filename}".format(
                                            subdir = sample_name,
                                            filename = filename_format ) )
            else:
                self.__f_subdir = False
                input_type_list.append( filename_format )

        for input_file in input_type_list:
            input_file_list = input_file_list + sorted( glob( Geno.dir[ 'data' ] + '/' + input_file ) )

        return input_file_list

    def subdir_exists( self ):
        return self.__f_subdir

