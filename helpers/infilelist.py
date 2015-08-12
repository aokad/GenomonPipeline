#  Copyright Human Genome Center, Institute of Medical Science, the University of Tokyo
#  @since 2015

"""
    Genomon input CSV/TSV/XLS file parser

"""

import sys
import os

class input_file_list:
    #
    # Interface
    #
    def __init__( self, filename ):

        try:
            if filename == None:
                raise Exception( 'input_file_list object: filename is not specified' )

            self.__input_file = filename
            self.input = {}
            self.compare = {}

            self.parse_file()

        except Exception as err:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print( "input_file_list.init: Unexpected error" )
            print( err )
            print ("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )
            raise

    def parse_file( self ):
        # file_extension is either csv, tsv, or xls

            file_extension = os.path.splitext( self.__input_file )[ 1 ]
            if file_extension.lower() == '.csv':
                return_value = self.parse_csv()
            elif file_extension.lower() == '.tsv':
                return_value = self.parse_tsv()
            elif file_extension.lower() == '.xlsx':
                return_value = self.parse_xlsx()

            return return_value

    def parse_tsv( self ):
        with open( self.__input_file ) as f:
            mode = ''
            for line in f:
                line = line.strip()
                if line:
                    if line[ 0:7 ] == '[Input]': # header
                        mode = 'input'
                    elif line[ 0:9 ] == '[Compare]': # header
                        mode = 'compare'
                    elif line[ 0 ] == '#': # comment
                        continue
                    elif line and mode == 'input':
                        line_split = line.replace( ';', ',' ).split( '\t' )
                        self.input[ line_split[ 0 ].strip() ] = [ line_split[ 1 ].strip(), line_split[ 2 ].strip() ]

                    elif line and mode == 'compare':
                        line_split = line.replace( ';', ',' ).split( '\t' )
                        self.compare[ line_split[ 0 ].strip() ] = line_split[ 1 ].strip()


    def parse_csv( self ):
        import csv
        with open( self.__input_file ) as f:
            csv_obj = csv.reader( f )
            mode = ''
            for line in csv_obj:
                if line:
                    data_type = line[ 0 ].strip()
                    if data_type == '[Input]': # header
                        mode = 'input'
                    elif data_type == '[Compare]': # header
                        mode = 'compare'
                    elif data_type[ 0 ] == '#': # header or comment
                        continue
                    elif mode == 'input':
                        data_type = line[ 0 ].strip().replace( ';', ',' )
                        fastq1 = line[ 1 ].strip().replace( ';', ',' )
                        fastq2 = line[ 1 ].strip().replace( ';', ',' )
                        self.input[ data_type ] = [ fastq1, fastq2 ]

                    elif mode == 'compare':
                        data_type = line[ 0 ].replace( ';', ',' )
                        self.compare[ data_type ] = line[ 1 ].strip().replace( ';', ',' )

    def parse_xlsx( self ):
        import xlrd
        mode = ''
        workbook = xlrd.open_workbook(  self.__input_file )
        worksheet = workbook.sheet_by_index( 0 )

        for rownum in xrange( worksheet.nrows ):
            x = worksheet.row_values( rownum )
            data = self.data2str( x[ 0 ].strip() ) if x else None
            if data:
                if data == '[Input]': # header
                    mode = 'input'
                elif data == '[Compare]': # header
                    mode = 'compare'
                elif data[ 0 ] == '#': # header or comment
                    continue
                elif mode == 'input':
                    data = data.replace( ';', ',' )
                    fastq1 = self.data2str( self.data2str( x[ 1 ].strip().replace( ';', ',' ) ) )
                    fastq2 = self.data2str( self.data2str( x[ 2 ].strip().replace( ';', ',' ) ) )
                    self.input[ data ] = [ fastq1, fastq2 ]

                elif mode == 'compare':
                    data = data.replace( ';', ',' )
                    self.compare[ data ] = self.data2str( x[ 1 ].strip().replace( ';', ',' ) )


    def data2str( self, x ):
        if isinstance( x, basestring ):
            x = x.replace( "\t", "" )
        if type( x ) == type ( u'' ):
            data =  x.encode( 'utf-8' )
        else:
            data = str( x )

        return data

    def get_sample_list( self ):
        return self.input

    def get_compare( self ):
        return self.compare

