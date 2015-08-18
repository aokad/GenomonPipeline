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
            self.compare = []
            self.control_panel = {}

            if not self.parse_file():
                raise Exception( 'Failed to parse input file.' )

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
            else:
                return_value = False

            return return_value

    def parse_tsv( self ):
        return_value = True
        with open( self.__input_file ) as f:
            mode = ''
            for line in f:
                line_l = line.strip().replace( ' ', '' ).lower()
                if line_l:
                    if line_l[ 0:7 ] == '[input]': # header
                        mode = 'input'
                        continue
                    elif line_l[ 0:9 ] == '[compare]': # header
                        mode = 'compare'
                        continue
                    elif line_l[ 0:15 ] == '[controlpanel]': # header
                        mode = 'controlpanel'
                        continue
                    elif line_l[ 0 ] == '#': # comment
                        continue

                if line_l and mode == 'input':
                    line_split = line.replace( ';', ',' ).split( '\t' )
                    self.input[ line_split[ 0 ].strip() ] = [ line_split[ 1 ].strip(), line_split[ 2 ].strip() ]

                elif line and mode == 'compare':
                    line_split = line.replace( ';', ',' ).split( '\t' )
                    line_split = map( str.strip, line_split )
                    normal = line_split[ 0 ].strip() if len( line_split ) >= 1 and line_split[ 0 ] else None
                    tumor = line_split[ 1 ].strip() if len( line_split ) >= 2 and line_split[ 1 ] else None
                    control_panel = line_split[ 2 ].strip() if len( line_split ) == 3  and line_split[ 2 ] else None

                    if normal or tumor or control_panel:
                        self.compare.append( ( normal, tumor, control_panel) )

                elif line_l and mode == 'controlpanel':
                    line_split = line.split( '\t' )
                    control_panel_key = line_split[ 0 ].strip()
                    self.control_panel[ control_panel_key ] = []
                    for item_tmp in line_split[ 1: ]:
                        for item in item_tmp.split( ',' ):
                            self.control_panel[ control_panel_key ].append( item.strip() )

        return return_value


    def parse_csv( self ):
        import csv
        return_value = True
        with open( self.__input_file ) as f:
            csv_obj = csv.reader( f )
            mode = ''
            for line in csv_obj:

                if line:
                    data_type = line[ 0 ].strip().replace( ' ', '' )
                    data_type_l = data_type.lower()
                    if data_type_l == '[input]': # header
                        mode = 'input'
                        continue
                    elif data_type_l == '[compare]': # header
                        mode = 'compare'
                        continue
                    elif data_type_l == '[controlpanel]': # header
                        mode = 'controlpanel'
                        continue
                    elif data_type_l and data_type_l[ 0 ] == '#': # header or comment
                        continue

                if line and mode == 'input':
                    data_type = line[ 0 ].strip().replace( ';', ',' )
                    fastq1 = line[ 1 ].strip().replace( ';', ',' )
                    fastq2 = line[ 1 ].strip().replace( ';', ',' )
                    self.input[ data_type ] = [ fastq1, fastq2 ]

                elif mode == 'compare':
                    normal = line[ 0 ].strip().replace( ';', ',' ) if len ( line ) >= 1 and line[ 0 ] else None
                    tumor = line[ 1 ].strip().replace( ';', ',' ) if len ( line ) >= 2 and line[ 1 ] else None
                    controlpanel = line[ 2 ].strip().replace( ';', ',' )  if len ( line ) == 3 and line[ 2 ] else None
                    if normal or tumor or controlpanel:
                        self.compare.append( ( normal, tumor, controlpanel) )
                                         

                elif line and mode == 'controlpanel':
                    control_panel_key = line[ 0 ].strip()
                    self.control_panel[ control_panel_key ] = []
                    for item in line[ 1: ]:
                        self.control_panel[ control_panel_key ].append( item.strip() )

        return return_value

    def parse_xlsx( self ):
        import xlrd

        return_value = True
        mode = ''

        workbook = xlrd.open_workbook(  self.__input_file )
        worksheet = workbook.sheet_by_index( 0 )

        for rownum in xrange( worksheet.nrows ):
            line = worksheet.row_values( rownum )
            data = self.data2str( line[ 0 ].strip().replace( ' ', '' ) ) if line else None
            data_l = data.lower()
            if data:
                if data_l == '[input]': # header
                    mode = 'input'
                    continue
                elif data_l == '[compare]': # header
                    mode = 'compare'
                    continue
                elif data_l == '[controlpanel]': # header
                    mode = 'controlpanel'
                    continue
                elif data[ 0 ] == '#': # header or comment
                    continue

            if data and mode == 'input':
                data = data.replace( ';', ',' )
                fastq1 = self.data2str( self.data2str( line[ 1 ].strip().replace( ';', ',' ) ) )
                fastq2 = self.data2str( self.data2str( line[ 2 ].strip().replace( ';', ',' ) ) )
                self.input[ data ] = [ fastq1, fastq2 ]

            elif mode == 'compare':
                data  = data.replace( ';', ',' ) if data else None
                data1 = self.data2str( line[ 1 ].strip().replace( ';', ',' ) ) if len( line ) >= 2 and line[ 1 ] else None
                data2 = self.data2str( line[ 2 ].strip().replace( ';', ',' ) ) if len( line ) == 3 and line[ 2 ]  else None
                if data or data1 or data2:
                    self.compare.append( ( data, data1, data2 ) )

            elif data and mode == 'controlpanel':
                self.control_panel[ data ] = []
                for item_tmp in line[ 1: ]:
                    for item in item_tmp.split( ',' ):
                        if item:
                            self.control_panel[ data ].append( self.data2str( item.strip() ) )


        return return_value


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

    def get_control( self ):
        return self.control_panel
