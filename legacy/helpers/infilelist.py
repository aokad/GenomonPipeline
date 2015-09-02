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
            self.__input_file = filename
            self.input = {}
            self.compare = []
            self.control_panel = {}
            self.err_msg = "Unknown"
            
            if filename == None:
                self.err_msg = 'input_file_list object: filename is not specified'
            
            else:
                result = self.parse_file()
                self.err_msg = result

        except Exception as err:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            
            err_text = "input_file_list.init: Unexpected error." \
                        + err \
                        + "{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno)

            self.err_msg = err_text

    #
    # each extension, call parse method
    # 
    def parse_file( self ):
        # file_extension is either csv, tsv, or xlsx
        err_msg = self.__input_file + ":"

        file_extension = os.path.splitext( self.__input_file )[ 1 ]

        # CSV
        if file_extension.lower() == '.csv':
            typ, msg = self.parse_data(self.parse_csv())
            if (typ == False):
                typ, msg = self.parse_data(self.parse_tsv())
                if (typ == False):
                    err_msg += "list file is invalid format.\n"
                elif (len(msg) > 0):
                    err_msg += "list file is invalid format.\n" + msg
                else:
                    # success type tsv.
                    err_msg = ""
            elif (len(msg) > 0):
                err_msg += "list file is invalid format.\n" + msg 
            else:
                # success type csv
                err_msg = ""

        # TSV
        elif file_extension.lower() == '.tsv':
            typ, msg = self.parse_data(self.parse_tsv())
            if (typ == False):
                typ, msg = self.parse_data(self.parse_csv())
                if (typ == False):
                    err_msg += "list file is invalid format.\n"
                elif (len(msg) > 0):
                    err_msg += "list file is invalid format.\n" + msg
                else:
                    # success type csv.
                    err_msg = ""
            elif (len(msg) > 0):
                err_msg += "list file is invalid format.\n" + msg 
            else:
                # success type tsv
                err_msg = ""

        # EXCEL
        elif file_extension.lower() == '.xlsx':
            typ, msg = self.parse_data(self.parse_xlsx())
            if (len(msg) > 0):
                err_msg += "list file is invalid xlsx format.\n" + msg
            else:
                err_msg = ""
                
        # old EXCEL
        elif file_extension.lower() == '.xls':
            err_msg += "This software does not support to xls format.\n"
        
        # OTHER
        else:
            typ, msg1 = self.parse_data(self.parse_csv())
            if (typ == False):
                typ, msg2 = self.parse_data(self.parse_tsv())
                if (typ == False):
                    err_msg += "list file is invalid format.\n" \
                              + "### try read csv format ###\n" \
                              + msg1 \
                              + "### try read tsv format ###\n" \
                              + msg2
                              
                elif (len(msg2) > 0):
                    err_msg += "list file is invalid format.\n" + msg2
                else:
                    # success type tsv.
                    err_msg = ""
            elif (len(msg1) > 0):
                err_msg += "list file is invalid format.\n" + msg1
            else:
                # success type csv
                err_msg = ""

        return err_msg

    #
    # read tsv file, convert to list object
    #    
    def parse_tsv( self ):

        with open( self.__input_file ) as f:
            data = []
            trim = []
            for line in f:
                cells = line.split('\t')
                
                temp = []
                for cell in cells:
                    temp.append(cell)
                
                data.append(temp)

            trim = self.trim_data(data)
        
        return trim

    #
    # read csv file, convert to list object
    #     
    def parse_csv( self ):
        
        import csv
        
        with open( self.__input_file ) as f:
            csv_obj = csv.reader( f )
            data = []
            trim = []
            for cells in csv_obj:
                temp = []
                for cell in cells:
                    temp.append(cell)

                data.append(temp)

            trim = self.trim_data(data)
            
        return trim

    #
    # read xlsx file, convert to list object
    # 
    def parse_xlsx( self ):
        
        import xlrd
        
        workbook = xlrd.open_workbook(  self.__input_file )
        worksheet = workbook.sheet_by_index( 0 )
        
        data = []
        
        for rownum in xrange( worksheet.nrows ):
            cells = worksheet.row_values( rownum )
 
            temp = []
            for cell in cells:
                text = self.data2str(cell)
                text_split = text.split(',')
                
                for t in text_split:
                    temp.append(t)

            data.append(temp)

        return self.trim_data(data)

    #
    # parse list data
    #
    def parse_data( self, data ):
    
        err_msg = ''
        mode = ''
        
        num_compare = 0
        num_controlpanel = 0
        
        # input file check flgs
        num_input = 0
        flg_input = False
        chk_fastq = False
        chk_bam = False
        chk_pair = False
        chk_single = False
        
        self.input = {}
        self.compare = []
        self.control_panel = {}
 
        for row in data:
            header = row[0].lower()
            # header
            if '[input]' in header:
                mode = 'input'
                flg_input = True
                continue
            elif '[compare]' in header:
                mode = 'compare'
                continue
            elif '[controlpanel]' in header:
                mode = 'controlpanel'
                continue

            # section data
            if mode == 'input':
                row = self.trim_row(row)

                if len(row) < 2:
                    err_msg += "[input] section, column is too few. " + self.get_rowtext(row) + "\n"
                if len(row) > 3:
                    err_msg += "[input] section, column is too many. " + self.get_rowtext(row) + "\n"
                else:
                    files = []
                    if len(row) == 2:
                        chk_single = True
                    if len(row) == 3:
                        chk_pair = True
                        
                    f_ckeck = []    
                    for i in range(1,len(row)):
                        files.append(row[i].replace(';', ','))
                    
                        f_split = row[i].split(";")
                        
                        for f in f_split:
                            if not os.path.exists( f ):
                                err_msg += "[input] section, No.{0} file is not exits. [{file}]\n".format(num_input+1, file=f)

                            file_extension = os.path.splitext(f)[1]
                            if file_extension.lower() == ".bam":
                                chk_bam = True
                            else:
                                chk_fastq = True
#                            else:
#                                err_msg += "[input] section, No.{0} file is not fastq/bam. [{file}]\n".format(num_input+1, file=f)
                            
                            # check duplex
                            if f in f_ckeck:
                                err_msg += "[input] section, No.{0} file is duplex. [{file}]\rn".format(num_input+1, file=f)
                                
                            f_ckeck.append(f)
                            
                    self.input[ row[0] ] = files
                    num_input += 1
                    
            elif mode == 'compare':
                if len(row) < 3:
                    err_msg += "[compare] section, No.{0} column is too few. ".format(num_compare+1) + self.get_rowtext(row) + "\n"

                else:
                    normal        = row[0] if len(row[0]) >= 1 and row[0] else None
                    tumor         = row[1] if len(row[1]) >= 1 and row[1] else None
                    controlpanel  = row[2] if len(row[2]) >= 1 and row[2] else None
                    
                    if normal or tumor or controlpanel:
                        self.compare.append( ( normal, tumor, controlpanel) )
                        num_compare += 1

            elif mode == 'controlpanel':
                row = self.trim_row(row)

                if len(row) < 2:
                    err_msg += "[controlpanel] section, list item is none. " + self.get_rowtext(row) + "\n"
                else:
                    control_panel_key = row[0]
                    self.control_panel[ control_panel_key ] = []
                    for item in row[ 1: ]:
                        self.control_panel[ control_panel_key ].append( item )

                    if len(row)-1 > 20:
                        err_msg += "[controlpanel] section, list item is over 20. " + control_panel_key + "\n"
                    
                    num_controlpanel += 1
                
        if (flg_input == True):
            if (num_input < 1):
                err_msg += "[input] section, data is none.\n"
            if (chk_single == True and chk_pair == True):
                err_msg += "[input] section, cannot mix single/pair.\n"
            if (chk_fastq == True and chk_bam == True):
                err_msg += "[input] section, cannot mix fastq/bam.\n"
            if (chk_bam == True and chk_pair == True):
                err_msg += "[input] section, cannot use bam to pair.\n"
        if num_compare < 1:
            err_msg += "[compare] section, data is none.\n"
#        if num_controlpanel < 1:
#            err_msg += "[controlpanel] section, data is none.\n"

        typ = True
        if num_input == 0 and num_compare == 0 and num_controlpanel == 0:
            typ = False
            
        return typ, err_msg

    #
    # remove blank row, comment row '#', white space, end code'\n'
    #
    # list 
    #  [[A, B, C], [D, E, ,F], [, , ,], [G, ,]]
    #  to
    #  [[A, B, C], [D, E, ,F], [G, ,]]
    #
    def trim_data( self, data ):
    
        trim = []
        
        for row in data:
            
            temp = []
            not_blank = False
            
            for text in row:
                # remove space
                text = text.replace( ' ', '' )
            
                # remove end code
                text = text.rstrip("\r\n")

                if (len(text) > 0):
                    # skip comment
                    if text[0] == '#':  # comment
                        break
                    
                    not_blank = True
                    
                temp.append(text)
            
            if (not_blank == True):
                trim.append(temp)
    
        return trim

    #
    # list [A, B, ,C] -> [A, B, C]
    #        
    def trim_row( self, row ):
        trim = []

        for cell in row:
            if (len(cell) > 0):
                trim.append(cell)
        
        return trim
    
    #
    # list [A, B, ,C] -> str "[A, B, C]"
    #
    def get_rowtext( self, row ):
        text = ""

        for cell in row:
            if (len(text) == 0):
                text = "["
            else:
                text += ", "
                
            text += "'" + cell + "'"
        
        text += "]"
        
        return text
        
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

    def get_err_msg( self ):
        return self.err_msg


