#! /usr/bin/python
"""

This is a Python parser module for annotation database


"""
################################################################################
#
# genome size file parser
#
################################################################################
class GenSize:

    def __init__( self, filename ):
        self.data = {}

        try:
            file = open( filename, 'r' )
            Error = self.parse_file( file )
            file.close()

        except ValueError:
            print "Error open/close/parse file"

        if False == Error:
            print "Error parse database file"

    def parse_file( self, file ):

        for line in file:
            line_split = line.split()
            self.data[ line_split[ 0 ] ] = line_split[ 1 ]

        return True

    def get_chr( self ):

        return self.data.keys()

    def get_size( self, chr ):

        if chr in self.data:
            return_value = self.data[ chr ]
        else:
            return_value = 0

        return return_value
