#! /usr/bin/python
"""

This is a Python parser module for annotation database


"""
import random
from genomesize import GenSize

################################################################################
#
# Bed file parser
#
################################################################################
class BedExtract:

    def __init__( self, filename, genome_size, bin_size = 1000, chr_str = False ):
        self.data = {}

        try:
            file = open( filename, 'r' )
            Error = self.parse_file( file )
            file.close()

            self.gen_size = GenSize( genome_size )
            self.bin_size = bin_size
            self.chr_str = chr_str

        except ValueError:
            print "Error open/close/parse file"

        if False == Error:
            print "Error parse database file"

    def parse_file( self, file ):

        for line in file:
            line_split = line.split()
            if not line_split[ 0 ] in self.data:
                self.data[ line_split[ 0 ] ] = {}
            self.data[ line_split[ 0 ] ][ line_split[ 1 ] ]= line_split[ 2 ]

        return True

    def get_chr( self ):

        return self.data.keys()

    def get_random_interval( self ):
        return_value = [ None, None ]
        while not return_value[ 0 ]:
            chr_tmp = random.randrange( 1, 24 )
            if chr_tmp == 23 or chr_tmp == 24:
                chr_tmp = ( 'X', 'Y' )[ chr_tmp - 23 ]

            if self.chr_str:
                chr_return = 'chr' + str( chr_tmp )
            else:
                chr_return = str( chr_tmp )

            pos = random.randrange( 0, int( self.gen_size.get_size( chr_return ) ) - self.bin_size )
            if self.get_pos( chr_return, pos, self.bin_size ):

                return_value = [ chr_return, pos ]

        return return_value

    def get_pos( self, chr, pos, size ):

        return_value = False

        if self.chr_str:
            chr = str( chr )
        else:
            chr = 'chr' + str( chr )

        if chr in self.data:
            for start in sorted( self.data[ chr ].keys() ):
                #print "L:{0}:{1}-{2}-{3}".format( chr, start, pos, self.data[ chr ] [ start ] )
                if int( start ) <= int( pos ) and int( pos ) + self.bin_size <= int( self.data[ chr ][ start ] ):
                    return_value = True      
                    #print "HIT!"
                    break

        return return_value

