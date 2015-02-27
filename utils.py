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

