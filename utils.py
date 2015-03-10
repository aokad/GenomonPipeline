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

########################################
def whoami():
    return inspect.stack()[1][3]

########################################
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
def make_script_file_name( function_name ):
    """
    Split files by line number

    """

    try:
        now = datetime.now()
        shell_script_name = res.file_timestamp_format.format(
                                        name=function_name,
                                        year=now.year,
                                        month=now.month,
                                        day=now.day,
                                        hour=now.hour,
                                        min=now.minute,
                                        msecond=now.microsecond )

        shell_script_full_path = "{script}/{file}.sh".format(
                                        script = Geno.dir[ 'script' ],
                                        file = shell_script_name )

    except IOError as ( errno, strerror ):
        log.error( "split_file failed." )
        log.error( "IOError {0}]{1}".format( errno, strerror ) )
    except:
        log.error( "split_file failed." )
        log.error( "Unexpected error: {1}".format( sys.exc_info()[0] ) )

    return shell_script_full_path

########################################
def make_sample_file_name( filename,
                           file_fmt,
                           dir = None,
                           ext = None ):
    """
    {dir}
    {base}
    {ext}
    """
    tmp_name = os.path.splitext( os.path.basename( filename ) )
    basename = tmp_name[ 0 ]
    if not ext:
        ext = tmp_name[ 1 ]
    if not dir:
        dir = os.path.dirname( filename )

    ret_name = file_fmt.format(
                    dir = dir,
                    base = basename,
                    ext = ext )
    return ret_name

