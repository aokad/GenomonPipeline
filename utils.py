"""
Utilty functions

"""
import os
import inspect
from datetime import datetime
from genomon_rc import file_timestamp_format as file_time
from __main__ import Geno

#####################################################################
#
# Private modules
#
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


########################################
def make_script_file_name( function_name ):
    """
    Split files by line number

    """

    now = datetime.now()
    shell_script_name = file_time.format(
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

    return shell_script_full_path

########################################
def make_sample_file_name( filename,
                           file_fmt,
                           dir = None,
                           subdir = None,
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

    if subdir:
        ret_name = file_fmt.format(
                        dir = dir,
                        subdir = subdir,
                        base = basename,
                        ext = ext )
    else:
        ret_name = file_fmt.format(
                        dir = dir,
                        base = basename,
                        ext = ext )
    return ret_name

