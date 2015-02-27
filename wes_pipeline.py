"""
wes_pipeline.py

"""

import sys
import os
import shutil
from datetime import datetime
from ruffus import *
from runtask import RunTask

#####################################################################
#
# Private modules
#
from __main__ import *
import genomon_rc as res
from utils import *

#####################################################################
#
#   STAGE 0 data preparation
#
starting_file_list = [
        [ './test/a_1_0.txt', './test/a_1_1.txt' ],
        [ './test/a_2_0.txt', './test/a_2_1.txt' ]
        ]


def check_file_exists(input_file, output_file):
    if not os.path.exists(output_file):
        return True, "Missing file %s" % output_file
    else:
        return False, "File %s exists" % output_file

#####################################################################
#
#   STAGE 1 split fastq
#
@parallel( starting_file_list )
@check_if_uptodate( check_file_exists )
def stage_1(
        input_file,
        output_file,
        ):
    """
        Stage 1

    """
    return_code = True

    try:
        log.info( "# {function} : {inputfile}, {outputfile}".format(
                    function = whoami(),
                    inputfile = input_file,
                    outputfile = output_file
                    ) )
        shutil.copyfile( input_file, output_file )

    except IOError as (errno, strerror):
        log.error( "{function}: I/O error({num}): {error}".format(
                        function = whoami(),
                        num = errno,
                        error = strerror)
                )
        return_code = False

    except ValueError:
        log.error( "{function}: ValueError".format(
                        function = whoami()
                    )
                )
        return_code = False

    except:
        log.error( "function: Unexpected error: {error} ".format(
                    function = whoami(),
                    error = sys.exc_info()[0]
                    )
                )
        return_code = False


    return return_code

#####################################################################
#
#   STAGE 2 fastq to bam
#
@transform( stage_1, suffix( "1.txt" ), "2.txt" )
@follows( stage_1 )
@check_if_uptodate( check_file_exists )
def stage_2(
        input_file,
        output_file
        ):
    """
        Stage 2

    """
    return_code = True

    try:
        log.info( "# {function} : {inputfile}, {outputfile}".format(
                    function = whoami(),
                    inputfile = input_file,
                    outputfile = output_file
                    ) )
        shutil.copyfile( input_file, output_file )

    except IOError as (errno, strerror):
        log.error( "{function}: I/O error({num}): {error}".format(
                        function = whoami(),
                        num = errno,
                        error = strerror)
                )
        return_code = False

    except ValueError:
        log.error( "{function}: ValueError".format(
                        function = whoami()
                    )
                )
        return_code = False

    except:
        log.error( "function: Unexpected error: {error} ".format(
                    function = whoami(),
                    error = sys.exc_info()[0]
                    )
                )
        return_code = False


    return return_code


#####################################################################
#
#   STAGE 3
#
@transform( stage_2, suffix( "2.txt" ), "3.txt" )
@follows( stage_2 )
@check_if_uptodate( check_file_exists )
def stage_3(
        input_file,
        output_file
        ):
    """
        Stage 3

    """
    return_code = True

    try:
        log.info( "# {function} : {inputfile}, {outputfile}".format(
                    function = whoami(),
                    inputfile = input_file,
                    outputfile = output_file
                    ) )
        shutil.copyfile( input_file, output_file )

    except IOError as (errno, strerror):
        log.error( "{function}: I/O error({num}): {error}".format(
                        function = whoami(),
                        num = errno,
                        error = strerror)
                )
        return_code = False

    except ValueError:
        log.error( "{function}: ValueError".format(
                        function = whoami()
                    )
                )
        return_code = False

    except:
        log.error( "function: Unexpected error: {error} ".format(
                    function = whoami(),
                    error = sys.exc_info()[0]
                    )
                )
        return_code = False


    return return_code

#####################################################################
#
#   STAGE 4
#
@transform( stage_3, suffix( "3.txt" ), "4.txt" )
@follows( stage_3 )
@check_if_uptodate( check_file_exists )
def stage_4(
        input_file,
        output_file
        ):
    """
        Stage 4

    """
    return_code = True

    try:
        log.info( "# {function} : {inputfile}, {outputfile}".format(
                    function = whoami(),
                    inputfile = input_file,
                    outputfile = output_file
                    ) )
        shutil.copyfile( input_file, output_file )

    except IOError as (errno, strerror):
        log.error( "{function}: I/O error({num}): {error}".format(
                        function = whoami(),
                        num = errno,
                        error = strerror)
                )
        return_code = False

    except ValueError:
        log.error( "{function}: ValueError".format(
                        function = whoami()
                    )
                )
        return_code = False

    except:
        log.error( "function: Unexpected error: {error} ".format(
                    function = whoami(),
                    error = sys.exc_info()[0]
                    )
                )
        return_code = False


    return return_code


#####################################################################
#
#   LAST STAGE 
#

@follows( stage_4 )
def last_function():
    log.info( "Genomon pipline has finished successflly!" )
    return True
