"""
test_pipeline.py

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
from resource import genomon_rc as res
from resource import star_resource as star_res
from utils import *
from sample import Sample

def check_file_exists(
    input_file,
    output_file
    ):

    if not os.path.exists( output_file ):
        return True, 'Missing {file}'.format( file = output_file )
    else:
        return False, '{file} exists.'.format( file = output_file )

def generate_params_for_stage_1():
    for input_file, output_file in ( ( '~/tmp/test1a.txt', '~/tmp/test1b.txt' ),
                                     ( '~/tmp/test2a.txt', '~/tmp/test2b.txt' ) ):
        yield input_file, output_file

def generate_params_for_stage_2():
    for input_file, output_file in ( ( '~/tmp/test1b.txt', '~/tmp/test1c.txt' ),
                                     ( '~/tmp/test2b.txt', '~/tmp/test2c.txt' ) ):
        yield input_file, output_file


@active_if( 'stage_1' in Geno.job.get_job( 'tasks' )[ 'TEST' ] )
@files( generate_params_for_stage_1 )
@check_if_uptodate( check_file_exists )
def stage_1(
        input_file,
        output_file,
        ):
    print "stage_1"
    shutil.copyfile( os.path.expanduser( input_file ), os.path.expanduser( output_file ) )
    return True

@follows( stage_1 )
@active_if( 'stage_2' in Geno.job.get_job( 'tasks' )[ 'TEST' ] )
@files( generate_params_for_stage_2 )
@check_if_uptodate( check_file_exists )
def stage_2(
        input_file,
        output_file
        ):
    print "stage_2"
    shutil.copyfile( os.path.expanduser( input_file ), os.path.expanduser( output_file ) )
    return True

@follows( stage_2 )
def last_function():
    log.info( "Genomon pipline has finished successflly!" )
    return True
