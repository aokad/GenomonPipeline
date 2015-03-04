#  Copyright Human Genome Center, Institute of Medical Science, the University of Tokyo
#  @since 2015

import yaml
from __main__ import *

"""
    Genomon job configuration file parse
    How to use it:
        ge_job = genomonjob( 'path to the job configuration file' )
        ge_job.get( 'bwa' )

"""
class genomon_job:
    def __init__( self, job_file = None, log = None ):

        self.__log = log
        if job_file != None:
            self.open_job( job_file )

    def open_job( self, job_file = None ):
        if job_file != None:
            self.__job_file = job_file

        try:
            if self.__job_file == None:
                self.__log.error( "genomon_job.get: job file is not loaded properly." )
                raise

            f = open( self.__job_file )
            self.__job = yaml.load( f )
            f.close()

        except IOError as (errno, stderror ):
            self.__log.error( "genomon_job.open_job: IOError: error number: {num}, std_error: {stderr}".format(
                        num = errno, stderr = stderror ) )
        except:
            self.__log.error( "genomon_job.open_job: unexpected error:", sys.exc_info()[0] )

    def get( self, item ):
        if self.__job != None:
            return_item = self.__job.get( item )
            if None == return_item:
                self.__log.error( "genomon_job.get: specified item \"{item}\" has not been found in yaml.".format(
                                    item = item ) )
            return return_item
        else:
            self.__log.error( "genomon_job.get: job file is not loaded properly." )
            return None
