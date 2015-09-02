#  Copyright Human Genome Center, Institute of Medical Science, the University of Tokyo
#  @since 2015
import sys
import os
import yaml
from __main__ import *
from resource.genomon_rc import default_param as default_param
from resource.genomon_rc import default_job as default_job
from resource.genomon_rc import default_task as default_task
from resource.genomon_rc import date_format as date_format
import job_check


"""
    Genomon job configuration file parse
    How to use it:
        ge_job = genomonjob( 'path to the job configuration file' )
        ge_job.get( 'bwa' )

"""

class genomon_job:
    #
    # Interface
    #
    def __init__( self,
                  job_file = None,
                  param_file = None,
                  Genomon_dir = None,
                  project_root = None,
                  compare_list = None,
                  task = None,
                  now = None,
                  log = None ):

        try:
            self.__log = log

            if not Genomon_dir:
                raise
            self.__Genomon_dir = Genomon_dir
            self.__project_root = project_root
            self.__compare_list = compare_list
            self.__task = task
            self.__now = date_format.format( year = now.year, month = now.month, day = now.day )
            self.__sample_group = 'analysis'

            self.open_job( job_file )
            self.open_param( param_file )

        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            self.__log.error( "genomon_job.init: unexpected error:", sys.exc_info()[0] )
            self.__log.error("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )
            raise Exception( 'genomon_job.__init__' )

    def open_param( self, param_file = None ):
        self.__param_file = param_file

        try:
            if self.__param_file != None:
                f = open( self.__param_file )
                self.__param = yaml.load( f )
                f.close()
            else:
                self.__param = None

            f = open( self.__Genomon_dir + '/' + default_param )
            self.__default_param = yaml.load( f )
            f.close()

        except IOError as (errno, stderror ):
            self.__log.error( "genomon_job.open_job: IOError: error number: {num}, std_error: {stderr}".format(
                        num = errno, stderr = stderror ) )
            raise Exception( "genomon_job.open_param" )

        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            self.__log.error( "genomon_job.open_job: unexpected error:", sys.exc_info()[0] )
            self.__log.error("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )
            raise Exception, "open_param"

    def compare_list_to_pairs( self, compare_list ):
        if compare_list:
            return_data = {}
            for normal, tumor, list in compare_list:
                if normal and tumor:
                    return_data[ normal ] = tumor
                elif normal:
                    if 'Normal' in return_data:
                        return_data[ 'Normal' ].append( normal )
                    else:
                        return_data[ 'Normal' ] = [ normal ]

                elif tumor:
                    if 'Disease' in return_data:
                        return_data[ 'Disease' ].append( tumor )
                    else:
                        return_data[ 'Disease' ] = [ tumor ]

            return return_data
                


    def open_job( self, job_file = None ):
        self.__job_file = job_file

        try:
            f = open( self.__Genomon_dir + '/' + default_job )
            job_yaml = self.__default_job = yaml.load( f )
            f.close()

            if self.__job_file:
                f = open( self.__job_file )
                job_yaml = self.__job = yaml.load( f )
                f.close()
            else:
                self.__job = None

            if self.__compare_list:
                job_yaml[ 'control_disease_pairs' ] = self.compare_list_to_pairs( self.__compare_list )

            if self.__project_root:
                job_yaml[ 'project_root' ] = self.__project_root

            if self.__task:
                job_yaml[ 'tasks' ] =  {}
                job_yaml[ 'tasks' ][ self.__task ] = yaml.load( default_task )[ self.__task ]
                f.close()

            if self.__sample_group:
                job_yaml[ 'sample_group' ] = self.__sample_group

            if not 'analysis_data' in job_yaml:
                job_yaml[ 'analysis_data' ] = self.__now

        except IOError as (errno, stderror ):
            self.__log.error( "genomon_job.open_job: IOError: error number: {num}, std_error: {stderr}".format(
                        num = errno, stderr = stderror ) )
            raise Exception( 'genomon_job.open_job' )

        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            self.__log.error( "genomon_job.open_job: unexpected error:", sys.exc_info()[0] )
            self.__log.error("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )
            raise Exception( 'genomon_job.open_job' )

    def get_param( self, task, item ):
        return_value = None

        if self.__param:
            param_dict = self.__param.get( task )
            if param_dict:
                if item in param_dict:
                    return_value = param_dict[ item ]

        if None == return_value:
            if item in self.__default_param[ task ]:
                return_value = self.__default_param[ task ][ item ]
            else:
                return_value = None

        return return_value


    def get_job( self, item ):
        return_item = None

        if self.__job != None:
            return_item = self.__job.get( item )

        if return_item == None:
            return_item = self.__default_job.get( item )

        return return_item

    def check_file( self, keyword_file ):
        try:
            f = open( keyword_file )
            f_yaml = yaml.load( f )
            if self.__job and job_check.Job_file_check( self.__job, f_yaml ):
                return_value = True
            else:
                return_value = False

            if ( return_value and job_check.Param_file_check( self.__job, 
                                                              self.__param if  hasattr(self, '__param') else self.__default_param,
                                                              f_yaml ) ):
                return_value = True
            else:
                return_value = False

            f.close()

            return return_value

        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            self.__log.error( "genomon_job.open_job: unexpected error:", sys.exc_info()[0] )
            self.__log.error("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) )
            raise Exception( 'genomon_job.check_file' )
