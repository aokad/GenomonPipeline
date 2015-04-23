#  Copyright Human Genome Center, Institute of Medical Science, the University of Tokyo
#  @since 2015

from __main__ import *
from resource.genomon_rc import job_config_default as default_values

"""
    Genomon job status management object

"""
class genomon_status:
    #
    # Interface
    #
    def __init__( self, status_file, datetime ):

        if status_file == None or datetime == None:
            raise

        self.__file_name = status_file
        self.__datetime = status_file

    def save_status( self, function, input, return_code ):
        input_tmp = os.path.basename( input ).split( '.' )[0]
        file_obj = open( "{fileprefix}_{function}_{inputfile}_r{returncode}".format(
                                        fileprefix = self.__file_name,
                                        function = function,
                                        inputfile = input_tmp,
                                        returncode = return_code ),
                        'w' )
        file_obj.write( str( return_code ) )
        file_obj.close()

    def task_failed( self ):
        for file_name in glob( self.__file_name + '*' ):
            if file_name[:-1] != '0':
                return file_name
        return '0'

