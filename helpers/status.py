#  Copyright Human Genome Center, Institute of Medical Science, the University of Tokyo
#  @since 2015

import re
from __main__ import *
from resource.genomon_rc import job_config_default as default_values

"""
    Genomon job status management object

"""
class genomon_status:
    #
    # Interface
    #
    def __init__( self, config_dir, header, timestamp ):

        if config_dir == None or header == None or timestamp == None:
            raise

        self.__config_dir = config_dir
        self.__header = header
        self.__timestamp = timestamp

    def save_status( self, function, input, return_code, use_subdir = False):

        input_tmp = os.path.basename( input ).split( '.' )[ 0 ]
        if use_subdir:
            subdir = os.path.basename( os.path.split( input )[ 0 ] )
            input_tmp = subdir + '_' + input_tmp

        file_obj = open( "{dir}/{header}_{timestamp}_{function}_{inputfile}_r{returncode}".format(
                                        dir = self.__config_dir,
                                        header = self.__header,
                                        timestamp = self.__timestamp,
                                        function = function,
                                        inputfile = input_tmp,
                                        returncode = return_code ),
                        'w' )
        file_obj.write( str( return_code ) )
        file_obj.close()

    def check_exit_status( self, function, input, use_subdir = False ):

        #
        # file name example
        # genomon_20150602_1626_491694_bam_stats_merge_ATL_data_markdup_r0
        #
        input_tmp = os.path.basename( input ).split( '.' )[ 0 ]
        if use_subdir:
            subdir = os.path.basename( os.path.split( input )[ 0 ] )
            input_tmp = subdir + '_' + input_tmp

        file_name = "{dir}/{header}_*_{function}_{inputfile}_r*".format(
                                        dir = self.__config_dir,
                                        header = self.__header,
                                        function = function,
                                        inputfile = input_tmp )
        file_name_reg = "{header}_.+_{function}_{inputfile}_r([0-9]+)$".format(
                                        header = self.__header,
                                        function = function,
                                        inputfile = input_tmp )

        pattern = re.compile( file_name_reg )
        for file in glob( file_name ):
            match = pattern.search( file )
            if match:
                match_list = match.groups()
                if len( match_list ) == 1:
                    return_code = int( match_list[ 0 ] )
            else:
                return_code = 1

            if  return_code == 0:
                break
        else:
            return_code = 1

        return return_code
