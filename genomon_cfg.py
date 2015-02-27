#  Copyright Human Genome Center, Institute of Medical Science, the University of Tokyo
#  @since 2015

import ConfigParser

"""
    Genomon system configuration file parser
    How to use it:
        ge_cfg = genomon_config( 'path to file' )
        ge_cfg.get( 'REFERENCE', 'hg19_fasta' )

"""
class genomon_config( object ):
    def __init__( self, config_file = None, log = None):

        self.__log = log
        if config_file != None:
            self.__config_file = config_file
            self.open_cfg( config_file )

    def open_cfg( self, config_file ):
        self.__config_file = config_file

        try:
            if self.__config_file != None:
                self.__conf = ConfigParser.SafeConfigParser()
                self.__conf.read( self.__config_file )
            else:
                self.__log.error( "genomon_config.open_cfg: configuration file is not loaded properly." )
                raise

        except IOError as (errno, stderror ):
            self.__log.error( "genomon_config.open_cfg: IOError: error number: {num}, std_error: {stderr}".format(
                        num = errno, stderr = stderror ) )
        except:
            self.__log.error( "genomon_config.open_cfg: Unexpected error:", sys.exc_info()[0] )


    def get( self, section, item ):
        if self.__conf != None:
            return_string =  self.__conf.get( section, item )
            if None == return_string:
                self.__log.error( "genomon_config.get: {sec}: {item} is not defined in system config file.".format(
                                sec = section,
                                item = item
                            ))
            return return_string
        else:
            self.__log.error( "genomon_config.get: configuration file is not loaded properly." )
            return None

