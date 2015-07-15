import os
from glob import glob

def get_dir ( dir_tree, cwd, dir_name ):
    """
    return the path to the specified directory by dir_tree

    """
    if isinstance( dir_tree, dict ):
        for dir_tmp in dir_tree.keys():
            cwd_tmp = cwd + '/' + dir_tmp
            if isinstance( dir_tmp, str) and dir_tmp  == dir_name:
                return cwd_tmp
            if ( isinstance( dir_tree[ dir_tmp ], dict ) or
                 isinstance( dir_tree[ dir_tmp ], list ) ):
                dir_returned =  get_dir( dir_tree[ dir_tmp ], cwd_tmp, dir_name )

                if None != dir_returned:
                    return dir_returned
            
    elif isinstance( dir_tree, list ):
        n = 0
        for dir_tmp in dir_tree:
            if isinstance( dir_tmp, str):
                cwd_tmp = cwd + '/' + dir_tmp
            elif isinstance( dir_tmp, dict):
                cwd_tmp = cwd + '/' + dir_tmp.keys()[ 0 ]

            if ( ( isinstance( dir_tmp, str) and dir_tmp == dir_name ) or
                 ( isinstance( dir_tmp, dict) and dir_tmp.keys()[0] == dir_name ) ):
                return cwd_tmp
            else:
                if ( isinstance( dir_tree[ n ], dict ) or
                     isinstance( dir_tree[ n ], list ) ):
                    dir_returned =  get_dir( dir_tree[ n ], cwd, dir_name )
                    if None != dir_returned:
                        return dir_returned
            n = n + 1
    else:
        if isinstance( dir_tmp, str) and dir_tmp  == dir_name:
            cwd_tmp = cwd + '/' + dir_tmp
            return cwd_tmp

    return None

#
# Job file checker
#
def Job_file_check( job_yaml, keywords ):
    """
        Job file checker

    """
    #
    # Mandatory difinitions
    for keyword in keywords[ 'mandatory' ]:
        if not ( keyword in job_yaml ):
            print( "Keyword '{keyword}' is not defined.".format(
                        keyword = keyword ) )
            return False

    #
    # Dependency
    #
    # 1) input_file_dir exists.
    # 2) project_root exists.
    # 3) input_file exists.
    # 4) pair_id check: pair_id is defined, if {pair_id} exists in qsub_cmd.
    # 5) file_ext check: check if the file_name contains file_ext.
    # 6) sample_subdir exits in input_file_dir.
    #

    #
    # 1), 2)
    data_dir_list = ( 'input_file_dir', 'project_root' )
    for dir_tmp in data_dir_list:
        if not glob( os.path.expanduser( job_yaml[ dir_tmp ] ) ):
            print( "The directory specified in '{dir}:' does not exist.".format( dir = dir_tmp ) )
            return False

    #
    # 3), 4), 5), 6), 7), 8)
    pair_id_list = None
    if 'pair_id' in job_yaml:
        pair_id_list = job_yaml[ 'pair_id' ]
    elif job_yaml[ 'input_file_type' ] == 'pair_fastq':
        print( "'pair_id' is not defined, though 'input_file_type' equals 'pair_fastq'."  )
        return False

    if isinstance( job_yaml[ 'file_name' ], list ):
        file_name_list = job_yaml[ 'file_name' ]
    else:
        file_name_list = [ job_yaml[ 'file_name' ] ]

    for file_name_str in file_name_list: 
        if pair_id_list != None and -1 == file_name_str.find( 'pair_id' ):
            print( "'file_name:' does not contain '{pair_id}'." )
            return False

        if 'file_ext' in job_yaml and -1 == file_name_str.find( job_yaml[ 'file_ext' ] ):
            print( "'file_name:' does not contain '{file_ext}'.".format( file_ext = job_yaml[ 'file_ext' ] ) )
            return False

        dir_list = []
        if 'input_file_dir' in job_yaml:
            if 'sample_subdir' in job_yaml:
                if pair_id_list != None:
                    for pair_id in pair_id_list:
                        dir_list.append( job_yaml[ 'input_file_dir' ] + '/' +
                                         job_yaml[ 'sample_subdir' ] + '/' +
                                         file_name_str.format( pair_id = pair_id ) )
                else:
                    dir_list.append( job_yaml[ 'input_file_dir' ] + '/' +
                                     job_yaml[ 'sample_subdir' ] + '/' +
                                     file_name_str )

            else:
                if pair_id_list != None:
                    for pair_id in pair_id_list:
                        dir_list.append( job_yaml[ 'input_file_dir' ] + '/' +
                                         file_name_str.format( pair_id = pair_id ) )
                else:
                    dir_list.append( job_yaml[ 'input_file_dir' ] + '/' +
                                     file_name_str )

        for dir_tmp in dir_list:
            if not glob( dir_tmp ):
                print( "Specified input file ( {file_name} ) do not exist.".format( file_name = dir_tmp ) )
                return False

    #
    # project_dir_tree output directory check
    #
    print( "" )
    if 'project_dir_tree' in job_yaml:
        dir_tree = job_yaml[ 'project_dir_tree' ]
        task_ids = job_yaml[ 'tasks' ]
        for task in task_ids:
            for nec_dir in keywords[ 'out_dir' ].keys():
                if ( nec_dir in task_ids[ task ] and
                     not get_dir( dir_tree, '', keywords[ 'out_dir' ][ nec_dir ] ) ):
                    print( "{dir} is not found in 'project_dir_tree:'.".format( dir = nec_dir ) )
                    return False

    return True

#
# Analysis parameter file check
#
def Param_file_check( job_yaml, param_yaml, keyword_file ):
    """
    Analsysi parameter file checker

    """

    #
    # Parameters
    #
    for task in job_yaml[ 'tasks' ].keys():
        for task_process in job_yaml[ 'tasks' ][ task ]:
            if task_process in keyword_file[ 'parameters' ].keys():
                for key_process in keyword_file[ 'parameters' ][ task_process ]:
                    if not( key_process in param_yaml[ task_process ].keys() ):
                        print( "Keyword '{keyword}' is not defined for {step}.".format(
                                keyword = word,
                                step = keyword[ word ] ) )
                        return False

    return True

#
# System configuration file check
#
def System_config_file_check( system_config, keyword_file ):

    system_data = keyword_file[ 'system' ]
    for type in system_data.keys():
        for data in system_data[ type ].keys():
            if not system_config.get( type, data ):
                print( "Keyword '{type}:{data}' is not defined in system configuration file.".format(
                            type = type, data = data ) )
                return False

    return True
