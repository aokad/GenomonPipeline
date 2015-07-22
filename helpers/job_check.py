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
    # 7) cmd_options need to be defined for specified process in 'tasks'.
    #

    #
    # 1), 2)
    data_dir_list = ( 'input_file_dir', 'project_root' )
    for dir_tmp in data_dir_list:
        if not glob( os.path.expanduser( job_yaml[ dir_tmp ] ) ):
            print( "The directory specified in '{dir}:' does not exist.".format( dir = dir_tmp ) )
            return False

    #
    # 3), 4), 5), 6)
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
    # 7) cmd_options
    #
    for task in job_yaml[ 'tasks' ].values()[0]:
        if not 'cmd_options' in job_yaml:
            print( "cmd_option is not defined." )
            return False

        elif not ( task in job_yaml[ 'cmd_options' ].keys() ):
            print( "cmd_option for {task} is not defined.".format( task = task ) )
            return False

    #
    # project_dir_tree output directory check
    #
    if 'project_dir_tree' in job_yaml:
        dir_tree = job_yaml[ 'project_dir_tree' ]
        if not dir_tree:
            print( "project_dir_tree is empty." )
            return False

        task_ids = job_yaml[ 'tasks' ]
        for task in task_ids.values()[0]:
            if ( task in keywords[ 'out_dir' ].keys() and
                not get_dir( dir_tree, '', keywords[ 'out_dir' ][ task ] ) ):
                print( "{dir} for {task} is not found in 'project_dir_tree:'.".format(
                            dir = keywords[ 'out_dir' ][ task ], task = task ) )
                return False

    #
    # If bam_read_group contains PI, the value needs to be a number.
    #
    if 'bam_read_group' in job_yaml:
        bam_read_group_str = job_yaml[ 'bam_read_group' ]
        PI_id = bam_read_group_str.find( 'PI:' )
        num_end_id = bam_read_group_str[ PI_id + 3: ].find( "\t" )
        if not bam_read_group_str[ PI_id + 3: num_end_id ].isdigit():
            print( "bam_read_group: 'PI:' needs to have number set." )
            return False


    return True

#
# Analysis parameter file check
#
def Param_file_check( job_yaml, param_yaml, keyword_file ):
    """
    Analysis parameter file checker

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
    """
    System configuration file checker

    """

    try:
        system_data = keyword_file[ 'system' ]
        for section in system_data.keys():
            for item in system_data[ section ].keys():
                data = system_config.get( section, item )
                if not data:
                    print( "Config file: [{section}] {item} is not defined in system configuration file.".format(
                                section = section, item = item ) )
                    return False

                elif system_data[ section ][ item ] == 'FILE' and not os.path.exists( data ):
                    print( "Config file: [{section}] {item} {data} does not exists.".format(
                                section = section, item = item, data = data ) )
                    return False

            if section != 'ENV':
                for item, data in system_config.items( section ):
                    if data != 'True' and data != 'False' and not os.path.exists( data ):
                        if len( glob( data + '*' ) ) == 0 :
                            print( "Config file: [{section}] {item} {data} does not exists.".format(
                                    section = section, item = item, data = data ) )
                            return False

        return True

    except:
        print( "Config file: get {type}:{item} failed.".format( type = type, item = item, data = data ) )
        return False
