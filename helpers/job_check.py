import os
from glob import glob

def get_dir ( dir_tree, cwd, dir_name ):
    """
    return the path to the specified directory by dir_tree

    """
    if isinstance( dir_tree, dict ):
        for dir_tmp in dir_tree.keys():
            cwd_tmp = cwd + '/' + dir_tmp

            if dir_tmp == dir_name:
                return cwd_tmp

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
        if isinstance( dir_tree, str) and dir_tree  == dir_name:
            cwd_tmp = cwd + '/' + dir_tree
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
    # 1) input_file exists.
    # 2) pair_id check: pair_id is defined, if {pair_id} exists in qsub_cmd.
    # 3) file_ext check: check if the file_name contains file_ext.
    # 4) sample_name exits in input_file_dir.
    # 5) cmd_options need to be defined for specified process in 'tasks'.
    #
    # All removed, because of a specification change.

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

        data_dir_list = get_dir( dir_tree, '', 'data' ).split( '/' )
        dir_tree_tmp = dir_tree
        for id in range( 1, len( data_dir_list ) ):
            dir_tree_tmp = dir_tree_tmp[ data_dir_list[ id ] ]

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

