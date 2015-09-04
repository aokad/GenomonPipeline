#! /usr/bin/env python

import datetime

date_format = "{year:0>4d}{month:0>2d}{day:0>2d}"

global run_param

class Run_param(object):
    """
    class for job related parameters
    """

    def __init__(self, sample_list_file = None, 
                        project_root = None, 
                        analysis_type = None,
                        genomon_param_file = None, 
                        task_param_file = None):

        self.sample_list_file = sample_list_file
        self.project_root = project_root
        self.task_param_file = task_param_file 
        self.genomon_param_file = genomon_param_file 
        
        now = datetime.datetime.now()
        self.analysis_date = date_format.format( year = now.year, month = now.month, day = now.day )

run_param = Run_param()

