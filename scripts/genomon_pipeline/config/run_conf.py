#! /usr/bin/env python

import datetime

date_format = "{year:0>4d}{month:0>2d}{day:0>2d}"

global run_conf

class Run_conf(object):
    """
    class for job related parameters
    """

    def __init__(self, sample_conf_file = None, 
                        project_root = None, 
                        analysis_type = None,
                        genomon_conf_file = None, 
                        task_conf_file = None):

        self.sample_conf_file = sample_conf_file
        self.project_root = project_root
        self.task_conf_file = task_conf_file 
        self.genomon_conf_file = genomon_conf_file 
        
        now = datetime.datetime.now()
        self.analysis_date = date_format.format( year = now.year, month = now.month, day = now.day )

run_conf = Run_conf()

