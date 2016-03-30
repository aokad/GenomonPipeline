#! /usr/bin/env python

import datetime

date_format = "{year:0>4d}{month:0>2d}{day:0>2d}"
timestamp_format = "{year:0>4d}{month:0>2d}{day:0>2d}_{hour:0>2d}{min:0>2d}_{msecond:0>2d}"

global run_conf

class Run_conf(object):
    """
    class for job related parameters
    """

    def __init__(self, sample_conf_file = None, 
                        project_root = None, 
                        analysis_type = None,
                        genomon_conf_file = None):

        self.sample_conf_file = sample_conf_file
        self.project_root = project_root
        self.genomon_conf_file = genomon_conf_file 
        
        now = datetime.datetime.now()
        self.analysis_date = date_format.format( year = now.year, month = now.month, day = now.day )
        self.analysis_timestamp = timestamp_format.format( year = now.year, month = now.month, day = now.day, hour=now.hour, min=now.minute, msecond=now.microsecond)

run_conf = Run_conf()

