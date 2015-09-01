#! /usr/bin/env python

date_format = {year:0>4d}{month:0>2d}{day:0>2d}

class genomon_job:
    """
    class for job related parameters
    """
    
    def __init__( self,
                  project_root = None,
                  sample_group = None,
                  job_file = None,
                  param_file = None,
                  config_file = None ):

    self.project_root = project_root
    self.analysis_date = date_format.format( year = now.year, month = now.month, day = now.day )
    self.sample_group = "analysis" 

    self.job_file = job_file
    self.param_file = param_file
    self.config_file = config_file
 
