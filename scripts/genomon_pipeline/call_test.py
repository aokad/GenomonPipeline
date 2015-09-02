#! /usr/bin/evn python

from genomon_pipeline.run_param import *
from genomon_pipeline.task_param import *
 
def call_test():

    # print run_param.sample_list_file
    print run_param.analysis_date

    print task_param.get("fisher", "man_AF_normal")

