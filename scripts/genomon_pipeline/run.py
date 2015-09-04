#! /usr/bin/env python

from genomon_pipeline.genomon_param import *
from genomon_pipeline.task_param import *
from genomon_pipeline.run_param import *
from genomon_pipeline.sample_list import *

import genomon_pipeline.call_test  

def main(args):

    ###
    # set run_param
    run_param.sample_list_file = args.sample_list_file
    run_param.project_root = args.project_root
    run_param.genomon_param_file = args.genomon_param_file
    run_param.task_param_file = args.task_param_file
    ###

    ###
    # read sample list file
    sample_list.parse_file(run_param.sample_list_file)
    ###

    ###
    # just test
    print run_param.sample_list_file    
    print run_param.analysis_date
    # genomon_pipeline.call_test.call_test()
    print sample_list.fastq
    print sample_list.compare
    print sample_list.control_panel
    ###

    ###
    # set and check genomon_param config data
    genomon_param.read("param.cfg")
    genomon_param_check()
    ###

    task_param.read(run_param.task_param_file)
    task_param_check()

    print task_param.get("bwa_mem", "qsub_option")
    genomon_pipeline.call_test.call_test() 

    print "shine"

# print genomon_param.get("bwa", "option")
    
