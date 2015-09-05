#! /usr/bin/env python

from genomon_pipeline.genomon_param import *
from genomon_pipeline.task_param import *
from genomon_pipeline.run_param import *
from genomon_pipeline.sample_list import *
from genomon_pipeline.rna_pipeline import rna_pipeline_run


def main(args):

    ###
    # set run_param
    run_param.sample_list_file = args.sample_list_file
    run_param.analysis_type = args.analysis_type
    run_param.project_root = args.project_root
    run_param.genomon_param_file = args.genomon_param_file
    run_param.task_param_file = args.task_param_file
    ###

    ###
    # read sample list file
    sample_list.parse_file(run_param.sample_list_file)
    ###

    ###
    # set and check genomon_param config data
    genomon_param.read(run_param.genomon_param_file)
    genomon_param_check()
    ###

    ###
    # set and check task parameter config data    
    task_param.read(run_param.task_param_file)
    task_param_check()
    ###

    if run_param.analysis_type == "dna":
        raise NotImplementedError("DNA pipeline is still in progress")
    elif run_param.analysis_type == "rna":
        rna_pipeline_run()
    else:
        raise NotImplementedError("Just DNA and RNA pipeline is prepared")    
