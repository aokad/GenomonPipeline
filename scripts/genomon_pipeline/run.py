#! /usr/bin/env python

from ruffus import *
from genomon_pipeline.config.genomon_conf import *
from genomon_pipeline.config.task_conf import *
from genomon_pipeline.config.run_conf import *
from genomon_pipeline.config.sample_conf import *


def main(args):

    ###
    # set run_conf
    run_conf.sample_conf_file = args.sample_conf_file
    run_conf.analysis_type = args.analysis_type
    run_conf.project_root = os.path.abspath(args.project_root)
    run_conf.genomon_conf_file = args.genomon_conf_file
    run_conf.task_conf_file = args.task_conf_file
    run_conf.drmaa = args.drmaa
    ###

    ###
    # read sample list file
    sample_conf.parse_file(run_conf.sample_conf_file)
    ###

    ###
    # set and check genomon_conf config data
    genomon_conf.read(run_conf.genomon_conf_file)
    genomon_conf_check()
    ###

    ###
    # set and check task parameter config data    
    task_conf.read(run_conf.task_conf_file)
    task_conf_check()
    ###

    if run_conf.analysis_type == "dna":
        import dna_pipeline
    elif run_conf.analysis_type == "rna":
        import rna_pipeline
    else:
        raise NotImplementedError("Just DNA and RNA pipeline is prepared")

    pipeline_run(
                 verbose = 3, 
                 multiprocess = 20,
                 multithread = 20
                )

    
