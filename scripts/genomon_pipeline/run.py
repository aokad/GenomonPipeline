#! /usr/bin/env python

from ruffus import *
from genomon_pipeline.config.genomon_conf import *
from genomon_pipeline.config.run_conf import *
from genomon_pipeline.config.sample_conf import *


def main(args):

    ###
    # set run_conf
    run_conf.sample_conf_file = args.sample_conf_file
    run_conf.analysis_type = args.analysis_type
    run_conf.project_root = os.path.abspath(args.project_root)
    run_conf.genomon_conf_file = args.genomon_conf_file
    run_conf.drmaa = False if args.disable_drmaa else True

    ###
    # read sample list file
    sample_conf.parse_file(run_conf.sample_conf_file)

    ###
    # set genomon_conf and task parameter config data
    genomon_conf.read(run_conf.genomon_conf_file)
    
    if run_conf.analysis_type == "dna":
        dna_genomon_conf_check()
        dna_software_version_set()
        import dna_pipeline
    elif run_conf.analysis_type == "rna":
        rna_genomon_conf_check()
        import rna_pipeline
    else:
        raise NotImplementedError("Just DNA and RNA pipeline is prepared")

    if not (args.param_check):
        pipeline_run(
                     verbose = args.verbose, 
                     multiprocess = args.multiprocess
                    )

        
