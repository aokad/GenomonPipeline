#! /usr/bin/env python

import sys
import os
import ConfigParser
from genomon_pipeline.config.task_conf import *

global genomon_conf

genomon_conf = ConfigParser.SafeConfigParser()

dna_reference_list = ["ref_fasta",
                      "interval_list",
                      "hg19_genome",
                      "gaptxt",
                      "bait_file",
                      "simple_repeat_tabix_db",
                      "HGVD_tabix_db",
                      "HGMD_tabix_db",
                      "inhouse_tumor_tabix_db",
                      "inhouse_normal_tabix_db"
                      ]
           
dna_software_list = ["blat",
                     "bwa",
                     "samtools",
                     "bedtools",
                     "biobambam",
                     "PCAP",
                     "genomon_sv",
                     "mutfilter",
                     "ebfilter",
                     "fisher",
                     "mutanno",
                     "annovar"
                     ]

rna_reference_list = ["star_genome"
                      ]
           
rna_software_list = ["samtools",
                     "tophat2",
                     "STAR",
                     "STAR-Fusion",
                     "fusionfusion"
                     ]

err_msg = 'No target File : \'%s\' for the %s key in the section of %s' 


def dna_genomon_conf_check():
    """
    function for checking the validity of genomon_conf for DNA analysis
    """

    section = "REFERENCE"
    for key in dna_reference_list:

        if key == "inhouse_normal_tabix_db":
            if task_conf.has_option("annotation", "active_inhouse_normal_flag"):
                flag = task_conf.get("annotation", "active_inhouse_normal_flag")
                if flag == "True":
                    value = genomon_conf.get(section, key)
                    if not os.path.exists(value):
                        raise ValueError(err_msg % (value, key, section))
            continue
        
        if key == "inhouse_tumor_tabix_db":
            if task_conf.has_option("annotation", "active_inhouse_tumor_flag"):
                flag = task_conf.get("annotation", "active_inhouse_tumor_flag")
                if flag == "True":
                    value = genomon_conf.get(section, key)
                    if not os.path.exists(value):
                        raise ValueError(err_msg % (value, key, section))
            continue
            
        if key == "HGVD_tabix_db":
            if task_conf.has_option("annotation", "active_HGVD_flag"):
                flag = task_conf.get("annotation", "active_HGVD_flag")
                if flag == "True":
                    value = genomon_conf.get(section, key)
                    if not os.path.exists(value):
                        raise ValueError(err_msg % (value, key, section))
            continue
            
        if key == "HGMD_tabix_db":
            if task_conf.has_option("annotation", "active_HGMD_flag"):
                flag = task_conf.get("annotation", "active_HGMD_flag")
                if flag == "True":
                    value = genomon_conf.get(section, key)
                    if not os.path.exists(value):
                        raise ValueError(err_msg % (value, key, section))
            continue
            
        value = genomon_conf.get(section, key)
        if not os.path.exists(value):
            raise ValueError(err_msg % (value, key, section))

    section = "SOFTWARE"
    for key in dna_software_list:
        
        if key == "annovar":
            if task_conf.has_option("annotation", "active_annovar_flag"):
                flag = task_conf.get("annotation", "active_annovar_flag")
                if flag == True:
                    value = genomon_conf.get(section, key)
                    if not os.path.exists(value):
                        raise ValueError(err_msg % (value, key, section))
            continue
            
        value = genomon_conf.get(section, key)
        if not os.path.exists(value):
            raise ValueError(err_msg % (value, key, section))

    pass

def rna_genomon_conf_check():
    """
    function for checking the validity of genomon_conf for RNA analysis
    """

    section = "REFERENCE"
    for key in rna_reference_list:
        value = genomon_conf.get(section, key)
        if not os.path.exists(value):
            raise ValueError(err_msg % (value, key, section))

    section = "SOFTWARE"
    for key in rna_software_list:
        value = genomon_conf.get(section, key)
        if not os.path.exists(value):
            raise ValueError(err_msg % (value, key, section))

    pass

