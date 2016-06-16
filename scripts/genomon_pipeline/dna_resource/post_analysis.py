#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Res_PostAnalysis(Stage_task):

    task_name = "post_analysis"

    script_template = """
#!/bin/bash
#
# Set SGE
#
#$ -S /bin/bash         # set shell in UGE
#$ -cwd                 # execute at the submitted dir
pwd                     # print current working directory
hostname                # print hostname
date                    # print date
set -xv

# set python environment
export PYTHONHOME={pythonhome}
export PATH=$PYTHONHOME/bin:$PATH
export PYTHONPATH={pythonpath}

{genomon_pa} run {mode} {output_dir} {genomon_root} {sample_sheet} \
--config_file {config_file} \
--samtools {samtools} --bedtools {bedtools} \
--input_file_case1 "{input_file_case1}" \
--input_file_case2 "{input_file_case2}" \
--input_file_case3 "{input_file_case3}" \
--input_file_case4 "{input_file_case4}"
"""

    def __init__(self, qsub_option, script_dir):
        super(Res_PostAnalysis, self).__init__(qsub_option, script_dir)
        
    def sample_split_case(self, li, unpair=True, unpanel=True):
        tmr_nrml_list = []
        tmr_nrml_none = []
        tmr_none_list = []
        tmr_none_none = []
    
        for item in li:
            if item[1]== None:
                if unpair == True:
                    if item[2] == None:
                        if unpanel == True: tmr_none_none.append(item[0])
                    else:
                        tmr_none_list.append(item[0])
            else:
                if item[2] == None:
                    if unpanel == True: tmr_nrml_none.append(item[0])
                else:
                    tmr_nrml_list.append(item[0])
        
        return {"case1": tmr_nrml_list, "case2": tmr_nrml_none, "case3": tmr_none_list, "case4": tmr_none_none}


    def output_files(self, mode, samples, base_dir, genomon_conf):
        import ConfigParser
        pa_conf = ConfigParser.RawConfigParser()
        pa_conf.read(genomon_conf.get("post_analysis", "config_file"))

        if len(samples) == 0:
            return {}
        
        if mode == "qc":
            return {"all": base_dir + '/' + pa_conf.get("merge_format_qc", "output_all")}
        
        section = ""
        if mode == "mutation":
            section = "merge_format_mutation"
        elif mode == "sv":
            section = "merge_format_sv"
        else:
            return {}
            
        li_output_files = {}
        case1=False
        case2=False
        case3=False
        case4=False
        include_unpair = pa_conf.getboolean(section, "include_unpair")
        include_unpanel = pa_conf.getboolean(section, "include_unpanel")
        for complist in samples:
            if (complist[1] != None and complist[2] != None): case1 = True
            if (complist[1] == None and complist[2] != None and include_unpair == True): case3 = True
            if (complist[1] != None and complist[2] == None and include_unpanel == True): case2 = True
            if (complist[1] == None and complist[2] == None and include_unpair == True and include_unpanel == True): case4 = True
    
        if pa_conf.getboolean(section, "all_in_one"):
            li_output_files["filt_all"] = base_dir + '/' + pa_conf.get(section, "output_filt_all")
            if pa_conf.getboolean(section, "include_unfilt"):
                li_output_files["all"] = base_dir + '/' + pa_conf.get(section, "output_all")
                
        if pa_conf.getboolean(section, "separate"):
            if case1 == True: 
                li_output_files["filt_case1"] = base_dir + '/' + pa_conf.get(section, "output_filt_case1")
                if pa_conf.getboolean(section, "include_unfilt"):
                    li_output_files["case1"] = base_dir + '/' + pa_conf.get(section, "output_case1")
            if case2 == True: 
                li_output_files["filt_case2"] = base_dir + '/' + pa_conf.get(section, "output_filt_case2")
                if pa_conf.getboolean(section, "include_unfilt"):
                    li_output_files["case2"] = base_dir + '/' + pa_conf.get(section, "output_case2")
            if case3 == True: 
                li_output_files["filt_case3"] = base_dir + '/' + pa_conf.get(section, "output_filt_case3")
                if pa_conf.getboolean(section, "include_unfilt"):
                    li_output_files["case3"] = base_dir + '/' + pa_conf.get(section, "output_case3")
            if case4 == True: 
                li_output_files["filt_case4"] = base_dir + '/' + pa_conf.get(section, "output_filt_case4")
                if pa_conf.getboolean(section, "include_unfilt"):
                    li_output_files["case4"] = base_dir + '/' + pa_conf.get(section, "output_case4")
        
        return li_output_files
        