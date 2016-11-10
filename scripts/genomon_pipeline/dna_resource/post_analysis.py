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

{genomon_pa} {mode} {output_dir} {genomon_root} {sample_sheet} \
--config_file {config_file} \
--samtools {samtools} --bedtools {bedtools} \
--input_file_case1 "{input_file_case1}" \
--input_file_case2 "{input_file_case2}" \
--input_file_case3 "{input_file_case3}" \
--input_file_case4 "{input_file_case4}" || exit $?
"""

    def __init__(self, qsub_option, script_dir):
        super(Res_PostAnalysis, self).__init__(qsub_option, script_dir)

    def output_files(self, mode, samples, genomon_root, sample_conf_name, genomon_conf):
        
        import os
        import ConfigParser
        
        pa_conf = ConfigParser.RawConfigParser()
        pa_conf.read(genomon_conf.get("post_analysis", "config_file"))

        analysis_dir = genomon_root + "/" + mode
        pa_dir = genomon_root + "/post_analysis/" + sample_conf_name
        
        if mode == "qc":
            di_outputs = { "outputs": [], "run_analysis":False, "run_pa":False}
        
            if len(samples) == 0:
                return di_outputs

            for sample in samples:
                if not os.path.exists(analysis_dir + "/" + sample + "/" + sample + pa_conf.get("result_format_qc", "suffix")):
                    di_outputs["run_analysis"] = True
                    break
            
            if di_outputs["run_analysis"] == True: di_outputs["run_pa"] = True
            
            output = pa_dir + '/' + pa_conf.get("merge_format_qc", "output_all")
            if not os.path.exists(output): di_outputs["run_pa"] = True
            
            di_outputs["outputs"].append(output)
            
            return di_outputs
            
        section = ""
        section_in = ""
        
        if mode == "mutation":
            section = "merge_format_mutation"
            section_in = "result_format_mutation"
        elif mode == "sv":
            section = "merge_format_sv"
            section_in = "result_format_sv"
        else:
            return {}
            
        di_outputs = { \
            "case1":{"output_filt":"", "output_unfilt":"", "run_analysis":False, "run_pa":False, "samples":[]}, \
            "case2":{"output_filt":"", "output_unfilt":"", "run_analysis":False, "run_pa":False, "samples":[]}, \
            "case3":{"output_filt":"", "output_unfilt":"", "run_analysis":False, "run_pa":False, "samples":[]}, \
            "case4":{"output_filt":"", "output_unfilt":"", "run_analysis":False, "run_pa":False, "samples":[]}, \
            "all":  {"output_filt":"", "output_unfilt":"", "run_analysis":False, "run_pa":False}, }
        li_outputs = []
        
        if len(samples) == 0:
            di_outputs["outputs"] = []
            return di_outputs
        
        include_unpair = pa_conf.getboolean(section, "include_unpair")
        include_unpanel = pa_conf.getboolean(section, "include_unpanel")
        for complist in samples:
            type = ""
            if (complist[1] != None and complist[2] != None): type = "case1"
            if (complist[1] != None and complist[2] == None and include_unpanel == True): type = "case2"
            if (complist[1] == None and complist[2] != None and include_unpair == True): type = "case3"
            if (complist[1] == None and complist[2] == None and include_unpair == True and include_unpanel == True): type = "case4"
            
            if type == "": continue
            
            di_outputs[type]["samples"].append(complist[0])
            if not os.path.exists(analysis_dir + "/" + complist[0] + "/" + complist[0] + pa_conf.get(section_in, "suffix")):
                di_outputs[type]["run_analysis"] = True
                di_outputs["all"]["run_analysis"] = True
        
        # each case
        for key in di_outputs:
            if key == "all": continue
            if pa_conf.getboolean(section, "separate") == False: continue
            if len(di_outputs[key]["samples"]) == 0: continue
            
            if di_outputs[key]["run_analysis"] == True: di_outputs[key]["run_pa"] = True
            
            output_filt = pa_dir + '/' + pa_conf.get(section, "output_filt_" + key)
            if not os.path.exists(output_filt): di_outputs[key]["run_pa"] = True
            
            output_unfilt = ""
            if pa_conf.getboolean(section, "include_unfilt"):
                output_unfilt = pa_dir + '/' + pa_conf.get(section, "output_" + key)
                if not os.path.exists(output_unfilt): di_outputs[key]["run_pa"] = True
            
            di_outputs[key]["output_filt"] = output_filt
            di_outputs[key]["output_unfilt"] = output_unfilt
            
            if di_outputs[key]["run_pa"] == True:
                li_outputs.append(output_filt)
                li_outputs.append(output_unfilt)
                di_outputs["all"]["run_pa"] = True
       
        # all
        if pa_conf.getboolean(section, "all_in_one") == True:
            
            if di_outputs["all"]["run_analysis"] == True: di_outputs["all"]["run_pa"] = True
            
            output_filt = pa_dir + '/' + pa_conf.get(section, "output_filt_all")
            if not os.path.exists(output_filt): di_outputs["all"]["run_pa"] = True
            
            output_unfilt = ""
            if pa_conf.getboolean(section, "include_unfilt"):
                output_unfilt = pa_dir + '/' + pa_conf.get(section, "output_all")
                if not os.path.exists(output_unfilt): di_outputs["all"]["run_pa"] = True
            
            di_outputs["all"]["output_filt"] = output_filt
            di_outputs["all"]["output_unfilt"] = output_unfilt
            
            if di_outputs["all"]["run_pa"] == True:
                li_outputs.append(output_filt)
                li_outputs.append(output_unfilt)
       
        di_outputs["outputs"] = []
        di_outputs["outputs"].extend(li_outputs)
        
        return di_outputs
