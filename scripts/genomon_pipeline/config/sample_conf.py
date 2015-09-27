#! /usr/bin/env python

import sys
import os
from genomon_pipeline.config.run_conf import *

class Sample_conf(object):

    def __init__(self):

        self.fastq = {}
        self.bam_tofastq = {}
        self.bam_import = {}
        self.compare = []
        self.control_panel = {}
        # 
        # should add the file exist check here ?
        #
    

    def parse_file(self, file_path):

        file_ext = os.path.splitext(file_path)[1]

        file_data = []
        if file_ext.lower() == '.csv':
            file_data = self.parse_csv(file_path)
        elif file_ext.lower() == '.txt' or file_ext.lower() == '.tsv':
            file_data = self.parse_tsv(file_path)
        # elif file_ext.lower() == '.xlsx':
            # file_data = self.parse_xlsx(file_path)
        else:
            # 
            # should treat other cases ??
            #
            raise NotImplementedError("currently, we can just accept tsv and csv formats")
 

        file_data_trimmed = []
        for line_data in file_data:
       
            # skip empty lines
            if len(line_data) == 0: continue
 
            # line starting with '#' is comment
            if line_data[0].startswith('#'): continue
             
            # remove spaces
            line_data = map(lambda x: x.strip(' '), line_data)

            # skip if all the elements are empty
            if len(line_data) == line_data.count(''): continue

            file_data_trimmed.append(line_data)


        self.parse_data(file_data_trimmed)


    def parse_csv(self, file_path):

        _file_data = []
        import csv
        with open(file_path, 'r') as hIN:
            csv_obj = csv.reader(hIN)
            for cells in csv_obj:
                tempdata = []
                for cell in cells:
                    tempdata.append(cell)
                _file_data.append(tempdata)
    
        return _file_data


    def parse_tsv(self, file_path):

        _file_data = []
        with open(file_path, 'r') as hIN:
            for line in hIN:
                F = line.rstrip('\n').split('\t')
                _file_data.append(F)

        return _file_data


    def parse_data(self, _data ):
    
        mode = ''
        
        sampleID_list = []
        for row in _data:
            if row[0].startswith('['):

                # header
                if row[0].lower() == '[fastq]':
                    mode = 'fastq'
                    continue
                elif row[0].lower() == '[bam_tofastq]':
                    mode = 'bam_tofastq'
                    continue
                elif row[0].lower() == '[bam_import]':
                    mode = 'bam_import'
                    continue
                elif row[0].lower() == '[compare]':
                    mode = 'compare'
                    continue
                elif row[0].lower() == '[controlpanel]':
                    mode = 'controlpanel'
                    continue
                else:
                    err_msg = "Section name should be either of [fastq], [bam_tofastq], [bam_import], " + \
                              "[compare] or [controlpanel]. " + \
                              "Also, sample name should not start with '['."
                    raise ValueError(err_msg)
            
            
            # section data
            if mode == 'fastq':

                sampleID = row[0]
                # 'None' is presereved for special string
                if sampleID == 'None':
                    err_msg = "None can not be used as sampleID"
                    raise ValueError(err_msg)

                if sampleID in sampleID_list:
                    err_msg = sampleID + " is duplicated."
                    raise ValueError(err_msg)

                sampleID_list.append(sampleID)

                if len(row) not in [2, 3]:
                    err_msg = sampleID + ": the path for read1 (and read2) should be provided"
                    raise ValueError(err_msg)

                sequence1 = row[1].split(';')
                sequence2 = row[2].split(';')

                for seq in sequence1 + sequence2:
                    if not os.path.exists(seq):
                        err_msg = sampleID + ": " + seq +  " does not exists" 
                        raise ValueError(err_msg)

                self.fastq[sampleID] = [sequence1, sequence2]

            elif mode == 'bam_tofastq':

                sampleID = row[0]
                # 'None' is presereved for special string
                if sampleID == 'None':
                    err_msg = "None can not be used as sampleID"
                    raise ValueError(err_msg)

                if sampleID in sampleID_list:
                    err_msg = sampleID + " is duplicated."
                    raise ValueError(err_msg)

                sampleID_list.append(sampleID)

                if len(row) not in [2, 3]:
                    err_msg = sampleID + ": only one bam file is allowed"
                    raise ValueError(err_msg)

                sequence = row[1]
                if not os.path.exists(sequence):
                    err_msg = sampleID + ": " + sequence +  " does not exists"
                    raise ValueError(err_msg)

                self.bam_tofastq[sampleID] = sequence
                
            elif mode == 'bam_import':

                sampleID = row[0]
                # 'None' is presereved for special string
                if sampleID == 'None':
                    err_msg = "None can not be used as sampleID"
                    raise ValueError(err_msg)

                if sampleID in sampleID_list:
                    err_msg = sampleID + " is duplicated."
                    raise ValueError(err_msg)

                sampleID_list.append(sampleID)

                if len(row) not in [2, 3]:
                    err_msg = sampleID + ": only one bam file is allowed"
                    raise ValueError(err_msg)

                sequence = row[1]
                if not os.path.exists(sequence):
                    err_msg = sampleID + ": " + sequence +  " does not exists"
                    raise ValueError(err_msg)
                
                sequence_prefix, ext = os.path.splitext(sequence)
                if (not os.path.exists(sequence + '.bai')) or (not os.path.exists(sequence_prefix + '.bai')):
                    err_msg = sampleID + ": " + sequence +  " index does not exists"
                    raise ValueError(err_msg)

                self.bam_import[sampleID] = sequence


            elif mode == 'compare':

                tumorID = row[0]
                if tumorID not in sampleID_list:
                    err_msg = "[compare] section, " + tumorID + " is not defined"
                    raise ValueError(err_msg)

                normalID = row[1] if len(row) >= 2 and row[1] not in ['', 'None'] else None
                controlpanelID = row[2] if len(row) >= 3 and row[2] not in ['', 'None'] else None

                if normalID is not None and normalID not in sampleID_list:
                    err_msg = "[compare] section, " + normalID + " is not defined"
                    raise ValueError(err_msg)

                self.compare.append((tumorID, normalID, controlpanelID))


            elif mode == 'controlpanel':

                if len(row) <= 1:
                    err_msg = "[controlpanel] section, list item is none for the row: " + ','.join(row)
                    raise ValueError(err_msg)

                controlpanelID = row[0]

                for sample in row[1:]:
                    if sample not in sampleID_list:
                        err_msg = "[controlpanel] section, " + sample + " is not defined in " + \
                                    "controlpanelID: " + controlpanelID
                        raise ValueError(err_msg)
 
                self.control_panel[controlpanelID] = row[1:]



        # check whether controlpanleID in compare section is defined
        # for comp in self.compare:
        #     if comp[2] is not None and comp[2] not in self.controlpanel:
        #         err_msg = "[compare] section, controlpanelID: " + comp[2] + " is not defined"
        #         raiseValueError(err_msg)


    # get the paths where bam files imported from another projects will be located"
    def get_linked_bam_import_path(self):
        linked_bam_import_path = []
        for sample in sample_conf.bam_import:
            linked_bam_import_path.append(run_conf.project_root + '/bam/' + sample + '/' + sample + '.bam')
        return linked_bam_import_path

    def sample2bam(self,sample):
        bam = ''
        if (self.bam_import.has_key(sample)):
            bam = run_conf.project_root + '/bam/' + sample + '/' + sample + '.bam'
        else:
            bam = run_conf.project_root + '/bam/' + sample + '/' + sample + '.markdup.bam'
        return bam 

    def get_control_panel_list(self,panel_name):
        control_panel_bam = []
        for sample in self.control_panel[panel_name]:
            control_panel_bam.append(self.sample2bam(sample))
        return control_panel_bam

    def get_disease_and_control_panel_bam(self):
        unique_bams = []
        for complist in sample_conf.compare:
            panel_name = complist[2]
            unique_bams.extend(self.get_control_panel_list(panel_name))
            disease_sample = complist[1]
            unique_bams.append(run_conf.project_root + '/bam/' + disease_sample + '/' + disease_sample + '.markdup.bam')
        result_list = list(set(unique_bams))       
        return result_list

global sample_conf 
sample_conf = Sample_conf()


