#! /usr/bin/env python

import sys
import os

class Sample_list(object):

    def __init__(self, file_path):

        self.fastq = {}
        self.compare = []
        self.control_panel = {}

        # 
        # should add the file exist check here ?
        #

        self.parse_file(file_path)



    def parse_file(self, file_path):

        file_ext = os.pat.splitext(file_path)[1]]

        file_data = []
        if file_ext.lower() == '.csv':
            file_data = self.parse_csv(file_path)
        elif file_ext.lower() == '.txt' or file_ext.lower() == '.tsv':
            file_data = self.parse_tsv(file_path)
        elif file_ext.lower() == '.xlsx':
            file_data = self.parse_xlsx(file_path)
            # self.parse_data(self.parse_xlsx(filepath))
        else:
            # 
            # should treat other cases ??
            #
            raise NotImplementedError("Currently, we can just accept tsv, csv and xlsx formats")
 

        file_data_trimmed = []
        for line_data in file_data:
        
            # line starting with '#' is comment
            if line_data[0].startswith('#'): continue
             
            # remove spaces
            line_data = map(lambda x: x.strip(' '), line_data)

            # skip if all the elements are empty
            if len(line_data) == line_data.count(''): continue

            file_data_trimmed.append(line_data)


    def parse_csv(self, file_path):

        _file_data = []
        import csv
        with open(file_path, 'r') as hIN:
            csv_obj = csv.reader(hIN):
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
              _file_data.append(F)

        return _file_data


