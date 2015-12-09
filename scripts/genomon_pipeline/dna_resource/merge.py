#! /usr/bin/env python

from genomon_summary.stage_task import *

class Res_Merge(Stage_task):

    task_name = "merge"

    script_template = """
#!/bin/bash
#
# Set SGE
#
#$ -S /bin/bash         # set shell in UGE
#$ -cwd                 # execute at the submitted dir
#$ -e {log}             # log file directory
#$ -o {log}             # log file directory
pwd                     # print current working directory
hostname                # print hostname
date                    # print date
set -xv

# cat {input1} {input2} {input3} {input4} > {output}
"""

    def __init__(self, qsub_option, script_dir):
        super(Res_Merge, self).__init__(qsub_option, script_dir)
        
    def mkxls(self, files, excel_file):
        import xlwt
        import re
        import numpy
        
        # No.1 bamstat
        bamstats_header = []
        bamstats_string = []    
        f=open(files[0])
        string=f.read()
        f.close()
        input_file = string.split("\n")
        for line in input_file:
            line = line.replace("\n", "")
            line_split = line.split("\t")
            
            if len(bamstats_header) == 0:
                bamstats_header += line_split
            else:
                if len(line_split) > 1:
                    bamstats_string.append(line_split)

        num_index = 6
        range_values = range(num_index, len(bamstats_header))
        bamstats_values = numpy.loadtxt(files[0], delimiter='\t', skiprows=1, 
                                        usecols = range_values)
        
        # No.5 coverage.py
        coverage_header = []
        coverage_values = []  
        f=open(files[1])
        string=f.read()
        f.close()
        input_file = string.split("\n")
        for line in input_file:
            line = line.replace("\n", "")
            line_split = line.split("\t")
            if len(coverage_header) == 0:
                coverage_header += line_split
            else:
                if len(line_split) > 1:
                    coverage_values += line_split
    
        #
        # Make Excel file
        #
        wb = xlwt.Workbook()
        ws = wb.add_sheet('data')
    
        # write bamstats data, total read group
        num_readgroup = len(bamstats_string)
        for i in range(0, len(bamstats_header)):
            ws.write(0, i, bamstats_header[i])
    
            if bamstats_header[i] == 'readgroup':
                ws.write(1, i, "*")
    
            elif bamstats_header[i] in [
                                      '#_mapped_bases',
                                      '#_mapped_bases_r1',
                                      '#_mapped_bases_r2',
                                      '#_divergent_bases',
                                      '#_divergent_bases_r1',
                                      '#_divergent_bases_r2',
                                      '#_total_reads',
                                      '#_total_reads_r1',
                                      '#_total_reads_r2',
                                      '#_mapped_reads',
                                      '#_mapped_reads_r1',
                                      '#_mapped_reads_r2',
                                      '#_mapped_reads_properly_paired',
                                      '#_gc_bases_r1',
                                      '#_gc_bases_r2',
                                      '#_duplicate_reads',
                                     ]:
                if len(bamstats_values) == len(range_values):
                    ws.write(1, i, numpy.sum(bamstats_values[i-num_index]))
                else:
                    ws.write(1, i, numpy.sum(bamstats_values[:,i-num_index]))

            elif bamstats_header[i] in [
                                      'mean_insert_size',
                                      'insert_size_sd',
                                      'median_insert_size',
                                      ]:
                if len(bamstats_values) == len(range_values):
                    ws.write(1, i, numpy.average(bamstats_values[i-num_index]))
                else:
                    ws.write(1, i, numpy.average(bamstats_values[:,i-num_index]))

            elif bamstats_header[i] in [
                                      'read_length_r1',
                                      'read_length_r2',
                                     ]:

                ws.write(1, i, int(float(bamstats_string[0][i])))
    
            else:
                ws.write(1, i, bamstats_string[0][i])
    
        # write bamstats data, each read group
        for i in range(0, num_readgroup):
            for j in range(0, len(bamstats_string[i])):
                if bamstats_header[j] in [
                                      'read_length_r1',
                                      'read_length_r2',
                                      '#_mapped_bases',
                                      '#_mapped_bases_r1',
                                      '#_mapped_bases_r2',
                                      '#_divergent_bases',
                                      '#_divergent_bases_r1',
                                      '#_divergent_bases_r2',
                                      '#_total_reads',
                                      '#_total_reads_r1',
                                      '#_total_reads_r2',
                                      '#_mapped_reads',
                                      '#_mapped_reads_r1',
                                      '#_mapped_reads_r2',
                                      '#_mapped_reads_properly_paired',
                                      '#_gc_bases_r1',
                                      '#_gc_bases_r2',
                                      '#_duplicate_reads'
                                      ]:
                    ws.write(i + 2, j, int(bamstats_string[i][j]))
    
                elif bamstats_header[j] in [
                                      'mean_insert_size',
                                      'insert_size_sd',
                                      'median_insert_size',
                                      ]:
                    ws.write(i + 2, j, float(bamstats_string[i][j]))
    
                else:
                    ws.write(i + 2, j, bamstats_string[i][j])
                
        # write coverage data
        x_pos = len(bamstats_header)

        for i in range(0, len(coverage_header)):
            ws.write(0, i + x_pos, coverage_header[i])
            search_ratio = re.search("[0-9]{1,10}x_ratio$", coverage_header[i])
            search_x = re.search("[0-9]{1,10}x$", coverage_header[i])

            if coverage_header[i] in ["non-N_total_depth", "non-N_bases"] \
                or search_x != None:
                ws.write(1, i + x_pos, int(float(coverage_values[i])))
    
            elif coverage_header[i] in ['average_depth', 'depth_stdev'] \
                or search_ratio != None:
                ws.write(1, i + x_pos, float(coverage_values[i]))
    
            else:
                ws.write(1, i + x_pos, coverage_values[i])
        
        # save
        wb.save(excel_file)
   
    
    def Excel2TSV(self, ExcelFile, TSVFile):
         import xlrd
         workbook = xlrd.open_workbook(ExcelFile)
         worksheet = workbook.sheet_by_name('data')
         tsvfile = open(TSVFile, 'wb')
    
         for rownum in xrange(worksheet.nrows):
             data_to_write = []
             for x in worksheet.row_values(rownum):
                 if isinstance(x, basestring):
                     x = x.replace("\t", "")
                 if type(x) == type(u''):
                     data_to_write.append( x.encode('utf-8'))
                 else:
                     data_to_write.append(str(x))
    
             tsvfile.write('\t'.join(data_to_write) + '\n')
    
         tsvfile.close()
    
    
