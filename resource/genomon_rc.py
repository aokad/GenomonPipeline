#!/usr/bin/python
"""

Genome Analysis Pipeline
resource file


"""
#
# General
#
date_format = "{year:0>4d}{month:0>2d}{day:0>2d}"
timestamp_format =  "{year:0>4d}{month:0>2d}{day:0>2d}_{hour:0>2d}{min:0>2d}_{msecond:0>6d}"
file_timestamp_format = "{name}_{year:0>4d}{month:0>2d}{day:0>2d}_{hour:0>2d}{min:0>2d}_{msecond:0>6d}"

#
# Misc
#

#
# Default directory three
#
dir_tree_resource = \
"""
project_directory:
    data:
        sample_date:
            sample_name
    results:
        - all_summary
        - analysis_date:
            sample_date_sample_name:
                - config
                - script:
                - log
                - out:
                    - fastq
                    - bam
                    - annotation
                    - mutation
                    - cufflinks
                    - cuffdiff
                    - cummeRbund
                    - cnv
                    - fusion
                    - sv
                    - summary
"""

end_dir_list = ( 'config',
                 'script',
                 'log',
                 'fastq',
                 'bam',
                 'annotation',
                 'mutation',
                 'cufflinks',
                 'cuffdiff',
                 'cummeRbund',
                 'cnv',
                 'fusion',
                 'sv',
                 'summary',
                 'all_summary'
                 )

subdir_list = ( 'fastq',
                'bam',
                'annotation',
                'mutation',
                'cnv',
                'fusion',
                'sv',
                'cufflinks',
                'cuffdiff',
                'cummeRbund',
                'summary',
                'all_summary' )

data_ext_list = { 'fastq':          'fastq',
                  'bam':            'bam',
                  'annotation':     'txt',
                  'summary':        'txt',
                  'all_summary':    'txt'
}

#
# script files to copy
#
script_files = ( 'shell/utility.sh',
                 'shell/sleep.sh',
                 'shell/interval.sh',
                 'shell/interval_list.sh',
                 'perl/fastqNPadding.pl',
                 'python/bamfilter.py',
                 'python/fisher.py',
                 'python/coverage.py',
                 'python/bedparse.py',
                 'python/genomesize.py',
                 'python/merge_cov.py',
                 'python/mkxls.py',
                 'python/xl2tsv.py',
                 'R/cummeRbund.R',
                 'R/report.R'
        )

#
# Job configuration file default values
#
job_config_default ={
    'split_fastq':
        { 'fastq_filter': False,
          'split_fastq_line_number': 160000000 },
    'cutadapt':
        { 'adaptor': [ 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN', 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN' ] },
    'bwa_mem':
        { 'min_score': 20 },
    'fisher_mutation_call':
        { 'max_indel': 2,
          'max_distance': 5,
          'map_quality': 30,
          'base_quality': 15,
          'mismatch_rate': 0.07,
          'min_depth': 9 },
    'use_biobambam': False,
    'bam_read_group': '@RG\\tID:Unknown\\tSM:Unknown\\tLB:Unknown\\tPL:Unknown\\tPU:Unknown\\tCN:unknown'
}

#
# Environment variables
#
env_list = {
    'libmaus_PATH': [ 'LD_LIBRARY_PATH' ],
    'drmaa_PATH':   [ 'DRMAA_LIBRARY_PATH' ]
}

