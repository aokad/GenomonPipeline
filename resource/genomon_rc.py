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
        sample_name
    results:
        sample_name:
            - config
            - script:
            - log
            - out:
                - all_summary
                - fastq
                - bam
                - summary
                - analysis_date:
                    - annotation
                    - mutation
                    - cufflinks
                    - cuffdiff
                    - cummeRbund
                    - sv
                    - itd
                    - star
                    - star_fusion
                    - fusionfusion
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
                 'fusionfusion',
                 'sv',
                 'itd',
                 'star',
                 'star_fusion',
                 'summary',
                 'all_summary'
                 )

subdir_list = ( 'fastq',
                'bam',
                'annotation',
                'mutation',
                'fusion',
                'sv',
                'itd',
                'cufflinks',
                'cuffdiff',
                'cummeRbund',
                'star',
                'star_fusion',
                'fusionfusion',
                'summary',
                'all_summary' )

dir_task_list = { 'fastq':          'split_fastq',
                  'bam':            'bwa_mem',
                  'bam':            'merge_bam',
                  'bam':            'markduplicates',
                  'summary':        'bam_stats',
                  'mutation':       'fisher_mutation_call',
                  'annotation':     'annotation',
                  'itd':            'itd_detection',
                  'star':           'star',
                  'star_fusion':    'star_fusion',
                  'fusionfusion':   'fusionfusion'
}
                
data_ext_list = { 'fastq':          'fastq',
                  'bam':            'bam',
                  'itd':            'txt',
                  'annotation':     'txt',
                  'fusionfusion':   'txt',
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
                 'python/ToVCF.py',
                 'python/ANNOVAR.py',
                 'R/cummeRbund.R',
                 'R/report.R'
        )

#
# Environment variables
#
env_list = {
    'libmaus_PATH': [ 'LD_LIBRARY_PATH' ],
    'drmaa_PATH':   [ 'DRMAA_LIBRARY_PATH' ]
}

#
# Default parameter file name
#
default_file_name = 'resource/default_param.yaml'

