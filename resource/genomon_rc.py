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
# Terminal directory list
#
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
                 )
#
# Directories that need to create subdirs.
dir_nontask_list = {
                'config',
                'script',
                'log'
}

dir_task_list = { 'fastq':          [ 'split_fastq' ],
                  'bam':            [ 'bwa_mem', 'merge_bam', 'markduplicates' ],
                  'summary':        [ 'bam_stats' ],

                  'mutation':       [ 'fisher_mutation_call' ],
                  'annotation':     [ 'annotation' ],

                  'itd':            [ 'itd_detection' ],

                  'sv':             [ 'sv_detection'],

                  'star':           [ 'star' ],
                  'star_fusion':    [ 'star_fusion' ],
                  'fusionfusion':   [ 'fusionfusion' ],
                  'cufflinks':      [ 'cufflinks' ],
                  'cuffdiff':       [ 'cuffdiff' ],
                  'cummeRbund':     [ 'cummeRbund' ],
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
default_param   = 'resource/default_param.yaml'
default_config  = 'resource/default_config.cfg'
default_job     = 'resource/default_job.yaml'
job_file_words  = 'resource/job_file_words.yaml'

default_task = """
DNA:
    - split_fastq
    - bwa_mem
    - merge_bam
    - markduplicates
    - bam_stats
    - fisher_mutation_call
    - itd_detection
    - annotation

RNA:
    - star
    - star_fusion
    - fusionfusion
"""
