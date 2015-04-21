#!/usr/bin/python
"""

Genome Analysis Pipeline
resource file


"""
#
# General
#
date_format = "{year:0>4d}{month:0>2d}{day:0>2d}"
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
        analysis_date:
            sample_date_sample_name:
                - config
                - script:
                - log
                - out:
                    - fastq
                    - bam
                    - annotation
                    - mutation
                    - cnv
                    - fusion
                    - sv
"""

end_dir_list = ( 'config', 'script', 'log', 'fastq', 'bam', 'annotation', 'mutation', 'cnv', 'fusion', 'sv' )
subdir_list = ( 'fastq', 'bam', 'annotation', 'mutation', 'cnv', 'fusion', 'sv' )
data_ext_list = { 'fastq':      'fastq',
                  'bam':        'bam',
                  'annotation': 'txt'
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
                 'python/fisher.py'
        )

#
# Job configuration file default values
#
job_config_default ={
    'max_indel': 2,
    'max_distance': 5,
    'map_quality': 30,
    'base_quality': 15,
    'mismatch_rate': 0.07,
    'min_score': 20,
    'min_depth': 9,
    'rg_id': 'Unknown',
    'sample_desc': 'Unknown',
    'platform': 'Unknown',
    'platform_unit': 'Unknown',
    'library': 'Unknown',
    'seq_center': 'Unknown',
    'pred_med_insert': '500',
    'fastq_filter': False,
    'split_fastq_line_number': 160000000,
    'use_biobambam': False
}

#
# Environment variables
#
env_list = { 'libmaus_PATH': [ 'LD_LIBRARY_PATH' ] }

