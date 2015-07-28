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
            - all_summary
            - config
            - script:
            - log
            - out:
                - fastq
                - bam
                - summary
                - analysis_date:
                    - annotation
                    - mutation
                    - cufflinks
                    - cuffdiff
                    - cummeRbund
                    - cnv
                    - sv
                    - itd
                    - star
                    - star_fusion
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
                'cnv',
                'fusion',
                'sv',
                'itd',
                'cufflinks',
                'cuffdiff',
                'cummeRbund',
                'star',
                'star_fusion',
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
                  'star_fusion':    'star_fusion'
}
                
data_ext_list = { 'fastq':          'fastq',
                  'bam':            'bam',
                  'itd':            'txt',
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
                 'python/ToVCF.py',
                 'python/ANNOVAR.py',
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
    'markduplicates':
        { 'java_memory': '6G' },
    'fisher_mutation_call':
        { 'max_indel': 2,
          'max_distance': 5,
          'map_quality': 30,
          'base_quality': 15,
          'mismatch_rate': 0.07,
          'min_depth': 9 },
    'star_genome': 
        { 'additional_params': '--sjdbOverhang 99' },
    'star': 
        { 'additional_params': '--outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
                                --outReadsUnmapped None --chimSegmentMin 15 --chimJunctionOverhangMin 15 \
                                --alignMatesGapMax 200000 --alignIntronMax 200000 --outSAMtype BAM SortedByCoordinate' },
    'star_fusion': 
        { 'additional_params': ' ' },
    'tophat2':
        { 'additional_params': '-p 8' },
    'cuffdiff':
        { 'additional_params': '-p 4' },
    'cufflinks':
        { 'additional_params': ' ' },
    'itd_detection':
        { 'create_ctrl_panel': False,
          'itd_ctrl_panel_files': '',
          'ctrl_panel_normal': '',
          'additional_params': '' },
    'annotattion':
        { 'additional_params': '--buildver hg19 --verdbsnp 131' },
    'use_biobambam': False,
    'bam_read_group': '@RG\\tID:Unknown\\tSM:Unknown\\tLB:Unknown\\tPL:Unknown\\tPU:Unknown\\tCN:unknown',
}

#
# Environment variables
#
env_list = {
    'libmaus_PATH': [ 'LD_LIBRARY_PATH' ],
    'drmaa_PATH':   [ 'DRMAA_LIBRARY_PATH' ]
}

