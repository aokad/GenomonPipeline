#!/bin/bash
#
# Set SGE
#
#$ -S /bin/bash # set shell in UGE
#$ -o /home/eigos/Log
#$ -e /home/eigos/Log
#$ -cwd         # execute at the submitted dir

pwd             # print current working directory
hostname        # print hostname
date            # print date
echo arg1=$1    # print 1st argument of shell script

python $1 genomon.py --config_file /home/eigos/Data/GenomonProj/genomon_test/genomon.cfg \
                     --job_file /home/eigos/Data/GenomonProj/genomon_test/genomon.job \
                     --jobs 6 \
                     --verbose 10
