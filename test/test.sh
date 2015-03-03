#!/bin/bash
#
# Set SGE
#
#$ -S /bin/bash # set shell in UGE
#$ -cwd         # execute at the submitted dir

pwd             # print current working directory
hostname        # print hostname
date            # print date
echo arg1=$1    # print 1st argument of shell script

python $1 genomon.py --config_file /home/[user]/GenomonProj/script/test/genomon.cfg \
                     --job_file /home/[user]/GenomonProj/script/test/genomon.job \
                     --jobs 6 \
                     --verbose 10
