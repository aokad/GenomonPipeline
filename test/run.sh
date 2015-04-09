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
set -xv
/usr/local/package/python2.7/2.7.2/bin/python $2 /home/w3varann/tools/Genomon/genomon.py \
           --config_file /home/w3varann/tools/Genomon/genomon.cfg \
           --job_file $1 \
           --jobs 8
