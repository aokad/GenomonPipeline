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
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/package/python2.7/current/lib
export PYTHONPATH=$PYTHONPATH:/home/w3varann/.local/lib/python2.7/site-packages
/usr/local/package/python2.7/2.7.2/bin/python $2 ./genomon.py \
           --config_file /home/eigos/Data/Genomon/db/genomon.cfg \
           --job_file $1 \
           --jobs 8 \
           --verbose 10 \
           --abpath
