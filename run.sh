#!/bin/bash
#
# Set SGE
#
#$ -S /bin/bash # set shell in UGE
#$ -cwd         # execute at the submitted dir
pwd             # print current working directory
hostname        # print hostname
date            # print date
set -xv
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/package/python2.7/current/lib
export PYTHONPATH=$PYTHONPATH:/home/w3varann/.local/lib/python2.7/site-packages
/usr/local/package/python2.7/2.7.2/bin/python $3 /home/w3varann/tools/Genomon/genomon.py \
           --config_file /home/w3varann/tools/Genomon/samples/genomon.cfg \
           --job_file $1 \
           --param_file $2 \
           --jobs 10 \
           --verbose 10
