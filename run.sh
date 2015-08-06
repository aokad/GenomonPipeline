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

export PYTHONHOME=/usr/local/package/python2.7/current
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PYTHONHOME}/lib:/home/w3varann/.local/lib
export R_LIBS=/home/w3varann/.R
export R_HOME=/home/w3varann/.R


/usr/local/package/python2.7/current/bin/python $3 /home/w3varann/tools/Genomon/genomon.py \
           --config_file /home/w3varann/tools/Genomon/samples/genomon.cfg \
           --job_file $1 \
           --param_file $2 \
           --verbose 1 \
           --drmaa

