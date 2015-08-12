#!/bin/bash
#
# Set SGE
#
#$ -S /bin/bash # set shell in UGE
#$ -cwd         # execute at the submitted dir
pwd             # print current working directory
hostname        # print hostname
date            # print date

PYTHONHOME=/usr/local/package/python2.7/current
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PYTHONHOME}/lib:/home/w3varann/.local/lib
R_LIBS=/home/w3varann/.R
R_HOME=/home/w3varann/.R

GENOMON_DIR=`dirname $0`
/usr/local/package/python2.7/current/bin/python ${GENOMON_DIR}/Genomon --verbose 1 $1 $2 $3 

