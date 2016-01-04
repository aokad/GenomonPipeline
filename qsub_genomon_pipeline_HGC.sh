#$ -S /bin/bash         # set shell in UGE
#$ -cwd                 # execute at the submitted dir

export PYTHONHOME=/usr/local/package/python/current2.7
export PYTHONPATH=~/.local/lib/python2.7/site-packages
export PATH=${PYTHONHOME}/bin:~/.local/bin:$PATH
export LD_LIBRARY_PATH=${PYTHONHOME}/lib:${LD_LIBRARY_PATH}
export DRMAA_LIBRARY_PATH=/geadmin/N1GE/lib/lx-amd64/libdrmaa.so.1.0

TARGET_PATH=$6
echo "${TARGET_PATH}/genomon_pipeline $1 $2 $3 $4 $5"
${TARGET_PATH}/genomon_pipeline $1 $2 $3 $4 $5

