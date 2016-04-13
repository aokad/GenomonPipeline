# /bin/bash 

export PYTHONHOME=/usr/local/package/python/current2.7
export PYTHONPATH=~/.local/lib/python2.7/site-packages
export PATH=${PYTHONHOME}/bin:~/.local/bin:$PATH
export LD_LIBRARY_PATH=${PYTHONHOME}/lib:${LD_LIBRARY_PATH}
export DRMAA_LIBRARY_PATH=/geadmin/N1GE/lib/lx-amd64/libdrmaa.so.1.0

ulimit -u 4096

TARGET_PATH=`dirname $0`

target_pipeline=$1
sample_conf=$2
project_dir=$3
genomon_conf=$4

echo "Genomon is checking parameters ..."
${TARGET_PATH}/genomon_pipeline --param_check $target_pipeline $sample_conf $project_dir $genomon_conf || exit $?
echo "Parameters check is complete."

mkdir -p ${project_dir}/log || exit $?
echo "Genomon created the '${project_dir}/log' directory"

/home/geadmin/N1GE/bin/lx-amd64/orgqsub -o ${project_dir}/log -e ${project_dir}/log -l s_vmem=64G,mem_req=128G -l ljob -r no ${TARGET_PATH}/qsub_genomon_pipeline_HGC.sh $target_pipeline $sample_conf $project_dir $genomon_conf ${TARGET_PATH}


