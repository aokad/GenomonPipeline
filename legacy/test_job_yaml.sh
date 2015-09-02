export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/package/python2.7/current/lib
export PYTHONPATH=/home/w3varann/.local/lib/python2.7/site-packages

/usr/local/package/python2.7/current/bin/python $4 /home/eigos/Data/Genomon/job_file_check \
        -k /home/w3varann/tools/Genomon/resource/job_file_words.yaml \
        -s $1 \
        -j $2 \
        -p $3
