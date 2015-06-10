export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/package/python2.7/current/lib
export PYTHONPATH=$PYTHONPATH:/home/w3varann/.local/lib/python2.7/site-packages
/usr/local/package/python2.7/2.7.2/bin/python $4 ./job_file_check.py \
        -k ./db/job_file_words.yaml \
        -s $1 \
        -j $2 \
        -p $3
