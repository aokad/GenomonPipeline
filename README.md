# Genomon
Genomon pipeline in 2015

Test version is released on 2015/03/03.


How to test it.

First you need to install python modules
```
  pip install --user ruffus
  pip install --user pyyaml
```

1. Download Genomnn.zip and extract it.
2. Copy Genomon-master/test/run.sh to your working directory.
3. Edit the run.sh and change the path to genomon.py and genomon.cfg in Genomon-master.
4. Copy genomon_job.yaml to your working directory.
5. Edit the yaml file so that the file describes your job properly. Especially, project_root, input_file_dir, input_file_type, file_name, file_ext, pair_id.
6. Run run.sh [your job yaml file].


Or if you are on HGC Shirokane2 or Shirokane3, copy one of job yaml files from /home/w3varann/samples/*.yaml, edit it so that the file describes your job configuration and run the following command
```
/home/w3varann/tools/Genomon/test/run.sh [your job yaml file]
```
