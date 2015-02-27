# Genomon
Genomon pipeline in 2015

Test version is released on 2015/02/27.


How to test it.

First you need to install python modules
  pip install --user ruffus
  pip install --user mpi4py
  pip install --user yaml


1. Download zip file and extract it.
2. Make directory structure like the following.


3. Edit genomon.cfg to specify tools and reference.
4. Edit genomon.yaml to specify project_root, project_name, sample_name, and sample_date
5. Edit test.sh and specify the the path to genomon.py, genomon.cfg, genomon.yaml.
