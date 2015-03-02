# Genomon
Genomon pipeline in 2015

Test version is released on 2015/02/27.


How to test it.

First you need to install python modules
```
  pip install --user ruffus
  pip install --user mpi4py
  pip install --user yaml
  and so on.
```

1. Download zip file and extract it.
2. Make directory structure like the following.
```
├── scripts --> copy *.py scripts here.
└── test  --> this is the path for project_root in genomon.job
    └── project_name  --> this is proect in genomon.job.
        └── data
            └── 20150227 --> this is sample_date in genomon.job
                └── test01 --> this is sample_name in genomon.job.


```
3. Copy fastq files to test/project_name/data/20150227/test01.
4. Copy *.py files to scripts.
5. Edit genomon.cfg to specify tools and reference.
6. Edit genomon.yaml to specify project_root, project_name, sample_name, and sample_date
7. Edit test.sh and specify the the path to genomon.py, genomon.cfg, genomon.yaml.
8. Run test.sh
