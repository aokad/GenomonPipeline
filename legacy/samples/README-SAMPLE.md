# Sample configure files

This directory contains some sample parameter files for genomon.

---

# 1. Description

Running on Shirokane3, this files are not necessary to edit.

files:

```
{Genomon root}
└── sample
    ├── README-SAMPLE.md................This file
    ├── genomon.cfg.....................Settings for genomon required tools
    ├── templete_job.yaml...............Templete job file
    ├── templete_param.yaml.............Templete param file
    ├── sample_xxx_job.yaml.............Example job file (paired with param file)
    ├── sample_xxx_param.yaml...........Example param file (paired with job file)
    ├── sample.csv......................Example sample list(csv file)
    ├── sample.tsv......................Example sample list(tsv file)
    └── sample.xlsx.....................Example sample list(Excel file)

```

---

# 2. Sample Files

## check for initial setting 

[DNA]
 - sample1_job.yaml
 - sample1_param.yaml

## demo

### DNA analysis demo 1.

Use bwa_mem, fisher_mutation_call.

 - sample_dna1_job.yaml
 - sample_dna1_param.yaml

### DNA analysis demo 2.

Use bwa_mem, bam_stats.

 - sample_dna2_job.yaml
 - sample_dna2_param.yaml

### RNA analysis demo

Use star.

 - sample_rna1_job.yaml
 - sample_rna1_param.yaml


