# Genomon
Genomon pipeline in 2015

Test version is released on 2015/06/05.


How to run.

First you need to install python modules
```
  pip install --user ruffus
  pip install --user pyyaml
  pip install --user drmaa
```
Then, install the following tools.
```
  bedtools-2
  biobambam
  bowtie2
  bwa
  cufflinks-2
  GATK-3
  Genomon
  libmaus
  PCAP-core-dev
  picard-tools
  samtools
  tophat-2
```
Install R modules.
```
  cummeRbund
```
Finally,

1. Download Genomon-master.zip and extract it.
2. Copy Genomon-master/test/run.sh to your working directory.
3. Edit the run.sh and change the path to genomon.py and genomon.cfg in Genomon-master.
4. Copy [job description yaml file (sample?_job.yaml at samples directory ] to your working directory.
5. Edit the copied yaml file so that the file describes your job configuration. Especially, project_root, input_file_dir, input_file_type, file_name, file_ext, pair_id.
6. Copy [analysis parameter file (sample?_param.yaml at samples directory ] to your working directory.
7. Edit the copied yaml file, if necessary.
8. Run run.sh [your job yaml file] [your analysis parameter yaml file].


Or if you are on HGC Shirokane2 or Shirokane3, 

1. Copy a pair of job yaml and parameter yaml files from /home/w3varann/samples/*.yaml.
2. Edit them so that the files describe your job configuration.
3. Run the following command.
```
/home/w3varann/tools/Genomon/run.sh [your job yaml file] [you analysis parameter yaml file]
```
