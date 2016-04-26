GenomonPipeline
===============
Genomon is a cancer genome data analysis pipeline.

AUTHOR/SUPPORT
==============
Kenichi Chiba, Ai Okada and Yuichi Shiraishi

MANUAL
======
http://genomon-project.github.io/GenomonPages/

DEPENDENCY
==========
[ruffus](http://www.ruffus.org.uk/)
[drmaa](https://www.drmaa.org/)
[PyYAML](http://pyyaml.org/)

INSTALL
=======

```bash

# get latest genomon-pipeline
wget https://github.com/Genomon-Project/GenomonPipeline/archive/v2.2.5.tar.gz
tar xzvf v2.2.0.tar.gz
cd GenomonPipeline-2.2.0

# install genomon pipeline
python setup.py install --user

# get ruffus
wget https://github.com/bunbun/ruffus/archive/v2.6.3.tar.gz
tar xzvf v2.6.3.tar.gz
cd ruffus-2.6.3

# install ruffus
python setup.py install --user

# get PyYAML
git clone https://github.com/ravenac95/PyYAML
cd PyYAML

# install PyYAML
python setup.py install --user

# install drmaa
pip install drmaa --user


