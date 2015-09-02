#!/usr/bin/env python
"""Setup file and install script for NextGen sequencing analysis scripts.
"""
import sys
import os
import inspect
from setuptools import setup, find_packages

version = "0.1.1"
exec_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) 
os.chdir( exec_dir )
with open( exec_dir + "/requirements.txt", "r") as f:
    install_requires = [x.strip() for x in f.readlines()]

setup(name="Genomon",
      version=version,
      author="Genomon team at HGC.jp.",
      author_email="eigos@hgc.jp",
      description="Genomon pipelines for sequencing analysis",
      long_description=(open( exec_dir + '/README.md').read()),
      license="HGC",
      url="https://github.com/Genomon-Project/Genomon",
      packages=find_packages(),
      scripts=['Genomon', 'job_file_check', 'setenv.sh', 'test_job_yaml.sh', 'run.sh' ],
      install_requires=install_requires)
