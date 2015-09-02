#!/usr/bin/env python

from distutils.core import setup

setup(name='genomon_pipeline',
      version='2.0.0',
      description='Python tools for running genomon pipeline for cancer genome and transcriptome sequencing analysis',
      author='Kenichi Chiba, Eigo Shimizu, Yuichi Shiraishi',
      author_email='genomon_team@gamil.com',
      url='https://github.com/Genomon-Project/Genomon.git',
      package_dir = {'': 'lib'},
      packages=['genomon_pipeline'],
      scripts=['genomon_pipeline'],
      license='GPL-3'
     )
