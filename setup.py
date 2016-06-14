#!/usr/bin/env python

from distutils.core import setup

setup(name='genomon_pipeline',
      version='2.3.0',
      description='Python tools for running genomon pipeline for cancer genome and transcriptome sequencing analysis',
      author='Kenichi Chiba, Ai Okada and Yuichi Shiraishi',
      author_email='genomon.devel@gmail.com',
      url='https://github.com/Genomon-Project/Genomon.git',
      package_dir = {'': 'scripts'},
      packages=['genomon_pipeline', 'genomon_pipeline.rna_resource', 'genomon_pipeline.dna_resource', 'genomon_pipeline.config'],
      scripts=['genomon_pipeline'],
      license='License of GenomonPipeline'
     )
