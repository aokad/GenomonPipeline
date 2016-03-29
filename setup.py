#!/usr/bin/env python

from distutils.core import setup
exec(open('scripts/genomon_pipeline/__init__.py').read())

setup(name='genomon_pipeline',
      version=__version__,
      description='Python tools for running genomon pipeline for cancer genome and transcriptome sequencing analysis',
      author='Kenichi Chiba, Eigo Shimizu, Yuichi Shiraishi',
      author_email='genomon_team@gamil.com',
      url='https://github.com/Genomon-Project/Genomon.git',
      package_dir = {'': 'scripts'},
      packages=['genomon_pipeline', 'genomon_pipeline.rna_resource', 'genomon_pipeline.dna_resource', 'genomon_pipeline.config'],
      scripts=['genomon_pipeline'],
      license='GPL-3'
     )
