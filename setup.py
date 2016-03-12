# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 11:27:39 2016

@author: dangeles
"""

from distutils.core import setup

version= 0.13

setup(
  name = 'tissue_enrichment_analysis',
  packages = ['tissue_enrichment_analysis'], # this must be the same as the name above
  version = version,
  description = 'This package contains all the software used to implement\
  TEA in WormBase and remotely',
  author = 'David Angeles-Albores',
  author_email = 'dangeles@caltech.edu',
  url = 'https://github.com/dangeles/TissueEnrichmentAnalysis', # use the URL to the github repo
  download_url = 'https://github.com/dangeles/TissueEnrichmentAnalysis/tarball/{0}'.format(version),
  keywords = ['tissue enrichment', 'celegans', 'biology'], # arbitrary keywords
  install_requires=[
          'matplotlib', 'scipy', 'numpy', 'urllib', 'contextlib'
      ],
  classifiers = [
  'Programming Language :: Python :: 3.5'
  ],
  licenses= 'MIT'
)

