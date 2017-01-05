"""
A script to setup the package.

Created on Fri Mar 11 11:27:39 2016

@author: dangeles
"""
# -*- coding: utf-8 -*-

from distutils.core import setup
from setuptools import setup, find_packages
import os
import sys

version = '0.15.00'

# just type in python setup.py publish and
# this takes care of publishing to pypi!

# tag with git
if sys.argv[1] == 'tag':

    if len(sys.argv) > 2:
        if sys.argv[2] == '-m':
            print(sys.argv[2])
            print('please input your commit message')
            message = sys.stdin.readline().strip()
            print("git tag -a %s -m '%s'" % (version, message))
            os.system("git tag -a %s -m '%s'" % (version, message))
    else:
        os.system("git tag -a %s -m 'version %s'" % (version, version))

    os.system("git push --tags")
    sys.exit()

# publish to pypi
if sys.argv[-1] == 'publish':
    os.system("python setup.py sdist upload")
    os.system("python setup.py bdist_wheel upload")
#    print("You probably want to also tag the version now:")
#    print("  git tag -a %s -m 'version %s'" % (version, version))
#    print("  git push --tags")
    sys.exit()


# for pytest.py purposes -- just type in python setup.py test!
if sys.argv[-1] == 'test':
    test_requirements = [
        'pytest',
        'flake8',
        'coverage'
    ]
    try:
        modules = map(__import__, test_requirements)
    except ImportError as e:
        err_msg = e.message.replace("No module named ", "")
        msg = "%s is not installed. Install your test requirements." % err_msg
        raise ImportError(msg)
    os.system('py.test')
    sys.exit()

# Below this point is the rest of the setup() function


def readme():
    """A function to open a readme.rst file."""
    with open('README.rst') as f:
        return f.read()

setup(
  name='tissue_enrichment_analysis',
  packages=find_packages(exclude=("tests",)),
  version=version,
  description='This package contains all the software used to implement\
  TEA in WormBase and remotely',
  author='David Angeles-Albores',
  author_email='dangeles@caltech.edu',
  url='https://github.com/dangeles/TissueEnrichmentAnalysis',  # github repo
  download_url='https://github.com/dangeles/TissueEnrichmentAnalysis/tarball/{0}'.format(version),
  keywords=['tissue enrichment analysis', 'TEA',
            'RNAseq', 'celegans', 'biology'],  # arbitrary keywords
  install_requires=[
          'matplotlib', 'scipy', 'numpy'
      ],
  classifiers=['License :: OSI Approved :: MIT License',
               'Programming Language :: Python :: 3.5'
               ],
  licenses='MIT',
  scripts=['bin/tea']
)
