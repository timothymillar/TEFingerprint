#! /usr/bin/env python

import os
from setuptools import setup


def read_file(file_name):
    return os.path.join(os.path.dirname(__file__), file_name)


# read readme as long description
long_description = read_file('README.rst')


# read version from single location
__version__ = 'undefined'
version_file = os.path.dirname(os.path.realpath(__file__)) + \
               '/tefingerprint/version.py'
exec(open(version_file).read())


# set download url based on version
download_url = ('https://github.com/PlantandFoodResearch/TEFingerprint/'
                'archive/'
                'v' + __version__ + '.tar.gz')


setup(name='tefingerprint',
      version=__version__,
      author='Tim Millar',
      author_email='tim.millar@plantandfood.co.nz',
      url='https://github.com/PlantandFoodResearch/TEFingerprint',
      download_url=download_url,
      description='Toolkit for identifying transposon movement',
      long_description=long_description,
      scripts=['applications/tefingerprint',
               'applications/tef-extract-informative',
               'applications/tef-filter-gff'],
      packages=['tefingerprint',
                'tefingerprint.util',
                'tefingerprint.util.numpy',
                'tefingerprint._applications'],
      keywords=['biology',
                'bioinformatics',
                'transposon',
                'transposable element',
                'cluster'],
      classifiers=['Development Status :: 3 - Alpha',
                   'Intended Audience :: Science/Research',
                   'Programming Language :: Python :: 3.4',
                   'Topic :: Scientific/Engineering',
                   'Operating System :: Unix']
      )
