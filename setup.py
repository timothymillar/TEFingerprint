#! /usr/bin/env python

import os
from setuptools import setup


def read_file(file_name):
    return os.path.join(os.path.dirname(__file__), file_name)


setup(name='tefingerprint',
      version='0.1.3',
      author='Tim Millar',
      author_email='tim.millar@plantandfood.co.nz',
      url='https://github.com/PlantandFoodResearch/TEFingerprint',
      description='Toolkit for identifying transposon movement',
      long_description=read_file('README.rst'),
      scripts=['applications/tefingerprint',
               'applications/tefingerprint-extract-informative-reads'],
      packages=['tefingerprint'],
      classifiers=['Development Status :: 3 - Alpha']
      )
