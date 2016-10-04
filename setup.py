#! /usr/bin/env python

import os
from setuptools import setup


def read_file(file_name):
    return os.path.join(os.path.dirname(__file__), file_name)


setup(name='tectoolkit',
      version='0.0.1',
      author='Tim Millar',
      author_email='tim.millar@plantandfood.co.nz',
      url='https://github.com/PlantandFoodResearch/TECtoolkit',
      description='Toolkit for identifying transposable element movement in regenerant clones',
      long_description=read_file('README.MD'),
      scripts=['applications/tec'],
      packages=['tectoolkit'],
      classifiers=['Development Status :: 2 - Pre-Alpha']
      )
