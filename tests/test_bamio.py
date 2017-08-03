#! /usr/bin/env python

from tefingerprint import bamio
import pytest
import os
DATA_PATH = os.path.dirname(os.path.realpath(__file__)) + '/data/'


def test_extract_bam_references():
    query = bamio.extract_bam_reference_strings(DATA_PATH + 'testA-2017-06-08.bam')
    answer = ['chr1:1-25000']
    assert query == answer


def test_extract_multiple_bam_references():
    query = bamio.extract_bam_reference_strings(DATA_PATH + 'testA-2017-06-08.bam',
                                                DATA_PATH + 'testB-2017-06-08.bam')
    answer = ['chr1:1-25000']
    assert query == answer
