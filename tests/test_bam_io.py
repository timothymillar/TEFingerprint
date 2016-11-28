#! /usr/bin/env python

from tectoolkit import bam_io
import numpy as np
import numpy.testing as npt


def test_parse_flag():
    """Test function parse_flag"""
    flag = 0
    answer = np.array([False, False, False, False, False, False, False, False, False, False, False, False, ])
    npt.assert_array_equal(bam_io.parse_flag(flag), answer)

    flag = 1365
    answer = np.array([True, False, True, False, True, False, True, False, True, False, True, False])
    npt.assert_array_equal(bam_io.parse_flag(flag), answer)

    flag = 4095
    answer = np.array([True, True, True, True, True, True, True, True, True, True, True, True])
    npt.assert_array_equal(bam_io.parse_flag(flag), answer)


def test_flag_orientation():
    """Test for function flag_orientation
    Returns '+' for forward, '-' for reverse, None for unmapped"""
    flag = 0
    answer = '+'
    assert bam_io.flag_orientation(flag) == answer

    flag = 16
    answer = '-'
    assert bam_io.flag_orientation(flag) == answer

    flag = 20
    answer = None
    assert bam_io.flag_orientation(flag) == answer


def test_flag_attributes():
    """Test for function flag_attributes"""
    flag = 0
    answer = {'first in pair': False,
              'mate reverse strand': False,
              'mate unmapped': False,
              'not primary alignment': False,
              'read fails platform / vendor quality checks': False,
              'read is PCR or optical duplicate': False,
              'read mapped in proper pair': False,
              'read paired': False,
              'read reverse strand': False,
              'read unmapped': False,
              'second in pair': False,
              'supplementary alignment': False}
    assert bam_io.flag_attributes(flag) == answer

    flag = 1
    answer = {'first in pair': False,
              'mate reverse strand': False,
              'mate unmapped': False,
              'not primary alignment': False,
              'read fails platform / vendor quality checks': False,
              'read is PCR or optical duplicate': False,
              'read mapped in proper pair': False,
              'read paired': True,
              'read reverse strand': False,
              'read unmapped': False,
              'second in pair': False,
              'supplementary alignment': False}
    assert bam_io.flag_attributes(flag) == answer

    flag = 20
    answer = {'first in pair': False,
              'mate reverse strand': False,
              'mate unmapped': False,
              'not primary alignment': False,
              'read fails platform / vendor quality checks': False,
              'read is PCR or optical duplicate': False,
              'read mapped in proper pair': False,
              'read paired': False,
              'read reverse strand': True,
              'read unmapped': True,
              'second in pair': False,
              'supplementary alignment': False}
    assert bam_io.flag_attributes(flag) == answer

    flag = 4095
    answer = {'first in pair': True,
              'mate reverse strand': True,
              'mate unmapped': True,
              'not primary alignment': True,
              'read fails platform / vendor quality checks': True,
              'read is PCR or optical duplicate': True,
              'read mapped in proper pair': True,
              'read paired': True,
              'read reverse strand': True,
              'read unmapped': True,
              'second in pair': True,
              'supplementary alignment': True}
    assert bam_io.flag_attributes(flag) == answer