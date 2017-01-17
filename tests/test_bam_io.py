#! /usr/bin/env python

from tectoolkit import bam_io
import numpy as np
import numpy.testing as npt
import pytest


@pytest.mark.parametrize("flag,answer",
                         [(0, np.array([False,
                                        False,
                                        False,
                                        False,
                                        False,
                                        False,
                                        False,
                                        False,
                                        False,
                                        False,
                                        False,
                                        False])),
                          (1365, np.array([True,
                                           False,
                                           True,
                                           False,
                                           True,
                                           False,
                                           True,
                                           False,
                                           True,
                                           False,
                                           True,
                                           False])),
                          (4095, np.array([True,
                                           True,
                                           True,
                                           True,
                                           True,
                                           True,
                                           True,
                                           True,
                                           True,
                                           True,
                                           True,
                                           True]))],
                         ids=['flag: 0', 'flag: 1365', 'flag: 4095'])
def test_parse_flag(flag, answer):
    """Test factory function parse_flag"""
    npt.assert_array_equal(bam_io.parse_flag(flag), answer)


@pytest.mark.parametrize("flag,orientation",
                         [(0, '+'),
                          (16, '-'),
                          (20, None)],
                         ids=['forward', 'reverse', 'unmapped'])
def test_flag_orientation(flag, orientation):
    """Test factory for function flag_orientation
    Returns '+' for forward, '-' for reverse, None for unmapped"""
    assert bam_io.flag_orientation(flag) == orientation


@pytest.mark.parametrize("flag,attributes",
                         [(0, {'first in pair': False,
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
                               'supplementary alignment': False}),
                          (1, {'first in pair': False,
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
                               'supplementary alignment': False}),
                          (20, {'first in pair': False,
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
                                'supplementary alignment': False}),
                          (4095, {'first in pair': True,
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
                                  'supplementary alignment': True})],
                         ids=['flag: 0', 'flag: 1', 'flag: 20', 'flag: 4095'])
def test_flag_attributes(flag, attributes):
    """Test factory for function flag_attributes"""
    assert bam_io.flag_attributes(flag) == attributes


@pytest.mark.parametrize("string,answer",
                         [('103S111M1D36M', 148)])
def test_cigar_mapped_length(string, answer):
    assert bam_io._cigar_mapped_length(string) == answer


@pytest.mark.parametrize("strings,answer",
                         [([
                               "readA\t0\tchr1\t524\t60\t103S111M1D36M\t*\t0\t0\tGAAATACGT\t1/F11<12G\tNM:i:5\tMD:Z:6T24C35A43^T5T30\tAS:i:120\tXS:i:0\tME:Z:Copia-7_VV",
                               "readB\t16\tchr1\t600\t0\t103S111M1D36M\t*\t0\t0\tGAAATACGT\t1/F11<12G\tNM:i:5\tMD:Z:6T24C35A43^T5T30\tAS:i:120\tXS:i:0\tME:Z:Copia-7_VV",
                               "readC\t0\tchr1\t5\t0\t3M3I2M5S\t*\t0\t0\tGAAATACGT\t1/F11<12G\tNM:i:5\tMD:Z:6T24C35A43^T5T30\tAS:i:120\tXS:i:0\tME:Z:Copia-7_VV"],
                           [(671, 524, '+', 'readA'), (600, 747, '-', 'readB'), (9, 5, '+', 'readC')])])
def test_parse_sam_strings(strings, answer):
    """
    Test factory for hidden method _parse_sam_strings.
    Tests for parsing sets of forward, reverse and mixed reads.
    """
    query = list(bam_io._parse_sam_strings(strings))
    assert query == answer
