#! /usr/bin/env python

from tefingerprint import bamio
import pytest


@pytest.mark.parametrize("string,answer",
                         [('103S111M1D36M', 148)])
def test_cigar_mapped_length(string, answer):
    assert bamio._cigar_mapped_length(string) == answer


@pytest.mark.parametrize("strings,answer",
                         [([
                               "readA\t0\tchr1\t524\t60\t103S111M1D36M\t*\t0\t0\tGAAATACGT\t1/F11<12G\tNM:i:5\tMD:Z:6T24C35A43^T5T30\tAS:i:120\tXS:i:0\tME:Z:Copia-7_VV",
                               "readB\t0\tchr1\t600\t0\t103S111M1D36M\t*\t0\t0\tGAAATACGT\t1/F11<12G\tNM:i:5\tMD:Z:6T24C35A43^T5T30\tAS:i:120\tXS:i:0\tME:Z:Copia-7_VV",
                               "readC\t0\tchr1\t5\t0\t3M3I2M5S\t*\t0\t0\tGAAATACGT\t1/F11<12G\tNM:i:5\tMD:Z:6T24C35A43^T5T30\tAS:i:120\tXS:i:0\tME:Z:Copia-7_VV"],
                           [(524, 671, 'readA'), (600, 747, 'readB'), (5, 9, 'readC')])])
def test_parse_read_loci(strings, answer):
    """
    Test factory for hidden method _parse_sam_strings.
    Tests for parsing sets of forward, reverse and mixed reads.
    """
    query = list(bamio._parse_read_loci(strings))
    assert query == answer
