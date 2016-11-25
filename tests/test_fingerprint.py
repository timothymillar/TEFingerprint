#! /usr/bin/env python

import numpy as np
import numpy.testing as npt
from tectoolkit.classes import ReadGroup
from tectoolkit.fingerprint import Fingerprint
from tectoolkit.compare import FingerprintComparison


class TestReadGroup:
    def test__init__(self):
        """
        Test for __init__ method.
        Read array should be copied.
        """
        input_reads = np.array([(8, 2, '+', 'a'), (5, 3, '+', 'b'), (7, 5, '+', 'c')], dtype=ReadGroup.DTYPE_READ)
        query = ReadGroup(input_reads)
        assert query.reads is not input_reads
        npt.assert_array_equal(query.reads, input_reads)

    def test__iter__(self):
        """
        Test for __iter__ method.
        Iterating over :class:`ReadGroup` should iterate of the the wrapped numpy array of reads.
        """
        input_reads = np.array([(8, 2, '+', 'a'), (5, 3, '+', 'b'), (7, 5, '+', 'c')], dtype=ReadGroup.DTYPE_READ)
        query = ReadGroup(input_reads)
        npt.assert_array_equal(np.array([r for r in query], dtype=ReadGroup.DTYPE_READ), input_reads)

    def test__getitem__(self):
        """
        Test for __getitem__ method.
        Indexing should pass through to the wrapped numpy array of reads.
        """
        input_reads = np.array([(8, 2, '+', 'a'), (5, 3, '+', 'b'), (7, 5, '+', 'c')], dtype=ReadGroup.DTYPE_READ)
        query = ReadGroup(input_reads)
        npt.assert_array_equal(query['tip'], np.array([8, 5, 7]))
        npt.assert_array_equal(query['name'], np.array(['a', 'b', 'c']))
        npt.assert_array_equal(query[0:2], input_reads[0:2])
        assert tuple(query[1]) == (5, 3, '+', 'b')

    def test__len__(self):
        """
        Test for __len__ method.
        Should return len wrapped numpy array of reads.
        """
        input_reads = np.array([(8, 2, '+', 'a'), (5, 3, '+', 'b'), (7, 5, '+', 'c')], dtype=ReadGroup.DTYPE_READ)
        query = ReadGroup(input_reads)
        assert len(query) == 3

    def test_subset_by_locus(self):
        """
        Test for method subset_by_locus.
        Test for subsetting by 'tip' or 'tail' end of reads.
        """
        # subset using default end='tip'
        query = ReadGroup(np.array([(8, 2, '+', 'a'), (5, 3, '+', 'b'), (7, 5, '+', 'c')], dtype=ReadGroup.DTYPE_READ))
        answer = ReadGroup(np.array([(5, 3, '+', 'b'), (7, 5, '+', 'c')], dtype=ReadGroup.DTYPE_READ))
        npt.assert_array_equal(query.subset_by_locus(5, 7).reads, answer.reads)

        # subset using end='tail'
        query = ReadGroup(np.array([(8, 2, '+', 'a'), (5, 3, '+', 'b'), (7, 5, '+', 'c')], dtype=ReadGroup.DTYPE_READ))
        answer = ReadGroup(np.array([(8, 2, '+', 'a'), (5, 3, '+', 'b')], dtype=ReadGroup.DTYPE_READ))
        npt.assert_array_equal(query.subset_by_locus(1, 3, end='tail').reads, answer.reads)

    def test_parse_sam_strings(self):
        """
        Test for hidden method _parse_sam_strings.
        Test for parsing both forwards and reverse reads.
        """
        # forwards reads
        input_strings = ["Gypsy27_a\t0\tchr1\t2\t66M19S\t*\t0\t0\tAACCCTA\tFFFIIII\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66",
                         "Gypsy27_b\t0\tchr1\t3\t66M19S\t*\t0\t0\tACC\tFFI\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66",
                         "Gypsy27_c\t0\tchr1\t5\t66M19S\t*\t0\t0\tAAC\tFFF\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66"]
        query = list(ReadGroup._parse_sam_strings(input_strings))
        query.sort()
        answer = [(8, 2, '+', 'Gypsy27_a'), (5, 3, '+', 'Gypsy27_b'), (7, 5, '+', 'Gypsy27_c')]
        answer.sort()
        npt.assert_array_equal(query, answer)

        # reverse reads
        input_strings = ["Gypsy27_a\t16\tchr1\t2\t66M19S\t*\t0\t0\tAACCCTA\tFFFIIII\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66",
                         "Gypsy27_b\t16\tchr1\t4\t66M19S\t*\t0\t0\tACC\tFFI\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66",
                         "Gypsy27_c\t16\tchr1\t7\t66M19S\t*\t0\t0\tAAC\tFFF\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66"]
        query = list(ReadGroup._parse_sam_strings(input_strings))
        query.sort()
        answer = [(2, 8, '-', 'Gypsy27_a'), (4, 6, '-', 'Gypsy27_b'), (7, 9, '-', 'Gypsy27_c')]
        answer.sort()
        npt.assert_array_equal(query, answer)

    def test_from_sam_strings(self):
        """
        Test for method from_sam_strings.
        Test for parsing both forwards and reverse reads.
        """
        # forwards reads
        input_strings = ["Gypsy27_a\t0\tchr1\t2\t66M19S\t*\t0\t0\tAACCCTA\tFFFIIII\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66",
                         "Gypsy27_b\t0\tchr1\t3\t66M19S\t*\t0\t0\tACC\tFFI\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66",
                         "Gypsy27_c\t0\tchr1\t5\t66M19S\t*\t0\t0\tAAC\tFFF\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66"]
        query = ReadGroup.from_sam_strings(input_strings)
        query.sort(order='name')
        answer = ReadGroup(np.array([(8, 2, '+', 'Gypsy27_a'),
                                     (5, 3, '+', 'Gypsy27_b'),
                                     (7, 5, '+', 'Gypsy27_c')], dtype=ReadGroup.DTYPE_READ))
        npt.assert_array_equal(query.reads, answer.reads)

        # reverse reads
        input_strings = ["Gypsy27_a\t16\tchr1\t2\t66M19S\t*\t0\t0\tAACCCTA\tFFFIIII\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66",
                         "Gypsy27_b\t16\tchr1\t4\t66M19S\t*\t0\t0\tACC\tFFI\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66",
                         "Gypsy27_c\t16\tchr1\t7\t66M19S\t*\t0\t0\tAAC\tFFF\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66"]
        query = ReadGroup.from_sam_strings(input_strings)
        query.sort(order='name')
        answer = ReadGroup(np.array([(2, 8, '-', 'Gypsy27_a'),
                                     (4, 6, '-', 'Gypsy27_b'),
                                     (7, 9, '-', 'Gypsy27_c')], dtype=ReadGroup.DTYPE_READ))
        npt.assert_array_equal(query.reads, answer.reads)
