#! /usr/bin/env python

import pytest
import numpy as np
import numpy.testing as npt
from tectoolkit.reads import ReadGroup, UnivariateLoci


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

    @pytest.mark.parametrize("iterable",
                             [[(8, 2, '+', 'a')],
                              [(8, 2, '+', 'a'), (5, 3, '+', 'b'), (7, 5, '+', 'c')],
                              [(2, 8, '-', 'Gypsy27_a'), (5, 3, '+', 'Gypsy27_b'), (7, 9, '-', 'Gypsy27_c')]],
                             ids=['single', 'forward', 'mixed'])
    def test_from_iter(self, iterable):
        """
        Test factory for method from_iter using generators.
        """
        generator = (item for item in iterable)
        query = ReadGroup.from_iter(generator)
        answer = ReadGroup(np.array(iterable, dtype=ReadGroup.DTYPE_READ))
        query.sort()
        answer.sort()
        npt.assert_array_equal(query, answer)

    @pytest.mark.parametrize("query,answer,locus,end",
                             # loci contain read tips
                             [([(8, 2, '+', 'a'), (5, 3, '+', 'b'), (7, 5, '+', 'c')],
                               [(5, 3, '+', 'b'), (7, 5, '+', 'c')],
                               (5, 7),
                               'tip'),
                              # loci contain read tails
                              ([(8, 2, '+', 'a'), (5, 3, '+', 'b'), (7, 5, '+', 'c')],
                               [(5, 3, '+', 'b'), (8, 2, '+', 'a')],
                               (1, 3),
                               'tail'),
                              # loci contain read tips and tails
                              ([(8, 3, '+', 'a'), (5, 3, '+', 'b'), (5, 2, '+', 'c')],
                               [(5, 3, '+', 'b')],
                               (3, 5),
                               'both')],
                             ids=['tips', 'tails', 'mixed'])
    def test_subset_by_locus(self, query, answer, locus, end):
        """
        Test factory for method subset_by_locus.
        Tests for subsetting by 'tip', 'tail' or 'both' end(s) of reads.
        """
        query = ReadGroup.from_iter(query)
        query.sort()
        answer = ReadGroup.from_iter(answer)
        answer.sort()
        npt.assert_array_equal(query.subset_by_locus(*locus, end=end), answer)

    @pytest.mark.parametrize("strings,answer",
                             # forward reads
                             [(["Gypsy27_a\t0\tchr1\t2\t66M19S\t*\t0\t0\tAACCCTA\tFFFIIII\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66",
                                "Gypsy27_b\t0\tchr1\t3\t66M19S\t*\t0\t0\tACC\tFFI\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66",
                                "Gypsy27_c\t0\tchr1\t5\t66M19S\t*\t0\t0\tAAC\tFFF\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66"],
                               [(8, 2, '+', 'Gypsy27_a'), (5, 3, '+', 'Gypsy27_b'), (7, 5, '+', 'Gypsy27_c')]),
                              # reverse reads
                              (["Gypsy27_a\t16\tchr1\t2\t66M19S\t*\t0\t0\tAACCCTA\tFFFIIII\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66",
                                "Gypsy27_b\t16\tchr1\t4\t66M19S\t*\t0\t0\tACC\tFFI\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66",
                                "Gypsy27_c\t16\tchr1\t7\t66M19S\t*\t0\t0\tAAC\tFFF\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66"],
                               [(2, 8, '-', 'Gypsy27_a'), (4, 6, '-', 'Gypsy27_b'), (7, 9, '-', 'Gypsy27_c')]),
                              # forward and reverse reads
                              (["Gypsy27_a\t0\tchr1\t2\t66M19S\t*\t0\t0\tAACCCTA\tFFFIIII\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66",
                                "Gypsy27_b\t0\tchr1\t3\t66M19S\t*\t0\t0\tACC\tFFI\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66",
                                "Gypsy27_c\t16\tchr1\t7\t66M19S\t*\t0\t0\tAAC\tFFF\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66"],
                               [(8, 2, '+', 'Gypsy27_a'), (5, 3, '+', 'Gypsy27_b'), (7, 9, '-', 'Gypsy27_c')]
                               )],
                             ids=['tips', 'tails', 'mixed'])
    def test_parse_sam_strings(self, strings, answer):
        """
        Test factory for hidden method _parse_sam_strings.
        Tests for parsing sets of forward, reverse and mixed reads.
        """
        query = list(ReadGroup._parse_sam_strings(strings))
        query.sort()
        answer.sort()
        assert query == answer

    @pytest.mark.parametrize("strings,answer",
                             # forward reads
                             [(["Gypsy27_a\t0\tchr1\t2\t66M19S\t*\t0\t0\tAACCCTA\tFFFIIII\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66",
                                "Gypsy27_b\t0\tchr1\t3\t66M19S\t*\t0\t0\tACC\tFFI\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66",
                                "Gypsy27_c\t0\tchr1\t5\t66M19S\t*\t0\t0\tAAC\tFFF\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66"],
                               ReadGroup.from_iter([(8, 2, '+', 'Gypsy27_a'),
                                                    (5, 3, '+', 'Gypsy27_b'),
                                                    (7, 5, '+', 'Gypsy27_c')])),
                              # reverse reads
                              (["Gypsy27_a\t16\tchr1\t2\t66M19S\t*\t0\t0\tAACCCTA\tFFFIIII\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66",
                                "Gypsy27_b\t16\tchr1\t4\t66M19S\t*\t0\t0\tACC\tFFI\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66",
                                "Gypsy27_c\t16\tchr1\t7\t66M19S\t*\t0\t0\tAAC\tFFF\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66"],
                               ReadGroup.from_iter([(2, 8, '-', 'Gypsy27_a'),
                                                    (4, 6, '-', 'Gypsy27_b'),
                                                    (7, 9, '-', 'Gypsy27_c')])),
                              # forward and reverse reads
                              (["Gypsy27_a\t16\tchr1\t2\t66M19S\t*\t0\t0\tAACCCTA\tFFFIIII\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66",
                                "Gypsy27_b\t0\tchr1\t3\t66M19S\t*\t0\t0\tACC\tFFI\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66",
                                "Gypsy27_c\t16\tchr1\t7\t66M19S\t*\t0\t0\tAAC\tFFF\tNM:i:0\tMD:Z:66\tAS:i:66\tXS:i:66"],
                               ReadGroup.from_iter([(2, 8, '-', 'Gypsy27_a'),
                                                    (5, 3, '+', 'Gypsy27_b'),
                                                    (7, 9, '-', 'Gypsy27_c')]))],
                              ids=['tips', 'tails', 'mixed'])
    def test_from_sam_strings(self, strings, answer):
        """
        Test for method from_sam_strings.
        Test for parsing both forwards and reverse reads.
        """
        query = ReadGroup._from_sam_strings(strings)
        query.sort(order='name')
        answer.sort(order='name')
        npt.assert_array_equal(query, answer)


class TestUnivariateLoci:
    """
    Tests for class ReadLoci
    """
    def test_sort(self):
        """
        Test for method sort.
        """
        input_loci = np.array([(2, 4),
                               (3, 4),
                               (3, 3),
                               (4, 4),
                               (3, 99),
                               (1, 1)], dtype=UnivariateLoci.DTYPE_ULOCUS)
        query = UnivariateLoci(input_loci)
        query.sort()
        answer = np.array([(1, 1),
                           (2, 4),
                           (3, 3),
                           (3, 4),
                           (3, 99),
                           (4, 4)], dtype=UnivariateLoci.DTYPE_ULOCUS)
        npt.assert_array_equal(query.loci, answer)

    def test_from_iter(self):
        """
        Test for method from_iter.
        """
        loci = [(3, 6), (6, 8), (7, 9), (10, 12), (13, 13), (15, 25), (16, 17), (19, 20)]
        iterable = (locus for locus in loci)
        query = UnivariateLoci.from_iter(iterable)
        query.sort()
        answer = UnivariateLoci(np.array(loci, dtype=UnivariateLoci.DTYPE_ULOCUS))
        npt.assert_array_equal(query, answer)

    @pytest.mark.parametrize("loci,melted_loci",
                             # single locus spanning single base
                             [([(13, 13)],
                               [(13, 13)]),
                              # nested loci
                              ([(15, 25), (16, 17), (19, 20)],
                               [(15, 25)]),
                              # adjacent loci
                              ([(7, 9), (10, 12)],
                               [(7, 9), (10, 12)]),
                              # combined
                              ([(3, 6), (6, 8), (7, 9), (10, 12), (13, 13), (15, 25), (16, 17), (19, 20)],
                               [(3, 9), (10, 12), (13, 13), (15, 25)])],
                             ids=['single', 'nested', 'adjacent', 'combined'])
    def test_melt(self, loci, melted_loci):
        """
        Test for method melt.
        Method modifies loci in place.
        Test includes following edge cases:
         * Long locus completely overlaps short loci: (15, 25) & (16, 17) & (19, 20) --> (15, 25)
         * Adjacent loci do not get merged: (7, 9) & (10, 12) -->  (*, 9) & (10, *)
         * Locus may span a single base: (13, 13) --> (13, 13)
        """
        query = UnivariateLoci.from_iter(loci)
        query.melt()  # melt should automatically sort loci
        answer = UnivariateLoci.from_iter(melted_loci)
        answer.sort()
        npt.assert_array_equal(query.loci, answer)

    @pytest.mark.parametrize("loci,subset,locus,end",
                             # check if both 'start' and 'stop' ends are in locus
                             [([(3, 6), (6, 8), (7, 9), (10, 12), (13, 13), (15, 25), (16, 17), (19, 20)],
                               [(7, 9), (10, 12), (13, 13), (16, 17)],
                               (7, 17),
                               'both'),
                              # check if only 'start' end is in locus
                              ([(3, 6), (6, 8), (7, 9), (10, 12), (13, 13), (15, 25), (16, 17), (19, 20)],
                               [(7, 9), (10, 12), (13, 13), (15, 25), (16, 17)],
                               (7, 17),
                               'start')],
                             ids=['both', 'start'])
    def test_subset_by_locus(self, loci, subset, locus, end):
        """
        Test factory for method subset_by_locus.
        """
        query = UnivariateLoci.from_iter(loci)
        query.sort()
        subset = UnivariateLoci.from_iter(subset)
        subset.sort()
        npt.assert_array_equal(query.subset_by_locus(*locus, end=end), subset)

    def test_append(self):
        """
        Test for method append.
        """
        x = UnivariateLoci.from_iter([(3, 9), (10, 12), (13, 13), (15, 25)])
        y = UnivariateLoci.from_iter([(4, 11), (7, 22), (23, 33), (25, 35)])
        query = UnivariateLoci.append(x, y)
        answer = UnivariateLoci.from_iter([(3, 9), (4, 11), (7, 22), (10, 12), (13, 13), (15, 25), (23, 33), (25, 35)])
        npt.assert_array_equal(query, answer)
