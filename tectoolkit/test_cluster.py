#! /usr/bin/env python

import numpy as np
import numpy.testing as npt
from tectoolkit.cluster import _UnivariateLoci, FlatUnivariateDensityCluster, HierarchicalUnivariateDensityCluster


class TestUL:
    """
    Tests for hidden class _UnivariateLoci.
    """
    def test_melt_uloci(self):
        """
        Test for hidden method _melt_uloci.
        Test includes following edge cases:
         * Long locus completely overlaps short locus: (15, 25) & (16, 17) --> (15, 25)
         * Adjacent loci do not get merged: (7, 9) & (10, 12) -->  (*, 9) & (10, *)
         * Locus may span a single base: (13, 13) --> (13, 13)
        """
        query = np.array([(3, 6),
                          (6, 8),
                          (7, 9),
                          (10, 12),
                          (13, 13),
                          (15, 25),
                          (16, 17),
                          (19, 20)], dtype=_UnivariateLoci._ulocus)
        answer = np.array([(3, 9),
                           (10, 12),
                           (13, 13),
                           (15, 25)], dtype=_UnivariateLoci._ulocus)
        query = _UnivariateLoci._melt_uloci(query)
        npt.assert_array_equal(query, answer)

    def test_sort_uloci(self):
        """
        Test for hidden method _sort_uloci.
        By default, loci should be sorted by lower bound then upper bound.
        """
        query = _UnivariateLoci()
        query.loci = np.array([(2, 4),
                               (3, 4),
                               (3, 3),
                               (4, 4),
                               (3, 99),
                               (1, 1)], dtype=_UnivariateLoci._ulocus)
        answer = _UnivariateLoci()
        answer.loci = np.array([(1, 1),
                                (2, 4),
                                (3, 3),
                                (3, 4),
                                (3, 99),
                                (4, 4)], dtype=_UnivariateLoci._ulocus)
        query._sort_uloci()
        npt.assert_array_equal(query.loci, answer.loci)

    def test_locus_points(self):
        """
        Test for hidden method _locus_points.
        Points inside the inclusive boundary of the locus should be returned.
        """
        query = np.array([5, 9, 4, 1, 6, 8, 6, 2], dtype=int)
        locus = (5, 8)
        answer = np.array([5, 6, 8, 6], dtype=int)
        query = _UnivariateLoci._locus_points(locus, query)
        npt.assert_array_equal(query, answer)


class TestFUDC:
    """"""
    def test_flat_subcluster(self):
        """"""
        answer = np.array([(1, 2), (2, 5), (5, 6), (5, 7)], dtype=FlatUnivariateDensityCluster._ulocus)
        query = FlatUnivariateDensityCluster(4, 3)
        query.points = np.array([5, 2, 1, 5, 15, 1, 1, 6, 13, 5, 7, 14], dtype=int)
        npt.assert_array_equal(query._flat_subcluster(query.points, query.min_pts, query.eps),
                               answer)

    def test_flat_cluster(self):
        """"""
        answer = np.array([(1, 2), (6, 14)], dtype=FlatUnivariateDensityCluster._ulocus)
        query = query = FlatUnivariateDensityCluster(5, 4)
        query.points = np.array([9, 2, 2, 1, 7, 19, 1, 1, 6, 13, 10, 7, 14, 11, 11, ], dtype=int)
        npt.assert_array_equal(query._flat_cluster(query.points, query.min_pts, query.eps),
                               answer)