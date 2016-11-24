#! /usr/bin/env python

import numpy as np
import numpy.testing as npt
from tectoolkit.cluster import _UnivariateLoci, FUDC, HUDC


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
                          (19, 20)], dtype=_UnivariateLoci.DTYPE_ULOCUS)
        answer = np.array([(3, 9),
                           (10, 12),
                           (13, 13),
                           (15, 25)], dtype=_UnivariateLoci.DTYPE_ULOCUS)
        npt.assert_array_equal(_UnivariateLoci._melt_uloci(query), answer)

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
                               (1, 1)], dtype=_UnivariateLoci.DTYPE_ULOCUS)
        answer = _UnivariateLoci()
        answer.loci = np.array([(1, 1),
                                (2, 4),
                                (3, 3),
                                (3, 4),
                                (3, 99),
                                (4, 4)], dtype=_UnivariateLoci.DTYPE_ULOCUS)
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
        npt.assert_array_equal(_UnivariateLoci._locus_points(locus, query), answer)


class TestFUDC:
    """
    Tests for class FlatUnivariateDensityCluster.
    """
    def test_flat_subcluster(self):
        """
        Test for hidden method _flat_subcluster.

        """
        query = np.array([5, 2, 1, 5, 15, 1, 1, 6, 13, 5, 7, 14], dtype=int)
        answer = np.array([(1, 2), (2, 5), (5, 6), (5, 7)], dtype=FUDC.DTYPE_ULOCUS)
        npt.assert_array_equal(HUDC._flat_subcluster(query, 4, 3), answer)

    def test_flat_cluster(self):
        """
        Test for hidden method _flat_cluster.
        Most edge cases should be caught in tests for component methods.
        """
        query = np.array([9, 2, 2, 1, 7, 19, 1, 1, 6, 13, 10, 7, 14, 11, 11], dtype=int)
        answer = np.array([(1, 2), (6, 14)], dtype=FUDC.DTYPE_ULOCUS)
        npt.assert_array_equal(FUDC.flat_cluster(query, 5, 4), answer)

    def test_fit(self):
        """
        Test to run class via the method fit.
        Points passed to the fit method should be copied into new array.
        New copy of points should be sorted.
        Most edge cases should be caught in tests for component methods.
        """
        query = FUDC(5, 4)
        input_points = np.array([9, 2, 2, 1, 7, 19, 1, 1, 6, 13, 10, 7, 14, 11, 11], dtype=int)
        answer_points = np.array([1, 1, 1, 2, 2, 6, 7, 7, 9, 10, 11, 11, 13, 14, 19], dtype=int)
        answer_loci = np.array([(1, 2), (6, 14)], dtype=FUDC.DTYPE_ULOCUS)
        query.fit(input_points)
        assert query.points is not input_points
        npt.assert_array_equal(query.points, answer_points)
        npt.assert_array_equal(query.loci, answer_loci)


class TestHUDC:
    """
    Tests for class HierarchicalUnivariateDensityCluster.
    """
    def test_grow_tree_single(self):
        """
        Test for hidden method _grow_tree using data with simple densities cluster.
        Produce a tree with one (root) node.
        Minimum eps of 3 means point 5 is still within eps of 3 and 7
        """
        query = np.array([1, 1, 2, 2, 2, 3, 3, 5, 7, 7, 8, 8, 9, 9, 9, 11])  # Point 5 within 3 of 3 and 7
        answer = {'area': 48,
                  'base_eps': 5,
                  'base_locus': (1, 11),
                  'child_area': 0,
                  'children': None,
                  'selected': False}
        assert HUDC._grow_tree(query, 3, 5, 3) == answer  # Minimum eps of 3

    def test_grow_tree_nested(self):
        """
        Test for hidden method _grow_tree using data with nested densities clusters.
        Produce a tree with three nodes (root node and two child nodes).
        Minimum eps of 2 allows cluster to be split when points are not adjacent.
        """
        query = np.array([1, 1, 2, 2, 2, 3, 3, 5, 7, 7, 8, 8, 9, 9, 9, 11])
        answer = {'area': 64,
                  'base_eps': 5,
                  'base_locus': (1, 11),
                  'child_area': 0,
                  'children': [{'area': 7,
                                'base_eps': 1,
                                'base_locus': (1, 3),
                                'child_area': 0,
                                'children': None,
                                'selected': False},
                               {'area': 7,
                                'base_eps': 1,
                                'base_locus': (7, 9),
                                'child_area': 0,
                                'children': None,
                                'selected': False}],
                  'selected': False}
        assert HUDC._grow_tree(query, 3, 5, 2) == answer  # Minimum eps of 2

    def test_child_area(self):
        """
        Test for hidden method _child_area.
        Method modifies tree/node 'child_area' value in place.
        """
        query = {'area': 42,
                 'base_eps': 5,
                 'base_locus': (1, 11),
                 'child_area': 0,
                 'children': [{'area': 7,
                               'base_eps': 2,
                               'base_locus': (1, 3),
                               'child_area': 0,
                               'children': None,
                               'selected': False},
                              {'area': 7,
                               'base_eps': 2,
                               'base_locus': (6, 11),
                               'child_area': 0,
                               'children': [{'area': 3,
                                             'base_eps': 1,
                                             'base_locus': (6, 7),
                                             'child_area': 0,
                                             'children': [],
                                             'selected': False},
                                            {'area': 3,
                                             'base_eps': 1,
                                             'base_locus': (9, 9),
                                             'child_area': 0,
                                             'children': None,
                                             'selected': False}],
                               'selected': False}],
                 'selected': False}
        answer = {'area': 42,
                  'base_eps': 5,
                  'base_locus': (1, 11),
                  'child_area': 20,
                  'children': [{'area': 7,
                                'base_eps': 2,
                                'base_locus': (1, 3),
                                'child_area': 0,
                                'children': None,
                                'selected': False},
                               {'area': 7,
                                'base_eps': 2,
                                'base_locus': (6, 11),
                                'child_area': 6,
                                'children': [{'area': 3,
                                              'base_eps': 1,
                                              'base_locus': (6, 7),
                                              'child_area': 0,
                                              'children': [],
                                              'selected': False},
                                             {'area': 3,
                                              'base_eps': 1,
                                              'base_locus': (9, 9),
                                              'child_area': 0,
                                              'children': None,
                                              'selected': False}],
                                'selected': False}],
                  'selected': False}
        HUDC._child_area(query)
        assert query == answer

    def test_select_nodes(self):
        """
        Test for hidden method _select_nodes.
        Method modifies tree/node 'selected' value in place.
        """
        query = {'area': 19,
                 'base_eps': 10,
                 'base_locus': (1, 11),
                 'child_area': 20,
                 'children': [{'area': 5,
                               'base_eps': 5,
                               'base_locus': (1, 3),
                               'child_area': 0,
                               'children': None,
                               'selected': False},
                              {'area': 9,
                               'base_eps': 5,
                               'base_locus': (6, 11),
                               'child_area': 6,
                               'children': [{'area': 3,
                                             'base_eps': 2,
                                             'base_locus': (6, 7),
                                             'child_area': 0,
                                             'children': [],
                                             'selected': False},
                                            {'area': 3,
                                             'base_eps': 2,
                                             'base_locus': (9, 9),
                                             'child_area': 0,
                                             'children': None,
                                             'selected': False}],
                               'selected': False}],
                 'selected': False}
        answer = {'area': 19,
                  'base_eps': 10,
                  'base_locus': (1, 11),
                  'child_area': 20,
                  'children': [{'area': 5,
                                'base_eps': 5,
                                'base_locus': (1, 3),
                                'child_area': 0,
                                'children': None,
                                'selected': True},
                               {'area': 9,
                                'base_eps': 5,
                                'base_locus': (6, 11),
                                'child_area': 6,
                                'children': [{'area': 3,
                                              'base_eps': 2,
                                              'base_locus': (6, 7),
                                              'child_area': 0,
                                              'children': [],
                                              'selected': False},
                                             {'area': 3,
                                              'base_eps': 2,
                                              'base_locus': (9, 9),
                                              'child_area': 0,
                                              'children': None,
                                              'selected': False}],
                                'selected': True}],
                  'selected': False}
        HUDC._select_nodes(query)
        assert query == answer

    def test_retrieve_selected_loci_nested(self):
        """
        Test for hidden method _retrieve_selected_loci using a nested data set.
        """
        query = {'area': 19,
                 'base_eps': 10,
                 'base_locus': (1, 11),
                 'child_area': 20,
                 'children': [{'area': 5,
                               'base_eps': 5,
                               'base_locus': (1, 3),
                               'child_area': 0,
                               'children': None,
                               'selected': True},
                              {'area': 9,
                               'base_eps': 5,
                               'base_locus': (6, 11),
                               'child_area': 6,
                               'children': [{'area': 3,
                                             'base_eps': 2,
                                             'base_locus': (6, 9),
                                             'child_area': 0,
                                             'children': [{'area': 5,
                                                           'base_eps': 5,
                                                           'base_locus': (9, 9),
                                                           'child_area': 0,
                                                           'children': None,
                                                           'selected': False}],
                                             'selected': True},
                                            {'area': 3,
                                             'base_eps': 2,
                                             'base_locus': (11, 11),
                                             'child_area': 0,
                                             'children': None,
                                             'selected': True}],
                               'selected': False}],
                 'selected': False}
        answer = [(1, 3), [(6, 9), (11, 11)]]
        assert HUDC._retrieve_selected_loci(query) == answer

    def test_retrieve_selected_loci_simple(self):
        """
        Test for hidden method _retrieve_selected_loci using a non-nested data set.
        """
        query = {'area': 19,
                 'base_eps': 10,
                 'base_locus': (1, 11),
                 'child_area': 20,
                 'children': None,
                 'selected': True}
        answer = (1, 11)
        assert HUDC._retrieve_selected_loci(query) == answer

    def test_flatten_list_nested(self):
        """
        Test for hidden method _flatten_list using a nested data set.
        Method _flatten_list returns a generator which should be coerced to a list.
        """
        query = [(1, 3), [(6, 9), (11, 11)]]
        answer = [(1, 3), (6, 9), (11, 11)]
        assert list(HUDC._flatten_list(query)) == answer

    def test_flatten_list_tuple(self):
        """
        Test for hidden method _flatten_list using a tuple.
        Method _flatten_list returns a generator which should be coerced to a list.
        """
        query = (1, 11)
        answer = [(1, 11)]
        assert list(HUDC._flatten_list(query)) == answer

    def test_single_hierarchical_cluster(self):
        """
        Test for hidden method _single_hierarchical_cluster.
        """
        query = np.array([1, 1, 2, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 13, 21, 21, 21, 22, 22, 22, 23, 31], dtype=int)
        answer = np.array([(1, 4), (13, 23)], dtype=_UnivariateLoci.DTYPE_ULOCUS)
        npt.assert_array_equal(HUDC._single_hierarchical_cluster(query, 3, 10, 2), answer)

    def test_hierarchical_cluster(self):
        """
        Test for hidden method _hierarchical_cluster.
        """
        query = np.array([1, 2, 21, 22, 22, 22, 24, 38, 54, 54, 55, 56, 65, 65, 66, 67, 68, 90], dtype=int)
        answer = np.array([(21, 24), (54, 56), (65, 68)], dtype=_UnivariateLoci.DTYPE_ULOCUS)
        npt.assert_array_equal(HUDC._single_hierarchical_cluster(query, 3, 10, 2), answer)

    def test_fit(self):
        """
        Test to run class via the method fit.
        Points passed to the fit method should be copied into new array.
        New copy of points should be sorted.
        Most edge cases should be caught in tests for component methods.
        """
        query = HUDC(3, 10, 2)
        input_points = np.array([22, 54, 24, 22, 2, 21, 54, 22, 90, 38, 65, 67, 68, 56, 55, 65, 66, 1], dtype=int)
        answer_points = np.array([1, 2, 21, 22, 22, 22, 24, 38, 54, 54, 55, 56, 65, 65, 66, 67, 68, 90], dtype=int)
        answer_loci = np.array([(21, 24), (54, 56), (65, 68)], dtype=FUDC.DTYPE_ULOCUS)
        query.fit(input_points)
        assert query.points is not input_points
        npt.assert_array_equal(query.points, answer_points)
        npt.assert_array_equal(query.loci, answer_loci)
