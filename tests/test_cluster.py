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
    """
    Tests for class FlatUnivariateDensityCluster.
    """
    def test_flat_subcluster(self):
        """
        Test for hidden method _flat_subcluster.

        """
        answer = np.array([(1, 2), (2, 5), (5, 6), (5, 7)], dtype=FlatUnivariateDensityCluster._ulocus)
        query = FlatUnivariateDensityCluster(4, 3)
        query.points = np.array([5, 2, 1, 5, 15, 1, 1, 6, 13, 5, 7, 14], dtype=int)
        npt.assert_array_equal(query._flat_subcluster(query.points, query.min_pts, query.eps),
                               answer)

    def test_flat_cluster(self):
        """
        Test for hidden method _flat_cluster.
        Most edge cases should be caught in tests for component methods.
        """
        answer = np.array([(1, 2), (6, 14)], dtype=FlatUnivariateDensityCluster._ulocus)
        query = FlatUnivariateDensityCluster(5, 4)
        query.points = np.array([9, 2, 2, 1, 7, 19, 1, 1, 6, 13, 10, 7, 14, 11, 11], dtype=int)
        npt.assert_array_equal(query._flat_cluster(query.points, query.min_pts, query.eps),
                               answer)

    def test_integration_fit(self):
        """
        Integration test to run class via the method fit.
        Points passed to the fit method should be copied into new array.
        New copy of points should be sorted.
        Most edge cases should be caught in tests for component methods.
        """
        query = FlatUnivariateDensityCluster(5, 4)
        input_points = np.array([9, 2, 2, 1, 7, 19, 1, 1, 6, 13, 10, 7, 14, 11, 11], dtype=int)
        answer_points = np.array([1, 1, 1, 2, 2, 6, 7, 7, 9, 10, 11, 11, 13, 14, 19], dtype=int)
        answer_loci = np.array([(1, 2), (6, 14)], dtype=FlatUnivariateDensityCluster._ulocus)
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
        query = HierarchicalUnivariateDensityCluster(3, 5, 3)  # Minimum eps of 3
        query.points = np.array([1, 1, 2, 2, 2, 3, 3, 5, 7, 7, 8, 8, 9, 9, 9, 11])  # Point 5 within 3 of 3 and 7
        answer = {'area': 48,
                  'base_eps': 5,
                  'base_locus': (1, 11),
                  'child_area': 0,
                  'children': None,
                  'selected': False}
        assert query._grow_tree(query.points, query.min_pts, query.max_eps, query.min_eps) == answer

    def test_grow_tree_nested(self):
        """
        Test for hidden method _grow_tree using data with nested densities clusters.
        Produce a tree with three nodes (root node and two child nodes).
        Minimum eps of 2 allows cluster to be split when points are not adjacent.
        """
        query = HierarchicalUnivariateDensityCluster(3, 5, 2)  # Minimum eps of 2
        query.points = np.array([1, 1, 2, 2, 2, 3, 3, 5, 7, 7, 8, 8, 9, 9, 9, 11])
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
        assert query._grow_tree(query.points, query.min_pts, query.max_eps, query.min_eps) == answer

    def test_child_area(self):
        """
        Test for hidden method _child_area.
        Method modifies tree/node 'child_area' value in place.
        """
        hudc = HierarchicalUnivariateDensityCluster(3, 5, 2)
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
        hudc._child_area(query)
        assert query == answer

    def test_select_nodes(self):
        """
        Test for hidden method _select_nodes.
        Method modifies tree/node 'selected' value in place.
        """
        hudc = HierarchicalUnivariateDensityCluster(3, 10, 2)
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
        hudc._select_nodes(query)
        assert query == answer

    def test_retrieve_selected_loci_nested(self):
        """
        Test for hidden method _retrieve_selected_loci using a nested data set.
        """
        hudc = HierarchicalUnivariateDensityCluster(3, 10, 2)
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
        assert hudc._retrieve_selected_loci(query) == answer

    def test_retrieve_selected_loci_simple(self):
        """
        Test for hidden method _retrieve_selected_loci using a non-nested data set.
        """
        hudc = HierarchicalUnivariateDensityCluster(3, 10, 2)
        query = {'area': 19,
                 'base_eps': 10,
                 'base_locus': (1, 11),
                 'child_area': 20,
                 'children': None,
                 'selected': True}
        answer = (1, 11)
        assert hudc._retrieve_selected_loci(query) == answer

    def test_flatten_nested_list(self):
        """
        Test for hidden method _flatten using a nested data set.
        Method _flatten returns a generator which should be coerced to a list.
        """
        hudc = HierarchicalUnivariateDensityCluster(3, 10, 2)
        query = [(1, 3), [(6, 9), (11, 11)]]
        answer = [(1, 3), (6, 9), (11, 11)]
        assert list(hudc._flatten(query)) == answer

    def test_flatten_tuple(self):
        """
        Test for hidden method _flatten using a tuple.
        Method _flatten returns a generator which should be coerced to a list.
        """
        hudc = HierarchicalUnivariateDensityCluster(3, 10, 2)
        query = (1, 11)
        answer = [(1, 11)]
        assert list(hudc._flatten(query)) == answer

    def test_single_hierarchical_cluster(self):
        """
        Test for hidden method _single_hierarchical_cluster.
        """
        hudc = HierarchicalUnivariateDensityCluster(3, 10, 2)
        query = np.array([1, 1, 2, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 13, 21, 21, 21, 22, 22,22, 23, 31])
        answer = np.array([(1, 4), (13, 23)], dtype=_UnivariateLoci._ulocus)
        npt.assert_array_equal(hudc._single_hierarchical_cluster(query, 3, 10, 2), answer)
