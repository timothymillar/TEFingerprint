#! /usr/bin/env python3

import pytest
import numpy as np
import numpy.testing as npt
from tefingerprint.cluster import DBICAN, SDBICAN

class TestDBICAN:
    """
    Tests for class DBICAN.
    """
    @pytest.mark.parametrize("array,answer",
                             [(np.array([0, 0, 0, 3, 4, 7, 99], dtype=int),
                               True),
                              (np.array([99, 55, 55, 55, 2, 0], dtype=int),
                               False),
                              (np.array([0, 0, 4, 3, 4, 7, 99], dtype=int),
                               False)],
                             ids=['ascending', 'non-ascending', 'descending'])
    def test_sorted_ascending(self, array, answer):
        """
        Test for hidden method _sorted_ascending.
        """
        assert DBICAN._sorted_ascending(array) == answer

    @pytest.mark.parametrize("slices,melted_slices",
                             # single slice spanning single base
                             [(np.array([(13, 14)],
                                        dtype=DBICAN._DTYPE_SLICE),
                               np.array([(13, 14)],
                                        dtype=DBICAN._DTYPE_SLICE)),
                              # nested slices (this shouldn't happen in practice but is handled correctly)
                              (np.array([(15, 25), (16, 17), (19, 20)],
                                        dtype=DBICAN._DTYPE_SLICE),
                               np.array([(15, 25)],
                                        dtype=DBICAN._DTYPE_SLICE)),
                              # adjacent loci (slices are half open intervals)
                              (np.array([(7, 9), (9, 12)],
                                        dtype=DBICAN._DTYPE_SLICE),
                               np.array([(7, 9), (9, 12)],
                                        dtype=DBICAN._DTYPE_SLICE)),
                              # combined
                              (np.array([(3, 6), (6, 8), (7, 9), (10, 12), (12, 13), (15, 25), (16, 17), (19, 20)],
                                        dtype=DBICAN._DTYPE_SLICE),
                               np.array([(3, 6), (6, 9), (10, 12), (12, 13), (15, 25)],
                                        dtype=DBICAN._DTYPE_SLICE))
                              ],
                             ids=['single', 'nested', 'adjacent', 'combined'])
    def test_melt_slices(self, slices, melted_slices):
        """
        Test for hidden method _melt_slices.
        Test includes following edge cases:
         * Long slice completely overlaps short loci: (15, 25) & (16, 17) & (19, 20) --> (15, 25)
         * Adjacent slices do not get merged: (7, 9) & (9, 12) -->  (*, 9) & (9, *)
         * Slice may span a single value: (13, 14) --> (13, 14)
        """
        npt.assert_array_equal(DBICAN._melt_slices(slices), melted_slices)

    def test_subcluster(self):
        """
        Test for hidden method _subcluster.
        """
        array = np.array([1, 2, 21, 22, 22, 22, 24, 38, 54, 54, 55, 56, 65, 65, 66, 67, 68, 90],
                         dtype=int)
        slices = np.array([(2, 6), (3, 7), (8, 12), (12, 16), (13, 17)],
                          dtype=DBICAN._DTYPE_SLICE)
        npt.assert_array_equal(DBICAN._subcluster(array, 4, 5), slices)

    def test_flat_cluster(self):
        """
        Test for hidden method _cluster.
        Most edge cases should be caught in tests for component methods.
        """
        sub_slices = np.array([1, 2, 21, 22, 22, 22, 24, 38, 54, 54, 55, 56, 65, 65, 66, 67, 68, 90],
                              dtype=int)
        slices = np.array([(2, 7), (8, 12), (12, 17)],
                          dtype=DBICAN._DTYPE_SLICE)
        npt.assert_array_equal(DBICAN._cluster(sub_slices, 4, 5), slices)

    @pytest.mark.parametrize("array,slices",
                             [(np.array([], dtype=int),
                               np.array([], dtype=DBICAN._DTYPE_SLICE)),
                              (np.array([1], dtype=int),
                               np.array([], dtype=DBICAN._DTYPE_SLICE)),
                              (np.array([1, 2, 3, 4], dtype=int),
                               np.array([], dtype=DBICAN._DTYPE_SLICE)),
                              (np.array([1, 2, 3, 4, 5], dtype=int),
                               np.array([(0, 5)], dtype=DBICAN._DTYPE_SLICE))
                              ],
                             ids=['0<n', '1<n', '4<n', '5==n'])
    def test_DBICAN_few_reads(self, array, slices):
        """
        Test for method DBICAN with small arrays.
        Method udc should correctly handle an array of length 0 or greater.
        If the length of the array is less than 'n' then an empty array will
        always be returned.
        """
        npt.assert_array_equal(DBICAN.DBICAN(array, 5, 5), slices)

    def test_fit(self):
        """
        Test to run class via the method fit.
        Points passed to the fit method should be copied into new array.
        New copy of points should be sorted.
        Most edge cases should be caught in tests for component methods.
        """
        model = DBICAN(4, 5)
        input_array = np.array([1, 2, 21, 22, 22, 22, 24, 38, 54, 54, 55, 56, 65, 65, 66, 67, 68, 90],
                               dtype=int)
        answer_array = np.array([1, 2, 21, 22, 22, 22, 24, 38, 54, 54, 55, 56, 65, 65, 66, 67, 68, 90],
                                dtype=int)
        answer_slices = np.array([(2, 7), (8, 12), (12, 17)],
                                 dtype=DBICAN._DTYPE_SLICE)
        model.fit(input_array)
        assert model.input_array is not input_array
        npt.assert_array_equal(model.input_array, answer_array)
        npt.assert_array_equal(model.slices, answer_slices)

    def test_clusters(self):
        """
        Test for method clusters.
        """
        model = DBICAN(4, 5)
        input_array = np.array([1, 2, 21, 22, 22, 22, 24, 38, 54, 54, 55, 56, 65, 65, 66, 67, 68, 90], dtype=int)
        cluster_arrays = [np.array([21, 22, 22, 22, 24], dtype=int),
                          np.array([54, 54, 55, 56], dtype=int),
                          np.array([65, 65, 66, 67, 68], dtype=int)]
        model.fit(input_array)
        for result, answer in zip(model.clusters(), cluster_arrays):
            npt.assert_array_equal(result, answer)

    def test_cluster_extremities(self):
        """Test for method cluster_extremities"""
        model = DBICAN(4, 5)
        input_array = np.array([1, 2, 21, 22, 22, 22, 24, 38, 54, 54, 55, 56, 65, 65, 66, 67, 68, 90], dtype=int)
        cluster_extremities = [(21, 24), (54, 56), (65, 68)]
        model.fit(input_array)
        for result, answer in zip(model.cluster_extremities(), cluster_extremities):
            assert result == answer

    def test_labels(self):
        """Test for method labels"""
        udc_object = DBICAN(4, 5)
        input_array = np.array([1, 2, 21, 22, 22, 22, 24, 38, 54, 54, 55, 56, 65, 65, 66, 67, 68, 90], dtype=int)
        labels = np.array([-1, -1, 0, 0, 0, 0, 0, -1, 1, 1, 1, 1, 2, 2, 2, 2, 2, -1], dtype=int)
        udc_object.fit(input_array)
        npt.assert_array_equal(udc_object.labels(), labels)


class TestSDBICAN:
    """
    Tests for class SDBICAN.
    """
    def test_core_distances(self):
        """
        Test for hidden method _core_distances.
        Data is derived from real edge case.
        """
        input_array = np.array([   0,    0,   60,   61,   61,   61,   76,   78,  122,  122,  141,
                                 183,  251,  260,  260,  263,  263,  267,  267,  288,  288,  295,
                                 300,  310,  310,  317,  317,  334,  334,  335,  338,  338,  338,
                                 338,  340,  342,  342,  344,  344,  358,  367,  370,  370,  377,
                                 387,  402,  403,  410,  410,  410,  418,  418,  424,  424,  577,
                                 857,  879,  921,  921, 1007, 1031, 1051, 1051, 1059, 1071, 1071,
                                1080, 1094, 1094, 1110, 1110, 1113, 1113, 1183, 1189, 1200, 1200,
                                1217, 1234, 1234, 1591, 1620, 1620, 1662, 1686, 1707, 1755, 1828,
                                1828, 1848, 1848, 1848, 1848, 1851, 1851, 1852, 1917], dtype=int)
        point_eps = np.array([122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 123, 105,  44,
                               40,  40,  40,  40,  40,  40,  40,  40,  40,  38,  28,  28,  23,
                               23,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   9,   9,
                               20,  29,  32,  32,  37,  37,  37,  37,  37,  37,  37,  37,  37,
                               37,  37, 175, 214, 192, 159, 159,  87,  79,  59,  59,  54,  54,
                               54,  54,  54,  54,  54,  54,  54,  54, 106, 106, 106, 106, 123,
                              124, 124, 257, 228, 228, 186, 165, 144,  97,  89,  89,  89,  89,
                               89,  89,  89,  89,  89,  89], dtype=int)
        npt.assert_array_equal(SDBICAN._core_distances(input_array, 10), point_eps)

    def test_point_eps_hidden_peak(self):
        """
        Test for hidden method _point_eps.
        Data contains a single "hidden" eps peak between values 6 and 26
        (values themselves have low eps but a large gap between them).
        """
        input_array = np.array([0, 0, 3, 4, 4, 6, 26, 28, 28, 29, 32, 32], dtype=int)
        point_eps = np.array([3, 3, 1, 1, 1, 2, 2, 1, 1, 1, 3, 3], dtype=int)
        npt.assert_array_equal(SDBICAN._core_distances(input_array, 3), point_eps)

    @pytest.mark.parametrize("array,min_points,answer",
                             [(np.array([], dtype=int), 3, None),
                              (np.array([1], dtype=int), 3, None),
                              (np.array([0, 4, 6, 7, 7, 7, 7, 7, 8, 10, 20], dtype=int), 3, None),
                              (np.array([0, 0, 3, 4, 4, 6, 26, 28, 28, 29, 32, 32], dtype=int), 3, 22),
                              (np.array([   0,    0,   60,   61,   61,   61,   76,   78,  122,  122,  141,
                                          183,  251,  260,  260,  263,  263,  267,  267,  288,  288,  295,
                                          300,  310,  310,  317,  317,  334,  334,  335,  338,  338,  338,
                                          338,  340,  342,  342,  344,  344,  358,  367,  370,  370,  377,
                                          387,  402,  403,  410,  410,  410,  418,  418,  424,  424,  577,
                                          857,  879,  921,  921, 1007, 1031, 1051, 1051, 1059, 1071, 1071,
                                         1080, 1094, 1094, 1110, 1110, 1113, 1113, 1183, 1189, 1200, 1200,
                                         1217, 1234, 1234, 1591, 1620, 1620, 1662, 1686, 1707, 1755, 1828,
                                         1828, 1848, 1848, 1848, 1848, 1851, 1851, 1852, 1917], dtype=int), 10, 454)],
                             ids=['empty', 'short', 'none', 'hidden', 'real'])
    def test_fork_epsilon(self, array, min_points, answer):
        """
        Test for hidden method _fork_epsilon.
        """
        assert SDBICAN._fork_epsilon(array, min_points) == answer

    @pytest.mark.parametrize("nested,flat",
                             [([(1, 3), [(6, 9), (11, 11)]],
                               [(1, 3), (6, 9), (11, 11)]),
                              ([(1, 3), (6, 9), (11, 11)],
                               [(1, 3), (6, 9), (11, 11)]),
                              ((1, 11),
                               [(1, 11)])],
                             ids=['nested', 'flat', 'tuple'])
    def test_flatten_list(self, nested, flat):
        """
        Tests for hidden method _flatten_list.
        Method _flatten_list returns a generator which should be coerced
        to a list.
        """
        assert list(SDBICAN._flatten_list(nested)) == flat

    @pytest.mark.parametrize("array,min_points,max_eps,answer",
                             [(np.array([], dtype=int), 5, 5, np.array([], dtype=DBICAN._DTYPE_SLICE)),
                              (np.array([1], dtype=int), 5, 5, np.array([], dtype=DBICAN._DTYPE_SLICE)),
                              (np.array([1, 2, 3, 4], dtype=int), 5, 5, np.array([], dtype=DBICAN._DTYPE_SLICE)),
                              (np.array([1, 2, 3, 4, 5], dtype=int), 5, 5, np.array([(0, 5)], dtype=DBICAN._DTYPE_SLICE)),
                              (np.array([0, 2, 3, 7, 8, 12, 13, 17], dtype=int), 3, 5, np.array([(0, 8)], dtype=DBICAN._DTYPE_SLICE)),
                              (np.array([0, 2, 3, 7, 8, 12, 13, 16], dtype=int), 3, 5, np.array([(0, 3), (5, 8)], dtype=DBICAN._DTYPE_SLICE))],
                             ids=['empty', 'single', 'p < min-points', 'p == min-points', 'parent', 'children'])
    def test_sDBICAN(self, array, min_points, max_eps, answer):
        """
        Test for method sDBICAN with small arrays.
        Method sDBICAN should correctly handle an array of length 0 or greater.
        If the length of the array is less than 'n' then an empty array will
        always be returned.
        """
        npt.assert_array_equal(SDBICAN.sDBICAN(array, min_points, max_eps), answer)

    def test_sDBICAN_variations(self):
        """
        Test for 'aggressive' and 'conservative' variations

        In this data set with epsilon = 6 the same clusters are detected but
        support of the first
        child cluster [2, 5, 6, 7, 10, 11, 13] is 7 with the aggressive
        variation and 14
        with the conservative variation. The support of it's children is 9 in
        both variations so
        it is selected by the conservative variation but not by the aggressive
        variation.
        """
        array = np.array([2, 5, 6, 7, 10, 11, 13, 16, 20, 21, 22], dtype=int)
        answer_aggressive = np.fromiter([(1, 4), (4, 7), (8, 11)],
                                        dtype=SDBICAN._DTYPE_SLICE)
        answer_conservative = np.fromiter([(0, 7), (8, 11)],
                                          dtype=SDBICAN._DTYPE_SLICE)

        npt.assert_array_equal(SDBICAN.sDBICAN(array, 3, 6,
                                               aggressive_method=True),
                               answer_aggressive)
        npt.assert_array_equal(SDBICAN.sDBICAN(array, 3, 6,
                                               aggressive_method=False),
                               answer_conservative)

    def test_sDBICAN_initial_parent_support_calculation(self):
        """
        Test for method sDBICAN.
        Tests that the support initial clusters is calculated correctly.

        Cluster support is a measure of how far epsilon can be reduced while
        still forming a cluster.
        More specifically the sum of points retained within the cluster at
        each value of epsilon
        bellow the maximum value of epsilon of that cluster.
        In the case of top level clusters (i.e. those found at the initial
        values of epsilon: max_eps)
        that do not exist at max_eps - 1 their calculates support should be 0.
        They are still a valid cluster (and would be found with
        non-hierarchical clustering) but if they
        have child clusters then the child clusters are always selected.
        In this test when (max) epsilon = 3 the parent cluster is detected by
        the non-hierarchical
        version and is also detected by the hierarchical version but because
        it has 0 support it
        is discarded in favour of the child clusters.
        """

        array = np.array([1, 3, 4, 5, 7, 8, 9, 11], dtype=int)
        answer_1 = np.fromiter([(0, 8)], dtype=SDBICAN._DTYPE_SLICE)
        answer_2 = np.fromiter([(1, 4), (4, 7)], dtype=SDBICAN._DTYPE_SLICE)

        npt.assert_array_equal(DBICAN.DBICAN(array, 3, 4), answer_1)
        npt.assert_array_equal(DBICAN.DBICAN(array, 3, 3), answer_1)
        npt.assert_array_equal(DBICAN.DBICAN(array, 3, 2), answer_2)

        npt.assert_array_equal(SDBICAN.sDBICAN(array, 3, 4), answer_1)
        npt.assert_array_equal(SDBICAN.sDBICAN(array, 3, 3), answer_2)
        npt.assert_array_equal(SDBICAN.sDBICAN(array, 3, 2), answer_2)

    def test_sDBICAN_child_support_calculation(self):
        """
        Test for method sDBICAN.
        Tests that the support for child clusters is calculated correctly.

        If child support is incorrectly calculated as the
        sum(epsilon_minimum - core_distances) then core_distances
        larger than epsilon_minimum (of the parent cluster) lead to a
        negative value which lowers the calculated
        child support and results in the selection of the parent
        cluster (0, 12).
        If child support is correctly calculated as the
        sum(maximum(0, (epsilon_minimum - core_distances))) then
        points with core_distance larger than epsilon_minimum
        (i.e. those not present in the child clusters) don't
        affect calculated support and the child clusters are selected.
        """
        array = np.array([0, 1, 2, 7, 8, 13, 14, 16, 21, 23, 28, 30], dtype=int)
        answer = np.fromiter([(0, 3), (5, 8)], dtype=SDBICAN._DTYPE_SLICE)
        npt.assert_array_equal(SDBICAN.sDBICAN(array, 3), answer)

    def test_sDBICAN_default_max_eps(self):
        """
        Test for method sDBICAN.
        Tests that the correct value of max_eps is selected if none
        is specified.
        By default it should exclude the root node of the cluster tree.

        """
        array = np.array([0, 1, 2, 3, 4, 5, 7, 9, 10, 11, 12, 13, 14], dtype=int)
        answer = np.fromiter([(0, 6), (7, 13)], dtype=SDBICAN._DTYPE_SLICE)
        npt.assert_array_equal(SDBICAN.sDBICAN(array, 3), answer)

    def test_sDBICAN_real(self):
        """Test for method sDBICAN"""
        input_array = np.array([   0,    0,   60,   61,   61,   61,   76,   78,  122,  122,  141,
                                 183,  251,  260,  260,  263,  263,  267,  267,  288,  288,  295,
                                 300,  310,  310,  317,  317,  334,  334,  335,  338,  338,  338,
                                 338,  340,  342,  342,  344,  344,  358,  367,  370,  370,  377,
                                 387,  402,  403,  410,  410,  410,  418,  418,  424,  424,  577,
                                 857,  879,  921,  921, 1007, 1031, 1051, 1051, 1059, 1071, 1071,
                                1080, 1094, 1094, 1110, 1110, 1113, 1113, 1183, 1189, 1200, 1200,
                                1217, 1234, 1234, 1591, 1620, 1620, 1662, 1686, 1707, 1755, 1828,
                                1828, 1848, 1848, 1848, 1848, 1851, 1851, 1852, 1917], dtype=int)
        answer_slices = np.fromiter([(0, 55), (56, 80), (83, 97)], dtype=SDBICAN._DTYPE_SLICE)
        npt.assert_array_equal(SDBICAN.sDBICAN(input_array, 10, max_eps=200, min_eps=10), answer_slices)

    def test_fit(self):
        """
        Test to run class via the method fit.
        Points passed to the fit method should be copied into new array.
        New copy of points should be sorted.
        Most edge cases should be caught in tests for component methods.
        """
        hudc_object = SDBICAN(10, max_eps=200, min_eps=10)
        input_array = np.array([   0,    0,   60,   61,   61,   61,   76,   78,  122,  122,  141,
                                 183,  251,  260,  260,  263,  263,  267,  267,  288,  288,  295,
                                 300,  310,  310,  317,  317,  334,  334,  335,  338,  338,  338,
                                 338,  340,  342,  342,  344,  344,  358,  367,  370,  370,  377,
                                 387,  402,  403,  410,  410,  410,  418,  418,  424,  424,  577,
                                 857,  879,  921,  921, 1007, 1031, 1051, 1051, 1059, 1071, 1071,
                                1080, 1094, 1094, 1110, 1110, 1113, 1113, 1183, 1189, 1200, 1200,
                                1217, 1234, 1234, 1591, 1620, 1620, 1662, 1686, 1707, 1755, 1828,
                                1828, 1848, 1848, 1848, 1848, 1851, 1851, 1852, 1917], dtype=int)
        answer_array = np.array(input_array, copy=True)

        answer_slices = np.fromiter([(0, 55), (56, 80), (83, 97)], dtype=SDBICAN._DTYPE_SLICE)
        hudc_object.fit(input_array)
        assert hudc_object.input_array is not input_array
        npt.assert_array_equal(hudc_object.input_array, answer_array)
        npt.assert_array_equal(hudc_object.slices, answer_slices)
