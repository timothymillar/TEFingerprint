#! /usr/bin/env python

import pytest
import numpy as np
import numpy.testing as npt
from tectoolkit.cluster import UDC, HUDC

class TestUDC:
    """
    Tests for class FlatUnivariateDensityCluster.
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
        assert UDC._sorted_ascending(array) == answer

    @pytest.mark.parametrize("slices,melted_slices",
                             # single slice spanning single base
                             [(np.array([(13, 14)],
                                        dtype=UDC._DTYPE_SLICE),
                               np.array([(13, 14)],
                                        dtype=UDC._DTYPE_SLICE)),
                              # nested slices (this shouldn't happen in practice but is handled correctly)
                              (np.array([(15, 25), (16, 17), (19, 20)],
                                        dtype=UDC._DTYPE_SLICE),
                               np.array([(15, 25)],
                                        dtype=UDC._DTYPE_SLICE)),
                              # adjacent loci (slices are half open intervals)
                              (np.array([(7, 9), (9, 12)],
                                        dtype=UDC._DTYPE_SLICE),
                               np.array([(7, 9), (9, 12)],
                                        dtype=UDC._DTYPE_SLICE)),
                              # combined
                              (np.array([(3, 6), (6, 8), (7, 9), (10, 12), (12, 13), (15, 25), (16, 17), (19, 20)],
                                        dtype=UDC._DTYPE_SLICE),
                               np.array([(3, 6), (6, 9), (10, 12), (12, 13), (15, 25)],
                                        dtype=UDC._DTYPE_SLICE))
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
        npt.assert_array_equal(UDC._melt_slices(slices), melted_slices)

    def test_subcluster(self):
        """
        Test for hidden method _subcluster.
        """
        array = np.array([1, 2, 21, 22, 22, 22, 24, 38, 54, 54, 55, 56, 65, 65, 66, 67, 68, 90],
                         dtype=int)
        slices = np.array([(2, 6), (3, 7), (8, 12), (12, 16), (13, 17)],
                          dtype=UDC._DTYPE_SLICE)
        npt.assert_array_equal(UDC._subcluster(array, 5, 4), slices)

    def test_flat_cluster(self):
        """
        Test for hidden method _cluster.
        Most edge cases should be caught in tests for component methods.
        """
        sub_slices = np.array([1, 2, 21, 22, 22, 22, 24, 38, 54, 54, 55, 56, 65, 65, 66, 67, 68, 90],
                              dtype=int)
        slices = np.array([(2, 7), (8, 12), (12, 17)],
                          dtype=UDC._DTYPE_SLICE)
        npt.assert_array_equal(UDC._cluster(sub_slices, 5, 4), slices)

    @pytest.mark.parametrize("array,slices",
                             [(np.array([], dtype=int),
                               np.array([], dtype=UDC._DTYPE_SLICE)),
                              (np.array([1], dtype=int),
                               np.array([], dtype=UDC._DTYPE_SLICE)),
                              (np.array([1, 2, 3, 4], dtype=int),
                               np.array([], dtype=UDC._DTYPE_SLICE)),
                              (np.array([1, 2, 3, 4, 5], dtype=int),
                               np.array([(0, 5)], dtype=UDC._DTYPE_SLICE))
                              ],
                             ids=['0<n', '1<n', '4<n', '5==n'])
    def test_udc_few_reads(self, array, slices):
        """
        Test for method udc with small arrays.
        Method udc should correctly handle an array of length 0 or greater.
        If the length of the array is less than 'n' then an empty array will always be returned.
        """
        npt.assert_array_equal(UDC.udc(array, 5, 5), slices)

    def test_fit(self):
        """
        Test to run class via the method fit.
        Points passed to the fit method should be copied into new array.
        New copy of points should be sorted.
        Most edge cases should be caught in tests for component methods.
        """
        udc_object = UDC(4, 5)
        input_array = np.array([1, 2, 21, 22, 22, 22, 24, 38, 54, 54, 55, 56, 65, 65, 66, 67, 68, 90],
                               dtype=int)
        answer_array = np.array([1, 2, 21, 22, 22, 22, 24, 38, 54, 54, 55, 56, 65, 65, 66, 67, 68, 90],
                                dtype=int)
        answer_slices = np.array([(2, 7), (8, 12), (12, 17)],
                                 dtype=UDC._DTYPE_SLICE)
        udc_object.fit(input_array)
        assert udc_object.input_array is not input_array
        npt.assert_array_equal(udc_object.input_array, answer_array)
        npt.assert_array_equal(udc_object.slices, answer_slices)

    def test_clusters(self):
        """
        Test for method clusters.
        """
        udc_object = UDC(4, 5)
        input_array = np.array([1, 2, 21, 22, 22, 22, 24, 38, 54, 54, 55, 56, 65, 65, 66, 67, 68, 90], dtype=int)
        cluster_arrays = [np.array([21, 22, 22, 22, 24], dtype=int),
                          np.array([54, 54, 55, 56], dtype=int),
                          np.array([65, 65, 66, 67, 68], dtype=int)]
        udc_object.fit(input_array)
        for result, answer in zip(udc_object.clusters(), cluster_arrays):
            npt.assert_array_equal(result, answer)

    def test_cluster_extremities(self):
        """Test for method cluster_extremities"""
        udc_object = UDC(4, 5)
        input_array = np.array([1, 2, 21, 22, 22, 22, 24, 38, 54, 54, 55, 56, 65, 65, 66, 67, 68, 90], dtype=int)
        cluster_extremities = [(21, 24), (54, 56), (65, 68)]
        udc_object.fit(input_array)
        for result, answer in zip(udc_object.cluster_extremities(), cluster_extremities):
            assert result == answer

    def test_labels(self):
        """Test for method labels"""
        udc_object = UDC(4, 5)
        input_array = np.array([1, 2, 21, 22, 22, 22, 24, 38, 54, 54, 55, 56, 65, 65, 66, 67, 68, 90], dtype=int)
        labels = np.array([-1, -1, 0, 0, 0, 0, 0, -1, 1, 1, 1, 1, 2, 2, 2, 2, 2, -1], dtype=int)
        udc_object.fit(input_array)
        npt.assert_array_equal(udc_object.labels(), labels)


class TestHUDC:
    """
    Tests for class HierarchicalUnivariateDensityCluster.
    """
    def test_point_eps(self):
        """
        Test for hidden method _point_eps.
        Data is derived from real edge case containing "hidden" eps peaks.
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
        npt.assert_array_equal(HUDC._point_eps(input_array, 10), point_eps)

    def test_point_eps_hidden_peak(self):
        """
        Test for hidden method _point_eps.
        Data contains a single "hidden" eps peak between values 6 and 26
        (values themselves have low eps but a large gap between them).
        """
        input_array = np.array([0, 0, 3, 4, 4, 6, 26, 28, 28, 29, 32, 32], dtype=int)
        point_eps = np.array([3, 3, 1, 1, 1, 2, 2, 1, 1, 1, 3, 3], dtype=int)
        npt.assert_array_equal(HUDC._point_eps(input_array, 3), point_eps)

    def test_eps_splits(self):
        """
        Test for hidden method _eps_splits.
        Data is derived from real edge case containing "hidden" eps peaks.
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
        split_eps = np.array([122,  42, 453, 436], dtype=int)
        npt.assert_array_equal(HUDC._eps_splits(input_array, 10), split_eps)

    def test_eps_splits_hidden_peak(self):
        """
        Test for hidden method _eps_splits.
        Data contains a single "hidden" eps peak between values 6 and 26
        (values themselves have low eps but a large gap between them).
        """
        input_array = np.array([0, 0, 3, 4, 4, 6, 26, 28, 28, 29, 32, 32], dtype=int)
        split_eps = np.array([21], dtype=int)
        npt.assert_array_equal(HUDC._eps_splits(input_array, 3), split_eps)

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
        Method _flatten_list returns a generator which should be coerced to a list.
        """
        assert list(HUDC._flatten_list(nested)) == flat

    def test_hudc(self):
        """Test for method hudc"""
        input_array = np.array([   0,    0,   60,   61,   61,   61,   76,   78,  122,  122,  141,
                                 183,  251,  260,  260,  263,  263,  267,  267,  288,  288,  295,
                                 300,  310,  310,  317,  317,  334,  334,  335,  338,  338,  338,
                                 338,  340,  342,  342,  344,  344,  358,  367,  370,  370,  377,
                                 387,  402,  403,  410,  410,  410,  418,  418,  424,  424,  577,
                                 857,  879,  921,  921, 1007, 1031, 1051, 1051, 1059, 1071, 1071,
                                1080, 1094, 1094, 1110, 1110, 1113, 1113, 1183, 1189, 1200, 1200,
                                1217, 1234, 1234, 1591, 1620, 1620, 1662, 1686, 1707, 1755, 1828,
                                1828, 1848, 1848, 1848, 1848, 1851, 1851, 1852, 1917], dtype=int)
        answer_slices = np.fromiter([(0, 55), (56, 80), (83, 97)], dtype=HUDC._DTYPE_SLICE)
        npt.assert_array_equal(HUDC.hudc(input_array, 10, max_eps=200, min_eps=10), answer_slices)

    @pytest.mark.parametrize("array,slices",
                             [(np.array([], dtype=int),
                               np.array([], dtype=UDC._DTYPE_SLICE)),
                              (np.array([1], dtype=int),
                               np.array([], dtype=UDC._DTYPE_SLICE)),
                              (np.array([1, 2, 3, 4], dtype=int),
                               np.array([], dtype=UDC._DTYPE_SLICE)),
                              (np.array([1, 2, 3, 4, 5], dtype=int),
                               np.array([(0, 5)], dtype=UDC._DTYPE_SLICE))
                              ],
                             ids=['0<n', '1<n', '4<n', '5==n'])
    def test_hudc_few_reads(self, array, slices):
        """
        Test for method hudc with small arrays.
        Method hudc should correctly handle an array of length 0 or greater.
        If the length of the array is less than 'n' then an empty array will always be returned.
        """
        npt.assert_array_equal(HUDC.hudc(array, 5, 5), slices)

    def test_fit(self):
        """
        Test to run class via the method fit.
        Points passed to the fit method should be copied into new array.
        New copy of points should be sorted.
        Most edge cases should be caught in tests for component methods.
        """
        hudc_object = HUDC(10, max_eps=200, min_eps=10)
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

        answer_slices = np.fromiter([(0, 55), (56, 80), (83, 97)], dtype=HUDC._DTYPE_SLICE)
        hudc_object.fit(input_array)
        assert hudc_object.input_array is not input_array
        npt.assert_array_equal(hudc_object.input_array, answer_array)
        npt.assert_array_equal(hudc_object.slices, answer_slices)
