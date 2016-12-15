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
        assert UDC._sorted_ascending(array) == answer

    @pytest.mark.parametrize("slices,melted_slices",
                             # single slice spanning single base
                             [(np.array([(13, 13)],
                                        dtype=UDC._DTYPE_SLICE),
                               np.array([(13, 13)],
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
                              (np.array([(3, 6), (6, 8), (7, 9), (10, 12), (13, 13), (15, 25), (16, 17), (19, 20)],
                                        dtype=UDC._DTYPE_SLICE),
                               np.array([(3, 6), (6, 9), (10, 12), (13, 13), (15, 25)],
                                        dtype=UDC._DTYPE_SLICE))
                              ],
                             ids=['single', 'nested', 'adjacent', 'combined'])
    def test_melt_slices(self, slices, melted_slices):
        """
        Test for hidden method _melt_slices.
        Test includes following edge cases:
         * Long locus completely overlaps short loci: (15, 25) & (16, 17) & (19, 20) --> (15, 25)
         * Adjacent loci do not get merged: (7, 9) & (9, 12) -->  (*, 9) & (9, *)
         * Locus may span a single base: (13, 13) --> (13, 13)
        """
        npt.assert_array_equal(UDC._melt_slices(slices), melted_slices)


    def test__subcluster(self):
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


class TestHUDC:
    """
    Tests for class HierarchicalUnivariateDensityCluster.
    """
    def test_point_eps(selfs):
        pass

    def test__eps_splits(self):
        pass

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
        Test for hidden method _flatten_list using a nested data set.
        Method _flatten_list returns a generator which should be coerced to a list.
        """
        assert list(HUDC._flatten_list(nested)) == flat

    def test_grow_hudc_tree(self):
        pass

    def test_hudc(self):
        pass

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
