#! /usr/bin/env python

import numpy as np


class UDC(object):
    """Univariate Density Clusters"""
    _DTYPE_SLICE = np.dtype([('start', np.int64), ('stop', np.int64)])

    def __init__(self, min_points, eps):
        """"""
        self.min_pts = min_points
        self.eps = eps
        self.slices = np.array([], dtype=UDC._DTYPE_SLICE)
        self.input_array = np.array([], dtype=int)

    @staticmethod
    def _sorted_ascending(array):
        return np.sum(array[1:] - array[:-1] < 0) == 0

    @staticmethod
    def _melt_slices(slices):
        starts, stops = slices['start'], slices['stop']
        starts.sort()
        stops.sort()
        splits = np.append(np.array([False]), stops[:-1] <= starts[1:])  # True means gap between

        def _melter(starts, stops, splits):
            start = starts[0]
            stop = stops[0]
            for i in range(1, len(starts)):
                if splits[i]:
                    # there is a gap so yield slice and start new one
                    yield start, stop
                    start, stop = slices['start'][i], slices['stop'][i]
                else:
                    # no gap so merge slices
                    stop = stops[i]
            # yield final slice
            yield start, stop

        return np.fromiter(_melter(starts, stops, splits), dtype=UDC._DTYPE_SLICE)


    @staticmethod
    def _subcluster(array, eps, n):
        """
        The input array must be sorted in ascending order

        :param array: a sorted array of ints
        :param eps:
        :param n:
        :return:
        """
        assert UDC._sorted_ascending(array)
        offset = n - 1
        upper = array[offset:]
        lower = array[:-offset]
        selected = upper - lower <= eps
        lower_index = np.arange(0, len(lower))[selected]
        upper_index = np.arange(offset, len(array))[selected] + 1
        return np.fromiter(zip(lower_index, upper_index), dtype=UDC._DTYPE_SLICE)

    @staticmethod
    def _cluster(array, eps, n):
        """
        The input array must be sorted in ascending order

        :param array: a sorted array of ints
        :param eps:
        :param n:
        :return:
        """
        # sorted-ascending checked in method _subcluster
        slices = UDC._subcluster(array, eps, n)
        if len(slices) > 1:
            slices = UDC._melt_slices(slices)
        return slices

    @staticmethod
    def _subset_buy_slices(array, slices):
        return [array[left:right] for left, right in slices]

    @staticmethod
    def udc(array, n, eps):
        assert UDC._sorted_ascending(array)
        slices = UDC._cluster(array, eps, n)
        return slices

    def fit(self, array):
        """

        :param array:
        :return:
        """
        self.input_array = np.array(array, copy=True)
        self.slices = UDC.udc(self.input_array, self.min_pts, self.eps)

    def clusters(self):
        """
        Return values from the input array grouped into clusters

        :return:
        """
        return (self.input_array[lower:upper] for lower, upper in self.slices)

    def cluster_extremities(self):
        """
        Return minimum and maximum values from the input array found in each cluster

        :return:
        """
        return ((self.input_array[lower], self.input_array[upper - 1]) for lower, upper in self.slices)

    def labels(self):
        """
        Return cluster labels for input values

        :return:
        """
        labels = np.full(len(self.input_array), -1, int)
        for i, (lower, upper) in enumerate(self.slices):
            labels[lower:upper] += (i + 1)
        return labels


class HUDC(UDC):
    """Univariate Density Clusters"""

    def __init__(self, min_points, max_eps=None, min_eps=None):
        """"""
        self.min_pts = min_points
        self.max_eps = max_eps
        self.min_eps = min_eps
        self.slices = np.array([], dtype=UDC._DTYPE_SLICE)
        self.input_array = np.array([], dtype=int)

    @staticmethod
    def _point_eps(array, n):
        assert n > 1  # groups must contain at least two points
        offset = n - 1  # offset for indexing
        length = len(array)
        lower = array[0:length - offset]
        upper = array[offset:length]
        eps_values = upper - lower
        eps_2d = np.full((n, length), np.max(eps_values), dtype=int)
        for i in range(n):
            eps_2d[i, i:length - (offset - i)] = eps_values
        return np.min(eps_2d, axis=0)

    @staticmethod
    def _eps_splits(array, n):
        if len(array) <= n:
            # no peaks possible because all points must have the same eps
            return np.array([], dtype=np.int)

        offset = n - 1

        # calculate split eps using the 2d method
        eps_values = array[offset:] - array[:-offset]
        eps_2d = np.full((offset, len(eps_values) + offset - 1), np.max(eps_values), dtype=int)
        for i in range(offset):
            eps_2d[i, i:len(eps_values) + i] = eps_values
        splits = np.min(eps_2d, axis=0)
        splits -= 1  # Convert values to thresholds

        # Remove plateaus
        gradients = splits[1:] - splits[:-1]
        splits = splits[np.append(np.array([True]), gradients != 0)]

        # Remove non-peaks
        is_peak = np.logical_and(np.append(np.array([False]), splits[1:] > splits[:-1]),
                                 np.append(splits[:-1] > splits[1:], np.array([False])))
        return splits[is_peak]

    @staticmethod
    def _flatten_list(item):
        if isinstance(item, list):
            for element in item:
                for item in HUDC._flatten_list(element):
                    yield item
        else:
            yield item

    @staticmethod
    def _grow_hudc_tree(points, base_eps, n):

        splits = HUDC._eps_splits(points['value'], n)

        if len(splits) == 0:
            # there are no children so return slice indices
            return points['index'][0], points['index'][-1] + 1

        # compare based on largest peak (least dense place)
        threshold_eps = np.max(splits)

        # compare areas
        total_area = np.sum(base_eps - points['eps'])
        child_area = np.sum(threshold_eps - points['eps'])
        parent_area = total_area - child_area

        if parent_area > child_area:
            # parent is lager so return slice indices
            return points['index'][0], points['index'][-1] + 1

        else:
            # combined area of children is larger so divide and repeat
            child_points = HUDC._subset_buy_slices(points, HUDC._cluster(points['value'], threshold_eps, n))
            return [HUDC._grow_hudc_tree(points, threshold_eps, n) for points in child_points]

    @staticmethod
    def hudc(array, n, max_eps=None, min_eps=None):
        assert HUDC._sorted_ascending(array)
        points = np.empty(len(array), dtype=np.dtype([('value', np.int64),
                                                      ('index', np.int64),
                                                      ('eps', np.int64)]))
        points['value'] = array
        points['index'] = np.arange(len(array), dtype=int)
        points['eps'] = HUDC._point_eps(array, n)
        if not max_eps:
            # start at highest observed eps
            max_eps = np.max(points['eps'])
        if min_eps:
            # overwrite lower eps
            points['eps'][points['eps'] < min_eps] = min_eps

        # initial splits
        child_points = HUDC._subset_buy_slices(points, HUDC._cluster(points['value'], max_eps, n))

        # run
        clusters = [HUDC._grow_hudc_tree(points, max_eps, n) for points in child_points]
        return np.fromiter(HUDC._flatten_list(clusters), dtype=UDC._DTYPE_SLICE)

    def fit(self, array):
        """

        :param array:
        :return:
        """
        self.input_array = np.array(array, copy=True)
        self.slices = HUDC.hudc(self.input_array, self.min_pts, max_eps=self.max_eps, min_eps=self.min_eps)


if __name__ == '__main__':
    pass
