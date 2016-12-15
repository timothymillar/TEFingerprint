#! /usr/bin/env python

import numpy as np
from functools import reduce


class UnivariateLoci(object):
    """
    A collection of univariate loci that relate to reads mapped to a reference genome.

    Each locus is represented as a single element in a numpy array with the following slots:
     - start: np.int64
     - stop: np.int64
    Where 'start' and 'stop' are respectively the lower and upper coordinates that define a segment of the reference.
    These coordinates are inclusive, 1 based integers (i.e, use the SAM coordinate system).
    """
    DTYPE_ULOCUS = np.dtype([('start', np.int64),
                             ('stop', np.int64)])

    def __init__(self, loci=None):
        """
        Init method for :class:`ReadLoci`.

        :param loci: A numpy array.
        :type loci: :class:`numpy.ndarray`[(int, int)]
        """
        if loci is None:
            loci = []
        self.loci = np.array(loci, dtype=UnivariateLoci.DTYPE_ULOCUS, copy=True)

    def __iter__(self):
        """
        Iter method for :class:`ReadLoci`.
        Passes through to wrapped numpy array.

        :return: an iterator of loci
        :rtype: generator[(int, int)]
        """
        for locus in self.loci:
            yield locus

    def __getitem__(self, item):
        """
        Getitem method for :class:`ReadLoci`.
        Passes through to wrapped numpy array.

        :param item: Index
        :type item: int | slice | str | numpy.ndarray[int] | numpy.ndarray[bool]

        :return: An numpy array with dtype = :class:`ReadLoci`.DTYPE_ULOCUS
        :rtype: :class:`numpy.ndarray`[(int, int)]
        """
        return self.loci[item]

    def __len__(self):
        """
        Len method for :class:`ReadLoci`.

        :return: The number of loci in the collection
        :rtype: int
        """
        return len(self.loci)

    def sort(self, order=('start', 'stop')):
        """
        Sort loci in place by field(s).

        :param order: A valid field or list of fields in :class:`ReadLoci`, defaults to ['start', 'stop']
        :type order: str | list[str]
        """
        self.loci.sort(order=order)

    def melt(self):
        """
        Merge overlapping loci into a single loci.
        Loci are sorted and modified in place.

        Example::
            loci = ReadLoci.from_iterable([(1, 4), (2, 6), (6, 7), (8, 10), (9, 12), (14, 15)])
            list(loci)
            [(1, 4), (2, 6), (6, 7), (8, 10), (9, 12), (14, 15)]
            loci.melt()
            list(loci)
            [(1, 7), (8, 12), (14, 15)]
        """
        def _melter(loci):
            start = loci['start'][0]
            stop = loci['stop'][0]
            for i in range(1, len(loci)):
                if loci['start'][i] <= stop:
                    if loci['stop'][i] > stop:
                        stop = loci['stop'][i]
                    else:
                        pass
                else:
                    yield start, stop
                    start = loci['start'][i]
                    stop = loci['stop'][i]
            yield start, stop
        self.sort()
        self.loci = np.fromiter(_melter(self.loci), dtype=UnivariateLoci.DTYPE_ULOCUS)

    def subset_by_locus(self, start, stop, margin=0, end='both'):
        """
        Returns a new ReadGroup object containing (the specified end of) all reads within specified (inclusive) bounds.

        :param start: Lower bound
        :type start: int
        :param stop: Upper bound
        :type stop: int
        :param margin: A value to extend both bounds by, defaults to 0
        :param end: The read end that must fall within the bounds, must be 'tip' or 'tail', defaults to 'tip'
        :type end: str

        :return: The subset of reads that fall within the specified bounds
        :rtype: :class:`ReadGroup`
        """
        assert end in {'start', 'stop', 'both'}
        start -= margin
        stop += margin
        if end == 'both':
            loci = self.loci[np.logical_and(self.loci['start'] >= start, self.loci['stop'] <= stop)]
        else:
            loci = self.loci[np.logical_and(self.loci[end] >= start, self.loci[end] <= stop)]
        return UnivariateLoci(loci)

    @classmethod
    def from_iter(cls, iterable):
        """
        Construct an instance of :class:`ReadLoci` form an iterable.

        :param iterable: Iterable of tuples containing loci bounds
        :type iterable: iterable[(int, int)]

        :return: Instance of :class:`ReadLoci`
        :rtype: :class:`UnivariateLoci`
        """
        loci = UnivariateLoci(np.fromiter(iterable, dtype=UnivariateLoci.DTYPE_ULOCUS))
        loci.sort()
        return loci

    @classmethod
    def append(cls, x, y):
        """
        Combine two :class:`ReadLoci` objects into a single object.

        :param x: Instance of :class:`ReadLoci`
        :typev x: :class:`ReadLoci`
        :param y: Instance of :class:`ReadLoci`
        :type y: :class:`UnivariateLoci`

        :return: Instance of :class:`ReadLoci`
        :rtype: :class:`UnivariateLoci`
        """
        loci = UnivariateLoci(np.append(x.loci, y.loci))
        loci.sort()
        return loci


class UDC(object):
    """Univariate Density Clusters"""
    _DTYPE_SLICE = np.dtype([('start', np.int64), ('stop', np.int64)])

    def __init__(self, min_points, max_eps=None, min_eps=None):
        """"""
        self.min_pts = min_points
        self.max_eps = max_eps
        self.min_eps = min_eps
        self.slices = np.array([], dtype=UDC._DTYPE_SLICE)
        self.input_array = np.array([], dtype=int)

    @staticmethod
    def _sorted_ascending(array):
        return np.sum(array[1:] - array[:-1] < 0) == 0

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
        array.sort()
        offset = n - 1
        upper = array[offset:]
        lower = array[:-offset]
        selected = upper - lower <= eps
        lower_index = np.arange(0, len(lower))[selected]
        upper_index = np.arange(offset, len(array))[selected] + 1
        return np.fromiter(zip(lower_index, upper_index), dtype=UDC._DTYPE_SLICE)

    @staticmethod
    def _cluster(array, eps, n):
        slices = UDC._subcluster(array, eps, n)
        if len(slices) > 1:
            slices = UDC._melt_slices(slices)
        return slices

    @staticmethod
    def _subset_buy_slices(array, slices):
        return [array[left:right] for left, right in slices]

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
                for item in UDC._flatten_list(element):
                    yield item
        else:
            yield item

    @staticmethod
    def _grow_hudc_tree(points, base_eps, n):

        splits = UDC._eps_splits(points['value'], n)

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
            child_points = UDC._subset_buy_slices(points, UDC._cluster(points['value'], threshold_eps, n))
            return [UDC._grow_hudc_tree(points, threshold_eps, n) for points in child_points]

    @staticmethod
    def hudc(array, n, max_eps=None, min_eps=None):
        assert UDC._sorted_ascending(array)
        points = np.empty(len(array), dtype=np.dtype([('value', np.int64),
                                                      ('index', np.int64),
                                                      ('eps', np.int64)]))
        points['value'] = array
        points['index'] = np.arange(len(array), dtype=int)
        points['eps'] = UDC._point_eps(array, n)
        if not max_eps:
            # start at highest observed eps
            max_eps = np.max(points['eps'])
        if min_eps:
            # overwrite lower eps
            points['eps'][points['eps'] < min_eps] = min_eps

        # initial splits
        child_points = UDC._subset_buy_slices(points, UDC._cluster(points['value'], max_eps, n))

        # run
        clusters = [UDC._grow_hudc_tree(points, max_eps, n) for points in child_points]
        return UnivariateLoci.from_iter(UDC._flatten_list(clusters))

    def fit(self, array):
        """

        :param array:
        :return:
        """
        self.input_array = np.array(array, copy=True)
        self.slices = UDC.hudc(self.input_array, self.min_pts, max_eps=self.max_eps, min_eps=self.min_eps)

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



class FUDC(object):
    """Flat Univariate Density Cluster"""

    DTYPE_SUBCLUSTER = np.dtype([('locus', UnivariateLoci.DTYPE_ULOCUS), ('eps', np.int64)])

    def __init__(self, min_points, eps):
        """"""
        self.min_pts = min_points
        self.eps = eps
        self.clusters = UnivariateLoci()
        self.points = np.empty_like
        self.labels = np.array([])

    @staticmethod
    def _subclusters(points, min_pts):
        points.sort()
        offset = min_pts - 1
        upper = points[offset:]
        lower = points[:-offset]
        return np.fromiter(zip(zip(lower, upper), upper - lower), dtype=FUDC.DTYPE_SUBCLUSTER)

    @staticmethod
    def _flat_subcluster( points, min_pts, eps):
        """

        :param points:
        :param min_pts:
        :param eps:
        :return:
        """
        subs = FUDC._subclusters(points, min_pts)
        return UnivariateLoci(subs['locus'][subs['eps'] <= eps])

    @staticmethod
    def flat_cluster(points, min_pts, eps):
        """

        :param points:
        :param min_pts:
        :param eps:
        :return:
        """
        clusters = FUDC._flat_subcluster(points, min_pts, eps)
        if len(clusters) > 1:
            clusters.melt()
        return clusters

    def fit(self, points):
        """

        :param points:
        :return:
        """
        self.points = np.array(points, copy=True)
        self.points.sort()
        self.clusters = FUDC.flat_cluster(self.points, self.min_pts, self.eps)


class HUDC(FUDC):
    """Hierarchical Univariate Density Cluster"""
    _node_template = {'base_eps': None,
                      'base_locus': (None, None),
                      'area': 0,
                      'child_area': 0,
                      'selected': False,
                      'children': None}

    def __init__(self, min_points, max_eps, min_eps):
        """"""
        assert max_eps > min_eps > 1
        self.min_pts = min_points
        self.max_eps = max_eps
        self.min_eps = min_eps
        self.clusters = UnivariateLoci()
        self.points = np.empty_like
        self.labels = np.array([])

    @staticmethod
    def _locus_points(locus, points):
        """

        :param locus:
        :param points:
        :return:
        """
        start, stop = locus
        return points[np.logical_and(points >= start, points <= stop)]

    @staticmethod
    def _grow_tree(points, min_pts, eps, min_eps, area=0, base_eps=None, base_locus=None):
        """

        :param points:
        :param min_pts:
        :param eps:
        :param min_eps:
        :param area:
        :param base_eps:
        :param base_locus:
        :return:
        """
        if base_eps is None:
            base_eps = eps
        if base_locus is None:
            base_locus = (min(points), max(points))
        area += len(points)
        child_clusters = HUDC.flat_cluster(points, min_pts, eps - 1)
        if len(child_clusters) == 1:
            # branch doesn't fork
            if eps > min_eps:
                # branch grows
                return HUDC._grow_tree(points,
                                       min_pts,
                                       eps - 1,
                                       min_eps=min_eps,
                                       area=area,
                                       base_eps=base_eps,
                                       base_locus=base_locus)
            else:
                # branch terminates
                node = HUDC._node_template.copy()
                node['area'] = area
                node['base_eps'] = base_eps
                node['base_locus'] = base_locus
                node['children'] = None
                return node
        else:
            # branch forks
            child_points = [HUDC._locus_points(loci, points) for loci in child_clusters]
            node = HUDC._node_template.copy()
            node['area'] = area
            node['base_eps'] = base_eps
            node['base_locus'] = base_locus
            node['children'] = [HUDC._grow_tree(pts, min_pts, eps - 1, min_eps=min_eps) for pts in child_points]
            return node

    @staticmethod
    def _child_area(node):
        """

        :param node:
        :return:
        """
        if node['children']:
            node['child_area'] = sum([HUDC._child_area(n) for n in node['children']])
        else:
            node['child_area'] = 0
        return node['child_area'] + node['area']

    @staticmethod
    def _select_nodes(node):
        """

        :param node:
        :return:
        """
        if node['area'] > node['child_area']:
            node['selected'] = True
        else:
            for node in node['children']:
                HUDC._select_nodes(node)

    @staticmethod
    def _retrieve_selected_loci(node):
        """

        :param node:
        :return:
        """
        if node['selected']:
            return node['base_locus']
        elif node['children']:
            return [HUDC._retrieve_selected_loci(node) for node in node['children']]

    @staticmethod
    def _flatten_list(item):
        """

        :param item:
        :return:
        """
        if isinstance(item, list):
            for element in item:
                for item in HUDC._flatten_list(element):
                    yield item
        else:
            yield item

    @staticmethod
    def _single_hierarchical_cluster(points, min_pts, max_eps, min_eps):
        """

        :param points:
        :param min_pts:
        :param max_eps:
        :param min_eps:
        :return:
        """
        tree = HUDC._grow_tree(points, min_pts, max_eps, min_eps)
        HUDC._child_area(tree)
        HUDC._select_nodes(tree)
        loci = HUDC._flatten_list(HUDC._retrieve_selected_loci(tree))
        return UnivariateLoci.from_iter(loci)

    @staticmethod
    def _hierarchical_cluster(points, min_pts, max_eps, min_eps):
        """

        :param points:
        :param min_pts:
        :param max_eps:
        :param min_eps:
        :return:
        """
        base_loci = HUDC.flat_cluster(points, min_pts, max_eps)
        base_points = (HUDC._locus_points(locus, points) for locus in base_loci)
        loci_generator = (HUDC._single_hierarchical_cluster(points, min_pts, max_eps, min_eps) for points in base_points)
        clusters = reduce(UnivariateLoci.append, loci_generator)
        return clusters

    def fit(self, points):
        """

        :param points:
        :return:
        """
        self.points = np.array(points, copy=True)
        self.points.sort()
        self.clusters = HUDC._hierarchical_cluster(self.points, self.min_pts, self.max_eps, self.min_eps)

if __name__ == '__main__':
    pass
