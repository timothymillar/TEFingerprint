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
    _DTYPE_UGRAD = np.dtype([('index', np.int64), ('grad', np.int64)])
    _DTYPE_UPEAK = np.dtype([('left', np.int64), ('right', np.int64), ('eps', np.int64)])
    _DTYPE_UPOOL = np.dtype([('left', np.int64), ('right', np.int64)])
    _DTYPE_UPOINT_EPS = np.dtype([('point', np.int64), ('eps', np.int64)])

    def __init__(self, min_points, max_eps=None, min_eps=None):
        """"""
        self.min_pts = min_points
        self.max_eps = max_eps
        self.min_eps = min_eps
        self.clusters = UnivariateLoci()
        self.points = np.empty_like
        self.labels = np.array([])

    @staticmethod
    def univariate_eps(points, n):
        assert n > 1  # groups must contain atleast two points
        offset = n - 1  # offset for indexing
        length = len(points)
        lower = points[0:length - offset]
        upper = points[offset:length]
        eps_groups = upper - lower
        eps_2d = np.full((n, length), np.max(eps_groups), dtype=int)
        for i in range(n):
            eps_2d[i, i:length - (offset - i)] = eps_groups
        return np.min(eps_2d, axis=0)

    @staticmethod
    def _univariate_gradient(array):
        length = len(array)
        gradients = np.empty(length - 1, dtype=UDC._DTYPE_UGRAD)
        gradients['grad'] = array[1:] - array[:-1]  # calculate gradients
        gradients['index'] = np.arange(1, length, dtype=np.int)  # store index between compared values
        return gradients

    @staticmethod
    def _find_univariate_peaks(gradients):
        gradients = gradients[gradients['grad'] != 0]  # drop plateaus
        gradients['grad'] = np.array(gradients['grad']>0, dtype=np.int)  # convert to boolean like
        is_peak = np.where(gradients['grad'][1:] < gradients['grad'][:-1])[0]
        peaks = np.empty(len(is_peak), dtype=UDC._DTYPE_UPEAK)
        peaks['left'] = gradients['index'][is_peak]
        peaks['right'] = gradients['index'][is_peak + 1]
        return peaks

    @staticmethod
    def univariate_peaks(array):
        gradients = UDC._univariate_gradient(array)
        peaks = UDC._find_univariate_peaks(gradients)
        peaks['eps'] = array[peaks['left']]  # recover initial values from array
        return peaks

    @staticmethod
    def univariate_split_peaks(points):
        splits = np.empty(len(points) - 1, dtype=UDC._DTYPE_UPEAK)
        splits['eps'] = points[1:] - points[:-1]
        splits['left'] = np.arange(1, 1 + len(splits))

        # Merge platues
        gradients = splits['eps'][1:] - splits['eps'][:-1]
        splits = splits[np.append(True, gradients != 0)]  # remove indices following platue
        splits['right'][-1] = len(points) - 1  # Will always be last index
        splits['right'][:-1] = splits['left'][1:] - 1  # right index based on neighbour to left

        # Remove non-peaks
        is_peak = np.logical_and(np.append(False, splits['eps'][1:] > splits['eps'][:-1]),
                                 np.append(splits['eps'][:-1] > splits['eps'][1:], False))
        return splits[is_peak]

    @staticmethod
    def _split_pools(points, peaks):
        # organise split indices for pools.
        # A pool is the inverse of a peak.
        # There is one more pool than there are peaks
        # The right bound of a peak becomes the left bound of a pool and vice versa
        pools = np.empty(len(peaks) + 1, dtype=UDC._DTYPE_UPOOL)
        pools[0]['left'] = 0  # first index of first pool
        pools[-1]['right'] = len(points)  # second index of last pool
        pools['left'][1:] = peaks['right']  # first index of remaining pools
        pools['right'][:-1] = peaks['left']  # second index of remaining pools

        # split up points into their pools
        pool_points = [points[left:right] for left, right in pools]

        # remove points with eps above that of lowest peak (must be done after split)
        pool_points = [points[points['eps'] < np.min(peaks['eps'])] for points in pool_points]

        return pool_points

    @staticmethod
    def _flatten_list(item):
        if isinstance(item, list):
            for element in item:
                for item in UDC._flatten_list(element):
                    yield item
        else:
            yield item

    @staticmethod
    def _hudc_tree(points, base_eps):
        peaks = UDC.univariate_split_peaks(points['point'])

        if len(peaks) == 0:
            # there are no children so return coordinates
            return points['point'][0], points['point'][-1]

        # Flatten over high peaks
        peaks['eps'][peaks['eps'] > base_eps] = base_eps

        threshold_eps = np.max(peaks['eps'])  # compare based on largest peak (least dense point)
        total_area = np.sum(base_eps - points['eps'])
        child_area = np.sum(threshold_eps - points['eps'])  # area above threshold
        parent_area = total_area - child_area  # area bellow threshold

        if parent_area > child_area:
            # parent is lager so return coordinates
            return points['point'][0], points['point'][-1]

        else:
            # combined area of children is larger so divide and repeat
            # identify peaks with eps == threshold
            threshold_peaks = peaks[np.where(peaks['eps'] == threshold_eps)[0]]

            # remove points with eps above threshold (must be done after split)
            pool_points = UDC._split_pools(points, threshold_peaks)

            # recursion with subsets of points and threshold_eps as new base_eps
            return [UDC._hudc_tree(points, threshold_eps) for points in pool_points]

    @staticmethod
    def hudc(array, n, max_eps=None, min_eps=None):
        points = np.empty(len(array), dtype=UDC._DTYPE_UPOINT_EPS)
        points['point'] = array
        points['eps'] = UDC.univariate_eps(array, n)
        if max_eps:
            # remove points with greater eps
            points = points[points['eps'] <= max_eps]
            # initial base_eps
            base_eps = max_eps
        else:
            # base_eps is set to range of values in points
            base_eps = np.max(points['point']) - np.min(points['point'])
        if min_eps:
            # overwrite lower eps
            points['eps'][points['eps'] < min_eps] = min_eps

        # run
        clusters = UDC._hudc_tree(points, base_eps)
        return UnivariateLoci.from_iter(UDC._flatten_list(clusters))

    def fit(self, points):
        """

        :param points:
        :return:
        """
        self.points = np.array(points, copy=True)
        self.points.sort()
        self.clusters = UDC.hudc(self.points, self.min_pts, max_eps=self.max_eps, min_eps=self.min_eps)


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
