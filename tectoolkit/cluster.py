#! /usr/bin/env python

import numpy as np
from functools import reduce
from tectoolkit.gff import GffFeature


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

    def __init__(self, loci):
        """
        Init method for :class:`ReadLoci`.

        :param loci: A numpy array.
        :type loci: :class:`numpy.ndarray`[(int, int)]
        """
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

    def _sort_uloci(self, order=('start', 'stop')):
        """

        :param order:
        :return:
        """
        self.loci.sort(order=order)

    def sort(self, order=('start', 'stop')):
        """
        Sort loci in place by field(s).

        :param order: A valid field or list of fields in :class:`ReadLoci`, defaults to ['start', 'stop']
        :type order: str | list[str]
        """
        self.loci.sort(order=order)

    @staticmethod
    def _melt_uloci(loci):
        """

        :param loci:
        :return:
        """

        def _merge(loci):
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

        loci.sort(order=('start', 'stop'))
        loci = np.fromiter(_merge(loci), dtype=UnivariateLoci.DTYPE_ULOCUS)
        loci.sort(order=('start', 'stop'))
        return loci

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

        self.sort()
        self.loci = self._melt_uloci(self.loci)

    def subset_by_locus(self, start, stop, margin=0, end='start'):
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
        assert end in {'start', 'stop'}
        start -= margin
        stop += margin
        loci = self.loci[np.logical_and(self.loci[end] >= start, self.loci[end] <= stop)]
        return UnivariateLoci(loci)

    def as_gff(self, *args, **kwargs):
        """
        Converts an instance of :class:`ReadLoci` to a list of :class:`GffFeature` objects.
        Start and end fields are automatically filled by loci bounds.

        :param kwargs: Optional gff attributes
        :return: list[:class:`GffFeature`]
        """
        return [GffFeature(*args, start=locus[0], end=locus[1], **kwargs) for locus in self.loci]

    @classmethod
    def from_iterable(cls, iterable):
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


class FUDC(UnivariateLoci):
    """Flat Univariate Density Cluster"""
    def __init__(self, min_points, eps):
        """"""
        self.min_pts = min_points
        self.eps = eps
        self.loci = np.array([], dtype=FUDC.DTYPE_ULOCUS)
        self.points = np.empty_like
        self.labels = np.array([])

    @staticmethod
    def _flat_subcluster(points, min_pts, eps):
        """

        :param points:
        :param min_pts:
        :param eps:
        :return:
        """
        points.sort()
        offset = min_pts - 1
        upper = points[offset:]
        lower = points[:-offset]
        diff = upper - lower
        dense = diff <= eps
        lower = lower[dense]
        upper = upper[dense]
        loci = ((lower[i], upper[i]) for i in range(len(lower)))
        loci = np.fromiter(loci, dtype=FUDC.DTYPE_ULOCUS)
        return loci

    @staticmethod
    def flat_cluster(points, min_pts, eps):
        """

        :param points:
        :param min_pts:
        :param eps:
        :return:
        """
        loci = FUDC._flat_subcluster(points, min_pts, eps)
        if len(loci) > 1:
            loci = FUDC._melt_uloci(loci)
        return loci

    def fit(self, points):
        """

        :param points:
        :return:
        """
        self.points = np.array(points, copy=True)
        self.points.sort()
        self.loci = FUDC.flat_cluster(self.points, self.min_pts, self.eps)


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
        self.loci = np.array([], dtype=FUDC.DTYPE_ULOCUS)
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
        child_loci = HUDC.flat_cluster(points, min_pts, eps - 1)
        if len(child_loci) == 1:
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
            child_points = [HUDC._locus_points(loci, points) for loci in child_loci]
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
        return np.fromiter(loci, dtype=HUDC.DTYPE_ULOCUS)

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
        loci = reduce(np.append, loci_generator)
        return loci

    def fit(self, points):
        """

        :param points:
        :return:
        """
        self.points = np.array(points, copy=True)
        self.points.sort()
        self.loci = HUDC._hierarchical_cluster(self.points, self.min_pts, self.max_eps, self.min_eps)

if __name__ == '__main__':
    pass
