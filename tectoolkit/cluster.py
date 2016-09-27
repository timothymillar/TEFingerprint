#! /usr/bin/env python

import numpy as np
from functools import reduce

_ulocus = np.dtype([('start', np.int64),
                   ('stop', np.int64)])


def _melt_uloci(loci):
    """

    :param loci:
    :return:
    """
    def merge(loci):
        start = loci['start'][0]
        stop = loci['stop'][0]
        for i in range(1, len(loci)):
            if loci['start'][i] <= stop:
                stop = loci['stop'][i]
            else:
                yield start, stop
                start = loci['start'][i]
                stop = loci['stop'][i]
        yield start, stop

    loci.sort(order=('start', 'stop'))
    loci = np.fromiter(merge(loci), dtype=_ulocus)
    loci.sort(order=('start', 'stop'))
    return loci


class UnivariateLoci(object):
    locus = _ulocus

    def __init__(self, loci):
        self.loci = loci

    def __iter__(self):
        for locus in self.loci:
            yield locus

    def __getitem__(self, item):
        return self.loci[item]

    def __len__(self):
        return len(self.loci)

    def sort(self, order=('start', 'stop')):
        """

        :param order:
        :return:
        """
        self.loci.sort(order=order)

    def melt(self):
        """

        :return:
        """

        self.sort()
        self.loci = _melt_uloci(self.loci)

    @classmethod
    def from_iterable(cls, iterable):
        """

        :param iterable:
        :return:
        """
        loci = UnivariateLoci(np.fromiter(iterable, dtype=UnivariateLoci.locus))
        loci.sort()
        return loci

    @classmethod
    def merge(cls, x, y):
        """

        :param x:
        :param y:
        :return:
        """
        loci = UnivariateLoci(np.append(x.loci, y.loci))
        loci.sort()
        return loci


class UnivariateDensityCluster(object):
    """"""
    _ulocus = np.dtype([('start', np.int64),
                        ('stop', np.int64)])

    _node_template = {'base_eps': None,
                      'base_locus': (None, None),
                      'area': 0,
                      'child_area': 0,
                      'selected': False,
                      'children': None}

    def __init__(self, min_points, eps, method):
        assert method in {'flat', 'hierarchical'}
        self.min_pts = min_points
        self.eps = eps
        self.method = method
        self.loci = np.array([], dtype=UnivariateDensityCluster._ulocus)
        self.points = np.empty_like
        self.labels = np.array([])

    @classmethod
    def _melt_uloci(cls, loci):
        """

        :param loci:
        :return:
        """

        def merge(loci):
            start = loci['start'][0]
            stop = loci['stop'][0]
            for i in range(1, len(loci)):
                if loci['start'][i] <= stop:
                    stop = loci['stop'][i]
                else:
                    yield start, stop
                    start = loci['start'][i]
                    stop = loci['stop'][i]
            yield start, stop

        loci.sort(order=('start', 'stop'))
        loci = np.fromiter(merge(loci), dtype=UnivariateDensityCluster._ulocus)
        loci.sort(order=('start', 'stop'))
        return loci

    @classmethod
    def _locus_points(cls, locus, points):
        start, stop = locus
        return points[np.logical_and(points >= start, points <= stop)]

    def _flat_subcluster(self, points, min_pts, eps):
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
        loci = np.fromiter(loci, dtype=UnivariateDensityCluster._ulocus)
        return loci

    def _flat_cluster(self, points, min_pts, eps):
        """

        :param points:
        :param min_pts:
        :param eps:
        :return:
        """
        loci = self._flat_subcluster(points, min_pts, eps)
        if len(loci) > 1:
            loci = UnivariateDensityCluster._melt_uloci(loci)
        return loci

    def _grow_tree(self, points, min_pts, eps, min_eps=2, area=0, base_eps=None, base_locus=None):
        if base_eps is None:
            base_eps = eps
        if base_locus is None:
            base_locus = (min(points), max(points))
        area += len(points)
        child_loci = self._flat_cluster(points, min_pts, eps - 1)
        if len(child_loci) == 1:
            # branch doesn't fork
            if eps > min_eps:
                # branch grows
                return self._grow_tree(points,
                                       min_pts,
                                       eps - 1,
                                       min_eps=min_eps,
                                       area=area,
                                       base_eps=base_eps,
                                       base_locus=base_locus)
            else:
                # branch terminates
                node = UnivariateDensityCluster._node_template.copy()
                node['area'] = area
                node['base_eps'] = base_eps
                node['base_locus'] = base_locus
                node['children'] = None
                return node
        else:
            # branch forks
            child_points = [self._locus_points(loci, points) for loci in child_loci]
            node = UnivariateDensityCluster._node_template.copy()
            node['area'] = area
            node['base_eps'] = base_eps
            node['base_locus'] = base_locus
            node['children'] = [self._grow_tree(pts, min_pts, eps -1, min_eps=min_eps) for pts in child_points]
            return node

    def _child_area(self, node):
        if node['children']:
            node['child_area'] = sum([self._child_area(n) for n in node['children']])
        else:
            node['child_area'] = 0
        return node['child_area'] + node['area']

    def _select_nodes(self, node):
        if node['area'] > node['child_area']:
            node['selected'] = True
        else:
            for node in node['children']:
                self._select_nodes(node)

    def _retrieve_selected_loci(self, node):
        if node['selected']:
            return node['base_locus']
        elif node['children']:
            return [self._retrieve_selected_loci(node) for node in node['children']]

    def _flatten(self, lst):
        for item in lst:
            if not isinstance(item, list):
                yield item
            else:
                for x in self._flatten(item):
                    yield x

    def _single_hierarchical_cluster(self, points, min_pts, eps):
        tree = self._grow_tree(points, min_pts, eps)
        self._child_area(tree)
        self._select_nodes(tree)
        loci = self._flatten(self._retrieve_selected_loci(tree))
        return np.fromiter(loci, dtype=UnivariateDensityCluster._ulocus)

    def _hierarchical_cluster(self, points, min_pts, eps):
        base_loci = self._flat_cluster(points, min_pts, eps)
        base_points = (UnivariateDensityCluster._locus_points(locus, points) for locus in base_loci)
        loci_generator = (self._single_hierarchical_cluster(points, min_pts, eps) for points in base_points)
        loci = reduce(np.append, loci_generator)
        return loci

    def fit(self, points):
        self.points = np.array(points, copy=True)
        self.points.sort()
        if self.method == 'flat':
            self.loci = self._flat_cluster(self.points, self.min_pts, self.eps)
        elif self.method == 'hierarchical':
            self.loci = self._hierarchical_cluster(self.points, self.min_pts, self.eps)

if __name__ == '__main__':
    pass
