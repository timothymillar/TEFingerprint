#! /usr/bin/env python3

import warnings
import numpy as np


def core_distances(array, min_points):
    """
    Identify the core distance (minimum value of epsilon) for each point
    in an array of integers.

    :param array: An array of integers sorted in ascending order
    :type array: :class:`numpy.ndarray`[int]
    :param min_points: number of points used to form a subcluster
    :type min_points: int

    :return: An array of minimum eps values
    :rtype: :class:`numpy.ndarray`[int]
    """
    # groups must contain at least two points
    assert min_points > 1

    # offset for indexing because the minimum points includes itself
    offset = min_points - 1
    length = len(array)
    lower = array[0:length - offset]
    upper = array[offset:length]
    eps_values = upper - lower
    eps_2d = np.full((min_points, length), np.max(eps_values), dtype=int)
    for i in range(min_points):
        eps_2d[i, i:length - (offset - i)] = eps_values
    return np.min(eps_2d, axis=0)


class DBICAN(object):
    """
    Univariate integer implimentation of DBICAN.

    Clusters are identified as suitably dense, continuous regions in
    the data where the threshold density is specified by parameters
    minimum points and epsilon. Points in the array that do not fall
    inside a cluster are classified as noise points.

    The range in values among every set of min-points in the array is
    calculated and compared to epsilon. Sets of n points with a range
    equal to, or less than epsilon are classified as subclusters.
    Overlapping subclusters are then merged together to form clusters.

    :param min_points: The minimum number of points allowed in each
        (sub)cluster
    :type min_points: int
    :param epsilon: The maximum distance allowed among each set of
        n points
    :type epsilon: int
    """

    _DTYPE_SLICE = np.dtype([('lower', np.int64), ('upper', np.int64)])

    def __init__(self, min_points, epsilon):
        self.min_points = min_points
        self.epsilon = epsilon
        self.slices = np.array([], dtype=DBICAN._DTYPE_SLICE)
        self.input_array = np.array([], dtype=int)

    @staticmethod
    def _sorted_ascending(array):
        """
        Checks if an array is sorted in ascending order.

        :param array: a 1 dimensional array of numeric values
        :type array: :class:`numpy.ndarray`[int]

        :return: True if array is sorted in ascending order
        :rtype: bool
        """
        return np.sum(array[1:] - array[:-1] < 0) == 0

    @staticmethod
    def _melt_slices(slices):
        """
        Combines overlapping slices assuming half-open intervals.

        :param slices: An array of paired integers representing the lower
            and upper values of a slice
        :type slices: :class:`numpy.ndarray`[(int, int)]

        :return: An array of paired integers representing the lower and
            upper values of a slice
        :rtype: :class:`numpy.ndarray`[(int, int)]
        """
        lowers, uppers = slices['lower'], slices['upper']
        lowers.sort()
        uppers.sort()
        splits = np.append(np.array([False]),
                           uppers[:-1] <= lowers[1:])  # True means gap between

        def _melter(lowers, uppers, splits):
            l = lowers[0]
            u = uppers[0]
            for i in range(1, len(lowers)):
                if splits[i]:
                    # there is a gap so yield slice and start new one
                    yield l, u
                    l, u = slices['lower'][i], slices['upper'][i]
                else:
                    # no gap so merge slices
                    u = uppers[i]
            # yield final slice
            yield l, u

        # return slices formatted as a numpy array
        return np.fromiter(_melter(lowers, uppers, splits),
                           dtype=DBICAN._DTYPE_SLICE)

    @staticmethod
    def _subcluster(array, min_points, epsilon):
        """
        Calculate sub-clusters of an array and return their slice indices.

        The input array must be sorted in ascending order.
        If n is greater than the length of the array, an empty array
        is returned.

        :param array: An array of ints sorted in ascending order
        :type array: :class:`numpy.ndarray`[int]
        :param min_points: The maximum distance allowed in among points
            of each subcluster
        :type min_points: int
        :param epsilon: The number of points in each subcluster
        :type epsilon: int

        :return: An array of paired lower and upper indices for each
            subcluster found in the array
        :rtype: :class:`numpy.ndarray`[(int, int)]
        """
        assert DBICAN._sorted_ascending(array)

        offset = min_points - 1
        upper = array[offset:]
        lower = array[:-offset]
        selected = upper - lower <= epsilon
        lower_index = np.arange(0, len(lower))[selected]
        upper_index = np.arange(offset, len(array))[selected] + 1
        return np.fromiter(zip(lower_index, upper_index),
                           dtype=DBICAN._DTYPE_SLICE)

    @staticmethod
    def _cluster(array, min_points, epsilon):
        """
        Calculate clusters of an array and return their slice indices.
        The input array must be sorted in ascending order.
        If n is greater than the length of the array, an empty array
        is returned.

        :param array: An array of ints sorted in ascending order
        :type array: :class:`numpy.ndarray`[int]
        :param min_points: The maximum distance allowed in among each set
            of n points
        :type min_points: int
        :param epsilon: The minimum number of points in each (sub)cluster
        :type epsilon: int

        :return: An array of paired lower and upper indices for each
            cluster found in the array
        :rtype: :class:`numpy.ndarray`[(int, int)]
        """
        # sorted-ascending checked in method _subcluster
        slices = DBICAN._subcluster(array, min_points, epsilon)
        if len(slices) > 1:
            slices = DBICAN._melt_slices(slices)
        return slices

    def fit(self, array):
        """
        Fit an array to an DBICAN model.
        The input array must be sorted in ascending order.

        :param array: An array of integers sorted in ascending order
        :type array: :class:`numpy.ndarray`[int]
        """
        assert self._sorted_ascending(array)
        self.input_array = np.array(array, copy=True)
        self.slices = self._cluster(self.input_array,
                                    self.min_points,
                                    self.epsilon)

    @staticmethod
    def dbican(array, min_points, epsilon):
        """
        Provides functional use of :class:`DBICAN`.

        See documentation for :class:`DBICAN`.

        :param array: An array of ints sorted in ascending order
        :type array: :class:`numpy.ndarray`[int]
        :param min_points: The minimum number of points in each
            (sub)cluster
        :type min_points: int
        :param epsilon: The maximum distance allowed in among each
            set of n points
        :type epsilon: int

        :return: An array of paired lower and upper indices for
            each cluster found in the array
        :rtype: :class:`numpy.ndarray`[int, int]
        """
        model = DBICAN(min_points, epsilon)
        model.fit(array)
        return model.slices

    def clusters(self):
        """
        Return values from the input array grouped into clusters.
        Points classified as noise are not returned.

        :return: A generator object of subset of the input array
        :rtype: generator[:class:`numpy.ndarray`[int]]
        """
        return (self.input_array[lower:upper]
                for lower, upper in self.slices)

    def cluster_extremities(self):
        """
        Return minimum and maximum values from the input array found
        in each cluster.

        :return: A generator object of pairs of lower and upper values
            found in each cluster
        :rtype: generator[(int, int)]
        """
        return ((self.input_array[lower], self.input_array[upper - 1])
                for lower, upper in self.slices)

    def labels(self):
        """
        Return cluster labels for all input values.
        Clusters are labeled in ascending order starting from 0
        Points classified as noise are labeled as -1

        :return: An array of integer labels
        :rtype: :class:`numpy.ndarray`[int]
        """
        labels = np.full(len(self.input_array), -1, int)
        for i, (lower, upper) in enumerate(self.slices):
            labels[lower:upper] += (i + 1)
        return labels


class SDBICAN(DBICAN):
    """
    Univariate integer implimentation of Splitting-DBICAN.

    Clusters are identified as suitably dense, continuous regions in
    the data where the threshold density is specified by parameters
    minimum points and epsilon. Points in the array that do not fall
    inside a cluster are classified as noise points.

    The range in values among every set of min-points in the array is
    calculated and compared to epsilon. Sets of n points with a range
    equal to, or less than epsilon are classified as subclusters.
    Overlapping subclusters are then merged together to form clusters.

    Once the initial set of clusters has been identified following DBICAN,
    the support of each cluster is assessed in comparision to its child
    clusters. If support for a cluster is weak it is split into it's
    constituent child clusters which are in turn assessed recursively.

    The search space may be constrained by setting a
    minimum allowable value of epsilon.
    If a minimum value for epsilon is selected, core distances bellow
    the specified minimum epsilon are set to the specified minimum epsilon.
    Setting a higher minimum epsilon will result in more conservative
    splitting of clusters by SDBICAN.
    This option is not typically used and may be deprecated in future
    versions.

    :param min_points: The minimum number of points allowed in each
        cluster
    :type min_points: int
    :param epsilon: The maximum distance allowed among
        each set of min-points
    :type epsilon: int
    :param min_epsilon: An optional value for the minimum value of epsilon to
        be used when calculating cluster depth
    :type min_epsilon: int
    :param aggressive_method: use old more aggressive splitting method
        (deprecated)
    :type aggressive_method: bool
    """

    def __init__(self,
                 min_points,
                 epsilon,
                 min_epsilon=None,
                 aggressive_method=False):
        self.min_points = min_points
        self.epsilon = epsilon
        self.min_epsilon = min_epsilon
        self.aggressive_method = aggressive_method
        if aggressive_method:
            warnings.warn("The aggressive splitting method is deprecated "
                          "and may be removed in future versions.",
                          DeprecationWarning)
        self.slices = np.array([], dtype=DBICAN._DTYPE_SLICE)
        self.input_array = np.array([], dtype=int)

    @staticmethod
    def _core_distances(array, min_points):
        """
        Identify the core distance (minimum value of epsilon) for every
        point in an array of integers.

        For each point, the minimum epsilon is calculated among all
        subclusters containing that point.

        :param array: An array of integers sorted in ascending order
        :type array: :class:`numpy.ndarray`[int]
        :param min_points: number of points used to form a subcluster
        :type min_points: int

        :return: An array of minimum eps values
        :rtype: :class:`numpy.ndarray`[int]
        """
        return core_distances(array, min_points)

    @staticmethod
    def _fork_epsilon(array, min_points):
        """
        Identify eps 'splits' in an array by calculating epsilon of the
        gaps between values in the array.

        Identifies the minimum epsilon of a cluster prior to it forking
        into child clusters.
        If the cluster does not fork into children then None is returned.
        Note that in this cas the minimum epsilon of the cluster is then
        equal to the minimum core distance of points in that cluster.

        Example when clustering the array of vales:
            [0, 0, 3, 4, 4, 6, 26, 28, 28, 29, 32, 32]
        with min-points = 3 The core distance calculated for values 6 and
        26 is 2 and 2 respectively.
        However the minimum eps value calculated for the cluster that
        include both values is 22 (26 - 4 or 28 - 6).

        Once epsilon is calculated for all gaps between values in an array,
        peaks are identified and classified as 'splits'.
        Splits are points at which a parent cluster are split into
        child clusters.

        Eps of splits is reduced by a constant of 1 in order to be used
        as a threshold value. This is because an eps
        value of 21 implies that at eps = 21 a valid subcluster is formed.
        Therefore eps = 20 is the value at which a subcluster cannot be
        formed ad there is a split.

        :param array: An array of integers sorted in ascending order
        :type array: :class:`numpy.ndarray`[int]
        :param min_points: number of points used to form a subcluster
        :type min_points: int

        :return: An epsilon value or None
        :rtype: int
        """
        if len(array) <= min_points:
            # no forks possible because all points must have the same eps
            return None

        offset = min_points - 1

        # calculate split eps using the 2d method
        eps_values = array[offset:] - array[:-offset]
        eps_2d = np.full((offset, len(eps_values) + offset - 1),
                         np.max(eps_values),
                         dtype=int)
        for i in range(offset):
            eps_2d[i, i:len(eps_values) + i] = eps_values
        splits = np.min(eps_2d, axis=0)

        # Remove plateaus
        gradients = splits[1:] - splits[:-1]
        splits = splits[np.append(np.array([True]), gradients != 0)]

        # Remove non-peaks
        is_peak = np.logical_and(np.append(np.array([False]),
                                           splits[1:] > splits[:-1]),
                                 np.append(splits[:-1] > splits[1:],
                                           np.array([False])))
        splits = splits[is_peak]

        # If this method calculates epsilon of 5 it means the child
        # cluster starts at epsilon 4.999...
        if len(splits) == 0:
            # The cluster does not fork into child clusters at all
            return None
        else:
            # We only require the largest fork value
            return np.max(splits)

    @staticmethod
    def _flatten_list(item):
        """
        Flatten a nested list.

        If item is not a list it will be returned inside a list
        of length one.

        :param item: a list or any other object
        :type item: list[any] | any

        :return: A flat list
        :rtype: list[any]
        """
        if isinstance(item, list):
            for element in item:
                for item in SDBICAN._flatten_list(element):
                    yield item
        else:
            yield item

    def _traverse_cluster_tree(self,
                               local_points,
                               local_max_eps):
        """
        Traverse a tree of nested density clusters and recursively
        identify clusters based on their area.

        :param local_points: An array of points with the slots 'value',
            'index' and 'eps
        :type local_points: :class:`numpy.ndarray`[(int, int, int)]
        :param local_max_eps: Maximum distance allowed in among
            each set of n points
        :type local_max_eps: int

        :return: A nested list of upper and (half-open) indices
            of selected clusters
        :rtype: list[list[]]
        """
        # Values of epsilon bellow which the cluster forks
        fork_epsilon = self._fork_epsilon(local_points['value'],
                                          self.min_points)

        if fork_epsilon is None:
            # The cluster doesn't fork so it has no children
            # Epsilon_minimum would equal the minimum of core
            # distances but it's not needed
            return local_points['index'][0], local_points['index'][-1] + 1

        # If a cluster forks into children then it's minimum epsilon
        # is the value at which forks
        local_min_eps = fork_epsilon

        # Compare support for cluster and its children
        if self.aggressive_method:
            support = np.sum(local_max_eps -
                             np.maximum(local_min_eps,
                                        local_points['core_dist']))
        else:
            support = np.sum(self.epsilon -
                             np.maximum(local_min_eps,
                                        local_points['core_dist']))

        support_children = np.sum(np.maximum(0,
                                             local_min_eps -
                                             local_points['core_dist']))

        if support >= support_children:
            # Parent is supported so return slice indices
            return local_points['index'][0], local_points['index'][-1] + 1

        else:
            # Combined support of children is larger so divide
            # and repeat recursively:
            # A minimum epsilon of 5 means the child clusters technically
            # starts at epsilon 4.999...
            # we calculate the child clusters using epsilon 4 which will
            # produce the same clusters as 4.999...
            child_cluster_bounds = self._cluster(local_points['value'],
                                                 self.min_points,
                                                 local_min_eps - 1)
            child_points = (local_points[left:right]
                            for left, right in child_cluster_bounds)
            # but then use epsilon 5 as the new maximum epsilon so that
            # support is calculated from epsilon 4.999...
            return [self._traverse_cluster_tree(points,
                                                local_min_eps)
                    for points in child_points]

    def _run(self, array):
        """
        See documentation for :class: `SDBICAN`.

        :param array: An array of integers sorted in ascending order
        :type array: :class:`numpy.ndarray`[int]

        :return: An array of paired lower and upper indices for each
            cluster found in the array
        :rtype: :class:`numpy.ndarray`[(int, int)]
        """
        assert self._sorted_ascending(array)

        if len(array) < self.min_points:
            # not enough points to form a cluster
            return np.array([], dtype=DBICAN._DTYPE_SLICE)

        else:
            # create a structured array for tracking value, index and
            # core distance of each individual
            points = np.empty(len(array),
                              dtype=np.dtype([('value', np.int64),
                                              ('index', np.int64),
                                              ('core_dist', np.int64)]))
            points['value'] = array
            points['index'] = np.arange(len(array), dtype=int)
            points['core_dist'] = self._core_distances(array, self.min_points)
            if self.epsilon is None:
                warnings.warn("Automatic setting of epsilon is deprecated "
                              "and may be removed in future versions. "
                              "The value of epsilon should be set manually.",
                              DeprecationWarning)
                # start at first fork, i.e. bellow root node
                # as in HDBSCAN* (Campello 2015)
                self.epsilon = self._fork_epsilon(points['value'],
                                                  self.min_points) - 1
            if self.min_epsilon:
                # overwrite lower eps
                points['core_dist'] = np.maximum(self.min_epsilon,
                                                 points['core_dist'])

            # initial splits based on the specified max_eps
            initial_cluster_bounds = self._cluster(points['value'],
                                                   self.min_points,
                                                   self.epsilon)
            child_points = (points[left:right]
                            for left, right in initial_cluster_bounds)

            # recursively run on all clusters
            clusters = [self._traverse_cluster_tree(points,
                                                    self.epsilon)
                        for points in child_points]
            return np.fromiter(self._flatten_list(clusters),
                               dtype=DBICAN._DTYPE_SLICE)

    def fit(self, array):
        """
        Fit an array to a SDBICAN model.
        The input array must be sorted in ascending order.

        :param array: An array of integers sorted in ascending order
        :type array: :class:`numpy.ndarray`[int]
        """
        self.input_array = np.array(array, copy=True)
        self.slices = self._run(self.input_array)

    @staticmethod
    def sdbican(array,
                min_points,
                epsilon=None,
                min_epsilon=None,
                aggressive_method=False):
        """
        Provides functional use of :class:`SDBICAN`.

        See documentation for :class:`SDBICAN`.

        :param array: An array of integers sorted in ascending order
        :type array: :class:`numpy.ndarray`[int]
        :param min_points: The minimum number of points allowed in
            each (sub)cluster
        :type min_points: int
        :param epsilon: An optional value for maximum distance allowed
            among each set of n points
        :type epsilon: int
        :param min_epsilon: An optional value for the minimum value of
            eps to be used when calculating cluster depth
        :type min_epsilon: int
        :param aggressive_method: use old more aggressive splitting method
        :type aggressive_method: bool

        :return: An array of paired lower and upper indices for each
            cluster found in the array
        :rtype: :class:`numpy.ndarray`[(int, int)]
        """
        model = SDBICAN(min_points,
                        epsilon=epsilon,
                        min_epsilon=min_epsilon,
                        aggressive_method=aggressive_method)
        model.fit(array)
        return model.slices


if __name__ == '__main__':
    pass
