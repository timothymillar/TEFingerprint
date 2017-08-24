#! /usr/bin/env python

import numpy as np


class UDBSCANx(object):
    """
    Univariate implimentation of DBSCAN*.

    This is an implementation of the DBSCAN* (Campello et al 2015) clustering method for use on
    a sorted (in ascending order) univariate array of integers.

    Clusters are identified as suitably dense, continuous regions in the data where the threshold density
    is specified by parameters minimum points and epsilon. Points in the array that do not fall
    inside a cluster are classified as noise points.

    The range in values among every set of min-points in the array is calculated and compared to epsilon.
    Sets of n points with a range equal to, or less than epsilon are classified as subclusters.
    Overlapping subclusters are then merged together to form clusters.

    :param min_points: The minimum number of points allowed in each (sub)cluster
    :type min_points: int
    :param epsilon: The maximum distance allowed among each set of n points
    :type epsilon: int
    """

    _DTYPE_SLICE = np.dtype([('lower', np.int64), ('upper', np.int64)])

    def __init__(self, min_points, epsilon):
        self.min_points = min_points
        self.epsilon = epsilon
        self.slices = np.array([], dtype=UDBSCANx._DTYPE_SLICE)
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

        :param slices: An array of paired integers representing the lower and upper values of a slice
        :type slices: :class:`numpy.ndarray`[(int, int)]

        :return: An array of paired integers representing the lower and upper values of a slice
        :rtype: :class:`numpy.ndarray`[(int, int)]
        """
        lowers, uppers = slices['lower'], slices['upper']
        lowers.sort()
        uppers.sort()
        splits = np.append(np.array([False]), uppers[:-1] <= lowers[1:])  # True means gap between

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
        return np.fromiter(_melter(lowers, uppers, splits), dtype=UDBSCANx._DTYPE_SLICE)

    @staticmethod
    def _subcluster(array, epsilon, min_points):
        """
        Calculate subclusters of an array and return their slice indices.
        The input array must be sorted in ascending order.
        If n is greater than the length of the array, an empty array is returned.

        :param array: An array of ints sorted in ascending order
        :param epsilon: The number of points in each subcluster
        :param min_points: The maximum distance allowed in among points of each subcluster

        :return: An array of paired lower and upper indices for each subcluster found in the array
        :rtype: :class:`numpy.ndarray`[(int, int)]
        """
        assert UDBSCANx._sorted_ascending(array)
        offset = min_points - 1
        upper = array[offset:]
        lower = array[:-offset]
        selected = upper - lower <= epsilon
        lower_index = np.arange(0, len(lower))[selected]
        upper_index = np.arange(offset, len(array))[selected] + 1
        return np.fromiter(zip(lower_index, upper_index), dtype=UDBSCANx._DTYPE_SLICE)

    @staticmethod
    def _cluster(array, epsilon, min_points):
        """
        Calculate clusters of an array and return their slice indices.
        The input array must be sorted in ascending order.
        If n is greater than the length of the array, an empty array is returned.

        :param array: An array of ints sorted in ascending order
        :param epsilon: The minimum number of points in each (sub)cluster
        :param min_points: The maximum distance allowed in among each set of n points

        :return: An array of paired lower and upper indices for each cluster found in the array
        :rtype: :class:`numpy.ndarray`[(int, int)]
        """
        # sorted-ascending checked in method _subcluster
        slices = UDBSCANx._subcluster(array, epsilon, min_points)
        if len(slices) > 1:
            slices = UDBSCANx._melt_slices(slices)
        return slices

    @staticmethod
    def udbscanx(array, min_points, epsilon):
        """
        See documentation for :class: `UDBSCANx`.

        :param array: An array of ints sorted in ascending order
        :type array: :class:`numpy.ndarray`[int]
        :param min_points: The minimum number of points in each (sub)cluster
        :type n: int
        :param epsilon: The maximum distance allowed in among each set of n points
        :type epsilon: int

        :return: An array of paired lower and upper indices for each cluster found in the array
        :rtype: :class:`numpy.ndarray`[int, int]
        """
        assert UDBSCANx._sorted_ascending(array)
        slices = UDBSCANx._cluster(array, epsilon, min_points)
        return slices

    def fit(self, array):
        """
        Fit an array to a Univariate Density Cluster model.
        The input array must be sorted in ascending order.

        :param array: An array of integers sorted in ascending order
        :type array: :class:`numpy.ndarray`[int]
        """
        self.input_array = np.array(array, copy=True)
        self.slices = UDBSCANx.udbscanx(self.input_array, self.min_points, self.epsilon)

    def clusters(self):
        """
        Return values from the input array grouped into clusters.
        Points classified as noise are not returned.

        :return: A generator object of subset of the input array
        :rtype: generator[:class:`numpy.ndarray`[int]]
        """
        return (self.input_array[lower:upper] for lower, upper in self.slices)

    def cluster_extremities(self):
        """
        Return minimum and maximum values from the input array found in each cluster.

        :return: A generator object of pairs of lower and upper values found in each cluster
        :rtype: generator[(int, int)]
        """
        return ((self.input_array[lower], self.input_array[upper - 1]) for lower, upper in self.slices)

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


class UDBSCANxH(UDBSCANx):
    """
    Univariate DBSCAN* with Hierarchical component.

    This is a variation of the DBSCAN*/HDBSCAN* (Campello et al 2015) clustering methods which allows some flexibility
    when calculating density based clusters. While it has a hierarchical component it differs significantly from
    HDBSCAN* in how cluster support is calculated and clusters are selected.

    This implementation is for use on a sorted (in ascending order) univariate array of integers.

    This method is useful when there is some existing knowledge about the expected density of informative clusters
    that are present in the data. It preferred over DBSCAN* because it allows *some* flexibility in the density of
    clusters identified. However its flexibility is much more limited than HDBSCAN* to ensure that clusters are of
    similar density to one another.

    Using the definitions of Campello et al 2015, this method differs from HDBSCAN* when calculating support for
    clusters in that support of a cluster is calculated as:

    .. math::
        (\textbf{C}_i) = \sum_{\textbf{x}_j \in \textbf{C}_i} \varepsilon_{\text{max}}(\textbf{C}_i) - \varepsilon_{\text{min}}(\textbf{x}_j, \textbf{C}_i)

    where:

    .. math::
        \varepsilon_{\text{min}}(\textbf{x}_j, \textbf{C}_i)

    is calculated as:

    .. math::
        \text{max}\{d_\text{core}(\textbf{x}_j) , \varepsilon_{\text{min}}(\textbf{C}_i) \}


    Clusters selection is performed in a top down manner from the root of the tree. A given cluster is selected if
    its support is greater than the sum of support for all of its decedent clusters. A cluster cannot be selected
    if one of its ancestor clusters has been selected. This process results in a set of flat (non-overlapping) clusters.

    The search space may be constrained by setting global maximum and minimum allowable values of epsilon.
    If a maximum value for epsilon is selected, the search algorithm will start at the specified value. If a
    minimum value for epsilon is selected, core distances bellow the specified level are increased to that level.

    :param min_points: The minimum number of points allowed in each cluster
    :type min_points: int
    :param max_eps: An optional value for maximum distance allowed among each set of min-points
    :type max_eps: int
    :param min_eps: An optional value for the minimum value of eps to be used when calculating cluster depth
    :type min_eps: int
    """

    def __init__(self, min_points, max_eps=None, min_eps=None):
        self.min_pts = min_points
        self.max_eps = max_eps
        self.min_eps = min_eps
        self.slices = np.array([], dtype=UDBSCANx._DTYPE_SLICE)
        self.input_array = np.array([], dtype=int)

    @staticmethod
    def _core_distances(array, min_points):
        """
        Identify the core distance (minimum value of epsilon) for every point in an array of integers.
        For each point, the minimum epsilon is calculated among all subclusters containing that point.

        :param array: An array of integers sorted in ascending order
        :type array: :class:`numpy.ndarray`[int]
        :param min_points: number of points used to form a subcluster
        :type n: int

        :return: An array of minimum eps values
        :rtype: :class:`numpy.ndarray`[int]
        """
        assert min_points > 1  # groups must contain at least two points
        offset = min_points - 1  # offset for indexing because the minimum points includes itself
        length = len(array)
        lower = array[0:length - offset]
        upper = array[offset:length]
        eps_values = upper - lower
        eps_2d = np.full((min_points, length), np.max(eps_values), dtype=int)
        for i in range(min_points):
            eps_2d[i, i:length - (offset - i)] = eps_values
        return np.min(eps_2d, axis=0)

    @staticmethod
    def _fork_epsilon(array, min_points):
        """
        Identify eps 'splits' in an array by calculating epsilon of the gaps between values in the array.

        Identifies the minimum epsilon of a cluster prior to it forking into child clusters.
        If the cluster does not fork into children then None is returned. Note that in this cas the minimum epsilon
        of the cluster is then equal to the minimum core distance of points in that cluster.

        Example when clustering the array of vales:
            [0, 0, 3, 4, 4, 6, 26, 28, 28, 29, 32, 32]
        with min-points = 3 The core distance calculated for values 6 and 26 is 2 and 2 respectively.
        However the minimum eps value calculated for the cluster that include both values is 22 (26 - 4 or 28 - 6).

        Once eps is calculated for all gaps between values in an array, peaks are identified and classified as 'splits'.
        Splits are points at which a parent cluster are split into child clusters.

        Eps of splits is reduced by a constant of 1 in order to be used as a threshold value. This is because an eps
        value of 21 implies that at eps = 21 a valid subcluster is formed. Therefore eps = 20 is the value at which
        a subcluster cannot be formed ad there is a split.

        :param array: An array of integers sorted in ascending order
        :type array: :class:`numpy.ndarray`[int]
        :param min_points: number of points used to form a subcluster
        :type n: int

        :return: An epsilon value or None
        :rtype: int
        """
        if len(array) <= min_points:
            # no forks possible because all points must have the same eps
            return None

        offset = min_points - 1

        # calculate split eps using the 2d method
        eps_values = array[offset:] - array[:-offset]
        eps_2d = np.full((offset, len(eps_values) + offset - 1), np.max(eps_values), dtype=int)
        for i in range(offset):
            eps_2d[i, i:len(eps_values) + i] = eps_values
        splits = np.min(eps_2d, axis=0)

        # Remove plateaus
        gradients = splits[1:] - splits[:-1]
        splits = splits[np.append(np.array([True]), gradients != 0)]

        # Remove non-peaks
        is_peak = np.logical_and(np.append(np.array([False]), splits[1:] > splits[:-1]),
                                 np.append(splits[:-1] > splits[1:], np.array([False])))
        splits = splits[is_peak]

        # If this method calculates epsilon of 5 it means the child cluster starts at epsilon 4.999...
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
        If item is not a list it will be returned inside a list of length one

        :param item: a list or any other object
        :type item: list[any] | any

        :return: A flat list
        :rtype: list[any]
        """
        if isinstance(item, list):
            for element in item:
                for item in UDBSCANxH._flatten_list(element):
                    yield item
        else:
            yield item

    @staticmethod
    def _traverse_tree(points, epsilon_maximum, min_points):
        """
        Traverse a tree of nested density clusters and recursively identify clusters based on their area.
        This method is intimately tied to method 'hudc'.

        :param points: An array of points with the slots 'value', 'index' and 'eps
        :type points: :class:`numpy.ndarray`[(int, int, int)]
        :param epsilon_maximum: Maximum distance allowed in among each set of n points
        :type epsilon_maximum: int
        :param min_points: The minimum number of points allowed in each (sub)cluster
        :type min_points: int

        :return: A nested list of upper and (half-open) indices of selected clusters
        :rtype: list[list[]]
        """
        # Values of epsilon bellow which the cluster forks
        fork_epsilon = UDBSCANxH._fork_epsilon(points['value'], min_points)

        if fork_epsilon is None:
            # The cluster doesn't fork so it has no children
            # Epsilon_minimum is equal the minimum of core distances
            epsilon_minimum = np.min(points['core_dist'])
        else:
            # If a cluster forks into children then it's minimum epsilon is the value at which forks
            epsilon_minimum = fork_epsilon

        # calculate support for cluster and total of children
        support = np.sum(epsilon_maximum - np.maximum(epsilon_minimum, points['core_dist']))
        support_children = np.sum(np.maximum(0, epsilon_minimum - points['core_dist']))

        # bounds of cluster are index of input array, not genome positions
        cluster = {'index': (points['index'][0], points['index'][-1] + 1),
                   'epsilon_maximum': epsilon_maximum,
                   'epsilon_minimum': epsilon_minimum,
                   'support': support,
                   'support_children': support_children}

        if fork_epsilon is None:
            # cluster has no children
            cluster['stability_hat'] = support
            cluster['selected'] = True
            cluster['children'] = []

        if support > support_children:
            # children cannot be selected
            cluster['stability_hat'] = support
            cluster['selected'] = True
            cluster['children'] = []

        else:
            # Recurse down to children:
            cluster['stability_hat'] = None
            cluster['selected'] = False
            # A minimum epsilon of 5 means the child clusters technically starts at epsilon 4.999...
            # we calculate the child clusters using epsilon 4 which will produce the same clusters as 4.999...
            child_cluster_bounds = UDBSCANxH._cluster(points['value'], epsilon_minimum - 1, min_points)
            child_points = (points[left:right] for left, right in child_cluster_bounds)
            # but then use epsilon 5 as the new maximum epsilon so that support is calculated from epsilon 4.999...
            cluster['children'] = [UDBSCANxH._traverse_tree(points, epsilon_minimum, min_points)
                                   for points in child_points]

        return cluster

    @staticmethod
    def _node_stability(node):
        if node['stability_hat']:
            return node['stability_hat']
        else:
            child_stability_hat = sum([UDBSCANxH._node_stability(child) for child in node['children']])
            if node['support'] >= child_stability_hat:
                node['selected'] = True
                node['stability_hat'] = node['support']
            else:
                node['stability_hat'] = child_stability_hat
        return node['stability_hat']

    @staticmethod
    def _select_tree_clusters(node):
        if node['selected']:
            return node['index']
        else:
            return [UDBSCANxH._select_tree_clusters(child) for child in node['children']]

    @staticmethod
    def _build_tree(array, min_points, max_eps=None, min_eps=None):
        """
        See documentation for :class: `UDBSCANxH`.

        :param array: An array of integers sorted in ascending order
        :type array: :class:`numpy.ndarray`[int]
        :param min_points: The minimum number of points allowed in each (sub)cluster
        :type min_points: int
        :param max_eps: An optional value for maximum distance allowed among each set of n points
        :type max_eps: int
        :param min_eps: An optional value for the minimum value of eps to be used when calculating cluster depth
        :type min_eps: int

        :return: An array of paired lower and upper indices for each cluster found in the array
        :rtype: :class:`numpy.ndarray`[(int, int)]
        """
        assert UDBSCANxH._sorted_ascending(array)

        if len(array) < min_points:
            # not enough points to form a cluster
            return np.array([], dtype=UDBSCANx._DTYPE_SLICE)

        else:
            # create a structured array for tracking value, index and core distance of each individual
            points = np.empty(len(array), dtype=np.dtype([('value', np.int64),
                                                          ('index', np.int64),
                                                          ('core_dist', np.int64)]))
            points['value'] = array
            points['index'] = np.arange(len(array), dtype=int)
            points['core_dist'] = UDBSCANxH._core_distances(array, min_points)
            if not max_eps:
                # start at first fork, i.e. bellow root node as in HDBSCAN* (Campello 2015)
                max_eps = UDBSCANxH._fork_epsilon(points['value'], min_points) - 1
            if min_eps:
                # overwrite lower eps
                points['core_dist'] = np.maximum(min_eps, points['core_dist'])

            # initial splits based on the specified max_eps
            initial_cluster_bounds = UDBSCANxH._cluster(points['value'], max_eps, min_points)
            child_points = (points[left:right] for left, right in initial_cluster_bounds)
            # recursively run on all clusters
            # initialise with max_eps + 1 to ensure points with core_distance == max_eps are counted
            trees = [UDBSCANxH._traverse_tree(points, max_eps + 1, min_points) for points in child_points]

            # calculate node stability
            for tree in trees:
                UDBSCANxH._node_stability(tree)

            return trees
            #return np.fromiter(UDBSCANxH._flatten_list(clusters), dtype=UDBSCANx._DTYPE_SLICE)

    @staticmethod
    def udbscanxh(array, min_points, max_eps=None, min_eps=None):
        tree = UDBSCANxH._build_tree(array, min_points, max_eps=max_eps, min_eps=min_eps)
        slices = [UDBSCANxH._select_tree_clusters(node) for node in tree]
        slices = np.fromiter(UDBSCANxH._flatten_list(slices), dtype=UDBSCANxH._DTYPE_SLICE)
        return slices

    def fit(self, array):
        """
        Fit an array to a Hierarchical Univariate Density Cluster model.
        The input array must be sorted in ascending order.

        :param array: An array of integers sorted in ascending order
        :type array: :class:`numpy.ndarray`[int]
        """
        self.input_array = np.array(array, copy=True)
        self.slices = UDBSCANxH.udbscanxh(self.input_array, self.min_pts, max_eps=self.max_eps, min_eps=self.min_eps)


if __name__ == '__main__':
    pass
