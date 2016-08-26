import numpy as np
from statsmodels.nonparametric.kde import KDEUnivariate

locus = np.dtype([('start', np.int64),
                  ('stop', np.int64)])

sam_read = np.dtype([('tip', np.int64),
                     ('tail', np.int64),
                     ('length', np.int64),
                     ('reverse', np.bool)])


def parse_sam_strings(sam_strings, strand):
    """

    :param sam_strings:
    :param strand:
    :return:
    """
    def parse_sam_string(sam_string, strand):
        attr = sam_string.split("\t")
        start = int(attr[3])
        length = len(attr[9])
        end = start+length
        if strand == '+':
            reverse = True
            tip = end
            tail = start
            return tip, tail, length, reverse
        elif strand == '-':
            reverse = False
            tip = start
            tail = end
            return tip, tail, length, reverse

    reads = (parse_sam_string(string, strand) for string in sam_strings)
    reads = np.fromiter(reads, dtype=sam_read)
    return reads


def simple_subcluster(reads, minpts, eps):
    """

    :param points:
    :param minpts:
    :param eps:
    :return:
    """
    offset = minpts - 1
    upper = reads['tip'][offset:]
    lower = reads['tip'][:-offset]
    diff = upper - lower
    dense = diff <= eps
    lower = lower[dense]
    upper = upper[dense]
    loci = ((lower[i], upper[i]) for i in range(len(lower)))
    loci = np.fromiter(loci, dtype=locus)
    return loci


def merge_clusters(clusters):
    """

    :param clusters:
    :return:
    """
    def merge(clusters):
        start = clusters['start'][0]
        stop = clusters['stop'][0]
        for i in range(1, len(clusters)):
            if clusters['start'][i] <= stop:
                stop = clusters['stop'][i]
            else:
                yield start, stop
                start = clusters['start'][i]
                stop = clusters['stop'][i]
        yield start, stop
    return np.fromiter(merge(clusters), dtype=locus)


def simple_cluster(points, minpts, eps):
    """

    :param points:
    :param minpts:
    :param eps:
    :return:
    """
    subclusters = simple_subcluster(points, minpts, eps)
    return merge_clusters(subclusters)


def resample_points(points, n=None, reps=1):
    """

    :param points:
    :param n:
    :param reps:
    :return:
    """
    length = len(points)
    if n is None:
        n = length
    def resample(sample, n, max_index):
        indexes = np.random.choice(max_index, n)
        new_sample = sample[indexes]
        new_sample.sort()
        return new_sample
    return (resample(points, n, length) for _ in range(reps))


def bootstrapped_simple_cluster(points, minpts, eps, reps, p=95):
    """

    :param points:
    :param minpts:
    :param eps:
    :param reps:
    :param support:
    :return:
    """
    sample_sets = resample_points(points, reps=reps)
    replicates = (simple_cluster(sample, minpts, eps) for sample in sample_sets)
    depth = np.zeros(max(points))
    for clusters in replicates:
        for start, stop in clusters:
            depth[start:(stop+1)] +=1
    depth = (depth/reps)*100
    supported_points = np.where(depth >= p)[0]
    supported_points.sort()
    supported_clusters = simple_cluster(supported_points, minpts, eps)
    return supported_clusters



def cluster_points(cluster, points):
    """

    :param cluster:
    :param points:
    :return:
    """
    start, stop = cluster
    subset = points[start:(stop+1)]
    return subset


def cluster_kde(points, cluster, binwidth, margin):
    """

    :param points:
    :param binwidth:
    :param margin:
    :return:
    """
    start, stop = cluster
    subset = cluster_points(cluster, points)
    kde = KDEUnivariate(subset)
    kde.fit(kernel="gau", bw=np.int(binwidth))
    smooth = kde.evaluate(np.arange(start - margin, stop + margin, 1))
    return smooth

