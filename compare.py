import pysam
import libtec
import numpy as np
from scipy.stats import f_oneway


def compare(input_bams, reference, family, strand, eps, min_reads):
    """

    :param input_bams:
    :param references:
    :param families:
    :param strands:
    :param eps:
    :param min_reads:
    :return:
    """
    flag = libtec.strand2flag(strand)
    setof_samstrings = [pysam.view(flag[0], flag[1], bam, reference, family) for bam in input_bams]
    setof_reads = [libtec.parse_sam_strings(strings, strand) for strings in setof_samstrings]
    setof_clusters = [libtec.simple_subcluster(reads, min_reads, eps) for reads in setof_reads]
    union_clusters = np.append(setof_clusters[0], setof_clusters[1:])
    union_clusters.sort(order='start')
    union_clusters = libtec.merge_clusters(union_clusters)
    for cluster in union_clusters:
        print(cluster)
        setof_cluster_reads = np.array([libtec.reads_in_locus(reads, cluster) for reads in setof_reads])
        print(f_oneway(*setof_cluster_reads))
        counts = np.fromiter(map(len, setof_cluster_reads), dtype=int)
        print((counts/np.max(counts)).std())
