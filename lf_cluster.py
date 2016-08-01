#!/usr/bin/env python

import sys
import numpy as np
import pysam
import argparse
from readcluster import ReadCluster
from sklearn.cluster import DBSCAN


def parse_args(args):
    parser = argparse.ArgumentParser('Identify transposon flanking regions')
    parser.add_argument('--reference', type=str)
    parser.add_argument('--strand', nargs='?', default=False)
    parser.add_argument('--read_group', nargs='?', default=False)
    parser.add_argument('--eps', type=int, default=100,
                        help=("When using the DBSCAN method to identify "
                              "read clusters, eps is the minimum distance "
                              "allowable between two points for inclusion "
                              "in the the same neighbourhood"))
    parser.add_argument('--min_reads', type=int, default=5,
                        help=("When using the DBSCAN method to identify "
                              "read clusters, min_reads is the minimum number "
                              "of read tips found in a single neighbourhood "
                              "in order to count as a cluster"))
    return parser.parse_args(args)


def tip(read):
    if read.is_reverse:
        return read.pos
    else:
        return read.pos + read.qlen


def tips(sam):
    for read in sam:
        yield tip(read)


def udbscan(points, eps, min_reads):
    points2d = np.column_stack([points, np.zeros(len(points))])
    dbscan = DBSCAN(eps=eps, min_samples=min_reads).fit(points2d)
    labels = dbscan.labels_.astype(np.int)
    return labels


def call_clusters(sam, reference, read_group, strand, eps, min_reads):
    sam_tips = np.fromiter(tips(sam), np.int)
    cluster_labels = udbscan(sam_tips, eps, min_reads)
    for label in np.unique(cluster_labels)[1:]:
        cluster_tips = sam_tips[np.where(cluster_labels == label)]
        read_cluster = ReadCluster(reference,
                                   read_group,
                                   strand,
                                   cluster_tips)
        yield read_cluster


def main():
    args = parse_args(sys.argv[1:])
    sam = pysam.AlignmentFile("-", "rb")
    for cluster in call_clusters(sam,
                                 args.reference,
                                 args.read_group,
                                 args.strand,
                                 args.eps,
                                 args.min_reads):
        print(str(cluster))


if __name__ == '__main__':
    main()
