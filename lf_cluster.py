#!/usr/bin/env python

import sys
import numpy as np
import pysam
import argparse
import collections
from sklearn.cluster import DBSCAN


def parse_args(args):
    parser = argparse.ArgumentParser('Identify transposon flanking regions')
    parser.add_argument('input_bam')
    parser.add_argument('--eps', type=int, default=100,
                        help=("When using the DBSCAN method to identify "
                              "read clusters, eps is the minimum distance "
                              "allowable between two points for inclusion "
                              "in the the same neighbourhood"))
    parser.add_argument('--min_tips', type=int, default=5,
                        help=("When using the DBSCAN method to identify "
                              "read clusters, min_tips is the minimum number "
                              "of read tips found in a single neighbourhood "
                              "in order to count as a cluster"))
    return parser.parse_args(args)


def split_references(sam):
    if sam.nreferences == 1:
        yield sam
    else:
        for name in sam.references:
            yield sam.fetch(name)


def tip(read):
    if read.is_reverse:
        return read.pos + read.qlen
    else:
        return read.pos


def tips(sam):
    for read in sam:
        yield tip(read)


def cluster(tips):
    tips = np.fromiter(tips, np.int)
    input_tips = np.column_stack([tips, np.zeros(len(tips))])
    dbscan = DBSCAN(eps=100, min_samples=5).fit(input_tips)
    labels = dbscan.labels_.astype(np.int)
    for cluster_label in np.unique(labels)[1:]:
        cluster_tips = tips[np.where(labels == cluster_label)]
        yield(np.min(cluster_tips),
              np.max(cluster_tips),
              len(cluster_tips),
              np.mean(cluster_tips),
              np.median(cluster_tips),
              np.bincount(cluster_tips).argmax())


def main():
    args = parse_args(sys.argv[1:])
    sam = pysam.AlignmentFile(args.input_bam, 'rb')
    sams = split_references(sam)
    for sam in sams:
        for clust in cluster(tips(sam)):
            print(str(clust))


if __name__ == '__main__':
    main()