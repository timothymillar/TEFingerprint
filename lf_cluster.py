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
    parser.add_argument('-s', '--stringency', type=int, default=5,
                        help=("Stringency to determine if a hit maps "
                              "uniquely. This score is the "
                              "minimal allowed difference between the "
                              "scores of the first and second "
                              "hit before a read is assumed to be "
                              "mapped to it's correct, unique, "
                              "location. (default=5)"))
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


def sub_cluster(parent_cluster, read_subset, **kwargs):
    """
    Returns a modified cluster with a subset of the original reads.
    Additional parameters can be added as **kwargs.
    Accepts a dictionary
    Returns a dictionary
    """
    child_cluster = {}

    # avoid passing reference to parent cluster or parent clusters reads
    for key in parent_cluster:
        if key == "reads":
            pass
        else:
            child_cluster[key] = parent_cluster[key]

    # add new attributes passed as kwargs to child cluster
    for key, value in kwargs.items():
        child_cluster[key] = value

    # add the explicitly passed reads to child cluster
    child_cluster["reads"] = read_subset

    return child_cluster


def split_families(cluster_generator):
    """
    Subdivides read-clusters based on read family.
    Accepts a dictionary generator.
    Returns a dictionary generator.
    """
    for parent_cluster in cluster_generator:
        families = collections.defaultdict(list)
        for read in parent_cluster["reads"]:
            try:
                family = read.qname.split('__')[0].rsplit('/', 1)[1]
            except IndexError:
                family = read.qname.split('__')[0]
            families[family].append(read)

        for family, reads in families.items():
            child_cluster = sub_cluster(parent_cluster, reads, family=family)
            yield child_cluster


def split_orientation(cluster_generator):
    """
    Subdivides read-clusters based on read orientation.
    Accepts a dictionary generator.
    Returns a dictionary generator.
    """
    for parent_cluster in cluster_generator:
        orientations = {"+": [],
                        "-": []}
        for read in parent_cluster["reads"]:
            if read.is_reverse:
                orientations["-"].append(read)
            else:
                orientations["+"].append(read)

        for orientation, reads in orientations.items():
            if len(reads) > 0:
                child_cluster = sub_cluster(parent_cluster, reads, orientation=orientation)
                yield child_cluster


def filter_unique(cluster_generator, args):
    """
    Filters read-clusters by the amount of uniquely mapped reads.
    threshold=args.min_diff
    Accepts a dictionary generator.
    Returns a dictionary generator.
    """
    for parent_cluster in cluster_generator:
        unique_reads = []
        for read in parent_cluster["reads"]:
            tag_as = -999
            tag_xs = -999
            for tag in read.tags:
                if tag[0] == 'AS':
                    tag_as = tag[1]
                if tag[0] == 'XS':
                    tag_xs = tag[1]

            score = tag_as - tag_xs
            if score >= args.stringency:
                unique_reads.append(read)
            else:
                pass

        if len(unique_reads) > 0:
            child_cluster = sub_cluster(parent_cluster, unique_reads, read_type='UNIQUE')
            yield child_cluster


def read_depth(cluster):
    """
    Calculates the read depth of a cluster.
    Accepts a dictionary.
    Returns an array.
    """
    depth = np.zeros((cluster["stop"] - cluster["start"]))
    for read in cluster["reads"]:
        depth[(read.pos - cluster["start"]):(read.pos + read.qlen - cluster["start"])] += 1
    return depth


def read_tips(cluster):#
    """
    Returns the read end positions of a cluster based on orientation.
    Returns right tip of forwards reads and left tip of reverse reads
    Accepts a dictionary.
    Returns an array.
    """
    tips = np.zeros(len(cluster["reads"]))
    count = 0
    for read in cluster["reads"]:
        if read.is_reverse:
            tips[count] = read.pos
        else:
            tips[count] = read.pos + read.qlen
        count += 1
    return tips.astype(np.int)


def identify_features_by_dbscan(cluster_generator, args):
    """
    Identifies features based on a DBSCAN clustering algorithm on tip positions
    :param args: command line arguments
    :param cluster_generator:  a dictionary generator
    :return: a dictionary generator
    """
    for cluster in cluster_generator:
        tips = read_tips(cluster).astype(np.int)
        input_tips = np.array(zip(tips, np.zeros(len(tips))), dtype=np.int)
        dbscan = DBSCAN(eps=args.eps, min_samples=args.min_tips).fit(input_tips)
        cluster["feature"] = np.zeros((cluster["stop"] - cluster["start"]), dtype=bool)
        labels = dbscan.labels_.astype(np.int)
        groups = np.unique(labels)
        groups = groups[groups >= 0]
        for group in groups:
            group_tips = tips[labels == group]
            cluster["feature"][min(group_tips) - 1: max(group_tips)] = True
        cluster["depth"] = read_depth(cluster)
        yield cluster


def extract_features(cluster_generator):
    """
    Extracts features for a gff file from clusters based of the feature attribute-array.
    Features are returned as dictionaries with start and stop attributes and inherit other attributes
    from their parent cluster.
    Accepts a dictionary generator.
    Returns a dictionary generator.
    """
    for cluster in cluster_generator:

        feature = {"start": 0, "stop": 0, "mean_depth": 0}
        for key in cluster:
            if key not in ["reads", "start", "stop", "feature", "depth"]:
                feature[key] = cluster[key]
            else:
                pass

        # determine position of features
        previously_in_feature = False
        for position, currently_in_feature in enumerate(cluster["feature"]):

            if not previously_in_feature and currently_in_feature:
                # start of a feature
                previously_in_feature = True
                feature["start"] = position + 1

            elif previously_in_feature and not currently_in_feature:
                # end of a feature
                previously_in_feature = False
                feature["stop"] = position
                feature["mean_depth"] = cluster["depth"][feature["start"] - 1: feature["stop"]].mean()
                feature["depth"] = cluster["depth"][feature["start"] - 1: feature["stop"]]
                yield feature


def format_features(feature_generator):
    """
    Formats feature dictionaries as strings for a gff file.
    Accepts a dictionary generator.
    Returns a string generator.
    """
    for feature in feature_generator:
        formated = "\t".join([feature["reference"],
                              "REFS",
                              "REFS." + feature["read_type"] + "." + feature["family"],
                              str(feature["start"]),
                              str(feature["stop"]),
                              str(feature["mean_depth"]),
                              feature["orientation"],
                              ".",
                              "ID=reps" + feature["reference"] + str(feature["start"]) +
                              feature["family"] + ";Name=" + feature["family"]])
        yield formated


def output_features(formatted_features):
    for feature in formatted_features:
        print(feature)


def main():
    sam = pysam.AlignmentFile(infile, "rb")
    sams = split_references(sam)
    for sam in sams:
        for clust in cluster(tips(sam)):
            print(str(clust))


if __name__ == '__main__':
    main()