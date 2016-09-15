#! /usr/bin/env python

import pysam
import libtec


def fingerprint(input_bam, reference, family, strand, eps, min_reads):
    """

    :param input_bam:
    :param references:
    :param families:
    :param strands:
    :param eps:
    :param min_reads:
    :return:
    """
    reads = libtec.read_sam_reads(input_bam, reference, family, strand)
    clusters = libtec.simple_cluster(reads, min_reads, eps)
    for cluster in clusters:
        gff = libtec.GffFeature(reference,
                                start=cluster['start'],
                                end=cluster['stop'],
                                strand=strand,
                                ID="{0}_{1}_{2}_{3}".format(family, reference, strand, cluster['start']),
                                Name=family)
        print(str(gff))


if __name__ == '__main__':
    pass