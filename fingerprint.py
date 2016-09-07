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
    flag = libtec.strand2flag(strand)
    sam_strings = pysam.view(flag[0], flag[1], input_bam, reference, family)
    reads = libtec.parse_sam_strings(sam_strings, strand)
    clusters = libtec.simple_cluster(reads, min_reads, eps)
    for cluster in clusters:
        gff = libtec.GffFeature(reference,
                                start=cluster['start'],
                                end=cluster['stop'],
                                strand=strand,
                                ID="{0}_{1}_{2}_{3}".format(family, reference, strand, cluster['start']),
                                Name=family)
        print(str(gff))
