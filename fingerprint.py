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
    clusters = libtec.simple_subcluster(reads, min_reads, eps)
    [print(cluster) for cluster in clusters]