#! /usr/bin/env python

import pysam


SAM_FLAGS = {'+': ('-F', '20'),
             '-': ('-f', '16'),
             '.': ('-F', '4')}


def read_bam_strings(input_bam, reference='', family='', strand='.'):
    """

    :param input_bam:
    :param reference:
    :param family:
    :param strand:
    :return:
    """
    flag = SAM_FLAGS[strand]
    sam_strings = (string for string in pysam.view(*flag, input_bam, reference, family).splitlines())
    return sam_strings


def read_bam_references(input_bam):
    """

    :param input_bam:
    :return:
    """
    with pysam.AlignmentFile(input_bam, 'rb') as bam:
        references = bam.references
        return references

if __name__ == '__main__':
    pass
