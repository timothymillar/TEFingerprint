#! /usr/bin/env python

import pysam
import numpy as np


SAM_FLAGS = {'+': ('-F', '20'),
             '-': ('-f', '16'),
             '.': ('-F', '4')}


def read_bam_strings(input_bam, reference='', family='', strand='.'):
    """
    Read a subset of strings from an indexed Bam file.

    :param input_bam: The path to an indexed sorted Bam file
    :type input_bam: str
    :param reference: Select reads from a single reference or slice of reference
    :type reference: str
    :param family: Select reads that have names beginning with 'family'
    :type family: str
    :param strand: Select reads that are found on a specific strand, '+' or '-'
    :type strand: str

    :return: A list of Sam formatted strings
    :rtype: list[str]
    """
    flag = SAM_FLAGS[strand]
    sam_strings = np.array(pysam.view(*flag, input_bam, reference).splitlines())
    in_family = np.array([string.startswith(family) for string in sam_strings])
    return sam_strings[in_family]


def read_bam_references(input_bam):
    """
    Read all reference names from a Bam file.

    :param input_bam: The path to an indexed sorted Bam file
    :type input_bam: str

    :return: A list of reference names
    :rtype: list[str]
    """
    with pysam.AlignmentFile(input_bam, 'rb') as bam:
        references = bam.references
        return references


def read_bam_reference_lengths(input_bam):
    """
    Read lengths of references from a Bam file

    :param input_bam: The path to an indexed sorted Bam file
    :type input_bam: str

    :return: A dictionary of the form {<reference>: <length>}
    :rtype: dict[str, int]
    """
    with pysam.AlignmentFile(input_bam, 'rb') as bam:
        references = bam.header['SQ']
        references = {r['SN']: r['LN'] for r in references}
        return references

if __name__ == '__main__':
    pass
