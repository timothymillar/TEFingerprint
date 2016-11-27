#! /usr/bin/env python

import pysam
import numpy as np


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
    SAM_FLAGS = {'+': ('-F', '20'),
                 '-': ('-f', '16'),
                 '.': ('-F', '4')}
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


def parse_flag(flag):
    """
    Parses a SAM flag into a boolean array.

    :param flag: SAM flag
    :type flag: int

    :return: Boolean array
    :rtype: :class:`numpy.array`[bool]
    """
    attributes = np.zeros(12, dtype=np.bool)
    bits = np.fromiter(map(int, tuple(bin(int(flag)))[:1:-1]), dtype=np.bool)
    attributes[:bits.shape[0]] = bits
    return attributes


def flag_orientation(flag):
    """
    Determines whether a SAM flag indicates that its read is on the forward or reverse strand.

    :param flag: SAM flag
    :type flag: int

    :return: '+', '-' or None
    :rtype: str | None
    """
    attributes = parse_flag(flag)
    if attributes[2]:  # read is unmapped
        return None
    elif attributes[4]:  # read is reversed
        return '-'
    else:  # read is forwards
        return '+'


def flag_attributes(flag):
    """
    Parses a SAM flag into a dictionary with attributes as keys and booleans as values.

    :param flag: SAM flag
    :type flag: int

    :return: Dictionary with attributes as keys and booleans as values
    :rtype: dict[str, bool]
    """
    attributes = ("read paired",
                  "read mapped in proper pair",
                  "read unmapped mate unmapped",
                  "read reverse strand",
                  "mate reverse strand",
                  "first in pair",
                  "second in pair",
                  "not primary alignment",
                  "read fails platform / vendor quality checks",
                  "read is PCR or optical duplicate",
                  "supplementary alignment")
    values = parse_flag(flag)
    return dict(zip(attributes, values))


if __name__ == '__main__':
    pass
