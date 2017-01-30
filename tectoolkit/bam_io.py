#! /usr/bin/env python

import os
import re
import pysam
import numpy as np
from tectoolkit.reads import ReadGroup


def read_bam_strings(input_bam, reference='', strand='.'):
    """
    Read  strings from an indexed Bam file.

    :param input_bam: The path to an indexed sorted Bam file
    :type input_bam: str
    :param reference: Select reads from a single reference or slice of reference
    :type reference: str
    :param strand: Select reads that are found on a specific strand, '+' or '-'
    :type strand: str

    :return: A list of Sam formatted strings
    :rtype: list[str]
    """
    SAM_FLAGS = {'+': ('-F', '20'),
                 '-': ('-f', '16', '-F', '4')}
    flag = SAM_FLAGS[strand]
    return np.array(pysam.view(*flag, input_bam, reference).splitlines())


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
                  "read unmapped",
                  "mate unmapped",
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


def _cigar_mapped_length(cigar):
    """
    Calculate length of the mapped section of a read from a sam CIGAR string.
    This length is calculated based on the length of the section of reference genome which the
    read is mapped to. Therefore, deletions are counted and insertions are not.
    Values for the symbols 'M', 'D', 'N', 'P', 'X' and '=' count towards length.

    :param cigar: A sam format CIGAR string thant may contain the symbols 'MIDNSHP=Q'
    :type cigar: str

    :return: length of the mapped section of read as it appears on the reference genome
    :rtype: int
    """
    return sum([int(i[0:-1]) for i in re.findall(r"[0-9]+[MIDNSHP=X]", cigar) if (i[-1] in 'MDNP=X')])


def _parse_sam_strings(strings, strand=None):
    """
    Parses a collection of SAM formatted strings into a tuple generator.

    :param strings: A collection of SAM formatted strings
    :type strings: iterable[str]
    :param strand: Strand ('+' or '-') of all reads (if known)
    :type strand: str

    :return: An iterable of mapped SAM read positions and names
    :rtype: generator[(int, int, str, str)]
    """

    def _parse_sam_string(string, strand):
        attr = string.split("\t")
        name = str(attr[0])
        start = int(attr[3])
        length = _cigar_mapped_length(attr[5])
        end = start + length - 1  # 1 based indexing used in SAM format
        if strand is None:
            strand = flag_orientation(int(attr[1]))
        if strand == '+':
            tip = end
            tail = start
            return tip, tail, strand, name
        elif strand == '-':
            tip = start
            tail = end
            return tip, tail, strand, name

    assert strand in ['+', '-', None]
    reads = (_parse_sam_string(string, strand) for string in strings)
    return reads


def read_bam_into_stranded_groups(bam, reference, strand, groups, group_tag='ME'):
    """
    Read a section of a bam file and return as a generator of :class:`ReadGroup`.
    One :class:`ReadGroup` is returned per group.
    Reads are sorted into any group for which there tag starts with the group name.

    :param bam: The path to an indexed sorted Bam file
    :type bam: str
    :param reference: Select reads from a single reference or slice of reference
    :type reference: str
    :param strand: Select reads that are found on a specific strand, '+' or '-'
    :type strand: str
    :param group_tag: Sam format tag which holds group information for each read
    :type group_tag: str
    :param groups: A list of group names
    :type groups: list[str]

    :return: A generator of :class:`ReadGroup`
    :rtype: generator[:class:`ReadGroup`]
    """
    assert strand in ['+', '-']
    strings = read_bam_strings(bam, reference=reference, strand=strand)
    group_tag = '\t' + group_tag + ':[Zi]:'
    tags = np.array([re.split(group_tag, s)[1].split('\t')[0] for s in strings])
    generator = (strings[tags.astype('U{0}'.format(len(group))) == group] for group in groups)
    generator = (_parse_sam_strings(strings, strand=strand) for strings in generator)
    return (ReadGroup.from_iter(reads,
                                reference=reference,
                                grouping=group,
                                source=os.path.basename(bam)) for group, reads in zip(groups, generator))


def read_bam_into_groups(bam, reference, groups, group_tag='ME'):
    """

    :param bam:
    :param reference:
    :param groups:
    :param group_tag:
    :return:
    """
    read_groups = zip(read_bam_into_stranded_groups(bam, reference, '+', groups, group_tag='ME'),
                      read_bam_into_stranded_groups(bam, reference, '-', groups, group_tag='ME'))
    return (ReadGroup.append(forward, reverse) for forward, reverse in read_groups)

if __name__ == '__main__':
    pass
