#! /usr/bin/env python

import os
import re
import pysam
import numpy as np


def _read_bam_strings(bam, reference, strand):
    """
    Read  strings from an indexed Bam file.

    :param bam: The path to an indexed sorted Bam file
    :type bam: str
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
    return np.array(pysam.view(*flag, bam, reference).splitlines())


def _read_bam_reference_lengths(bam):
    """
    Read lengths of references from a Bam file

    :param bam: The path to an indexed sorted Bam file
    :type bam: str

    :return: A dictionary of the form {<reference>: <length>}
    :rtype: dict[str, int]
    """
    with pysam.AlignmentFile(bam, 'rb') as bam:
        references = bam.header['SQ']
        references = {r['SN']: r['LN'] for r in references}
        return references


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


def _parse_read_loci(strings):
    """
    Parses a collection of SAM formatted strings into a tuple generator.

    :param strings: A collection of SAM formatted strings
    :type strings: iterable[str]
    :param strand: Strand ('+' or '-') of all reads
    :type strand: str

    :return: An iterable of mapped SAM read positions and names
    :rtype: generator[(int, int, str, str)]
    """
    for string in strings:
        attr = string.split("\t")
        name = str(attr[0])
        start = int(attr[3])
        length = _cigar_mapped_length(attr[5])
        stop = start + length - 1  # 1 based indexing used in SAM format
        yield start, stop, name



def _bam_strand_read_loci(bam, reference, strand, categories, tag='ME'):
    source = os.path.basename(bam)
    strings = _read_bam_strings(bam, reference, strand)
    tag = '\t' + tag + ':[Zi]:'
    tags = np.array([re.split(tag, s)[1].split('\t')[0] for s in strings])
    for category in categories:
        category_strings = strings[tags.astype('U{0}'.format(len(category))) == category]
        for start, stop, name in _parse_read_loci(category_strings):
            yield reference, strand, start, stop, name, category, source


def bam_read_loci(bam, reference, categories, tag='ME'):
    for locus in _bam_strand_read_loci(bam, reference, '+', categories, tag=tag):
        yield locus
    for locus in _bam_strand_read_loci(bam, reference, '-', categories, tag=tag):
        yield locus
