#! /usr/bin/env python

import os
import re
import pysam
import numpy as np
from itertools import product


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


def extract_bam_references(*args):
    """
    Returns a list of references (including their range) from the header(s) of one or more bam file(s).
    If multiple bam file are specified, they must all have identical references (including their range).

    :param args: The path(s) to one or more indexed sorted Bam file(s)
    :type args: iterable[str]

    :return: A list of strings with the structure 'name:minimum-maximum'
    :rtype: list[str]
    """
    bam_refs = [_read_bam_reference_lengths(bam) for bam in args]
    assert all(ref == bam_refs[0] for ref in bam_refs)
    references = ['{0}:0-{1}'.format(k, v) for k, v in bam_refs[0].items()]
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
    :rtype: generator[(int, int, str)]
    """
    for string in strings:
        attr = string.split("\t")
        name = str(attr[0])
        start = int(attr[3])
        length = _cigar_mapped_length(attr[5])
        stop = start + length - 1  # 1 based indexing used in SAM format
        yield start, stop, name


def _read_bam_ref_strand_loci(bam, reference, strand, categories, tag='ME'):
    """
    Read a single strand of a single reference from a bam file and categories by their associated transposon.
    If a reference is specified by name only, its range will be extracted from the bam and appended.

    :param bam: path to a bam file
    :type bam: str
    :param reference: target reference of format 'name' or 'name:minimum-maximum'
    :type reference: str
    :param strand: target strand '+' or '-'
    :type strand: str
    :param categories: a list of transposon categories
    :type categories: list[str]
    :param tag: the two letter sam tag containing the transposon associated with each read
    :type tag: str

    :return: a generator of category tuples and loci tuple generators
    :rtype: generator[((str, str, str, str), generator[((int, int, str))])]
    """
    source = os.path.basename(bam)
    strings = _read_bam_strings(bam, reference, strand)
    if ':' in reference:
        pass
    else:
        reference += ':0-{0}'.format(_read_bam_reference_lengths(bam)[reference])
    tag = '\t' + tag + ':[Zi]:'
    tags = np.array([re.split(tag, s)[1].split('\t')[0] for s in strings])
    for category in categories:
        category_strings = strings[tags.astype('U{0}'.format(len(category))) == category]
        loci = _parse_read_loci(category_strings)
        yield (reference, strand, category, source), loci


def _read_bam_ref_loci(bam, reference, categories, tag='ME'):
    """
    Read both strands of a single reference from a bam file and categories by their associated transposon.
    If a reference is specified by name only, its range will be extracted from the bam and appended.

    :param bam: path to a bam file
    :type bam: str
    :param reference: target reference of format 'name' or 'name:minimum-maximum'
    :type reference: str
    :param categories: a list of transposon categories
    :type categories: list[str]
    :param tag: the two letter sam tag containing the transposon associated with each read
    :type tag: str

    :return: a generator of category tuples and loci tuple generators
    :rtype: generator[((str, str, str, str), generator[((int, int, str))])]
    """
    for block in _read_bam_ref_strand_loci(bam, reference, '+', categories, tag=tag):
        yield block
    for block in _read_bam_ref_strand_loci(bam, reference, '-', categories, tag=tag):
        yield block


def extract_bam_reads(bams, categories, references=None, tag='ME'):
    """
    Extract reads from one or more bam file(s) and categories reads by their reference, strand, associated transposon,
    and source bam file.
    If a reference is specified by name only, its range will be extracted from the bam and appended.

    :param bams: path(s) to a one or more bam file(s)
    :type bams: str | list[str]
    :param references: target reference(s) of format 'name' or 'name:minimum-maximum'
    :type references: str | list[str]
    :param categories: target transposon categories
    :type categories: str | list[str]
    :param tag: the two letter sam tag containing the transposon associated with each read
    :type tag: str

    :return: a generator of category tuples and loci tuple generators
    :rtype: generator[((str, str, str, str), generator[((int, int, str))])]
    """
    if isinstance(bams, str):
        bams = [bams]

    if isinstance(categories, str):
        categories = [categories]

    if references is None:
        references = extract_bam_references(*bams)
    if isinstance(references, str):
        references = [references]

    # run jobs
    jobs = product(bams, references)
    for bam, reference in jobs:
        for block in _read_bam_ref_loci(bam, reference, categories, tag=tag):
            yield block
