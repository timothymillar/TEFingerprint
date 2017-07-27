#! /usr/bin/env python

import os
import pysam
from itertools import product
from collections import deque


def _get_reference_lengths(bam):
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


def _get_references(*args):
    """
    Return references names and intervals for a set of indexed bam files.
    Bam files must have the same reverence names and lengths.

    :param args: paths to one or more bam files
    :type args: [str]

    :return: reference names and 1-based intervals
    :rtype: generator[str, int, int]
    """
    bam_references = [_get_reference_lengths(bam) for bam in args]  # get references of all bams
    assert all(ref == bam_references[0] for ref in bam_references)  # assert that references are the same for all bams
    for name, length in bam_references[0].items():  # yield name and interval of references
        yield name, 1, length


def extract_bam_reference_strings(*args):
    """
    Return references names and intervals for a set of indexed bam files.
    Bam files must have the same reverence names and lengths.

    :param args: paths to one or more bam files
    :type args: [str]

    :return: reference names and 1-based intervals with the form 'name:start-stop'
    :rtype: list[str]
    """
    references = list(_get_references(*args))
    strings = ['{0}:{1}-{2}'.format(*reference) for reference in references]
    return strings


def _parse_reference_strings(bams, reference_strings=None):
    """
    Take a set of bams and a set of reference strings and return reference tuples.
    Reference strings may be a reference name or name and interval with the form 'name:start-stop'
    If no references are specified then all references recorded are returned.

    :param bams: paths to one or more bam files
    :type bams: iterable[str]
    :param reference_strings:
    :type reference_strings: iterable[str]

    :return: reference names and intervals suitable for use with pysam.AlignmentFile.fetch()
    :rtype: generator[(str, int, int)]
    """
    valid_references = {i[0]: (i[1], i[2]) for i in _get_references(*bams)}

    # return all valid references if none specified
    if reference_strings is None:
        for name, interval in valid_references.items():
            yield name, interval[0], interval[1]

    # else check and parse strings
    else:
        for reference in reference_strings:
            if ':' in reference:
                # assume reference contains a string representation of an interval
                name, interval = reference.split(':')
                assert name in valid_references.keys()
                start, stop = interval.split('-')
                yield name, int(start), int(stop)
            else:
                # assume reference name only so use the full interval from the dict
                assert reference in valid_references.keys()
                interval = valid_references[reference]
                yield reference, interval[0], interval[1]


def _create_loci_keys(bams, references, categories):
    """
    Returns all possible loci keys for a selection of bam files, references, categories

    :param bams: list of bam files
    :type bams: list[str]
    :param references: list of references
    :type references: list[(str, int, int)]
    :param categories: list of categories
    :type categories: list[str]

    :return: all possible loci keys
    :rtype: generator[(str, str, str, str)]
    """
    strands = ['+', '-']
    for bam, reference, category, strand in product(bams, references, categories, strands):
        yield '{0}:{1}-{2}'.format(*reference), strand, category, os.path.basename(bam)


def _parse_reads(bam, references, categories, quality=0, tag='ME'):
    """
    Parses reads from an indexed bam file into a format suitable for TEFingerprint.
    Reads are excluded if they are not in one of the specified reference-slices, categories or bellow the
    minimum mapping quality

    :param bam: path to an indexed bam file
    :type bam: str
    :param references: list of references with intervals
    :type references: list[str]
    :param categories: list of categories of reads i.e. transposon super families
    :type categories: list[str]
    :param quality: minimum mapping quality for reads
    :type quality: int
    :param tag: sam-tag containing the category information about each read (defaults to 'ME')
    :type tag: str

    :return: key-locus pairs that match the above criteria
    :rtype: generator[((str, str, str, str), (int, int, str))]
    """
    alignment = pysam.AlignmentFile(bam, 'r')
    for reference in references:
        for read in alignment.fetch(*reference):
            # filter out reads below specified mapping quality
            if read.mapping_quality >= quality:
                # filter out reads that are not in a specified category
                category_matches = tuple(filter(lambda x: read.get_tag(tag).startswith(x), categories))
                if category_matches:
                    key = ('{0}:{1}-{2}'.format(*reference),
                           '-' if read.is_reverse else '+',
                           max(category_matches, key=len),
                           os.path.basename(bam))
                    # TODO: use pysam 0-based half-open indices or sam 1-based closed indices?
                    locus = (read.blocks[0][0] + 1,  # adjust pysam indexing to sam indexing
                             read.blocks[-1][-1],
                             read.qname)
                    yield key, locus
                else:
                    pass
            else:
                pass


def _group_reads_by_keys(keys, reads):
    """
    Groups key-locus pairs by their keys

    :param keys: keys of the form (reference, strand, category, source)
    :type keys: iterable[str, str, str, str]
    :param reads: key-locus pairs of the form ((reference, strand, category, source), (start, stop, name))
    :type reads: iterable[(iterable[(str, str, str, str)], iterable[(int, int, str)])]

    :return: key-loci pairs
    :rtype: generator[((str, str, str, str), collections.Deque[(int, int, str)])]
    """
    queues = {key: deque() for key in keys}
    for key, locus in reads:
        queues[key].append(locus)

    for key, loci in queues.items():
        yield key, loci


def _process_bams(bams, references, categories, quality=0, tag='ME'):
    """
    Processes a set of bams to extract and group reads into a format suitable for TEFingerprint.

    :param bams: list of paths to bam files
    :type bams: list[str]
    :param references: list of references with intervals
    :type references: list[(str, int, int)]
    :param categories: list of categories of reads i.e. transposon super families
    :type categories: list[str]
    :param quality: minimum mapping quality
    :type quality: int
    :param tag: sam-tag containing the category information about each read (defaults to 'ME')
    :type tag: str

    :return: key-loci pairs
    :rtype: generator[((str, str, str, str), collections.Deque[(int, int, str)])]
    """
    keys = _create_loci_keys(bams, references, categories)
    for bam in bams:
        reads = _parse_reads(bam, references, categories, quality=quality, tag=tag)
        blocks = _group_reads_by_keys(keys, reads)
        for block in blocks:
            yield block


def extract_bam_reads(bams, categories, references=None, quality=0, tag='ME'):
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
    :param quality: minimum mapping quality
    :type quality: int
    :param tag: the two letter sam tag containing the transposon associated with each read
    :type tag: str

    :return: a generator of category tuples and loci tuple generators
    :rtype: generator[((str, str, str, str), generator[((int, int, str))])]
    """
    if isinstance(bams, str):
        bams = [bams]

    if isinstance(categories, str):
        categories = [categories]

    if isinstance(references, str):
        references = [references]

    references = list(_parse_reference_strings(bams, reference_strings=references))

    return _process_bams(bams, references, categories, quality=quality, tag=tag)

if __name__ == '__main__':
    pass
