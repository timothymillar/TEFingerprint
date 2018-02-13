#! /usr/bin/env python3

import os
import pysam
import numpy as np
from itertools import product
from collections import deque
from tefingerprint import loci
from tefingerprint import interval
from tefingerprint.gff import decode_column, decode_attribute
from tefingerprint.compression import zopen


def extract_references_from_bam(bam):
    """
    Extract names of references from a bam file.

    :param bam: path to a bam file
    :type bam: str

    :return: list of bam reference names
    :rtype: list[str]
    """
    with pysam.AlignmentFile(bam, 'rb') as bam:
        references = bam.header['SQ']
        return [r['SN'] for r in references]


def extract_references_from_bams(*args):
    """
    Extract names of references from a multiple bam files.

    :param args: paths to bam files
    :type args: list[str]

    :return: list of bam reference names
    :rtype: list[str]
    """
    reference_lists = [extract_references_from_bam(bam) for bam in args]
    for l in reference_lists:
        assert l == reference_lists[0]
    return reference_lists[0]


def _extract_bam_read_data(bam, reference, quality=0, tags=None):
    """
    Extract and parse read data from a single reference of a
    single bam file.

    Read positions are returned in sam format (one-based inclusive).
    Only extracts reads aligned to the specified reference.


    :param bam: Path to a bam file
    :type bam: str
    :param reference: Name of reference
    :type reference: str
    :param quality: Minimum mapping quality of reads
    :type quality: int
    :param tags: List of sam tags to include in output
    :type tags: list[str]

    :return: iterable of dictionaries containing read metadata
    :rtype: iterable[dict]
    """
    gen = pysam.AlignmentFile(bam, 'rb')
    for read in gen.fetch(region=reference):
        if read.mapping_quality >= quality:
            data = {
                'reference':reference.split(':')[0],
                'strand': '-' if read.is_reverse else '+',
                'start': read.blocks[0][0] + 1,  # adjust for pysam indexing
                'stop': read.blocks[-1][-1],
                'source': os.path.basename(bam)
            }
            if tags is not None:
                data.update({tag: read.get_tag(tag) for tag in tags})
            yield data


def extract_informative_read_tips(bams,
                                  references,
                                  categories,
                                  quality=0,
                                  tag='ME'):
    """
    Extract the tips of 'informative' reads from one or more bam files.

    Informative reads are those that flank potential transposon
    insertions.
    The specific element (mate element) that each read is linked to
    should be stored using a sam tag which is 'ME' by default.
    Reads are categorised by transposon (super-)families by matching
    family names to the start of each reads mate-element name.

    :param bams: Path(s) to one or more bam files
    :type bams: str | list[str]
    :param references: Name(s) of one or more bam references
    :type references: str | list[str]
    :param categories: Name(s) of one or more transposon (super-)families
    :type categories: str | list[str]
    :param quality: Minimum mapping quality of reads
    :type quality: int
    :param tag: Sam tag containing each reads mate element name
    :type tag: str

    :return: A set of contigs of read tips categorised by reference,
        strand, category (family), and source (bam file name)
    :rtype: :class:`loci2.ContigSet`
    """
    if isinstance(bams, str):
        bams = [bams]

    if isinstance(references, str):
        references = [references]

    if isinstance(categories, str):
        categories = [categories]

    keys = product([ref.split(':')[0] for ref in references],
                   ['+', '-'],
                   categories,
                   [os.path.basename(bam) for bam in bams])
    dictionary = {loci.Header(*key): deque() for key in keys}

    for bam in bams:
        for reference in references:
            for read in _extract_bam_read_data(bam,
                                               reference,
                                               quality=quality,
                                               tags=[tag]):

                # match to a category
                category_matches = tuple(filter(lambda x:
                                                read[tag].startswith(x),
                                                categories))

                # only include reads for specified categories
                if category_matches:

                    # longest matching category is the best category
                    category = max(category_matches, key=len)

                    # read header
                    header = loci.Header(reference=read['reference'],
                                         strand=read['strand'],
                                         category=category,
                                         source=read['source'])

                    # append loci data to que
                    tip = read['start'] if \
                        read['strand'] == '-' else \
                        read['stop']
                    dictionary[header].append((tip, read[tag]))

    dtype = np.dtype([('tip', np.int64), ('element', 'O')])
    return loci.ContigSet(*(loci.Contig(header, np.array(data,
                                                         dtype=dtype))
                            for header, data in dictionary.items()))


def _extract_bam_anchor_insert_data(bam, reference, quality=0):
    """
    Extract 'anchor' read inserts from a single reference of a
    single bam file.

    :param bam: Path to a bam file
    :type bam: str
    :param reference: Name of reference
    :type reference: str
    :param quality: Minimum mapping quality of reads
    :type quality: int

    :return: iterable of dictionaries containing read metadata
    :rtype: iterable[dict]
    """
    with pysam.AlignmentFile(bam, 'rb') as bam:
        for r in bam.fetch(region=reference):
            if not r.is_reverse \
                    and not r.is_unmapped \
                    and r.is_paired \
                    and r.is_read1 \
                    and r.mate_is_reverse \
                    and not r.mate_is_unmapped \
                    and r.reference_id == r.next_reference_id \
                    and r.mapping_quality >= quality \
                    and r.mpos > r.blocks[-1][-1]:
                yield r.blocks[-1][-1], r.mpos + 1


def extract_anchor_intervals(bams,
                             references,
                             known_transposons,
                             insert_size,
                             quality=0):
    """
    Extract 'anchor' read inserts from one or more bam files

    Anchor reads are paired reads in which neither has been mapped
    to a known transposon.
    The pair has then been mapped to a reference genome.
    Assuming that the insert size of the pair is smaller than the
    length of a transposon, the insert can be used to indicate a
    section of the samples genome in which there are no transposons
    on at least one allele.
    This can be used to infere heterozygousity of transposons
    insertions.

    Known transposon inserts from the reference genome are required for
    checking that anchor inserts overlapping these transposon are of
    a sensible length.

    Anchor reads are compressed to their interval unions for efficiency.

    :param bams: Path(s) to one or more bam files
    :type bams: str | list[str]
    :param references: Name(s) of one or more bam references
    :type references: str | list[str]
    :param known_transposons: Transposons known from the reference genome
    :type known_transposons: :class:`loci2.ContigSet`
    :param insert_size: Read pair insert size
    :type insert_size: int
    :param quality: Minimum maping quality of anchor reads
    :type quality: int

    :return: A set of contigs of unions of anchor inserts categorised
        by reference, strand, and source (bam file name)
    :rtype: :class:`loci2.ContigSet`
    """
    if isinstance(bams, str):
        bams = [bams]

    if isinstance(references, str):
        references = [references]

    # simplify known transposon headers for comparison
    known_transposons = known_transposons.map(lambda x:
                                              loci.mutate_header(x,
                                                                 strand='.',
                                                                 category=None,
                                                                 source=None),
                                              append_duplicate_headers=True)

    jobs = product(bams, references)
    dtype = np.dtype([('start', np.int64), ('stop', np.int64)])
    intervals = loci.ContigSet()

    for bam, reference in jobs:
        header = loci.Header(reference=reference.split(':')[0],
                             source=os.path.basename(bam),
                             strand='.')
        anchors = np.fromiter(_extract_bam_anchor_insert_data(bam,
                                                              reference,
                                                              quality=quality),
                              dtype=dtype)
        anchor_lengths = interval.lengths(anchors)

        # calculate lengths on known tranposons within each anchor interval
        reference_name = reference.split(':')[0]
        local_tes_header = loci.Header(reference=reference_name,
                                       strand='.')
        local_tes = known_transposons[local_tes_header]
        contained_te_lengths = interval.length_of_contains(anchors,
                                                           local_tes.loci)

        # filter anchors based on insert size
        adjusted_anchor_lengths = anchor_lengths - contained_te_lengths
        anchors = anchors[adjusted_anchor_lengths <= insert_size]

        # use unions of filtered anchors as loci
        intervals.add(loci.unions(loci.Contig(header=header, array=anchors)))

    return intervals


def extract_gff_intervals(gff, references, categories):
    """
    Extract known transposon intervals from a gff anotation
    file.

    :param gff: Path to a gff file of transposon anotations
    :type gff: str
    :param references: Name(s) of one or more bam references
    :type references: str | list[str]
    :param categories: Name(s) of one or more transposon (super-)families
    :type categories: str | list[str]

    :return: A set of contigs of read tips categorised by reference,
        strand, category (family), and source (bam file name)
    :rtype: :class:`loci2.ContigSet`
    """
    if isinstance(references, str):
        references = [references]

    if isinstance(categories, str):
        categories = [categories]

    source = os.path.basename(gff)
    references = [reference.split(':')[0] for reference in references]

    keys = product(references, categories)
    dictionary = {loci.Header(reference=key[0],
                              category=key[1],
                              source=source): deque() for key in keys}

    with zopen(gff, 'rb') as infile:
        for line in infile:
            line = line.decode().split('\t')

            # match to reference:
            reference = decode_column(line[0])
            if reference in references:

                # match to a category
                feature_type = decode_column(line[2])
                category_matches = tuple(filter(lambda x:
                                                feature_type.startswith(x),
                                                categories))

                # only include reads for specified categories
                if category_matches:

                    # longest matching category is the best category
                    category = max(category_matches, key=len)

                    header = loci.Header(reference=reference,
                                         category=category,
                                         source=source)

                    dictionary[header].append((int(line[3]),
                                               int(line[4]),
                                               feature_type))

    dtype = np.dtype([('start', np.int64),
                      ('stop', np.int64),
                      ('element', '<O')])
    return loci.ContigSet(*(loci.Contig(header, np.array(data, dtype=dtype))
                            for header, data in dictionary.items()))


if __name__ == "__main__":
    pass
