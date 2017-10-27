#! /usr/bin/env python

import os
import pysam
import numpy as np
from itertools import product
from collections import deque
from tefingerprint import loci2


def extract_bam_references(bam):
    with pysam.AlignmentFile(bam, 'rb') as bam:
        references = bam.header['SQ']
        return (r['SN'] for r in references)


def _parse_bam_read_data(bam, reference, quality=0, tags=None):
    gen = pysam.AlignmentFile(bam, 'rb')
    for read in gen.fetch(region=reference):
        if read.mapping_quality >= quality:
            data = {
                'reference':reference.split(':')[0],
                'strand': '-' if read.is_reverse else '+',
                'start': read.blocks[0][0] + 1,  # adjust pysam indexing to samtools indexing
                'stop': read.blocks[-1][-1],
                'source': os.path.basename(bam)
            }
            if tags is not None:
                data.update({tag: read.get_tag(tag) for tag in tags})
            yield data


def extract_informative_read_tips(bams, references, categories, quality=0, tag='ME'):
    """"""
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
    dictionary = {loci2.Header(*key): deque() for key in keys}

    for bam in bams:
        for reference in references:
            for read in _parse_bam_read_data(bam, reference, quality=quality, tags=[tag]):

                # match to a category
                category_matches = tuple(filter(lambda x: read[tag].startswith(x), categories))

                # only include reads for specified categories
                if category_matches:

                    # longest matching category is the best category
                    category = max(category_matches, key=len)

                    # read header
                    header = loci2.Header(reference=read['reference'],
                                          strand=read['strand'],
                                          category=category,
                                          source=read['source'])

                    # append loci data to que
                    tip = read['start'] if read['strand'] == '-' else read['stop']
                    dictionary[header].append((tip, read[tag]))

    dtype = np.dtype([('tip', np.int64), ('element', 'O')])
    return loci2.ContigSet(*(loci2.Contig(header, np.array(data, dtype=dtype))
                             for header, data in dictionary.items()))


def _parse_bam_anchor_insert_data(bam, reference, max_size, quality=0):
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
                    and r.mpos > r.blocks[-1][-1] \
                    and r.mpos - r.blocks[-1][-1] < max_size:
                yield r.blocks[-1][-1], r.mpos + 1


def extract_anchor_intervals(bams, references, max_size, quality=0):
    if isinstance(bams, str):
        bams = [bams]

    if isinstance(references, str):
        references = [references]

    jobs = product(bams, references)
    dtype = np.dtype([('start', np.int64), ('stop', np.int64)])
    intervals = loci2.ContigSet()

    for bam, reference in jobs:
        header = loci2.Header(reference=reference.split(':')[0], source=os.path.basename(bam))
        data = _parse_bam_anchor_insert_data(bam, reference, max_size, quality=quality)
        intervals.add(loci2.unions(loci2.Contig(header=header, array=np.fromiter(data, dtype=dtype))))

    return intervals


def extract_gff_intervals(gff, references, categories):
    if isinstance(references, str):
        references = [references]

    if isinstance(categories, str):
        categories = [categories]

    source = os.path.basename(gff)
    references = [reference.split(':')[0] for reference in references]

    keys = product(references, categories)
    dictionary = {loci2.Header(reference=key[0],
                               category=key[1],
                               source=source): deque() for key in keys}

    with open(gff) as infile:
        for line in infile:
            line = line.split('\t')

            # match to reference:
            if line[0] in references:

                # match to a category
                category_matches = tuple(filter(lambda x: line[2].startswith(x), categories))

                # only include reads for specified categories
                if category_matches:

                    # longest matching category is the best category
                    category = max(category_matches, key=len)

                    header = loci2.Header(reference=line[0],
                                          category=category,
                                          source=source)

                    dictionary[header].append((line[3], line[4], line[2]))

    dtype = np.dtype([('start', np.int64), ('stop', np.int64), ('element', '<O')])
    return loci2.ContigSet(*(loci2.Contig(header, np.array(data, dtype=dtype))
                             for header, data in dictionary.items()))








