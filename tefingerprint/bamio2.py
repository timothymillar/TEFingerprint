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
    dictionary = {contig.Header(*key): deque() for key in keys}

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











