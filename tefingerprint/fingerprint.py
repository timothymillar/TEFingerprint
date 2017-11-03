#! /usr/bin/env python

import numpy as np
from collections import Counter
from functools import reduce
from tefingerprint import utils
from tefingerprint import bamio2
from tefingerprint import loci2
from tefingerprint.cluster import core_distances as _core_distances


def count(bins, reads, trim=True, n_common_elements=0):
    sources = np.array(list({ctg.header.source for ctg in list(reads.contigs())}))
    sources.sort()

    # munging dtypes for counts
    dtype_element_count = np.dtype([('name', 'O'),
                                    ('count', np.int64)])
    dtype_elements = np.dtype([(str(i), dtype_element_count) for i in range(n_common_elements)])
    dtype_sample_count = np.dtype([('name', 'O'),
                                   ('count', np.int64),
                                   ('element', dtype_elements)])
    if n_common_elements == 0:
        dtype_sample_count = utils.remove_dtype_field(dtype_sample_count, 'element')
    dtype_samples = np.dtype([(str(i), dtype_sample_count) for i, _ in enumerate(sources)])
    dtype_new = np.dtype([('median', np.int64),
                          ('edge', np.int64),
                          ('sample', dtype_samples)])

    # new bins based on previous with additional slots for counts
    bins = bins.map(lambda x: loci2.add_field(x, dtype_new))

    # loop through bins
    for bin_contig in bins.contigs():

        # headers of read contigs to count
        read_headers = [bin_contig.header.mutate(source=source) for source in sources]

        # names of samples
        for i, header in enumerate(read_headers):
            bin_contig.loci['sample'][str(i)]['name'] = header.source

        # iterate through bins
        for locus in bin_contig.loci:

            # boolean masks of contained reads for each read set
            read_masks = [np.logical_and(reads[header].loci['tip'] >= locus['start'],
                                         reads[header].loci['tip'] <= locus['stop']) for header in read_headers]

            # read loci data of contained reads for each read set
            read_loci = [reads[header].loci[mask] for header, mask in zip(read_headers, read_masks)]

            # counts of reads
            for i, r in enumerate(read_loci):

                # total count of reads from each sample
                locus['sample'][str(i)]['count'] = len(r)

                # most common elements per sample per locus
                if n_common_elements > 0:
                    for j, pair in enumerate(Counter(r['element']).most_common(n_common_elements)):
                        locus['sample'][str(i)]['element'][str(j)] = pair

            # find median of cluster
            combined_read_tips = reduce(np.append, (tips['tip'] for tips in read_loci))
            combined_read_tips.sort()
            locus['median'] = np.median(combined_read_tips)

            # find edge of 90% core
            core_dists = _core_distances(combined_read_tips, int(len(combined_read_tips) * 0.99))
            if bin_contig.header.strand == '+':
                locus['edge'] = combined_read_tips[np.where(core_dists == np.min(core_dists))[0][-1]]
            elif bin_contig.header.strand == '-':
                locus['edge'] = combined_read_tips[np.where(core_dists == np.min(core_dists))[0][0]]

            # trim the potentially buffered locus to the first and last read tips
            if trim:
                locus['start'] = np.min(combined_read_tips)
                locus['stop'] = np.max(combined_read_tips)
    return bins


def _join_sorter(forward, reverse):
    dtype_sort = np.dtype([('value', np.int64), ('strand', np.int8), ('index', np.int64)])
    f = np.empty(len(forward), dtype=dtype_sort)
    f['value'] = forward
    f['strand'] = 0
    f['index'] = np.arange(0, len(forward))

    r = np.empty(len(reverse), dtype=dtype_sort)
    r['value'] = reverse
    r['strand'] = 1
    r['index'] = np.arange(0, len(reverse))

    clusters = np.append(f, r)
    clusters.sort(order=('value', 'strand'))

    prev = None
    for clust in clusters:
        if prev is None:
            if clust['strand'] == 1:
                # reverse cluster can't be paired
                yield (None, clust['index'])
            else:
                # store forward cluster as prev
                prev = clust
        else:
            if clust['strand'] == 0:
                # both forward so store cluster as prev
                yield (prev['index'], None)
                prev = clust
            else:
                # cluster is reverse and prev is forward
                yield (prev['index'], clust['index'])
                prev = None
    if prev is not None:
        # final cluster is un paired
        yield (prev['index'], None)


def join_clusters(clusters):
    """


    :param clusters: a set of contigs that are not un-stranded (i.e. are on the + or - strand)
    :return:
    """
    joint_clusters = loci2.ContigSet()

    dtype_joint = np.dtype([("start", np.int64),
                            ("stop", np.int64),
                            ("strand", '<U1'),
                            ("paired", np.int8),
                            ("forward", clusters.dtype_loci),
                            ("reverse", clusters.dtype_loci)])

    # new headers based on old but un stranded
    new_headers = {h.mutate(strand='.') for h in clusters.headers}

    for header in new_headers:

        # get forward and reverse loci for this key
        forward = clusters[header.mutate(strand='+')]
        reverse = clusters[header.mutate(strand='-')]

        # sort them into pairs based on median
        pairs = list(_join_sorter(forward['median'], reverse['median']))

        # create an empty array based on the pairs
        joint_loci = np.empty(len(pairs), dtype=dtype_joint)

        for i, (f, r) in enumerate(pairs):
            if f is not None:
                joint_loci[i]['forward'] = forward[f]
                joint_loci[i]['start'] = joint_loci[i]['forward']['edge']
            if r is not None:
                joint_loci[i]['reverse'] = reverse[r]
                joint_loci[i]['stop'] = joint_loci[i]['reverse']['edge']
            if f is not None and r is not None:
                joint_loci[i]['paired'] = 1

        # fill in blanks based on present data
        missing_forward = joint_loci['start'] == 0
        missing_reverse = joint_loci['stop'] == 0

        joint_loci['start'][missing_forward] = joint_loci['stop'][missing_forward] - 1
        joint_loci['stop'][missing_reverse] = joint_loci['start'][missing_reverse] + 1

        # strand based on 1 or 2 clusters
        joint_loci['strand'] = '.'
        joint_loci['strand'][missing_forward] = '-'
        joint_loci['strand'][missing_reverse] = '+'

        joint_clusters.add(loci2.Contig(header, joint_loci))

    return joint_clusters

