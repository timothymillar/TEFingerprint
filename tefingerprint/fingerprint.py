#! /usr/bin/env python

import numpy as np
from collections import Counter
from functools import reduce
from tefingerprint import utils
from tefingerprint import interval
from tefingerprint import bamio2
from tefingerprint import loci2
from tefingerprint.cluster import core_distances as _core_distances


def count_reads(clusters, reads, trim=True, n_common_elements=0):
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
    clusters = clusters.map(lambda x: loci2.add_field(x, dtype_new))

    # loop through bins
    for bin_contig in clusters.contigs():

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
    return clusters


def _known_insertion_matcher(contig, known_insertions, buffer=0):
    if contig.header.strand == '+':
        for cluster in contig.loci:
            mask = interval.singular_intersects_lower_bound(cluster['median'],
                                                            cluster['stop'] + buffer,
                                                            known_insertions)
            if np.any(mask):
                yield known_insertions['element'][mask][0]
            else:
                yield ''
    elif contig.header.strand == '-':
        for cluster in contig.loci:
            mask = interval.singular_intersects_upper_bound(cluster['start'] - buffer,
                                                            cluster['median'],
                                                            known_insertions)
            if np.any(mask):
                yield known_insertions['element'][mask][-1]
            else:
                yield ''


def match_known_insertions(clusters, known_insertions, buffer=0):

    matched = loci2.ContigSet()

    # make known insertion headers un-stranded and drop origin file
    known_insertions = known_insertions.map(lambda x: loci2.mutate_header(x, strand='.', source=None))

    # loop through contigs
    for contig in clusters.contigs:

        # get relevant known insertions
        known = known_insertions[contig.header.mutate(strand='.')]

        matches = np.array(list(_known_insertion_matcher(contig, known, buffer=buffer)),
                           dtype=np.dtype([('known', '<O')]))

        matched.add(loci2.Contig(contig.header, utils.bind_arrays(contig.loci, matches)))

    return matched


def _pair_clusters(forward, reverse, buffer=0):
    """sorts clusters into pairs with knowledge of known insertions or distances"""
    dtype_sort = np.dtype([('start', np.int64),
                           ('stop', np.int64),
                           ('median', np.int64),
                           ('known', '<O'),
                           ('strand', '<U1'),
                           ('index', np.int64)])

    f = np.empty(len(forward), dtype=dtype_sort)
    f['stop'] = forward.loci['stop']
    f['median'] = forward.loci['median']
    f['known'] = forward.loci['known']
    f['strand'] = '+'
    f['index'] = np.arange(0, len(forward))

    r = np.empty(len(reverse), dtype=dtype_sort)
    r['start'] = reverse.loci['start']
    r['median'] = reverse.loci['median']
    r['known'] = reverse.loci['known']
    r['strand'] = '-'
    r['index'] = np.arange(0, len(reverse))

    clusters = np.append(f, r)
    clusters.sort(order=('median', 'strand'))

    prev = None
    for clust in clusters:
        if prev is None:
            if clust['strand'] == '-':
                # reverse cluster can't be paired
                yield (None, clust['index'])
            else:
                # store forward cluster as prev
                prev = clust
        else:
            if clust['strand'] == '+':
                # both forward so store cluster as prev
                yield (prev['index'], None)
                prev = clust
            else:
                # cluster is reverse and prev is forward
                if prev['known'] != '' and prev['known'] == clust['known']:
                    # clusters match same known element
                    yield (prev['index'], clust['index'])
                elif prev['stop'] + buffer >= clust['start'] - buffer:
                    # they are within reasonable range of one another (2 * buffer)
                    yield (prev['index'], clust['index'])
                else:
                    # they are not a good match
                    yield (prev['index'], None)
                    yield (None, clust['index'])
                prev = None
    if prev is not None:
        # final cluster is un-paired
        yield (prev['index'], None)



def join_clusters(clusters, known_insertions):
    """


    :param clusters: a set of contigs that are not un-stranded (i.e. are on the + or - strand)
    :return:
    """
    joint_clusters = loci2.ContigSet()

    dtype_join = utils.append_dtypes(clusters.dtype_loci,
                                     np.dtype([("ID", "<O"),
                                               ("paired", np.int8),
                                               ("known", '<O')]))

    # new headers based on old but un-stranded
    new_headers = {h.mutate(strand='.') for h in clusters.headers}

    # make known insertion headers un-stranded and drop origin file
    known_insertions = known_insertions.map(lambda x: loci2.mutate_header(x, strand='.', source=None))

    # template for creating insertion IDs
    insertion_id_template = '{reference}_{start}_{stop}_{category}'

    for header in new_headers:

        # get forward and reverse loci for this key
        forward = clusters[header.mutate(strand='+')]
        reverse = clusters[header.mutate(strand='-')]

        # get known insertions of the same reference and family
        known = known_insertions[header.mutate(source=None)]

        # sort them into pairs based on median
        pairs = _join_sorter(forward['median'], reverse['median'])

        # create arrays for the new data
        forward_join_data = np.empty(len(forward), dtype=dtype_join)
        reverse_join_data = np.empty(len(reverse), dtype=dtype_join)
        for f, r in pairs:
            if f is not None and r is not None:
                insertion_id = insertion_id_template.format(header.reference, f['value'], r['value'], header.category)
                forward_join_data[f['index']]["ID"] = insertion_id
                reverse_join_data[r['index']]["ID"] = insertion_id
            elif f is not None:
                insertion_id = insertion_id_template.format(header.reference, f['value'], '.', header.category)
                forward_join_data[f['index']]["ID"] = insertion_id
            elif r is not None:
                insertion_id = insertion_id_template.format(header.reference, '.', r['value'], header.category)
                reverse_join_data[r['index']]["ID"] = insertion_id

        # combine existing data with join data and add to new contig set
        joint_clusters.add(loci2.Contig(header.mutate(strand='+'),
                                        utils.bind_arrays(forward, forward_join_data)))
        joint_clusters.add(loci2.Contig(header.mutate(strand='-'),
                                        utils.bind_arrays(reverse, reverse_join_data)))

    return joint_clusters

