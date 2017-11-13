#! /usr/bin/env python

import numpy as np
from itertools import product
from functools import reduce
from collections import Counter
from multiprocessing import Pool
from tefingerprint import util
from tefingerprint import interval
from tefingerprint import bamio2
from tefingerprint import loci2


def fingerprint(bams,
                categories,
                references,
                minimum_reads,
                epsilon,
                minimum_epsilon=0,
                n_common_elements=0,
                hierarchical=True,
                fingerprint_buffer=0,
                join_distance=0,
                quality=0,
                transposon_tag='ME',
                annotation=None,
                cores=1):
    """
    Create a transposon fingerprint of one or more bam files.

    :param bams:
    :param categories:
    :param references:
    :param minimum_reads:
    :param epsilon:
    :param minimum_epsilon:
    :param n_common_elements:
    :param hierarchical:
    :param fingerprint_buffer:
    :param join_distance:
    :param quality:
    :param transposon_tag:
    :param annotation:
    :param cores:
    :return:
    """

    if isinstance(bams, str):
        bams = [bams]

    if isinstance(references, str):
        references = [references]

    if isinstance(categories, str):
        categories = [categories]

    if references is [None]:
        references = bamio2.extract_references_from_bams(*bams)

    jobs = product([bams],
                   [annotation],
                   [categories],
                   references,  # job per reference
                   [quality],
                   [transposon_tag],
                   [minimum_reads],
                   [epsilon],
                   [minimum_epsilon],
                   [n_common_elements],
                   [hierarchical],
                   [fingerprint_buffer],
                   [join_distance])

    result = loci2.ContigSet()

    if cores == 1:
        # run on a single process
        for job in jobs:
            result.update(_fingerprint_dispatch(*job).contigs())
    else:
        # create a pool of processes
        with Pool(cores) as pool:
            parts = pool.starmap(_fingerprint_dispatch, jobs)
            for part in parts:
                result.update(part.contigs())

    return result


def _fingerprint_dispatch(bams,
                          annotation,
                          categories,
                          reference,
                          quality,
                          transposon_tag,
                          minimum_reads,
                          epsilon,
                          minimum_epsilon,
                          n_common_elements,
                          hierarchical,
                          fingerprint_buffer,
                          join_distance):
    """dispatch a single job of a fingerprint"""

    # read informative reads
    reads = bamio2.extract_informative_read_tips(bams,
                                                 reference,
                                                 categories,
                                                 quality=quality,
                                                 tag=transposon_tag)

    # read known transposons if used
    if annotation is None:
        pass
    else:
        known = bamio2.extract_gff_intervals(annotation,
                                             reference,
                                             categories)

    # sort reads
    reads = reads.map(lambda x: loci2.sort(x, order='tip'))

    # cluster reads
    clusters = reads.map(lambda x:
                         loci2.cluster(x,
                                       'tip',
                                       minimum_reads,
                                       epsilon,
                                       minimum_epsilon=minimum_epsilon,
                                       hierarchical=hierarchical))

    # drop origin files and append
    clusters = clusters.map(lambda x:
                            loci2.mutate_header(x, source=None),
                            append_duplicate_headers=True)

    # union of clusters
    clusters = clusters.map(lambda x:
                            loci2.unions_buffered(x, fingerprint_buffer))

    # count reads in bins
    clusters = count_reads(clusters,
                           reads,
                           n_common_elements=n_common_elements)

    # match clusters to known elements if used
    if annotation is None:
        pass
    else:
        clusters = match_known_insertions(clusters,
                                          known,
                                          distance=join_distance)

    # join cluster pairs
    use_known_elements = False if annotation is None else True

    clusters = pair_clusters(clusters,
                             distance=join_distance,
                             use_known_elements=use_known_elements)
    return clusters


def count_reads(clusters, reads, trim=True, n_common_elements=0):
    """
    Counts the read tips contained in each cluster (locus).

    A total count of reads across all source (bam) files is calculated.
    The median read position amoung all reads is calculated.
    The n_common_elements parameter may be used to include counts of the 'n'
    most common elements found in each sample.
    By default, clusters are trimmed to the most distant reads they contain.

    Fields required in 'clusters':
        'start': int, 'stop': int, 'median': int

    Fields required in 'reads':
        'tip': int, 'element': str, 'source': str

    Fields appended to return value:
        'median': int, 'sample': ('name': str, 'count': int,
        'element': ( 'name': str, 'count': int))

    :param clusters: a collection of cluster loci (intervals)
    :type clusters: :class:`loci2.ContigSet`
    :param reads: a collection of read tips from one or more bam files
    :type reads: :class:`loci2.ContigSet`
    :param trim: if true the boundaries of each cluster are trimmed to
        the range of the reads they contain
    :type trim: bool
    :param n_common_elements: specifies how many of the most common
        elements to count for each sample
    :type n_common_elements: int

    :return: a collection of cluster loci (intervals) with read counts
    :rtype: :class:`loci2.ContigSet`
    """
    sources = np.array(list({ctg.header.source
                             for ctg in list(reads.contigs())}))
    sources.sort()

    # munging dtypes for counts
    dtype_element_count = np.dtype([('name', 'O'),
                                    ('count', np.int64)])
    dtype_elements = np.dtype([(str(i), dtype_element_count)
                               for i in range(n_common_elements)])
    dtype_sample_count = np.dtype([('name', 'O'),
                                   ('count', np.int64),
                                   ('element', dtype_elements)])
    if n_common_elements == 0:
        dtype_sample_count = util.remove_dtype_field(dtype_sample_count,
                                                     'element')
    dtype_samples = np.dtype([(str(i), dtype_sample_count)
                              for i, _ in enumerate(sources)])
    dtype_new = np.dtype([('median', np.int64),
                          ('sample', dtype_samples)])

    # new bins based on previous with additional slots for counts
    clusters = clusters.map(lambda x: loci2.add_field(x, dtype_new))

    # loop through bins
    for bin_contig in clusters.contigs():

        # headers of read contigs to count
        read_headers = [bin_contig.header.mutate(source=source)
                        for source in sources]

        # names of samples
        for i, header in enumerate(read_headers):
            bin_contig.loci['sample'][str(i)]['name'] = header.source

        # iterate through bins
        for locus in bin_contig.loci:

            # boolean masks of contained reads for each read set
            read_masks = [np.logical_and(reads[header].loci['tip'] >=
                                         locus['start'],
                                         reads[header].loci['tip'] <=
                                         locus['stop'])
                          for header in read_headers]

            # read loci data of contained reads for each read set
            read_loci = [reads[header].loci[mask]
                         for header, mask in zip(read_headers, read_masks)]

            # counts of reads
            for i, r in enumerate(read_loci):

                # total count of reads from each sample
                locus['sample'][str(i)]['count'] = len(r)

                # most common elements per sample per locus
                if n_common_elements > 0:
                    counts = Counter(r['element'])
                    common = counts.most_common(n_common_elements)
                    for j, pair in enumerate(common):
                        locus['sample'][str(i)]['element'][str(j)] = pair

            # find median of cluster
            combined_read_tips = reduce(np.append,
                                        (tips['tip'] for tips in read_loci))
            combined_read_tips.sort()
            locus['median'] = np.median(combined_read_tips)

            # trim potentially buffered locus to the first and last read tips
            if trim:
                locus['start'] = np.min(combined_read_tips)
                locus['stop'] = np.max(combined_read_tips)
    return clusters


def _known_insertion_matcher(contig, known_insertions, distance=0):
    """matches clusters on a contig to known elements"""
    if contig.header.strand == '+':
        for cluster in contig.loci:
            mask = interval.singular_intersects_lower(cluster['median'],
                                                      cluster['stop'] +
                                                      distance,
                                                      known_insertions.loci)
            if np.any(mask):
                yield known_insertions.loci['element'][mask][0]
            else:
                yield ''
    elif contig.header.strand == '-':
        for cluster in contig.loci:
            mask = interval.singular_intersects_upper(cluster['start'] -
                                                      distance,
                                                      cluster['median'],
                                                      known_insertions.loci)
            if np.any(mask):
                yield known_insertions.loci['element'][mask][-1]
            else:
                yield ''


def match_known_insertions(clusters, known_insertions, distance=0):
    """
    Match clusters to known insertions annotated in the genome.

    Known insertions are represented as an object of
    :class:`loci2.ContigSet` created from a gff file.
    Clusters are matched to a known insertion if they are for the
    same category and are within the specified distance of the
    insertions end.

    Fields required in 'clusters':
        'start': int, 'stop': int, 'median': int

    Fields required in 'known_insertions':
        'start': int, 'stop': int, 'element': str

    Fields appended to return value:
        'known_element': str

    :param clusters: a collection of cluster loci (intervals)
    :type clusters: :class:`loci2.ContigSet`
    :param known_insertions: a collection of cluster loci (intervals)
    :type known_insertions: :class:`loci2.ContigSet`
    :param distance: maximum distance for connecting a cluster to a
        known insertion
    :type distance: int

    :return: a collection of cluster loci (intervals) tagged with
        known insertions
    :rtype: :class:`loci2.ContigSet`
    """
    matched = loci2.ContigSet()

    # make known insertion headers un-stranded and drop origin file
    known_insertions = known_insertions.map(lambda x:
                                            loci2.mutate_header(x, strand='.',
                                                                source=None))

    # loop through contigs
    for contig in clusters.contigs():

        # get relevant known insertions
        known = known_insertions[contig.header.mutate(strand='.')]

        matches = np.array(list(_known_insertion_matcher(contig,
                                                         known,
                                                         distance=distance)))
        matches = np.array(matches,
                           dtype=np.dtype([('known_element', '<O')]))

        matched.add(loci2.Contig(contig.header, util.bind_arrays(contig.loci,
                                                                 matches)))

    return matched


def _cluster_pairer(forward, reverse, distance=0, use_known_elements=True):
    """Sorts clusters into pairs using known insertions and distances"""
    dtype_sort = np.dtype([('start', np.int64),
                           ('stop', np.int64),
                           ('median', np.int64),
                           ('known_element', '<O'),
                           ('strand', '<U1'),
                           ('index', np.int64)])

    f = np.empty(len(forward), dtype=dtype_sort)
    f['stop'] = forward.loci['stop']
    f['median'] = forward.loci['median']
    if use_known_elements:
        f['known_element'] = forward.loci['known_element']
    else:
        f['known_element'] = ''
    f['strand'] = '+'
    f['index'] = np.arange(0, len(forward))

    r = np.empty(len(reverse), dtype=dtype_sort)
    r['start'] = reverse.loci['start']
    r['median'] = reverse.loci['median']
    if use_known_elements:
        r['known_element'] = reverse.loci['known_element']
    else:
        r['known_element'] = ''
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
                if prev['known_element'] != '' and \
                   prev['known_element'] == clust['known_element']:
                    # clusters match same known element
                    yield (prev['index'], clust['index'])
                elif prev['stop'] + distance >= clust['start'] - distance:
                    # they are within reasonable range of one another
                    yield (prev['index'], clust['index'])
                else:
                    # they are not a good match
                    yield (prev['index'], None)
                    yield (None, clust['index'])
                prev = None
    if prev is not None:
        # final cluster is un-paired
        yield (prev['index'], None)


def pair_clusters(clusters, distance=0, use_known_elements=True):
    """
    Join matching clusters on opposite strands.

    Clusters of the same calgary are joined if they are within
    2 * distance of one another.
    Clusters may also be joined if they have both been matched
    to the same known element.

    Fields required in 'clusters':
        'start': int, 'stop': int, 'median': int, 'known_element': str

    Fields appended to return value:
        'ID': str, 'paired' int


    :param clusters: a collection of cluster loci (intervals)
    :type clusters: :class:`loci2.ContigSet`
    :param distance: the distance to search out from each cluster
    :type distance: int
    :param use_known_elements: specify whether to join pairs based on a
        common known element (default: True)
    :type use_known_elements: bool

    :return: a collection of cluster loci (intervals) with 'ID' and
        'paired' fields
    :rtype: :class:`loci2.ContigSet`
    """
    joint_clusters = loci2.ContigSet()

    dtype_join_data = np.dtype([("ID", "<O"), ("paired", np.int8)])

    # new headers based on old but un-stranded
    new_headers = {h.mutate(strand='.') for h in clusters.headers()}

    # template for creating insertion IDs
    template_id = '{0}_{1}_{2}_{3}'

    for header in new_headers:
        # get forward and reverse loci for this key
        forward = clusters[header.mutate(strand='+')]
        reverse = clusters[header.mutate(strand='-')]

        # sort them into pairs based on median
        pairs = _cluster_pairer(forward,
                                reverse,
                                distance=distance,
                                use_known_elements=use_known_elements)

        # create arrays for the new data
        forward_join_data = np.empty(len(forward), dtype=dtype_join_data)
        reverse_join_data = np.empty(len(reverse), dtype=dtype_join_data)
        for f, r in pairs:
            if f is not None and r is not None:
                insertion_id = template_id.format(header.reference,
                                                  forward.loci[f]['median'],
                                                  reverse.loci[r]['median'],
                                                  header.category)
                forward_join_data[f]["ID"] = insertion_id
                reverse_join_data[r]["ID"] = insertion_id
                forward_join_data[f]["paired"] = 1
                reverse_join_data[r]["paired"] = 1
            elif f is not None:
                insertion_id = template_id.format(header.reference,
                                                  forward.loci[f]['median'],
                                                  '.',
                                                  header.category)
                forward_join_data[f]["ID"] = insertion_id
            elif r is not None:
                insertion_id = template_id.format(header.reference,
                                                  '.',
                                                  reverse.loci[r]['median'],
                                                  header.category)
                reverse_join_data[r]["ID"] = insertion_id

        # combine existing data with join data and add to new contig set
        joint_clusters.add(loci2.Contig(header.mutate(strand='+'),
                                        util.bind_arrays(forward.loci,
                                                         forward_join_data)))
        joint_clusters.add(loci2.Contig(header.mutate(strand='-'),
                                        util.bind_arrays(reverse.loci,
                                                         reverse_join_data)))

    return joint_clusters


if __name__ == "__main__":
    pass
