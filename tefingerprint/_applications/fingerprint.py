#! /usr/bin/env python3

import numpy as np
from itertools import product
from functools import reduce
from collections import Counter
from multiprocessing import Pool
from tefingerprint import util
from tefingerprint.util.numpy import interval
from tefingerprint._applications import fingerprintio
from tefingerprint import loci


def fingerprint(bams,
                categories,
                references,
                minimum_reads,
                epsilon,
                minimum_epsilon=0,
                n_common_elements=0,
                splitting_method='conservative',
                fingerprint_buffer=0,
                join_distance=0,
                quality=0,
                transposon_tag='ME',
                annotation=None,
                colourise_output=True,
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
    :param splitting_method:
    :param fingerprint_buffer:
    :param join_distance:
    :param quality:
    :param transposon_tag:
    :param annotation:
    :param colourise_output:
    :param cores:
    :return:
    """

    if isinstance(bams, str):
        bams = [bams]

    if isinstance(references, str):
        references = [references]

    if isinstance(categories, str):
        categories = [categories]

    if references == [None]:
        references = fingerprintio.extract_references_from_bams(*bams)

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
                   [splitting_method],
                   [fingerprint_buffer],
                   [join_distance],
                   [colourise_output])

    result = loci.ContigSet()

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
                          splitting_method,
                          fingerprint_buffer,
                          join_distance,
                          colourise_output):
    """dispatch a single job of a fingerprint"""

    # read informative reads
    reads = fingerprintio.extract_informative_read_tips(bams,
                                                        reference,
                                                        categories,
                                                        quality=quality,
                                                        tag=transposon_tag)

    # read known transposons if used
    if annotation is None:
        pass
    else:
        known = fingerprintio.extract_gff_intervals(annotation,
                                                    reference,
                                                    categories)

    # sort reads
    reads = reads.map(lambda x: loci.sort(x, order='tip'))

    # cluster reads
    clusters = reads.map(lambda x:
                         loci.clusters(x,
                                      'tip',
                                       minimum_reads,
                                       epsilon,
                                       minimum_epsilon=minimum_epsilon,
                                       method=splitting_method))

    # buffered union of clusters
    clusters = clusters.map(lambda x:
                            loci.unions_buffered(x, fingerprint_buffer))

    # drop origin files and append
    clusters = clusters.map(lambda x:
                            loci.mutate_header(x, source=None),
                            append_duplicate_headers=True)

    # union of clusters
    clusters = clusters.map(loci.unions)

    # count reads in bins
    clusters = count_reads(clusters,
                           reads,
                           n_common_elements=n_common_elements)

    # colour features based on read count proportions
    if colourise_output:
        if len(bams) > 1:
            clusters = clusters.map(colourise_read_counts)

    # match clusters to known elements if used
    if annotation is None:
        pass
    else:
        clusters = match_known_insertions(clusters,
                                          known,
                                          distance=join_distance)

    # create cluster IDs
    clusters = clusters.map(create_contig_ids)

    # join cluster pairs
    use_known_elements = False if annotation is None else True

    clusters = pair_clusters(clusters,
                             distance=join_distance,
                             use_known_elements=use_known_elements)
    return clusters


def create_contig_ids(contig):
    """
    Create ids for contig based on category, reference, strand and position.

    :param contig: a collection of cluster loci (intervals)
    :type contig: :class:`loci.Contig`

    :return: a collection of cluster loci (intervals) with 'ID' field
    :rtype: :class:`loci.Contig`
    """
    template = '{0}_{1}_{2}_'.format(contig.header.category,
                                     contig.header.reference,
                                     contig.header.strand)
    if contig.header.strand == '+':
        position = 'stop'
    elif contig.header.strand == '-':
        position = 'start'
    else:
        assert False

    ids = [template + str(element[position]) for element in contig.loci]
    ids = np.array(ids)
    ids = np.array(ids, dtype=np.dtype([('ID', '<O')]))

    loci_data = util.numpy.array.bind(contig.loci, ids)

    return loci.Contig(contig.header, loci_data)


def count_reads(clusters, reads, trim=True, n_common_elements=0):
    """
    Counts the read tips contained in each cluster (locus).

    A total count of reads across all source (bam) files is calculated.
    The median read position amoung all reads is calculated.
    The n_common_elements parameter may be used to include counts of the 'n'
    most common elements found in each sample.
    By default, clusters are trimmed to the most distant reads they contain.

    Header fields required in 'clusters':
        reference, strand, category

    Loci fields required in 'clusters':
        'start': int, 'stop': int

    Header fields required in 'reads':
        reference, strand, category, source

    Loci fields required in 'reads':
        'tip': int, 'element': str

    Fields appended to return value:
        'median': int, 'sample': ('name': str, 'count': int,
        'element': ( 'name': str, 'count': int))

    :param clusters: a collection of cluster loci (intervals)
    :type clusters: :class:`loci.ContigSet`
    :param reads: a collection of read tips from one or more bam files
    :type reads: :class:`loci.ContigSet`
    :param trim: if true the boundaries of each cluster are trimmed to
        the range of the reads they contain
    :type trim: bool
    :param n_common_elements: specifies how many of the most common
        elements to count for each sample
    :type n_common_elements: int

    :return: a collection of cluster loci (intervals) with read counts
    :rtype: :class:`loci.ContigSet`
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
        dtype_sample_count = util.numpy.dtype.remove_field(dtype_sample_count,
                                                           'element')
    dtype_samples = np.dtype([(str(i), dtype_sample_count)
                              for i, _ in enumerate(sources)])
    dtype_new = np.dtype([('median', np.int64),
                          ('sample', dtype_samples)])

    # new bins based on previous with additional slots for counts
    clusters = clusters.map(lambda x: loci.add_field(x, dtype_new))

    # loop through bins
    for bin_contig in clusters.contigs():

        # headers of read contigs to count
        read_headers = [bin_contig.header.mutate(source=source)
                        for source in sources]

        # names of samples
        for i, header in enumerate(read_headers):
            bin_contig.loci['sample'][str(i)]['name'] = header.source

            # default name and value for count of common elements
            if n_common_elements > 0:
                for j in range(n_common_elements):
                    bin_contig.loci['sample'][str(i)]['element'][str(j)] = ('.', 0)

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


# Colours from R scheme "blues9"
_GFF_PROPORTION_COLOURS = {0: '#C6DBEF',
                           1: '#C6DBEF',
                           2: '#C6DBEF',
                           3: '#C6DBEF',
                           4: '#C6DBEF',
                           5: '#C6DBEF',
                           6: '#9ECAE1',
                           7: '#6BAED6',
                           8: '#2171B5',
                           9: '#08306B',
                           10: '#08306B'}


def colourise_read_counts(contig):
    """
    Add html colour field based on read count proportions.

    This field is rendered in genome browsers (e.g. IGV) when
    viewing the gff output.
    A darker colour is used for features where a single sample
    contains a greater proportion of the total read count.
    This funtion requires contig loci fields created by the
    "count_reads" function.

    :param contig: a contig of loci with read counts
    :type contig: :class:`Contig`
    
    :return: a contig of loci with a 'color' field
    :rtype: :class:`loci.Contig`
    """
    contig = loci.add_field(contig, np.dtype([('color', 'O')]))
    if len(contig) > 0:
        samples = [str(i)
                   for i in range(len(contig.loci[0]['sample'].dtype.descr))]
        if len(samples) > 1:
            for i, counts in enumerate(zip(*[contig.loci['sample'][s]['count']
                                             for s in samples])):
                proportion = int(10 * (max(counts) / sum(counts)))
                contig.loci['color'][i] = _GFF_PROPORTION_COLOURS[proportion]
        else:
            contig.loci['color'] = '#0112b6'  # default igv colour
    return contig


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
                yield '.'
    elif contig.header.strand == '-':
        for cluster in contig.loci:
            mask = interval.singular_intersects_upper(cluster['start'] -
                                                      distance,
                                                      cluster['median'],
                                                      known_insertions.loci)
            if np.any(mask):
                yield known_insertions.loci['element'][mask][-1]
            else:
                yield '.'


def match_known_insertions(clusters, known_insertions, distance=0):
    """
    Match clusters to known insertions annotated in the genome.

    Known insertions are represented as an object of
    :class:`loci.ContigSet` created from a gff file.
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
    :type clusters: :class:`loci.ContigSet`
    :param known_insertions: a collection of cluster loci (intervals)
    :type known_insertions: :class:`loci.ContigSet`
    :param distance: maximum distance for connecting a cluster to a
        known insertion
    :type distance: int

    :return: a collection of cluster loci (intervals) tagged with
        known insertions
    :rtype: :class:`loci.ContigSet`
    """
    matched = loci.ContigSet()

    # make known insertion headers un-stranded and drop origin file
    known_insertions = known_insertions.map(lambda x:
                                            loci.mutate_header(x, strand='.',
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

        matched.add(loci.Contig(contig.header,
                                util.numpy.array.bind(contig.loci,
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
        'start': int, 'stop': int, 'median': int,
        'known_element': str, 'ID': str

    Fields appended to return value:
        'pair' str


    :param clusters: a collection of cluster loci (intervals)
    :type clusters: :class:`loci.ContigSet`
    :param distance: the distance to search out from each cluster
    :type distance: int
    :param use_known_elements: specify whether to join pairs based on a
        common known element (default: True)
    :type use_known_elements: bool

    :return: a collection of cluster loci (intervals) with 'pair' field
    :rtype: :class:`loci.ContigSet`
    """
    joint_clusters = loci.ContigSet()

    dtype_join_data = np.dtype([("pair", "<O")])

    # new headers based on old but un-stranded
    new_headers = {h.mutate(strand='.') for h in clusters.headers()}

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
        forward_join_data["pair"] = '.'
        reverse_join_data = np.empty(len(reverse), dtype=dtype_join_data)
        reverse_join_data["pair"] = '.'
        for f, r in pairs:
            if f is not None and r is not None:
                forward_join_data[f]["pair"] = reverse.loci[r]["ID"]
                reverse_join_data[r]["pair"] = forward.loci[f]["ID"]
            else:
                pass

        # combine existing data with join data and add to new contig set
        joint_clusters.add(loci.Contig(header.mutate(strand='+'),
                                       util.numpy.array.bind(forward.loci,
                                                             forward_join_data)))
        joint_clusters.add(loci.Contig(header.mutate(strand='-'),
                                       util.numpy.array.bind(reverse.loci,
                                                             reverse_join_data)))

    return joint_clusters


if __name__ == "__main__":
    pass
