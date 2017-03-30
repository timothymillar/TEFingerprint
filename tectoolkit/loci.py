#! /usr/bin/env python

import numpy as np
from tectoolkit.bamio import bam_read_loci as _bam_read_loci
from tectoolkit.cluster import UDC, HUDC


def _loci_melter(array):
    """
    Combine overlapping loci into a single locus.

    Example::
        loci = numpy.array([(1, 4), (2, 6), (6, 7), (8, 10), (9, 12), (14, 15)],
                             dtype = np.dtype([('start', np.int64), ('stop', np.int64)]))
        _loci_melter(loci)
        list(loci)
        [(1, 7), (8, 12), (14, 15)]

    :param array: a structured numpy array with integer slots 'start' and 'stop'
    :type array: :class:`numpy.ndarray`[(int, int, ...)]

    :return: a generator of melted loci
    :rtype: generator[(int, int)]
    """
    array_starts = np.array(array['start'])
    array_stops = np.array(array['stop'])
    array_starts.sort()
    array_stops.sort()

    start = array_starts[0]
    stop = array_stops[0]

    for i in range(1, len(array_starts)):
        if array_starts[i] <= stop:
            if array_stops[i] > stop:
                stop = array_stops[i]
            else:
                pass
        else:
            yield start, stop
            start = array_starts[i]
            stop = array_stops[i]
    yield start, stop


def merge(*args):
    """
    Combine multiple objects with superclass :class:`_Loci`.
    Every object must be of the same type.
    Loci with the same group are overwritten.

    :param args: Objects to be merged
    :type args: iterable[:class:`_Loci`]

    :return: Merged object
    :rtype: :class:`_Loci`
    """
    assert len(set(map(type, args))) == 1
    merged = type(args[0])()
    for arg in args:
        merged._update_dict(arg._dict)
    return merged


class _Loci(object):
    """
    A collection of categorised univariate loci.

    Loci are categorised into groups based on a combination of states relating to their location and or origin:
    - reference: str
    - strand: str
    - category: str
    Where 'reference' takes the form 'name:min-max', 'strand' is '+' or '-', and 'category' identifies which
    category of transpson the loci represent.

    Each locus is represented as a single element in a numpy array with the following slots:
     - start: np.int64
     - stop: np.int64
    Where 'start' and 'stop' are respectively the lower and upper coordinates that define a segment of the reference.
    Coordinates are inclusive, 1 based integers (i.e, use the SAM coordinate system).
    """
    _DTYPE_LOCI = np.dtype([('start', np.int64), ('stop', np.int64)])

    _LOCI_DEFAULT_VALUES = (0, 0)

    _DTYPE_ARRAY = np.dtype([('reference', np.str_, 256),
                             ('strand', np.str_, 1),
                             ('category', np.str_, 256),
                             ('start', np.int64),
                             ('stop', np.int64)])

    def __init__(self):
        self._dict = {}

    def __len__(self):
        return sum([len(loci) for loci in self.loci()])

    def __repr__(self):
        return repr(self._dict)

    def groups(self):
        """
        Return the definition (key) for each group of loci.

        :return: a collection of tuples containing the character states of each group of loci
        :rtype: :class:`dict_keys`[(str, str, ...)]
        """
        return self._dict.keys()

    def loci(self):
        """
        Return the loci array (value) for each group of loci.

        :return: a collection of structured :class:`numpy.ndarray` containing the character states of each locus
        :rtype: :class:`dict_values`[:class:`numpy.ndarray`[(int, int, ...)]]
        """
        return self._dict.values()

    def items(self):
        """
        Return paired groups and loci where groups are tuples and loci are structured :class:`numpy.ndarray`.

        :return: a collection of key-value pairs
        :rtype: :class:`dict_values`[(str, str, ...): :class:`numpy.ndarray`[(int, int, ...)]]
        """
        return self._dict.items()

    def split(self):
        """
        Split object of :class:`_Loci` by its constituent groups.

        :return: An object of :class:`_Loci` for each grouping of loci
        :rtype: generator[:class:`_Loci`]
        """
        for group, loci in self.items():
            child = type(self)()
            child._dict[group] = loci
            yield child

    def _update_dict(self, dictionary):
        """
        Update the wrapped dictionary.
        """
        self._dict.update(dictionary)

    def as_array(self):
        """
        Convert all loci to a structured array sorted by location.

        :return: a structured array of loci
        :rtype: :class:`numpy.ndarray`
        """
        array = np.empty(0, type(self)._DTYPE_ARRAY)
        for key, loci in self.items():
            sub_array = np.empty(len(loci), type(self)._DTYPE_ARRAY)
            sub_array.fill((*key, *type(self)._LOCI_DEFAULT_VALUES))
            for slot in list(type(self)._DTYPE_LOCI.fields.keys()):
                sub_array[slot] = loci[slot]
            array = np.append(array, sub_array)
        array.sort(order=('reference', 'start', 'stop'))
        return array

    @staticmethod
    def _format_gff_feature(record):
        """
        Worker method to format a single array entry as a single gff formatted string.

        :param record: a single entry from a structured :class:`numpy.ndarray`
        :type record: :class:`numpy.void`

        :return: gff formatted string
        :rtype: str
        """
        identifier = '_'.join([str(record[field]) for field in ('reference', 'strand', 'category', 'start')])
        attributes = '{0}={1}'.format('category', record['category'])
        attributes = 'ID=' + identifier + ';' + attributes
        template = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}"
        return template.format(record['reference'].split(':')[0],
                               '.',
                               '.',
                               record['start'],
                               record['stop'],
                               '.',
                               record['strand'],
                               '.',
                               attributes)

    def as_gff(self):
        """
        Convert all loci to a GFF3 formatted string, sorted by location.

        :return: A multi-line GFF3 formatted string
        :rtype: str
        """
        array = self.as_array()
        return '\n'.join((self._format_gff_feature(record) for record in array))


class ReadLoci(_Loci):
    _DTYPE_LOCI = np.dtype([('start', np.int64), ('stop', np.int64), ('name', np.str_, 254)])

    _LOCI_DEFAULT_VALUES = (0, 0, '')

    _DTYPE_ARRAY = np.dtype([('reference', np.str_, 256),
                             ('strand', np.str_, 1),
                             ('category', np.str_, 256),
                             ('source', np.str_, 256),
                             ('start', np.int64),
                             ('stop', np.int64),
                             ('name', np.str_, 254)])

    @classmethod
    def from_bam(cls, bam, reference, categories, tag='ME'):
        """
        Create an object of :class:`ReadLoci` from a bam file.
        The bam file should contain aligned reads without pairs where each read is tagged with the category of
        transposon it corresponds to.
        The tag 'ME' is used by default.

        :param bam: The path to a bam file
        :type bam: str
        :param reference: The name or slice of a reference chromosome to obtain reads from
        :type reference: str
        :param categories: A list of transposon categories to group reads by
        :type categories: list[str]
        :param tag: The sam tag containing the transposon each read corresponds to
        :type tag: str

        :return: An instance of :class:`ReadLoci`
        :rtype: :class:`ReadLoci`
        """
        reads = ReadLoci()
        reads._update_dict({group: np.fromiter(loci, dtype=cls._DTYPE_LOCI)
                            for group, loci in _bam_read_loci(bam, reference, categories, tag=tag)})
        return reads

    def tips(self):
        """
        Returns a dictionary group-tips pairs where groups are the same those used in :class:`ReadLoci`
        and tips are an array of integers representing the front tip of each read based on its orientation.

        :return: A dictionary group-tips pairs
        :rtype: dict[(int, int, str): :class:`numpy.ndarray`[int]]
        """
        d = {}
        for group, loci in self.items():
            if group[1] == '+':
                d[group] = loci["stop"]
            else:
                d[group] = loci["start"]
        return d

    def fingerprint(self, min_reads, eps, min_eps=0, hierarchical=True):
        """
        Fingerprint a :class:`ReadLoci` object based on density of read tips.

        The algorithm used is determined by the value of 'hierarchical' which is True by default).
        If hierarchical is True, the :class:`HUDC` algorithm is used.
        If hierarchical is False, the :class:`UDC` algorithm is used.

        :param min_reads: Minimum number of reads required to form cluster in cluster analysis
        :type min_reads: int
        :param eps: The eps value to be used in the cluster analysis
        :type eps: int
        :param min_eps: The minimum eps value to be used in if the :class:`HUDC` algorithm is used
        :type min_eps: int
        :param hierarchical: Determines which clustering algorithm is used
        :type hierarchical: bool

        :return: The density based fingerprints of read loci
        :rtype: :class:`FingerPrint`
        """
        dictionary = {}

        for group, tips in self.tips().items():
            tips.sort()

            # fit model
            if hierarchical:
                model = HUDC(min_reads, max_eps=eps, min_eps=min_eps)
            else:
                model = UDC(min_reads, eps)
            model.fit(tips)

            # get new loci
            positions = np.fromiter(model.cluster_extremities(),
                                    dtype=FingerPrint._DTYPE_LOCI)

            # add to fingerprint
            dictionary[group] = positions

        fprint = FingerPrint()
        fprint._update_dict(dictionary)
        return fprint

    @staticmethod
    def _format_gff_feature(record):
        """
        Worker method to format a single array entry as a single gff formatted string.

        :param record: a single entry from a structured :class:`numpy.ndarray`
        :type record: :class:`numpy.void`

        :return: gff formatted string
        :rtype: str
        """
        attributes = ';'.join(['{0}={1}'.format(slot, record[slot]) for slot in ('category', 'source')])
        attributes = 'ID=' + record['name'] + ';' + attributes
        template = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}"
        return template.format(record['reference'].split(':')[0],
                               '.',
                               '.',
                               record['start'],
                               record['stop'],
                               '.',
                               record['strand'],
                               '.',
                               attributes)


class FingerPrint(_Loci):

    _DTYPE_ARRAY = np.dtype([('reference', np.str_, 256),
                             ('strand', np.str_, 1),
                             ('category', np.str_, 256),
                             ('source', np.str_, 256),
                             ('start', np.int64),
                             ('stop', np.int64)])

    @staticmethod
    def _format_gff_feature(record):
        """
        Worker method to format a single array entry as a single gff formatted string.

        :param record: a single entry from a structured :class:`numpy.ndarray`
        :type record: :class:`numpy.void`

        :return: gff formatted string
        :rtype: str
        """
        identifier = '_'.join([str(record[field]) for field in ('reference', 'strand', 'category', 'start')])
        attributes = ';'.join(['{0}={1}'.format(slot, record[slot]) for slot in ('category', 'source')])
        attributes = 'ID=' + identifier + ';' + attributes
        template = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}"
        return template.format(record['reference'].split(':')[0],
                               '.',
                               '.',
                               record['start'],
                               record['stop'],
                               '.',
                               record['strand'],
                               '.',
                               attributes)


class ComparativeBins(_Loci):

    @classmethod
    def from_union(cls, *args):
        """
        Constructs a set of comparative bins form one or more inheriting from :class:`_Loci`.

        :param args: Objects to be merged
        :type args: iterable[:class:`_Loci`]

        :return: Bins to be used for comparisons
        :rtype: :class:`ComparativeBins`
        """
        for arg in args:
            assert isinstance(arg, _Loci)
        groups = list(set([group[0:3] for arg in args for group in arg.groups()]))
        dictionary = {group: np.empty(0, dtype=ComparativeBins._DTYPE_LOCI) for group in groups}
        for arg in args:
            for group, loci in arg.items():
                dictionary[group[0:3]] = np.append(dictionary[group[0:3]], loci)

        for group, loci in dictionary.items():
            dictionary[group] = np.fromiter(_loci_melter(loci), dtype=ComparativeBins._DTYPE_LOCI)

        bins = ComparativeBins()
        bins._update_dict(dictionary)
        return bins

    def buffer(self, value):
        """
        Expand each locus by a set amount in both directions (mutates data in place).
        Loci can not be expanded beyond the bounds of their reference and will not overlap neighbouring loci within
        their group.

        :param value: the (maximum) amount to expand loci in both directions
        :type value: int
        """
        if value <= 0:
            pass
        else:
            for key, loci in self.items():
                reference = key[0]
                minimum, maximum = tuple(map(int, reference.split(':')[1].split('-')))
                difs = ((loci['start'][1:] - loci['stop'][:-1]) - 1) / 2
                difs[difs > value] = value
                loci['start'][1:] = loci['start'][1:] - np.floor(difs)
                loci['stop'][:-1] = loci['stop'][:-1] + np.ceil(difs)
                loci['start'][0] = max(loci['start'][0] - value, minimum)
                loci['stop'][-1] = min(loci['stop'][-1] + value, maximum)

    def compare(self, reads):
        """
        Compare the read counts within each bin using reads from multiple sources (bam files).

        :param reads: Reads to be compared
        :type reads: :class:`ReadLoci`

        :return: A comparison of read counts by source within each bin
        :rtype: :class:`Comparison`
        """
        assert isinstance(reads, ReadLoci)
        sources = np.array(list({group[3] for group in list(reads.groups())}))
        sources.sort()
        tips_dict = reads.tips()
        results = {}
        for group, bins in self.items():
            group_results = np.empty(len(bins), dtype=Comparison._DTYPE_LOCI)
            group_results['start'] = bins['start']
            group_results['stop'] = bins['stop']
            group_results['sources'] = [sources for _ in bins]

            sample_tips = [tips_dict[(*group, sample)] for sample in sources]
            group_results['counts'] = [np.array([np.sum(np.logical_and(tips >= start, tips <= stop)) for tips in sample_tips]) for start, stop in bins]

            results[group] = group_results

        comparison = Comparison()
        comparison._update_dict(results)
        return comparison


class Comparison(_Loci):
    _DTYPE_LOCI = np.dtype([('start', np.int64),
                            ('stop', np.int64),
                            ('sources', np.object),
                            ('counts', np.object)])

    _LOCI_DEFAULT_VALUES = (0, 0, None, None)

    _DTYPE_ARRAY = np.dtype([('reference', np.str_, 256),
                             ('strand', np.str_, 1),
                             ('category', np.str_, 256),
                             ('start', np.int64),
                             ('stop', np.int64),
                             ('sources', np.object),
                             ('counts', np.object)])

    _DTYPE_FLAT_ARRAY = np.dtype([('reference', np.str_, 256),
                                  ('strand', np.str_, 1),
                                  ('category', np.str_, 256),
                                  ('start', np.int64),
                                  ('stop', np.int64),
                                  ('source', np.str_, 256),
                                  ('count', np.int64)])

    _LOCI_FLAT_DEFAULT_VALUES = (0, 0, '', 0)

    @staticmethod
    def _format_gff_feature(record):
        """
        Worker method to format a single array entry as a single gff formatted string.

        :param record: a single entry from a structured :class:`numpy.ndarray`
        :type record: :class:`numpy.void`

        :return: gff formatted string
        :rtype: str
        """
        identifier = '_'.join([str(record[field]) for field in ('reference', 'strand', 'category', 'start')])
        attributes = '{0}={1}'.format('category', record['category'])
        attributes += ';' + ';'.join(['{0}={1}'.format(slot, ','.join(map(str, record[slot])))
                                      for slot in ('sources', 'counts')])
        attributes = 'ID=' + identifier + ';' + attributes
        template = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}"
        return template.format(record['reference'].split(':')[0],
                               '.',
                               '.',
                               record['start'],
                               record['stop'],
                               '.',
                               record['strand'],
                               '.',
                               attributes)

    def as_flat_array(self):
        """
        Convert all loci to a structured array sorted by location.
        This method creates one entry per sample per locus to avoid nested structures.

        :return: a structured array of loci
        :rtype: :class:`numpy.ndarray`
        """
        array = np.empty(0, self._DTYPE_FLAT_ARRAY)
        for key, loci in self.items():
            for locus in loci:
                sub_array = np.empty(len(locus['sources']), self._DTYPE_FLAT_ARRAY)
                sub_array.fill((*key, locus['start'], locus['stop'], '', 0))
                sub_array['source'] = locus['sources']
                sub_array['count'] = locus['counts']
                array = np.append(array, sub_array)
        array.sort(order=('reference', 'start', 'stop'))
        return array

    @staticmethod
    def _format_flat_gff_feature(record):
        """
        Worker method to format a single array entry as a single gff formatted string.

        :param record: a single entry from a structured :class:`numpy.ndarray`
        :type record: :class:`numpy.void`

        :return: gff formatted string
        :rtype: str
        """
        identifier = '_'.join([str(record[field]) for field in ('reference', 'strand', 'category', 'start')])
        attributes = ';'.join(['{0}={1}'.format(slot, record[slot]) for slot in ('category', 'source', 'count')])
        attributes = 'ID=' + identifier + ';' + attributes
        template = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}"
        return template.format(record['reference'].split(':')[0],
                               '.',
                               '.',
                               record['start'],
                               record['stop'],
                               '.',
                               record['strand'],
                               '.',
                               attributes)

    def as_flat_gff(self):
        """
        Convert all loci to a GFF3 formatted string, sorted by location.
        This method creates one feature per sample per locus to avoid attributes containing lists.

        :return: A multi-line GFF3 formatted string
        :rtype: str
        """
        array = self.as_flat_array()
        return '\n'.join((self._format_flat_gff_feature(record) for record in array))