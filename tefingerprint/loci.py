#! /usr/bin/env python

import numpy as np
from collections import Counter
from functools import reduce
from tefingerprint import bamio
from tefingerprint import cluster
from tefingerprint.utils import append_dtypes, flatten_numpy_element, flatten_dtype, flatten_dtype_fields, quote_str, drop_dtype_field


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
    if len(array) == 1:
        yield array['start'], array['stop']

    else:
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


def append(*args):
    """
    Combine multiple objects of superclass :class:`_Loci`.
    Every object must be of the same type.

    :param args: Objects to be merged
    :type args: iterable[:class:`_Loci`]

    :return: Merged object
    :rtype: :class:`GenomeLoci`
    """
    assert len(set(map(type, args))) == 1
    merged = type(args[0])(dtype_key=args[0].dtype_key,
                           dtype_loci=args[0].dtype_loci)
    for arg in args:
        merged.append(arg, inplace=True)
    return merged


class LociKey(object):
    """
    A Loci Key is used to specify a specific 'reference', 'strand', 'category' of a genome
    Where category is a category of transposon (e.g. a super-family).
    An optional forth slot 'source' is used to record the data source (i.e. bam file) when appropriate.

    This class is used as a dictionary key within instances of :class:`GenomeLoci` to label groups
    of similar loci

    :param reference: name of a reference chromosome with the form 'name:min-max'
    :type reference: str
    :param strand: strand of the reference chromosome '+' or '-'
    :type strand: str
    :param category: name of transposon category/family
    :type category: str
    :param source: optional name of source bam file
    :type source: str
    """
    __slots__ = ['reference', 'strand', 'category', 'source']

    def __init__(self, reference=None, strand=None, category=None, source=None):
        self.reference = reference
        self.strand = strand
        self.category = category
        self.source = source

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__slots__ == other.__slots__
        return False

    def __ne__(self, other):
        """Define a non-equality test"""
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.reference, self.strand, self.category, self.source))

    def __iter__(self):
        yield self.reference
        yield self.strand
        yield self.category
        if self.source:
            yield self.source

    def __repr__(self):
        if self.source:
            return repr((self.reference, self.strand, self.category, self.source))
        else:
            return repr((self.reference, self.strand, self.category))


class GenomeLoci(object):
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
    _LOCI_DEFAULT_VALUES = (0, 0)

    def __init__(self, dtype_key, dtype_loci):
        self.dtype_key = dtype_key
        self.dtype_loci = dtype_loci
        self.dtype_array = append_dtypes(dtype_key, dtype_loci)
        self._dict = {}

    def __len__(self):
        return sum([len(loci) for loci in self.loci()])

    def __repr__(self):
        return repr(self._dict)

    def __getitem__(self, item):
        if isinstance(item, LociKey):
            return self._dict[item]
        else:
            return self._dict[LociKey(*item)]

    def keys(self):
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
        Split object of :class:`GenomeLoci` by its constituent groups.

        :return: An object of :class:`GenomeLoci` for each grouping of loci
        :rtype: generator[:class:`GenomeLoci`]
        """
        for group, loci in self.items():
            child = type(self)()
            child._dict[group] = loci
            yield child

    def append(self, x, inplace=False):
        """
        Append a Loci object.

        :param x: the Loci object to append
        :param inplace: if true, loci are appended in place rather than returning a new object
        :type inplace: bool
        """
        assert type(self) == type(x)
        keys = self.keys()
        if inplace:
            for key, loci in x.items():
                if key in keys:
                    self._dict[key] = np.append(self._dict[key], loci)
                else:
                    self._dict[key] = loci
        else:
            result = type(self)()
            result._dict = self._dict.copy()
            for key, loci in x.items():
                if key in keys:
                    result._dict[key] = np.append(result._dict[key], loci)
                else:
                    result._dict[key] = loci
            return result

    def features(self):
        """
        Yields one tuple per locus.

        :return: tuple of values per locus
        :rtype: generator
        """
        for key, loci in self.items():
            for locus in loci:
                yield (*key, *locus)

    @classmethod
    def from_dict(cls, dictionary):  # TODO: setting key dtype
        """
        Creat from a dictionary with the correct keys and dtype

        :param dictionary: dictionary of key-loci pairs
        :type dictionary: dict[tuple(...), :class:`np.array`[...]

        :return: instance of :class:`GenomeLoci`
        """
        d = {}
        for key, loci in dictionary.items():
            key = LociKey(*key)  # use named tuple for keys
            d[key] = loci
        result = cls()
        result._dict = d
        return result

    @classmethod
    def from_array(cls, array, key_fields):  # TODO: setting dtypes
        """
        Create from an array with the correct dtype

        :param array: numpy array
        :type array: :class:`np.array`[...]

        :return: instance of :class:`GenomeLoci`
        """
        assert array.dtype == cls._DTYPE_ARRAY
        key_instances = array[list(cls._DTYPE_KEY.names)]
        keys = np.unique(key_instances)
        d = {}
        for key in keys:
            loci = array[list(cls._DTYPE_LOCI.names)][key_instances == key]
            key = LociKey(*key)
            d[key] = loci
        result = cls()
        result._dict = d
        return result

    def as_array(self, order=False):
        """
        Convert all loci to a structured array sorted by location.

        :return: a structured array of loci
        :rtype: :class:`numpy.ndarray`
        """
        array = np.fromiter(self.features(), dtype=self.dtype_array, count=len(self))

        if order:
            if isinstance(order, str):
                order = [order]
            index = np.argsort(array[order],
                               order=(tuple(order)))
            array = array[index]
        return array

    def as_flat_array(self, order=False):
        data = map(tuple, map(flatten_numpy_element, self.features()))
        array = np.fromiter(data, flatten_dtype(dtype=self.dtype_array), count=len(self))

        if order:
            if isinstance(order, str):
                order = [order]
            index = np.argsort(array[order],
                               order=(tuple(order)))
            array = array[index]

        return array

    def as_tabular_lines(self, sep=','):
        yield sep.join(map(quote_str, flatten_dtype_fields(self.dtype_array.descr))) + '\n'
        for f in self.features():
            yield sep.join(map(quote_str, flatten_numpy_element(f))) + '\n'

    def as_gff_lines(self, order=False,
                     reference='reference',
                     start='start',
                     stop='stop',
                     strand='strand',
                     category='category'):

        array = self.as_flat_array(order=order)

        attribute_fields = list(array.dtype.names)
        for field in (reference, start, stop, strand):
            attribute_fields.remove(field)

        template = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n"

        for record in array:
            if category in attribute_fields:
                identifier = "{0}:{1}-{2}_{3}_{4}".format(record[reference].split(':')[0],
                                                          record[start],
                                                          record[stop],
                                                          record[strand],
                                                          record[category])
            else:
                identifier = "{0}:{1}-{2}_{3}".format(record[reference].split(':')[0],
                                                      record[start],
                                                      record[stop],
                                                      record[strand])

            attributes = ('{0}={1}'.format(field, record[field]) for field in attribute_fields)
            attributes = 'ID=' + identifier + ';' + ';'.join(attributes)
            yield template.format(record['reference'].split(':')[0],
                                  '.',
                                  '.',
                                  record['start'],
                                  record['stop'],
                                  '.',
                                  record['strand'],
                                  '.',
                                  attributes)



class InformativeReadLoci(GenomeLoci):
    """
    A collection of named sam reads, categorised by origin.

    Loci are categorised into groups based on a combination of states relating to their location and origin:
     - reference: str
     - strand: str
     - category: str
     - source: str
    Where 'reference' takes the form 'name:min-max', 'strand' is '+' or '-', 'category' identifies which
    category of transpson the loci represent, and source is the name of the source bam file.

    Each locus is represented as a single element in a numpy array with the following slots:
     - start: np.int64
     - stop: np.int64
     - name: np.str_, 254
    Where 'start' and 'stop' are respectively the lower and upper coordinates that define a segment of the reference,
    and name is the reads name in the source bam file.
    Coordinates are inclusive, 1 based integers (i.e, use the SAM coordinate system).
    """

    @classmethod
    def from_bam(cls, bams, categories, references=None, quality=30, tag='ME'):
        """
        Create an object of :class:`ReadLoci` from a bam file.
        The bam file should contain aligned reads without pairs where each read is tagged with the category of
        transposon it corresponds to.
        The tag 'ME' is used by default.

        :param bams: The path to a bam file or a list of paths to multiple bam files
        :type bams: str | list[str]
        :param references: The name or slice of a reference chromosome to obtain reads from
        :type references: str
        :param categories: A list of transposon categories to group reads by
        :type categories: list[str]
        :param quality: minimum mapping quality
        :type quality: int
        :param tag: The sam tag containing the transposon each read corresponds to
        :type tag: str

        :return: An instance of :class:`ReadLoci`
        :rtype: :class:`InformativeReadLoci`
        """
        reads = InformativeReadLoci(dtype_key=np.dtype([('reference', np.str_, 256),
                                                        ('strand', np.str_, 1),
                                                        ('category', np.str_, 256),
                                                        ('source', np.str_, 256)]),
                                    dtype_loci=np.dtype([('start', np.int64),
                                                         ('stop', np.int64),
                                                         ('name', np.str_, 254)]))
        reads._dict = {LociKey(*group): np.fromiter(loci, dtype=reads.dtype_loci)
                       for group, loci in bamio.extract_bam_reads(bams,
                                                                  categories,
                                                                  references=references,
                                                                  quality=quality,
                                                                  tag=tag)}
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
            if group.strand == '+':
                d[group] = loci["stop"]
            else:
                d[group] = loci["start"]
        return d

    def named_tips(self):
        """
        """
        d = {}
        for group, loci in self.items():
            if group.strand == '+':
                d[group] = loci[['name', "stop"]]
            else:
                d[group] = loci[['name', "start"]]
            d[group].dtype.names = ('name', 'tip')
        return d

    def cluster(self, min_reads, eps, min_eps=0, hierarchical=True):
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
        :rtype: :class:`GenomicBins`
        """
        dictionary = {}

        fprint = GenomicBins(dtype_key=self.dtype_key,
                             dtype_loci=np.dtype([('start', np.int64),
                                                  ('stop', np.int64)]))

        for key, tips in self.tips().items():
            tips.sort()

            # fit model
            if hierarchical:
                model = cluster.UDBSCANxH(min_reads, max_eps=eps, min_eps=min_eps)
            else:
                model = cluster.UDBSCANx(min_reads, eps)
            model.fit(tips)

            # get new loci
            positions = np.fromiter(model.cluster_extremities(),
                                    dtype=fprint.dtype_loci)

            # add to fingerprint
            dictionary[key] = positions

        fprint._dict = dictionary
        return fprint


class GenomicBins(GenomeLoci):
    """
    A collection of categorised univariate loci comprising a density based fingerprint.

    Loci are categorised into groups based on a combination of states relating to their location and origin:
     - reference: str
     - strand: str
     - category: str
     - source: str
    Where 'reference' takes the form 'name:min-max', 'strand' is '+' or '-', and 'category' identifies which
    category of transpson the loci represent.

    Each locus is represented as a single element in a numpy array with the following slots:
     - start: np.int64
     - stop: np.int64
    Where 'start' and 'stop' are respectively the lower and upper coordinates that define a segment of the reference.
    Coordinates are inclusive, 1 based integers (i.e, use the SAM coordinate system).
    """

    def merge_sources(self):
        """"""
        new_bins = GenomicBins(dtype_key=drop_dtype_field(self.dtype_key, 'source'),
                               dtype_loci=self.dtype_loci)

        for key, loci in self.items():

            # Use the Fingerprint key to make a ComparativeBins key
            new_key = LociKey(reference=key.reference,
                              strand=key.strand,
                              category=key.category)
            # Populate the dictionary with loci from FingerPrint object
            if new_key not in new_bins._dict.keys():
                new_bins._dict[new_key] = loci
            else:
                new_bins._dict[new_key] = np.append(new_bins._dict[new_key], loci)

        return new_bins

    def melt(self, *args):

        new_bins = GenomicBins(dtype_key=self.dtype_key,
                               dtype_loci=np.dtype([('start', np.int64),
                                                    ('stop', np.int64)]))

        for key, loci in self._dict.items():
            if len(loci) == 0:
                new_bins._dict[key] = np.empty(0, dtype=new_bins.dtype_loci)
            else:
                new_bins._dict[key] = np.fromiter(_loci_melter(loci), dtype=new_bins.dtype_loci)

        return new_bins

    def buffer(self, value):
        """
        Expand each locus by a set amount in both directions (mutates object in place.
        Loci can not be expanded beyond the bounds of their reference and will not overlap neighbouring loci within
        their group.

        :param value: the (maximum) amount to expand loci in both directions
        :type value: int
        """
        if value <= 0:
            pass
        else:
            for key, loci in self.items():
                if len(loci) == 0:
                    pass
                elif len(loci) == 1:
                    reference = key.reference
                    minimum, maximum = tuple(map(int, reference.split(':')[1].split('-')))
                    loci['start'][0] = max(loci['start'][0] - value, minimum)
                    loci['stop'][-1] = min(loci['stop'][-1] + value, maximum)
                else:
                    reference = key.reference
                    minimum, maximum = tuple(map(int, reference.split(':')[1].split('-')))
                    difs = ((loci['start'][1:] - loci['stop'][:-1]) - 1) / 2
                    difs[difs > value] = value
                    loci['start'][1:] = loci['start'][1:] - np.floor(difs)
                    loci['stop'][:-1] = loci['stop'][:-1] + np.ceil(difs)
                    loci['start'][0] = max(loci['start'][0] - value, minimum)
                    loci['stop'][-1] = min(loci['stop'][-1] + value, maximum)

    def count_reads(self, reads, trim=True, n_common_elements=0):
        assert isinstance(reads, InformativeReadLoci)
        sources = np.array(list({key.source for key in list(reads.keys())}))
        sources.sort()

        dtype_element_count = np.dtype([('name', np.str_, 256),
                                        ('count', np.int64)])
        dtype_elements = np.dtype([(str(i), dtype_element_count) for i in range(n_common_elements)])

        dtype_sample_count = np.dtype([('name', np.str_, 256),
                                       ('count', np.int64),
                                       ('element', dtype_elements)])

        if n_common_elements == 0:
            dtype_sample_count = drop_dtype_field(dtype_sample_count, 'element')

        dtype_samples = np.dtype([(str(i), dtype_sample_count) for i, _ in enumerate(sources)])

        dtype_loci = np.dtype([('start', np.int64),
                               ('stop', np.int64),
                               ('median', np.int64),
                               ('sample', dtype_samples)])

        # BinCounts to collect result
        comparison = BinCounts(dtype_key=self.dtype_key,
                               dtype_loci=dtype_loci)

        # dictionary of tips from all reads
        all_tips = reads.named_tips()

        for group, bins in self.items():
            new_loci = np.empty(len(bins), dtype=dtype_loci)

            # if the loci aren't trimmed then they are the same as the bins
            new_loci['start'] = bins['start']
            new_loci['stop'] = bins['stop']

            # sources are always the same
            for i, name in enumerate(sources):
                new_loci['sample'][str(i)]['name'] = name

            # list of dictionary keys comparable to this group of bins
            comparable_keys = [LociKey(*group, sample) for sample in sources]

            # list of read tips comparable to this group of bins
            comparable_tips = [all_tips[key] for key in comparable_keys]

            # iterate through the new loci
            for locus in new_loci:

                # named reads tips within this locus
                local_tips = [tips[np.logical_and(tips['tip'] >= locus['start'], tips['tip'] <= locus['stop'])] for tips in
                              comparable_tips]

                # counts of reads
                for i, r in enumerate(local_tips):

                    # total count of reads from each sample
                    locus['sample'][str(i)]['count'] = len(r)

                    # most common elements per sample per locus
                    if n_common_elements > 0:
                        for j, pair in enumerate(Counter(r['name']).most_common(n_common_elements)):
                            locus['sample'][str(i)]['element'][str(j)] = pair

                # find median of cluster
                combined_tips = reduce(np.append, (tips['tip'] for tips in local_tips))
                locus['median'] = np.median(combined_tips)

                # trim the potentially buffered locus to the first and last read tips
                if trim:
                    locus['start'] = np.min(combined_tips)
                    locus['stop'] = np.max(combined_tips)

            # add loci to new dict
            comparison._dict[group] = new_loci

        return comparison


class BinCounts(GenomeLoci):
    """"""


if __name__ == '__main__':
    pass
