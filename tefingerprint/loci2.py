#! /usr/bin/env python

import numpy as np
from tefingerprint import utils
from tefingerprint.cluster import UDBSCANx as _UDBSCANx, UDBSCANxH as _UDBSCANxH

class Header(object):
    """
    A Loci Key is used to specify a specific 'reference', 'strand', 'category' of a genome
    Where category is a category of transposon (e.g. a super-family).
    An optional forth slot 'source' is used to record the data source (i.e. bam file) when appropriate.

    This class is used as a dictionary key within instances of :class:`GenomeLoci` to label groups
    of similar loci

    :param reference: name of a reference chromosome with the form 'name:min-max'
    :type reference: str | None
    :param strand: strand of the reference chromosome '+' or '-'
    :type strand: str | None
    :param category: name of transposon category/family
    :type category: str | None
    :param source: optional name of source bam file
    :type source: str | None
    """
    __slots__ = ['_reference', '_strand', '_category', '_source']

    def __init__(self, reference=None, strand=None, category=None, source=None):
        self._reference = reference
        self._strand = strand
        self._category = category
        self._source = source

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return hash(self) == hash(other)
        return False

    def __ne__(self, other):
        """Define a non-equality test"""
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self._reference, self._strand, self._category, self._source))

    def __iter__(self):
        for i in (self._reference, self._strand, self._category, self._source):
            if i:
                yield i

    def __repr__(self):
        return repr((self._reference, self._strand, self._category, self._source))

    @property
    def reference(self):
        return self._reference

    @property
    def strand(self):
        return self._strand

    @property
    def category(self):
        return self._category

    @property
    def source(self):
        return self._source

    @property
    def tuple(self):
        return tuple(i for i in self)

    @property
    def dtype(self):
        descr = []
        if self._reference:
            descr.append(('reference', 'O'))
        if self._strand:
            descr.append(('strand', '<U1'))
        if self._category:
            descr.append(('category', 'O'))
        if self._source:
            descr.append(('source', 'O'))
        return np.dtype(descr)

    def mutate(self, **kwargs):
        data = {'reference': self._reference,
                'strand': self._strand,
                'category': self._category,
                'source': self._source}
        data.update(kwargs)
        return Header(**data)


class ContigSet(object):

    def __init__(self, *args, append_duplicate_headers=False):
        self._dict = {}
        if len(args) == 0:
            pass
        else:
            assert len({ctg.header.dtype for ctg in args}) == 1
            assert len({ctg.loci.dtype for ctg in args}) == 1
            for arg in args:
                if arg.header in self._dict.keys():
                    if append_duplicate_headers:
                        self._dict[arg.header] = append(self._dict[arg.header], arg)
                    else:
                        raise ValueError('More than one contig with header {0}'.format(arg.header))
                else:
                    self._dict[arg.header] = arg

    def __len__(self):
        return sum(map(len, self._dict.values()))

    def __repr__(self):
        return repr(self._dict)

    def __getitem__(self, header):
        return self._dict.__getitem__(header)

    def _dtype_headers(self):
        consensus = {ctg.header.dtype for ctg in self._dict.values()}
        assert len(consensus) == 1
        return list(consensus)[0]

    def _dtype_loci(self):
        consensus = {ctg.loci.dtype for ctg in self._dict.values()}
        assert len(consensus) == 1
        return list(consensus)[0]

    def add(self, contig, append_duplicate_headers=False):
        if len(self._dict.keys()) == 0:
            pass
        else:
            assert contig.header.dtype == self._dtype_headers()
            assert contig.loci.dtype == self._dtype_loci()
        if contig.header in self._dict.keys():
            if append_duplicate_headers:
                self._dict[contig.header] = append(self._dict[contig.header], contig)
            else:
                raise ValueError('More than one contig with header {0}'.format(contig.header))
        else:
            self._dict[contig.header] = contig

    def update(self, contigs, append_duplicate_headers=False):
        for contig in contigs:
            self.add(contig, append_duplicate_headers=append_duplicate_headers)

    def map(self, function, append_duplicate_headers=False):
        return ContigSet(*map(function, self._dict.values()), append_duplicate_headers=append_duplicate_headers)

    def contigs(self):
        for contig in self._dict.values():
            yield contig

    def iter_values(self):
        for ctg in self._dict.values():
            for val in iter_values(ctg):
                yield val

    def as_array(self):
        data = map(tuple, map(utils.flatten_numpy_element, self.iter_values()))
        dtype = utils.flatten_dtype(utils.append_dtypes(self._dtype_headers(), self._dtype_loci()))
        return np.array([i for i in data], dtype=dtype)


class Contig(object):

    def __init__(self, header, array):
        assert isinstance(header, Header)
        self.header = header
        self.loci = array

    def __repr__(self):
        return repr((self.header, self.loci))

    def __len__(self):
        return len(self.loci)


def sort(contig, order=None):
    return Contig(contig.header, np.sort(contig.loci, order=order))


def iter_values(contig):
    for locus in contig.loci:
        yield (*contig.header.tuple, *locus)


def as_array(contig):
    data = map(tuple, map(utils.flatten_numpy_element, iter_values(contig)))
    dtype = utils.flatten_dtype(utils.append_dtypes(contig.header.dtype, contig.loci.dtype))
    return np.array([i for i in data], dtype=dtype)


def mutate_header(contig, **kwargs):
    return Contig(contig.header.mutate(**kwargs), contig.loci)


def append(contig_x, contig_y):
    assert contig_x.header == contig_y.header
    return Contig(contig_x.header, np.append(contig_x.loci, contig_y.loci))


def drop_field(contig, field):
    return Contig(contig.header, utils.remove_array_field(contig.loci, field))


def add_field(contig, field_dtype, value=None):
    array = np.empty(len(contig.loci), dtype=field_dtype)
    if value is None:
        pass
    else:
        array.fill(value)
    return Contig(contig.header, array=utils.bind_arrays(contig.loci, array))


def cluster(contig,
            field,
            minimum_reads,
            epsilon,
            minimum_epsilon=0,
            hierarchical=True,
            lower_bound='start',
            upper_bound='stop'):

    if hierarchical:
        model = _UDBSCANxH(minimum_reads, max_eps=epsilon, min_eps=minimum_epsilon)
    else:
        model = _UDBSCANx(minimum_reads, epsilon)
    model.fit(contig.loci[field])

    dtype = np.dtype([(lower_bound, np.int64), (upper_bound, np.int64)])
    return Contig(contig.header, np.fromiter(model.cluster_extremities(), dtype=dtype))


def _unionise(array, lower_bound, upper_bound):
    if len(array) == 1:
        yield array[lower_bound], array[upper_bound]

    else:
        array_starts = np.array(array[lower_bound])
        array_stops = np.array(array[upper_bound])
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


def unions(contig, lower_bound='start', upper_bound='stop'):

    dtype = np.dtype([(lower_bound, np.int64), (upper_bound, np.int64)])
    if len(contig.loci) == 0:
        return Contig(contig.header, np.empty(0, dtype=dtype))
    else:
        generator = _unionise(contig.loci, lower_bound, upper_bound)
        return Contig(contig.header, np.fromiter(generator, dtype=dtype))


def unions_buffered(contig, value, lower_bound='start', upper_bound='stop'):
    contig = unions(contig, lower_bound=lower_bound, upper_bound=upper_bound)

    if value <= 0:
        pass
    else:
        if len(contig.loci) == 0:
            pass
        elif len(contig.loci) == 1:
            contig.loci[lower_bound][0] -= value
            contig.loci[upper_bound][-1] += value
        else:
            difs = ((contig.loci[lower_bound][1:] - contig.loci[upper_bound][:-1]) - 1) / 2
            difs[difs > value] = value
            contig.loci[lower_bound][1:] = contig.loci[lower_bound][1:] - np.floor(difs)
            contig.loci[upper_bound][:-1] = contig.loci[upper_bound][:-1] + np.ceil(difs)
            contig.loci[lower_bound][0] -= value
            contig.loci[upper_bound][-1] += value
    return contig

