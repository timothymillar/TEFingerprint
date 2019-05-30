#! /usr/bin/env python3

import numpy as np
import tefingerprint.util as util
from tefingerprint.cluster import DBICAN as _DBICAN
from tefingerprint.cluster import SDBICAN as _SDBICAN


class Header(object):
    """
    Meta data for a collection (contig) of genomic loci.

    :param reference: name of a reference chromosome
    :type reference: str | None
    :param strand: strand of the reference chromosome '+' or '-'
    :type strand: str | None
    :param category: name of transposon category/family
    :type category: str | None
    :param source: optional name of source bam file
    :type source: str | None
    """
    __slots__ = ['_reference', '_strand', '_category', '_source']

    def __init__(self,
                 reference=None,
                 strand=None,
                 category=None,
                 source=None):
        self._reference = reference
        self._strand = strand
        self._category = category
        self._source = source

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return hash(self) == hash(other)
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self._reference,
                     self._strand,
                     self._category,
                     self._source))

    def __iter__(self):
        for i in (self._reference,
                  self._strand,
                  self._category,
                  self._source):
            if i:
                yield i

    def __repr__(self):
        return repr((self._reference,
                     self._strand,
                     self._category,
                     self._source))

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
        """
        A numpy dtype suitable for containing the header data

        :return: a numpy dtype
        :rtype: :class:`numpy.dtype`
        """
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
        """
        Create a copy of the header with altered values.

        Valid argument names are 'reference', 'strand',
        'category' and 'source'

        :param kwargs: values to alter in the new copy of the header
        :return: :class:`Header`
        """
        data = {'reference': self._reference,
                'strand': self._strand,
                'category': self._category,
                'source': self._source}
        data.update(kwargs)
        return Header(**data)


class Contig(object):
    """
    A contig is a collection of loci on a contiguous section of genome.

    The contig header defines the combination of 'reference', 'strand'
    and 'category', 'source' that distinguish this group of loci from
    others.
    The loci are a structured numpy array defining the position of each
    locus along the contig with additional metadata.
    Contigs will generally contain points or intervals.

    :param header: a header that defines this contig
    :type header: :class:`Header`
    :param loci: loci contained in this contig
    :type loci: :class:`numpy.array`
    """
    def __init__(self, header, loci):
        assert isinstance(header, Header)
        self.header = header
        self.loci = loci

    def __repr__(self):
        return repr((self.header, self.loci))

    def __len__(self):
        return len(self.loci)

    def __eq__(self, other):
        return self.header == other.header and np.array_equal(self.loci,
                                                              other.loci)

    def __ne__(self, other):
        return not self.__eq__(other)


def sort(contig, order=None):
    """
    Creates a copy of a contig with sorted loci

    :param contig: a contig of loci
    :type contig: :class:`Contig`
    :param order: field names to sort loci by
    :type order: str | tuple[str]

    :return: a contig of sorted loci
    :rtype: :class:`Contig`
    """
    return Contig(contig.header, np.sort(contig.loci, order=order))


def iter_values(contig):
    """
    Iterates values of the contig including the header values.

    For each element in the contigs array of loci, the values of header
    fields are returned with the values of that element.

    :param contig: a contig of loci
    :type contig: :class:`Contig`

    :return: an iterable of values
    :rtype generator[tuple[any]]
    """
    for locus in contig.loci:
        yield (*contig.header.tuple, *locus)


def as_array(contig):
    """
    Converts a contig to a flattened numpy array.

    Header values are included for every element in the array.
    Nested loci values are flattened by appending nested field names.

    :param contig: a contig of loci
    :type contig: :class:`Contig`

    :return: a numpy array with no nested values
    :rtype: :class:`numpy.array`
    """
    data = map(tuple, map(util.numpy.element.flatten, iter_values(contig)))
    dtype = util.numpy.dtype.flatten_fields(util.numpy.dtype.append(contig.header.dtype,
                                                                    contig.loci.dtype))
    return np.array([i for i in data], dtype=dtype)


def mutate_header(contig, **kwargs):
    """
    Return a copy of a contig with altered header data.

    Valid argument names are 'reference', 'strand', 'category'
    and 'source'.

    :param contig: a contig of loci
    :type contig: :class:`Contig`
    :param kwargs: values to alter in the new copy of the header
    :type: kwargs: list[str]

    :return: a contig of loci with modified header
    :rtype: :class:`Contig`
    """
    return Contig(contig.header.mutate(**kwargs), contig.loci)


def append(contig_x, contig_y):
    """
    Append the values of two contigs that have identical headers.

    :param contig_x: a contig of loci
    :type contig_x: :class:`Contig`
    :param contig_y: a contig of loci
    :type contig_y: :class:`Contig`

    :return: a contig containing the loci of contig_x and contig_y
    :rtype: :class:`Contig`
    """
    try:
        assert contig_x.header == contig_y.header
    except AssertionError:
        raise ValueError('Contigs must have identical headers')
    else:
        return Contig(contig_x.header, np.append(contig_x.loci, contig_y.loci))


def drop_field(contig, field):
    """
    Return a copy of a contig without data for the specified loci field.

    :param contig: a contig of loci
    :type contig: :class:`Contig`
    :param field: name of a field in the contigs loci
    :type field: str
    :return: a contig of loci without the specified field
    :rtype: :class:`Contig`
    """
    return Contig(contig.header,
                  util.numpy.array.remove_field(contig.loci, field))


def add_field(contig, field_dtype):
    """
    Add a new empty field to the the loci of a contig.

    :param contig: a contig of loci
    :type contig: :class:`Contig`
    :param field_dtype: a named numpy dtype for the added field
    :type field_dtype: :class:`numpy.dtype`

    :return: a contig of loci
    :rtype: :class:`Contig`
    """
    array = np.zeros(len(contig.loci), dtype=field_dtype)
    return Contig(contig.header,
                  loci=util.numpy.array.bind(contig.loci, array))


def clusters(contig,
             field,
             minimum_points,
             epsilon,
             minimum_epsilon=0,
             method='SDBICAN',
             lower_bound='start',
             upper_bound='stop'):
    """
    Return a contig of cluster intervals based on the loci of a contig.

    The new cluster intervals will contain only fields of the upper and
    lower bounds of clusters which may be named with the lower_bound
    and upper_bound parameters.
    The field specified for cluster analysis must contain integer values.
    By default a hierarchical density based clustering algorithm is used
    (see documentation for :class:`_UDBSCANxH`).
    Alternitively a non-hierarchical density based clustering algorithm
    may be used (see documentation for :class:`_UDBSCAN).

    :param contig: a contig of loci
    :type contig: :class:`Contig`
    :param field: the name of the loci field to perform cluster analysis on
    :type field: str
    :param minimum_points: minimum number of points required to form a cluster
    :type minimum_points: int
    :param epsilon: the maximum distance among point of a cluster
    :type epsilon: int
    :param minimum_epsilon: minimum value when calculating hierarchical splits
    :type minimum_epsilon: int
    :param method: use hierarchical or non-hierarchical version of
        the algorithm. one of 'DBICAN', 'SDBICAN', 'SDBICAN-aggressive'
    :type method: str
    :param lower_bound: field name to use for lower bound of clusters
    :type lower_bound: str
    :param upper_bound: field name to use for upper bound of clusters
    :type upper_bound: str

    :return: a contig of cluster interval loci
    :rtype: :class:`Contig`
    """
    assert method in {'DBICAN', 'SDBICAN', 'SDBICAN-aggressive'}

    if method in {'SDBICAN', 'SDBICAN-aggressive'}:
        if method == 'SDBICAN-aggressive':
            aggressive_method = True
        else:
            aggressive_method = False
        model = _SDBICAN(minimum_points,
                         epsilon,
                         min_epsilon=minimum_epsilon,
                         aggressive_method=aggressive_method)
    else:
        model = _DBICAN(minimum_points, epsilon)
    model.fit(contig.loci[field])

    dtype = np.dtype([(lower_bound, np.int64), (upper_bound, np.int64)])
    return Contig(contig.header, np.fromiter(model.cluster_extremities(),
                                             dtype=dtype))


def _unionise(array, lower_bound, upper_bound):
    """worker for unions function"""
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
    """
    Calculate interval unions among loci of a contig.

    The new cluster intervals will contain only fields of the upper and
    lower bounds of clusters which may be named with the lower_bound
    and upper_bound parameters.

    :param contig: a contig of loci
    :type contig: :class:`Contig`
    :param lower_bound: field name to use for lower bound of unions
    :type lower_bound: str
    :param upper_bound: field name to use for upper bound of unions
    :type upper_bound: str

    :return: a contig of non-overlapping union loci
    :rtype: :class:`Contig`
    """
    dtype = np.dtype([(lower_bound, np.int64), (upper_bound, np.int64)])
    if len(contig.loci) == 0:
        return Contig(contig.header, np.empty(0, dtype=dtype))
    else:
        generator = _unionise(contig.loci, lower_bound, upper_bound)
        return Contig(contig.header, np.fromiter(generator, dtype=dtype))


def unions_buffered(contig, buffer, lower_bound='start', upper_bound='stop'):
    """
    Calculate interval unions among loci of a contig and extend them
    with a non-overlapping buffer.

    The new cluster intervals will contain only fields of the upper
    and lower bounds of clusters which may be named with the
    lower_bound and upper_bound parameters.

    :param contig: a contig of loci
    :type contig: :class:`Contig`
    :param buffer: the amount to extend unions by
    :type buffer: int
    :param lower_bound: field name to use for lower bound of unions
    :type lower_bound: str
    :param upper_bound: field name to use for upper bound of unions
    :type upper_bound: str

    :return: a contig of non-overlapping union loci
    :rtype: :class:`Contig`
    """
    contig = unions(contig, lower_bound=lower_bound, upper_bound=upper_bound)

    if buffer <= 0:
        pass
    else:
        if len(contig.loci) == 0:
            pass
        elif len(contig.loci) == 1:
            contig.loci[lower_bound][0] -= buffer
            contig.loci[upper_bound][-1] += buffer
        else:
            difs = ((contig.loci[lower_bound][1:] -
                     contig.loci[upper_bound][:-1]) - 1) / 2
            difs[difs > buffer] = buffer
            contig.loci[lower_bound][1:] = (contig.loci[lower_bound][1:] -
                                            np.floor(difs))
            contig.loci[upper_bound][:-1] = (contig.loci[upper_bound][:-1] +
                                             np.ceil(difs))
            contig.loci[lower_bound][0] -= buffer
            contig.loci[upper_bound][-1] += buffer
    return contig


class ContigSet(object):
    """
    A set of contigs of genomic loci.

    A ContigSet is an unordered collection of unique contigs (each
    having a different header).
    Individual contigs may be copied from the set by indexing the
    contig set with it's header (i.e. each contigs header is a
    dictionary key to extract the contig) but contigs cannot be
    updated in place with this syntax.
    New contigs may be added to the ContigSet by the 'add' or
    'update' methods. If a contigs header clashes with that
    of a contig already contained in the set an error will be raised.
    Alternatively the loci of both contigs
    may be appended if the correct parameter is set.

    """
    def __init__(self, *args, append_duplicate_headers=False):
        self._dict = {}
        self._cached_array = None
        self._cached_array_order = None
        if len(args) == 0:
            pass
        else:
            assert len({ctg.header.dtype for ctg in args}) == 1
            assert len({ctg.loci.dtype for ctg in args}) == 1
            for arg in args:
                if arg.header in self._dict.keys():
                    if append_duplicate_headers:
                        self._dict[arg.header] = append(self._dict[arg.header],
                                                        arg)
                    else:
                        message = 'More than one contig with header {0}'
                        raise ValueError(message.format(arg.header))
                else:
                    self._dict[arg.header] = arg

    def __len__(self):
        return sum(map(len, self._dict.values()))

    def __repr__(self):
        return repr(self._dict)

    def __getitem__(self, header):
        return self._dict.__getitem__(header)

    def __eq__(self, other):
        headers_self = set(self.headers())
        headers_other = set(other.headers())

        if headers_self == headers_other:
            return all([self[h] == other[h] for h in headers_self])
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def _reset_array_cache(self):
        self._cached_array = None
        self._cached_array_order = None

    def dtype_headers(self):
        """
        Generate a dtype suitable for storing data from the headers of
        contained contigs.

        :return: a numpy dytpe
        :rtype: :class:`numpy.dtype`
        """
        consensus = {ctg.header.dtype for ctg in self._dict.values()}
        assert len(consensus) == 1
        return list(consensus)[0]

    def dtype_loci(self):
        """
        Generate a dtype matching the dtypes of contained contig loci.

        :return: a numpy dytpe
        :rtype: :class:`numpy.dtype`
        """
        consensus = {ctg.loci.dtype for ctg in self._dict.values()}
        assert len(consensus) == 1
        return list(consensus)[0]

    def headers(self):
        """
        Return the set of headers of all the contained contigs

        :return: set of contig headers
        :rtype: :class:`header`
        """
        return set(self._dict.keys())

    def add(self, contig, append_duplicate_headers=False):
        """
        Add a single contig to the contig set.

        The header and loci fields of the contig must match those already
        contained in the set.
        If the contigs header is identical to one already contained in the
        set an error is raised, this may be avoided using
        append_duplicate_headers to append the loci values of the
        matching contigs.

        :param contig: a contig of loci
        :type contig: :class:`Contig`
        :param append_duplicate_headers: specify whether to append the
            loci of contigs with duplicate headers
        :type append_duplicate_headers: bool
        """
        self._reset_array_cache()

        if len(self._dict.keys()) == 0:
            pass
        else:
            assert contig.header.dtype == self.dtype_headers()
            assert contig.loci.dtype == self.dtype_loci()
        if contig.header in self._dict.keys():
            if append_duplicate_headers:
                self._dict[contig.header] = append(self._dict[contig.header],
                                                   contig)
            else:
                message = 'More than one contig with header {0}'
                raise ValueError(message.format(contig.header))
        else:
            self._dict[contig.header] = contig

    def update(self, contigs, append_duplicate_headers=False):
        """
        Add a collection of contigs to the contig set.

        The header and loci fields of the contigs must match those
        already contained in the set. If any of the new contigs
        headers are identical to one already contained in the set an
        error is raised, this may be avoided using append_duplicate_headers
        to append the loci values of the matching contigs.

        :param contigs: a collection of contigs
        :type contigs: iterable[:class:`Contig`]
        :param append_duplicate_headers: specify whether to append the
            loci of contigs with duplicate headers
        :type append_duplicate_headers: bool
        """
        self._reset_array_cache()

        for contig in contigs:
            self.add(contig, append_duplicate_headers=append_duplicate_headers)

    def map(self, function, append_duplicate_headers=False):
        """
        Map a function to every contig within the set.

        The mapped function should return :class:`Contig` to be collected
        into a new :class:`ContigSet`.
        If the mapped function alters the contig headers and leads to a
        header clash an error is raised.
        This may be avoided by the 'append_duplicate_headers' parameter
        which appends the loci of contigs with identical headers.

        :param function: function to apply to every contig
        :type function: callable
        :param append_duplicate_headers: specify whether to append the
            loci of contigs with duplicate headers
        :type append_duplicate_headers: bool

        :return: a new ContigSet
        :rtype: :class:`ContigSet`
        """
        return ContigSet(*map(function, self._dict.values()),
                         append_duplicate_headers=append_duplicate_headers)

    def contigs(self):
        """
        Yield all contigs in the set.

        :return: a generator of all contigs
        :rtype: iterable[:class:`Contig`]
        """
        for contig in self._dict.values():
            yield contig

    def iter_values(self):
        """
        Iterates values of the contigs including the header values.

        For each element in each contigs array of loci, the values of
        header fields are returned with the values of that element.

        :return: an iterable of values
        :rtype generator[tuple[any]]
        """
        for ctg in self._dict.values():
            for val in iter_values(ctg):
                yield val

    def _rebuild_array_cache(self, order=False):
        if order == self._cached_array_order:
            pass
        else:
            self._reset_array_cache()

            data = map(tuple,
                       map(util.numpy.element.flatten, self.iter_values()))
            dtype = util.numpy.dtype.flatten_fields(
                util.numpy.dtype.append(self.dtype_headers(),
                                        self.dtype_loci()))
            self._cached_array = np.array([i for i in data], dtype=dtype)

            if order:
                if isinstance(order, str):
                    order = [order]
                index = np.argsort(self._cached_array[order],
                                   order=(tuple(order)))
                self._cached_array = self._cached_array[index]

    def as_array(self, order=False):
        """
        Converts a ContigSet to a flattened numpy array.

        Header values are included for every element in each contig.
        Nested loci values are flattened by appending nested field names.

        :param order: fields to sort rows by
        :type order: tuple[str]

        :return: a numpy array with no nested values
        :rtype: :class:`numpy.array`
        """
        self._rebuild_array_cache(order=order)
        return np.copy(self._cached_array)

    def as_tabular_lines(self, sep='\t', quote=False, order=False):
        """
        Converts a ContigSet to an iterable of strings.

        Header values are included for every element in each contig.
        Nested loci values are flattened by appending nested field names.

        :param sep: separator value (default: '\t')
        :type sep: str
        :param quote: logical to surround string values in quotations
        :type quote: bool

        :return: an iterable of string suitable for writing to a tabular
            plain text file
        :rtype: iterable[str]
        """
        self._rebuild_array_cache(order=order)

        if quote:
            string = util.misc.quote_str
        else:
            string = str

        dtype = util.numpy.dtype.append(self.dtype_headers(), self.dtype_loci())
        columns = util.numpy.dtype.flatten_field_names(dtype)

        yield sep.join(map(string, columns))

        for element in self._cached_array:
            yield sep.join(map(string, element))

    def as_gff_lines(self,
                     order=False,
                     reference_field='reference',
                     start_field='start',
                     stop_field='stop',
                     strand_field='strand',
                     type_field=None,
                     program_name='TEFingerprint'):
        """
        Converts a ContigSet to an iterable of strings.

        Header values are included for every element in each contig.
        Nested loci values are flattened by appending nested field names

        :param order: fields to sort the gff by
        :type order: tuple[str]
        :param reference_field: name of field containing reference genome
            contig name
        :type reference_field: str
        :param start_field: name of field containing start of feature
        :type start_field: str
        :param stop_field: name of field containing end of feature
        :type stop_field: str
        :param strand_field: name of field containing strand of feature
        :type strand_field: str
        :param type_field: name of field containing type of feature
        :type type_field: str
        :param program_name: name of program used to generate gff file
        :type program_name: str

        :return: an iterable of string suitable for writing to a gff
            plain text file
        :rtype: iterable[str]
        """
        self._rebuild_array_cache(order=order)

        attribute_fields = list(self._cached_array.dtype.names)
        for field in (reference_field, start_field, stop_field, strand_field):
            if field in attribute_fields:
                attribute_fields.remove(field)
        if type_field:
            attribute_fields.remove(type_field)

        template = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}"

        for record in self._cached_array:
            if type_field:
                type_value = record[type_field]
            else:
                type_value = '.'

            attributes = ';'.join(('{0}={1}'.format(
                util.gff3.encode_attribute(field),
                util.gff3.encode_attribute(str(record[field])))
                                   for field in attribute_fields))
            yield template.format(util.gff3.encode_column(record[reference_field]),
                                  util.gff3.encode_column(program_name),
                                  util.gff3.encode_column(type_value),
                                  record[start_field],
                                  record[stop_field],
                                  '.',
                                  record[strand_field],
                                  '.',
                                  attributes)


if __name__ == "__main__":
    pass
