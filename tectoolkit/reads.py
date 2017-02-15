#! /usr/bin/env python

import numpy as np


class Reads(object):
    """
    A pair (forward and reverse) of collections of mapped SAM read positions.


    Each read is represented as a single element in a numpy array with the following fields:
     - tip: np.int64
     - tail: np.int64
     - name: np.str_, 256
    Where 'tip' and 'tail' are the coordinates (1 based indexing) of the read mapped to its reference,
    'strand' is a single character ('-' or '+') indicating the strand the read mapped to
    and 'name' is a string of up to 256 characters containing the read name.
    """
    def __init__(self, forward=None, reverse=None):
        """
        Requires a :class:`StrandReads` object containing forward reads and/or a :class:`StrandReads`
        object containing reverse reads.

        If both strands are specified they must be instances of :class:`StrandReads` with identical
        attributes values (reference, source, grouping) with the exception of strand which must be '+'
        for the forward strand and '-' for the reverse strand.

        If only one strand is specified the other strand will be treated as an empty instance of
        :class:`StrandReads` with attributes inherited from the specified strand.

        :param forward: An instance of :class:`StrandReads` mapped to the forward strand
        :type forward: :class:`StrandReads` | None
        :param reverse: An instance of :class:`StrandReads` mapped to the reverse strand
        :type reverse: :class:`StrandReads` | None
        """
        assert forward is not None or reverse is not None  # at least one strand must be specified
        if forward is None:
            forward = StrandReads(np.array([], dtype=StrandReads.DTYPE_READ),
                                  reference=reverse.reference,
                                  grouping=reverse.grouping,
                                  source=reverse.source,
                                  strand='+')
        if reverse is None:
            reverse = StrandReads(np.array([], dtype=StrandReads.DTYPE_READ),
                                  reference=forward.reference,
                                  grouping=forward.grouping,
                                  source=forward.source,
                                  strand='-')
        assert isinstance(forward, StrandReads)
        assert isinstance(reverse, StrandReads)
        assert forward.strand == '+'
        assert reverse.strand == '-'
        assert forward.reference == reverse.reference
        assert forward.source == reverse.source
        assert forward.grouping == reverse.grouping
        self.strand = {'+': forward, '-': reverse}
        self.source = forward.source
        self.reference = forward.reference
        self.grouping = forward.grouping


class StrandReads(object):
    """
    A collection of mapped SAM read positions on a single strand.
    This class does not contain read sequences.

    Each read is represented as a single element in a numpy array with the following fields:
     - tip: np.int64
     - tail: np.int64
     - name: np.str_, 256
    Where 'tip' and 'tail' are the coordinates (1 based indexing) of the read mapped to its reference,
    'strand' is a single character ('-' or '+') indicating the strand the read mapped to
    and 'name' is a string of up to 256 characters containing the read name.
    """
    DTYPE_READ = np.dtype([('tip', np.int64),
                           ('tail', np.int64),
                           ('name', np.str_, 254)])

    def __init__(self, reads=None, reference=None, strand=None, grouping=None, source=None):
        """
        Init method for :class:`StrandReads`.

        :param reads: A numpy array.
        :type reads: :class:`numpy.ndarray`[(int, int, str)]
        :param reference: The (optional) name of the reference/chromosome reads are aligned to
        :type reference: str
        :param strand: `+`, or `-` indicating which strand reads are mapped to
        :type strand: str
        :param grouping: The (optional) group name/type of the reads
        :type grouping: str
        :param source: The (optional) name of the source file that reads were imported from
        :type source: str
        """
        assert strand in ['+', '-']
        if reads is None:
            reads = []
        self.reads = np.array(reads, dtype=StrandReads.DTYPE_READ, copy=True)
        self.reference = reference
        self.grouping = grouping
        self.source = source
        self.strand = strand

    def __iter__(self):
        """
        Iter method for :class:`StrandReads`.
        Passes through to wrapped numpy array.

        :return: An iterable of mapped SAM read positions and names
        :rtype: generator[(int, int, str)]
        """
        for read in self.reads:
            yield read

    def __getitem__(self, item):
        """
        Getitem method for :class:`StrandReads`.
        Passes through to wrapped numpy array.

        :param item:
        :type item: int | slice | str | numpy.ndarray[int] | numpy.ndarray[bool]

        :return: An numpy array with dtype = :class:`StrandReads`.DTYPE_READ
        :rtype: :class:`numpy.ndarray`[(int, int, str)]
        """
        return self.reads[item]

    def __len__(self):
        """
        Len method for :class:`StrandReads`.
        Passes through to wrapped numpy array.

        :return: Number of reads in group
        :rtype: int
        """
        return len(self.reads)

    def sort(self, order='tip'):
        """
        Sort reads in place by field(s).

        :param order: A valid field or list of fields in :class:`StrandReads`, defaults to 'tip'
        :type order: str | list[str]
        """
        self.reads.sort(order=order)

    def within_locus(self, start, stop, margin=0, end='tip'):
        """
        Returns an array of boolean values indicating whether each read falls  within specified (inclusive) bounds.

        :param start: Lower bound
        :type start: int
        :param stop: Upper bound
        :type stop: int
        :param margin: A value to extend both bounds by, defaults to 0
        :param end: The read end that must fall within the bounds, must be 'tip', 'tail' or 'both', defaults to 'tip'
        :type end: str

        :return: An array of boolean values of equal length to the :class:`StrandReads` instance
        :rtype: :class:`numpy.ndarray`[bool]
        """
        assert end in {'tip', 'tail', 'both'}
        start -= margin
        stop += margin
        if end == 'both':
            # 'tip' may be smaller or larger than 'tail'
            tip_within = np.logical_and(self.reads['tip'] >= start, self.reads['tip'] <= stop)
            tail_within = np.logical_and(self.reads['tail'] >= start, self.reads['tail'] <= stop)
            return np.logical_and(tip_within, tail_within)
        else:
            return np.logical_and(self.reads[end] >= start, self.reads[end] <= stop)

    def subset_by_locus(self, start, stop, margin=0, end='tip'):
        """
        Returns a new :class:`StrandReads` object containing all reads within specified (inclusive) bounds.

        :param start: Lower bound
        :type start: int
        :param stop: Upper bound
        :type stop: int
        :param margin: A value to extend both bounds by, defaults to 0
        :param end: The read end that must fall within the bounds, must be 'tip' or 'tail', defaults to 'tip'
        :type end: str

        :return: The subset of reads that fall within the specified bounds
        :rtype: :class:`StrandReads`
        """
        bindex = self.within_locus(start, stop, margin=margin, end=end)
        return StrandReads(self.reads[bindex], strand=self.strand)

    @staticmethod
    def append(x, y):
        """
        Combine two instances of :class:`StrandReads` with the same attributes (reference, grouping, source and strand).

        :param x: An instance of :class:`StrandReads`
        :type x: :class:`StrandReads`
        :param y: An instance of :class:`StrandReads`
        :type y: :class:`StrandReads`

        :return: An instance of :class:`StrandReads`
        :rtype: :class:`StrandReads`
        """
        assert x.reference == y.reference
        assert x.grouping == y.grouping
        assert x.source == y.source
        assert x.strand == y.strand
        return StrandReads(np.append(x.reads, y.reads),
                           reference=x.reference,
                           grouping=x.grouping,
                           source=x.source,
                           strand=x.strand)

    @staticmethod
    def from_iter(iterable, **kwargs):
        """
        Create an instance of :class:`StrandReads` from an iterable.

        :param iterable: an iterable of read positions
        :type iterable: iter[(int, int, str, str)]

        :return: An instance of :class:`StrandReads`
        :rtype: :class:`StrandReads`
        """
        reads = np.fromiter(iterable, dtype=StrandReads.DTYPE_READ)
        reads.sort(order=('tip', 'tail'))
        return StrandReads(reads, **kwargs)


class StrandLoci(object):
    """
    A collection of univariate loci that relate to reads mapped to a single strand of a reference genome.

    Each locus is represented as a single element in a numpy array with the following slots:
     - start: np.int64
     - stop: np.int64
    Where 'start' and 'stop' are respectively the lower and upper coordinates that define a segment of the reference.
    These coordinates are inclusive, 1 based integers (i.e, use the SAM coordinate system).
    """
    DTYPE_ULOCUS = np.dtype([('start', np.int64),
                             ('stop', np.int64)])

    def __init__(self, loci=None, strand=None):
        """
        Init method for :class:`StrandLoci`.

        :param loci: A numpy array.
        :type loci: :class:`numpy.ndarray`[(int, int)]
        """
        if loci is None:
            loci = []
        self.loci = np.array(loci, dtype=StrandLoci.DTYPE_ULOCUS, copy=True)
        assert strand in {'+', '-'}
        self.strand = strand

    def __iter__(self):
        """
        Iter method for :class:`StrandLoci`.
        Passes through to wrapped numpy array.

        :return: an iterator of loci
        :rtype: generator[(int, int)]
        """
        for locus in self.loci:
            yield locus

    def __getitem__(self, item):
        """
        Getitem method for :class:`StrandLoci`.
        Passes through to wrapped numpy array.

        :param item: Index
        :type item: int | slice | str | numpy.ndarray[int] | numpy.ndarray[bool]

        :return: An numpy array with dtype = :class:`StrandLoci`.DTYPE_ULOCUS
        :rtype: :class:`numpy.ndarray`[(int, int)]
        """
        return self.loci[item]

    def __len__(self):
        """
        Len method for :class:`StrandLoci`.

        :return: The number of loci in the collection
        :rtype: int
        """
        return len(self.loci)

    def sort(self, order=('start', 'stop')):
        """
        Sort loci in place by field(s).

        :param order: A valid field or list of fields in :class:`StrandLoci`, defaults to ['start', 'stop']
        :type order: str | list[str]
        """
        self.loci.sort(order=order)

    def melt(self):
        """
        Merge overlapping loci into a single loci.
        Loci are sorted and modified in place.

        Example::
            loci = StrandLoci.from_iterable([(1, 4), (2, 6), (6, 7), (8, 10), (9, 12), (14, 15)])
            list(loci)
            [(1, 4), (2, 6), (6, 7), (8, 10), (9, 12), (14, 15)]
            loci.melt()
            list(loci)
            [(1, 7), (8, 12), (14, 15)]
        """
        def _melter(loci):
            start = loci['start'][0]
            stop = loci['stop'][0]
            for i in range(1, len(loci)):
                if loci['start'][i] <= stop:
                    if loci['stop'][i] > stop:
                        stop = loci['stop'][i]
                    else:
                        pass
                else:
                    yield start, stop
                    start = loci['start'][i]
                    stop = loci['stop'][i]
            yield start, stop
        self.sort()

        if len(self.loci) == 0:
            pass
        else:
            self.loci = np.fromiter(_melter(self.loci), dtype=StrandLoci.DTYPE_ULOCUS)

    def within_locus(self, start, stop, margin=0, end='both'):
        """
        Returns an array of booleans indicating whether each of the loci fall within specified (inclusive) bounds.

        :param start: Lower bound
        :type start: int
        :param stop: Upper bound
        :type stop: int
        :param margin: A value to extend both bounds by, defaults to 0
        :param end: Locus end that must fall within the bounds, must be 'start', 'stop' or 'both', defaults to 'both'
        :type end: str

        :returns: An array of boolean values with length equal to that of the :class:`StrandLoci` instance
        :rtype: numpy.ndarray[bool]
        """
        assert end in {'start', 'stop', 'both'}
        start -= margin
        stop += margin
        if end == 'both':
            return np.logical_and(self.loci['start'] >= start, self.loci['stop'] <= stop)
        else:
            return np.logical_and(self.loci[end] >= start, self.loci[end] <= stop)

    def subset_by_locus(self, start, stop, margin=0, end='both'):
        """
        Returns a new :class:`StrandLoci` object containing (the specified end of) all reads within
        specified (inclusive) bounds.

        :param start: Lower bound
        :type start: int
        :param stop: Upper bound
        :type stop: int
        :param margin: A value to extend both bounds by, defaults to 0
        :param end: Locus end that must fall within the bounds, must be 'start', 'stop' or 'both', defaults to 'both'
        :type end: str

        :return: The subset of reads that fall within the specified bounds
        :rtype: :class:`StrandLoci`
        """
        bindex = self.within_locus(start, stop, margin=margin, end=end)
        return StrandLoci(self.loci[bindex], strand=self.strand)

    @classmethod
    def from_iter(cls, iterable, strand=None):
        """
        Construct an instance of :class:`StrandLoci` form an iterable.

        :param iterable: Iterable of tuples containing loci bounds
        :type iterable: iterable[(int, int)]

        :return: Instance of :class:`StrandLoci`
        :rtype: :class:`StrandLoci`
        """
        loci = StrandLoci(np.fromiter(iterable, dtype=StrandLoci.DTYPE_ULOCUS), strand=strand)
        loci.sort()
        return loci

    @classmethod
    def append(cls, x, y):
        """
        Combine two :class:`StrandLoci` objects into a single object.

        :param x: Instance of :class:`StrandLoci`
        :type x: :class:`StrandLoci`
        :param y: Instance of :class:`StrandLoci`
        :type y: :class:`StrandLoci`

        :return: Instance of :class:`StrandLoci`
        :rtype: :class:`StrandLoci`
        """
        if len(x) == 0:
            return StrandLoci(y.loci, strand=y.strand)
        elif len(y) == 0:
            return StrandLoci(x.loci, strand=x.strand)
        else:
            assert x.strand == y.strand
            loci = StrandLoci(np.append(x.loci, y.loci), strand=x.strand)
            loci.sort()
            return loci

if __name__ == '__main__':
    pass
