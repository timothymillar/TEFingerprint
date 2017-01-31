#! /usr/bin/env python

import numpy as np


class Reads(object):

    def __init__(self, forward=None, reverse=None):
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
    A collection of mapped SAM read positions. This class does not contain read sequences.

    Each read is represented as a single element in a numpy array with the following fields:
     - tip: np.int64
     - tail: np.int64
     - strand: np.str_, 1
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
        Init method for :class:`ReadGroup`.

        :param reads: A numpy array.
        :type reads: :class:`numpy.ndarray`[(int, int, str, str)]
        :param reference: The (optional) name of the reference/chromosome reads are aligned to
        :type reference: str
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
        Iter method for :class:`ReadGroup`.
        Passes through to wrapped numpy array.

        :return: An iterable of mapped SAM read positions and names
        :rtype: generator[(int, int, str, str)]
        """
        for read in self.reads:
            yield read

    def __getitem__(self, item):
        """
        Getitem method for :class:`ReadGroup`.
        Passes through to wrapped numpy array.

        :param item:
        :type item: int | slice | str | numpy.ndarray[int] | numpy.ndarray[bool]

        :return: An numpy array with dtype = :class:`ReadGroup`.DTYPE_READ
        :rtype: :class:`numpy.ndarray`[(int, int, str, str)]
        """
        return self.reads[item]

    def __len__(self):
        """
        Len method for :class:`ReadGroup`.
        Passes through to wrapped numpy array.

        :return: Number of reads in group
        :rtype: int
        """
        return len(self.reads)

    def sort(self, order='tip'):
        """
        Sort reads in place by field(s).

        :param order: A valid field or list of fields in :class:`ReadGroup`, defaults to 'tip'
        :type order: str | list[str]
        """
        self.reads.sort(order=order)

    def within_locus(self, start, stop, margin=0, end='tip'):
        """
        Returns a new ReadGroup object containing (the specified end of) all reads within specified (inclusive) bounds.

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
        bindex = self.within_locus(start, stop, margin=margin, end=end)
        return StrandReads(self.reads[bindex], strand=self.strand)

    @staticmethod
    def append(x, y):
        """

        :param x: ReadGroup
        :param y: ReadGroup
        :return: ReadGroup
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
        Create an instance of :class:`ReadGroup` from an iterable.

        :param iterable: an iterable of read positions
        :type iterable: iter[(int, int, str, str)]

        :return: An instance of :class:`ReadGroup`
        :rtype: :class:`StrandReads`
        """
        reads = np.fromiter(iterable, dtype=StrandReads.DTYPE_READ)
        reads.sort(order=('tip', 'tail'))
        return StrandReads(reads, **kwargs)


class UnivariateLoci(object):
    """
    A collection of univariate loci that relate to reads mapped to a reference genome.

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
        Init method for :class:`ReadLoci`.

        :param loci: A numpy array.
        :type loci: :class:`numpy.ndarray`[(int, int)]
        """
        if loci is None:
            loci = []
        self.loci = np.array(loci, dtype=UnivariateLoci.DTYPE_ULOCUS, copy=True)
        assert strand in {'+', '-'}
        self.strand = strand

    def __iter__(self):
        """
        Iter method for :class:`ReadLoci`.
        Passes through to wrapped numpy array.

        :return: an iterator of loci
        :rtype: generator[(int, int)]
        """
        for locus in self.loci:
            yield locus

    def __getitem__(self, item):
        """
        Getitem method for :class:`ReadLoci`.
        Passes through to wrapped numpy array.

        :param item: Index
        :type item: int | slice | str | numpy.ndarray[int] | numpy.ndarray[bool]

        :return: An numpy array with dtype = :class:`ReadLoci`.DTYPE_ULOCUS
        :rtype: :class:`numpy.ndarray`[(int, int)]
        """
        return self.loci[item]

    def __len__(self):
        """
        Len method for :class:`ReadLoci`.

        :return: The number of loci in the collection
        :rtype: int
        """
        return len(self.loci)

    def sort(self, order=('start', 'stop')):
        """
        Sort loci in place by field(s).

        :param order: A valid field or list of fields in :class:`ReadLoci`, defaults to ['start', 'stop']
        :type order: str | list[str]
        """
        self.loci.sort(order=order)

    def melt(self):
        """
        Merge overlapping loci into a single loci.
        Loci are sorted and modified in place.

        Example::
            loci = ReadLoci.from_iterable([(1, 4), (2, 6), (6, 7), (8, 10), (9, 12), (14, 15)])
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
            self.loci = np.fromiter(_melter(self.loci), dtype=UnivariateLoci.DTYPE_ULOCUS)

    def within_locus(self, start, stop, margin=0, end='both'):
        """
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
        Returns a new ReadGroup object containing (the specified end of) all reads within specified (inclusive) bounds.

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
        return UnivariateLoci(self.loci[bindex], strand=self.strand)

    @classmethod
    def from_iter(cls, iterable, strand=None):
        """
        Construct an instance of :class:`ReadLoci` form an iterable.

        :param iterable: Iterable of tuples containing loci bounds
        :type iterable: iterable[(int, int)]

        :return: Instance of :class:`ReadLoci`
        :rtype: :class:`UnivariateLoci`
        """
        loci = UnivariateLoci(np.fromiter(iterable, dtype=UnivariateLoci.DTYPE_ULOCUS), strand=strand)
        loci.sort()
        return loci

    @classmethod
    def append(cls, x, y):
        """
        Combine two :class:`ReadLoci` objects into a single object.

        :param x: Instance of :class:`ReadLoci`
        :type x: :class:`ReadLoci`
        :param y: Instance of :class:`ReadLoci`
        :type y: :class:`UnivariateLoci`

        :return: Instance of :class:`ReadLoci`
        :rtype: :class:`UnivariateLoci`
        """
        if len(x) == 0:
            return UnivariateLoci(y.loci, strand=y.strand)
        elif len(y) == 0:
            return UnivariateLoci(x.loci, strand=x.strand)
        else:
            assert x.strand == y.strand
            loci = UnivariateLoci(np.append(x.loci, y.loci), strand=x.strand)
            loci.sort()
            return loci

if __name__ == '__main__':
    pass
