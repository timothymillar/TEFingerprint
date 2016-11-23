#! /usr/bin/env python

import numpy as np
from tectoolkit.cluster import _UnivariateLoci
from tectoolkit.gff import GffFeature


class ReadGroup(object):
    """
    A collection of mapped SAM read positions. This class does not contain read sequences.

    Each read is represented as a single element in a numpy array with the following slots:
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
                           ('strand', np.str_, 1),
                           ('name', np.str_, 254)])

    def __init__(self, reads):
        """
        Init method for :class:`ReadGroup`.

        :param reads: A numpy array.
        :type reads: :class:`numpy.ndarray`[(int, int, str, str)]
        """
        self.reads = np.array(reads, dtype=ReadGroup.DTYPE_READ, copy=True)

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

        :param order:
        :return:
        """
        self.reads.sort(order=order)

    def strand(self):
        """
        Identifies the consensus strand of all reads in the group.

        :return: '+', '-' or '.' respectively if all reads are on forwards, reverse or combination of strands
        :rtype: str
        """
        variation = set(self.reads['strand'])
        if len(variation) > 1:
            return '.'
        else:
            return variation.pop()

    def subset_by_locus(self, start, stop, margin=0, end='tip'):
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
        :rtype: :class:`ReadGroup`
        """
        assert end in {'tip', 'tail'}
        start -= margin
        stop += margin
        reads = self.reads[np.logical_and(self.reads[end] >= start, self.reads[end] <= stop)]
        return ReadGroup(reads)

    @classmethod
    def _parse_sam_flag(cls, flag):
        """

        :param flag:
        :return:
        """
        attributes = np.zeros(12, dtype=np.bool)
        bits = np.fromiter(map(int, tuple(bin(int(flag)))[:1:-1]), dtype=np.bool)
        attributes[:bits.shape[0]] = bits
        return attributes

    @classmethod
    def _flag_orientation(cls, flag):
        """

        :param flag:
        :return:
        """
        attributes = cls._parse_sam_flag(flag)
        if attributes[2]:  # read is unmapped
            return None
        elif attributes[4]:  # read is reversed
            return '-'
        else:  # read is forwards
            return '+'

    @classmethod
    def _flag_attributes(cls, flag):
        attributes = ("read paired",
                      "read mapped in proper pair",
                      "read unmapped mate unmapped",
                      "read reverse strand",
                      "mate reverse strand",
                      "first in pair",
                      "second in pair",
                      "not primary alignment",
                      "read fails platform / vendor quality checks",
                      "read is PCR or optical duplicate",
                      "supplementary alignment")
        values = cls._parse_sam_flag(flag)
        return dict(zip(attributes, values))

    @classmethod
    def _parse_sam_strings(cls, strings, single_strand=None):
        """

        :param strings:
        :param single_strand:
        :return:
        """

        def _parse_sam_string(string, strand):
            """

            :param string:
            :param strand:
            :return:
            """
            attr = string.split("\t")
            name = str(attr[0])
            start = int(attr[3])
            length = len(attr[9])
            end = start + length
            if strand is None:
                strand = cls._flag_orientation(int(attr[1]))
            if strand == '+':
                tip = end
                tail = start
                return tip, tail, strand, name
            elif strand == '-':
                tip = start
                tail = end
                return tip, tail, strand, name
            elif strand is None:
                pass

        assert single_strand in ['+', '-', None]
        reads = (_parse_sam_string(string, single_strand) for string in strings)
        return reads

    @classmethod
    def from_sam_strings(cls, strings, strand=None):
        """

        :param strings:
        :param strand:
        :return:
        """
        reads = cls._parse_sam_strings(strings, single_strand=strand)
        reads = np.fromiter(reads, dtype=ReadGroup.DTYPE_READ)
        reads.sort(order=('tip', 'tail'))
        return ReadGroup(reads)


class ReadLoci(_UnivariateLoci):
    """"""
    def __init__(self, loci):
        """"""
        self.loci = np.array(loci, dtype=ReadLoci._ulocus, copy=True)

    def __iter__(self):
        """"""
        for locus in self.loci:
            yield locus

    def __getitem__(self, item):
        """"""
        return self.loci[item]

    def __len__(self):
        """"""
        return len(self.loci)

    def sort(self, order=('start', 'stop')):
        """

        :param order:
        :return:
        """
        self.loci.sort(order=order)

    def melt(self):
        """

        :return:
        """

        self.sort()
        self.loci = self._melt_uloci(self.loci)

    def subset_by_locus(self, start, stop, margin=0, end='start'):
        """

        :param start:
        :param stop:
        :param margin:
        :param end:
        :return:
        """
        assert end in {'start', 'stop'}
        start -= margin
        stop += margin
        loci = self.loci[np.logical_and(self.loci[end] >= start, self.loci[end] <= stop)]
        return ReadLoci(loci)

    def as_gff(self, *args, **kwargs):
        return [GffFeature(*args, start=locus[0], end=locus[1], **kwargs) for locus in self.loci]

    @classmethod
    def from_iterable(cls, iterable):
        """

        :param iterable:
        :return:
        """
        loci = ReadLoci(np.fromiter(iterable, dtype=ReadLoci._ulocus))
        loci.sort()
        return loci

    @classmethod
    def append(cls, x, y):
        """

        :param x:
        :param y:
        :return:
        """
        loci = ReadLoci(np.append(x.loci, y.loci))
        loci.sort()
        return loci

if __name__ == '__main__':
    pass
