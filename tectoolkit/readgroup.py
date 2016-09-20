#! /usr/bin/env python

import numpy as np


class ReadGroup(object):
    """"""
    read = np.dtype([('tip', np.int64),
                     ('tail', np.int64),
                     ('strand', np.str_, 1),
                     ('name', np.str_, 254)])

    def __init__(self, reads):
        self.reads = reads

    def __iter__(self):
        for read in self.reads:
            yield read

    def __getitem__(self, item):
        return self.reads[item]

    def sort(self, order='tip'):
        """

        :param order:
        :return:
        """
        self.reads.sort(order=order)

    def strand(self):
        """

        :return:
        """
        variation = set(self.reads['strand'])
        if len(variation) > 1:
            return '.'
        else:
            return variation.pop()

    def sub_group_by_locus(self, start, stop, end='tip'):
        """

        :param start:
        :param stop:
        :param end:
        :return:
        """
        assert end in {'tip', 'tail'}
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
        reads = np.fromiter(reads, dtype=ReadGroup.read)
        reads.sort(order=('tip', 'tail'))
        return ReadGroup(reads)


class ReferenceLoci(object):
    """"""
    locus = np.dtype([('start', np.int64),
                      ('stop', np.int64)])

    def __init__(self, loci):
        self.loci = loci

    def __iter__(self):
        for locus in self.loci:
            yield locus

    def __getitem__(self, item):
        return self.loci[item]

    def sort(self, order=('start', 'stop')):
        """

        :param order:
        :return:
        """
        self.loci.sort(order=order)

    def merge_overlapping(self):
        """

        :return:
        """

        def merge(loci):
            self.sort()
            start = loci['start'][0]
            stop = loci['stop'][0]
            for i in range(1, len(loci)):
                if loci['start'][i] <= stop:
                    stop = loci['stop'][i]
                else:
                    yield start, stop
                    start = loci['start'][i]
                    stop = loci['stop'][i]
            yield start, stop

        self.sort()
        self.loci = np.fromiter(merge(self.loci), dtype=ReferenceLoci.locus)

    @classmethod
    def from_simple_subcluster(cls, points, minpts, eps):
        """

        :param points:
        :param minpts:
        :param eps:
        :return:
        """
        points.sort()
        offset = minpts - 1
        upper = points[offset:]
        lower = points[:-offset]
        diff = upper - lower
        dense = diff <= eps
        lower = lower[dense]
        upper = upper[dense]
        loci = ((lower[i], upper[i]) for i in range(len(lower)))
        loci = np.fromiter(loci, dtype=ReferenceLoci.locus)
        return ReferenceLoci(loci)

    @classmethod
    def from_simple_cluster(cls, points, minpts, eps):
        """

        :param points:
        :param minpts:
        :param eps:
        :return:
        """
        reference_loci = cls.from_simple_subcluster(points, minpts, eps)
        reference_loci.merge_overlapping()
        return reference_loci
