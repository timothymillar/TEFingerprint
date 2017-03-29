#! /usr/bin/env python

import numpy as np
from tectoolkit.bamio import bam_read_loci as _bam_read_loci
from tectoolkit.cluster import UDC, HUDC


def _loci_melter(array):
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
    assert len(set(map(type, args))) == 1
    merged = type(args[0])()
    for arg in args:
        merged._update_dict(arg._dict)
    return merged


class _Loci(object):
    _DTYPE_LOCI = np.dtype([('start', np.int64), ('stop', np.int64)])

    _LOCI_DEFAULT_VALUES = (0, 0)

    _DTYPE_ARRAY = np.dtype([('reference', np.str_, 256),
                             ('strand', np.str_, 1),
                             ('category', np.str_, 256),
                             ('start', np.int64),
                             ('stop', np.int64)])

    def __init__(self):
        self._dict = {}

    def groups(self):
        return self._dict.keys()

    def loci(self):
        return self._dict.values()

    def items(self):
        return self._dict.items()

    def split(self):
        for group, loci in self.items():
            child = type(self)()
            child._dict[group] = loci
            yield child

    def _update_dict(self, dictionary):
        self._dict.update(dictionary)

    def as_dict(self):
        return self._dict.copy()

    def as_array(self):
        array = np.empty(0, type(self)._DTYPE_ARRAY)
        for key, loci in self.items():
            sub_array = np.empty(len(loci), type(self)._DTYPE_ARRAY)
            sub_array.fill((*key, *type(self)._LOCI_DEFAULT_VALUES))
            for slot in list(type(self)._DTYPE_LOCI.fields.keys()):
                sub_array[slot] = loci[slot]
            array = np.append(array, sub_array)
        return array

    @staticmethod
    def _format_gff_feature(record):
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
        array = self.as_array()
        array.sort(order=('reference', 'start', 'stop'))
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
        reads = ReadLoci()
        reads._update_dict({group: np.fromiter(loci, dtype=cls._DTYPE_LOCI)
                            for group, loci in _bam_read_loci(bam, reference, categories, tag=tag)})
        return reads

    def tips(self):
        d = {}
        for group, loci in self.items():
            if group[1] == '+':
                d[group] = loci["stop"]
            else:
                d[group] = loci["start"]
        return d

    def fingerprint(self, min_reads, eps, min_eps=0, hierarchical=True):

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
        pass

    def compare(self, reads):
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
        array = np.empty(0, self._DTYPE_FLAT_ARRAY)
        for key, loci in self.items():
            for locus in loci:
                sub_array = np.empty(len(locus['sources']), self._DTYPE_FLAT_ARRAY)
                sub_array.fill((*key, locus['start'], locus['stop'], '', 0))
                sub_array['source'] = locus['sources']
                sub_array['count'] = locus['counts']
                array = np.append(array, sub_array)
        return array

    @staticmethod
    def _format_flat_gff_feature(record):
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
        array = self.as_flat_array()
        array.sort(order=('reference', 'start', 'stop', 'source'))
        return '\n'.join((self._format_flat_gff_feature(record) for record in array))