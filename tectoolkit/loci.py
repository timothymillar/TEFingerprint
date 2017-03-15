#! /usr/bin/env python

import numpy as np
from tectoolkit.bamio import bam_read_loci as _bam_read_loci

_LOCI_SLOTS = [
    ('reference', np.str_, 254),
    ('strand', np.str_, 1),
    ('start', np.int64),
    ('stop', np.int64)]

_ID_SLOT = [
    ('name', np.str_, 254)]

_CATEGORY_SLOT = [
    ('category', np.str_, 254)]

_SOURCE_SLOT = [
    ('sample', np.str_, 254)]


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


class Loci(object):
    _DTYPE = np.dtype(_LOCI_SLOTS)

    def __init__(self, loci=None):
        if loci is None:
            loci = []
        self.array = np.array(loci, dtype=self._DTYPE, copy=True)

    def __len__(self):
        return len(self.array)

    @classmethod
    def from_iter(cls, iterable):
        return cls(np.fromiter(iterable, dtype=cls._DTYPE))

    def attributes(self):
        return list(self.array.dtype.fields.keys())

    def unique_categories(self):
        attributes = self.attributes()
        attributes.remove('start')
        attributes.remove('stop')
        return np.unique((self.array[attributes]))

    def sub_group(self, attributes):
        return (type(self)(self.array[self.array[attributes] == group])
                for group in np.unique((self.array[attributes])))

    def sort(self, order=('start', 'stop')):
        self.array.sort(order=order)

    def melt(self):
        melted = np.empty(0, dtype=self.array.dtype)

        attributes = self.attributes()
        attributes.remove('start')
        attributes.remove('stop')

        for group in self.sub_group(attributes):
            positions = np.fromiter(_loci_melter(group.array),
                                    dtype=np.dtype([('start', np.int64),
                                                    ('stop', np.int64)]))
            loci = np.empty(len(positions), dtype=melted.dtype)
            loci.fill(group.array[0])  # default values
            loci['start'] = positions['start']
            loci['stop'] = positions['stop']
            melted = np.append(melted, loci)

        self.array = melted


class ReadLoci(Loci):
    _DTYPE = np.dtype(_LOCI_SLOTS + _ID_SLOT + _CATEGORY_SLOT + _SOURCE_SLOT)

    @classmethod
    def from_bam(cls, reference, strand, tag='ME'):
        cls.from_iter(_bam_read_loci(reference, strand, tag=tag))


class FingerPrint(Loci):
    _DTYPE = np.dtype(_LOCI_SLOTS + _CATEGORY_SLOT + _SOURCE_SLOT)
