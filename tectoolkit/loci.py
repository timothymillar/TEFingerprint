#! /usr/bin/env python

import numpy as np
from functools import reduce
from tectoolkit.bamio import bam_read_loci as _bam_read_loci
from tectoolkit.cluster import UDC, HUDC

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


def _partial_fill(array, template):
    for slot in template.dtype.fields.keys():
        array[slot] = template[slot]


class Loci(object):
    DTYPE = np.dtype(_LOCI_SLOTS)

    def __init__(self, loci=None):
        if loci is None:
            loci = []
        self.array = np.array(loci, dtype=self.DTYPE, copy=True)

    def __len__(self):
        return len(self.array)

    @classmethod
    def from_iter(cls, iterable):
        return cls(np.fromiter(iterable, dtype=cls.DTYPE))

    def attributes(self):
        return list(self.array.dtype.fields.keys())

    def _tips(self):
        forward_tips = self.array['stop'][self.array['strand'] == '+']
        reverse_tips = self.array['start'][self.array['strand'] == '-']
        tips = np.append(forward_tips, reverse_tips)
        tips.sort()
        return tips

    def sub_structures(self, by=None, ignore=None):
        if not by:
            by = self.attributes()
        if ignore:
            for attribute in ignore:
                by.remove(attribute)
        return np.unique((self.array[by]))

    def sub_group(self, structure):
        attributes = list(structure.dtype.fields.keys())
        return type(self)(self.array[self.array[attributes] == structure])

    def split(self, by=None, ignore=None):
        for structure in self.sub_structures(by=by, ignore=ignore):
            yield structure, self.sub_group(structure)

    def sort(self, order=('start', 'stop')):
        self.array.sort(order=order)

    def melt(self, by=None, ignore=None):
        ignore = ignore + ['start', 'stop'] if ignore else ['start', 'stop']

        melted = np.empty(0, dtype=self.array.dtype)

        for case, loci in self.split(by=by, ignore=ignore):
            positions = np.fromiter(_loci_melter(loci.array),
                                    dtype=np.dtype([('start', np.int64),
                                                    ('stop', np.int64)]))
            array = np.empty(len(positions), dtype=melted.dtype)
            _partial_fill(array, case)
            array['start'] = positions['start']
            array['stop'] = positions['stop']
            melted = np.append(melted, array)

        self.array = melted


class ReadLoci(object):
    _DTYPE = np.dtype([('start', np.int64), ('stop', np.int64), ('name', np.str_, 254)])

    def __init__(self):
        self._hash = {}

    @classmethod
    def from_bam(cls, bam, reference, categories, tag='ME'):
        reads = ReadLoci()
        reads._hash = {key: np.fromiter(loci, dtype=cls._DTYPE)
                       for key, loci in _bam_read_loci(bam, reference, categories, tag=tag)}
        return reads

    def fingerprint(self, min_reads, eps, min_eps=0, hierarchical=True):

        fprint = FingerPrint()

        for key, loci in self._hash.items():

            # get tips
            if key[1] == '+':
                tips = loci["stop"]
            else:
                tips = loci["start"]
            tips.sort()

            # fit model
            if hierarchical:
                model = HUDC(min_reads, max_eps=eps, min_eps=min_eps)
            else:
                model = UDC(min_reads, eps)
            model.fit(tips)

            # get new loci
            positions = np.fromiter(model.cluster_extremities(),
                                    dtype=FingerPrint._DTYPE)

            # add to fingerprint
            fprint._hash[key] = positions

        return fprint


class FingerPrint(object):
    _DTYPE = np.dtype([('start', np.int64), ('stop', np.int64)])

    def __init__(self):
        self._hash = {}


class ComparativeBins(object):
    _DTYPE = np.dtype([('start', np.int64), ('stop', np.int64)])

    def __init__(self):
        self._hash = {}

    @classmethod
    def from_union(cls, *args):
        keys = list(set([key[0:3] for arg in args for key in arg._hash.keys()]))
        array_hash = {key: np.empty(0, dtype=ComparativeBins._DTYPE) for key in keys}
        for arg in args:
            for key, loci in arg._hash.items():
                array_hash[key[0:3]] = np.append(array_hash[key[0:3]], loci)

        for key, loci in array_hash.items():
            array_hash[key] = np.fromiter(_loci_melter(loci), dtype=ComparativeBins._DTYPE)

        bins = ComparativeBins()
        bins._hash = array_hash
        return bins

    def buffer(self, value):
        pass

    def compare(self, reads):
        pass
