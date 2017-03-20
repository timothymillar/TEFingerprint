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


class ReadLoci(Loci):
    DTYPE = np.dtype(_LOCI_SLOTS + _ID_SLOT + _CATEGORY_SLOT + _SOURCE_SLOT)

    @classmethod
    def from_bam(cls, bam, reference, categories, tag='ME'):
        return cls.from_iter(_bam_read_loci(bam, reference, categories, tag=tag))

    def fingerprint(self, min_reads, eps, min_eps=0, hierarchical=True):

        fingerprint = np.empty(0, dtype=FingerPrint.DTYPE)

        for case, reads in self.split(ignore=['start', 'stop', 'name']):
            tips = reads._tips()
            if hierarchical:
                model = HUDC(min_reads, max_eps=eps, min_eps=min_eps)
            else:
                model = UDC(min_reads, eps)
            model.fit(tips)
            positions = np.fromiter(model.cluster_extremities(),
                                    dtype=np.dtype([('start', np.int64),
                                                    ('stop', np.int64)]))
            array = np.empty(len(positions), dtype=FingerPrint.DTYPE)
            _partial_fill(array, case)
            array['start'] = positions['start']
            array['stop'] = positions['stop']
            fingerprint = np.append(fingerprint, array)

        return FingerPrint(fingerprint)


class FingerPrint(Loci):
    DTYPE = np.dtype(_LOCI_SLOTS + _CATEGORY_SLOT + _SOURCE_SLOT)


class ComparativeBin(Loci):
    DTYPE = np.dtype(_LOCI_SLOTS + _CATEGORY_SLOT)

    @classmethod
    def from_union(cls, loci):
        assert isinstance(loci, Loci)
        bins = cls(loci.array[list(cls.DTYPE.fields.keys())])
        bins.melt()
        return bins

    def buffer(self, value):
        pass

    def compare(self, reads):

        samples = np.unique(reads.array['sample'])
        count_dtype = list(zip(samples, [np.int64] * len(samples)))
        comparison_dtype = np.dtype(_LOCI_SLOTS + _CATEGORY_SLOT + _SOURCE_SLOT + [('counts', count_dtype)])

        for case, bins in self.split(ignore=['start', 'stop']):
            case_reads = reads.sub_group(case)
            sample_tips = [case_reads.tips[case_reads['sample'] == sample] for sample in samples]
            for start, stop, in bins.array['start', 'stop']:
                [np.sum(np.logical_and(tips >= start, tips <= stop)) for tips in sample_tips]







class Comparison(Loci):
    DTYPE = np.dtype(_LOCI_SLOTS + _CATEGORY_SLOT + _SOURCE_SLOT)