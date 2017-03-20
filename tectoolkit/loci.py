#! /usr/bin/env python

import numpy as np
from functools import reduce
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
        merged.dict.update(arg.dict)
    return merged


class ReadLoci(object):
    _DTYPE = np.dtype([('start', np.int64), ('stop', np.int64), ('name', np.str_, 254)])

    _DTYPE_FLAT = np.dtype([('reference', np.str_, 256),
                            ('strand', np.str_, 1),
                            ('category', np.str_, 256),
                            ('source', np.str_, 256),
                            ('start', np.int64),
                            ('stop', np.int64),
                            ('name', np.str_, 254)])

    def __init__(self):
        self.dict = {}

    @classmethod
    def from_bam(cls, bam, reference, categories, tag='ME'):
        reads = ReadLoci()
        reads.dict = {key: np.fromiter(loci, dtype=cls._DTYPE)
                      for key, loci in _bam_read_loci(bam, reference, categories, tag=tag)}
        return reads

    def as_array(self):
        array = np.empty(0, ReadLoci._DTYPE_FLAT)
        for key, loci in self.dict.items():
            sub_array = np.empty(len(loci), ReadLoci._DTYPE_FLAT)
            sub_array.fill((*key, 0, 0, ''))
            sub_array['start'] = loci['start']
            sub_array['stop'] = loci['stop']
            sub_array['name'] = loci['name']
            array = np.append(array, sub_array)
        return array

    def fingerprint(self, min_reads, eps, min_eps=0, hierarchical=True):

        fprint = FingerPrint()

        for key, loci in self.dict.items():

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
            fprint.dict[key] = positions

        return fprint


class FingerPrint(object):
    _DTYPE = np.dtype([('start', np.int64),
                       ('stop', np.int64)])

    _DTYPE_FLAT = np.dtype([('reference', np.str_, 256),
                            ('strand', np.str_, 1),
                            ('category', np.str_, 256),
                            ('source', np.str_, 256),
                            ('start', np.int64),
                            ('stop', np.int64)])

    def __init__(self):
        self.dict = {}

    def as_array(self):
        array = np.empty(0, FingerPrint._DTYPE_FLAT)
        for key, loci in self.dict.items():
            sub_array = np.empty(len(loci), FingerPrint._DTYPE_FLAT)
            sub_array.fill((*key, 0, 0))
            sub_array['start'] = loci['start']
            sub_array['stop'] = loci['stop']
            array = np.append(array, sub_array)
        return array


class ComparativeBins(object):
    _DTYPE = np.dtype([('start', np.int64), ('stop', np.int64)])

    _DTYPE_FLAT = np.dtype([('reference', np.str_, 256),
                            ('strand', np.str_, 1),
                            ('category', np.str_, 256),
                            ('start', np.int64),
                            ('stop', np.int64)])

    def __init__(self):
        self.dict = {}

    @classmethod
    def from_union(cls, *args):
        keys = list(set([key[0:3] for arg in args for key in arg.dict.keys()]))
        array_dict = {key: np.empty(0, dtype=ComparativeBins._DTYPE) for key in keys}
        for arg in args:
            for key, loci in arg.dict.items():
                array_dict[key[0:3]] = np.append(array_dict[key[0:3]], loci)

        for key, loci in array_dict.items():
            array_dict[key] = np.fromiter(_loci_melter(loci), dtype=ComparativeBins._DTYPE)

        bins = ComparativeBins()
        bins.dict = array_dict
        return bins

    def buffer(self, value):
        pass

    def compare(self, reads):
        pass

    def as_array(self):
        array = np.empty(0, ComparativeBins._DTYPE_FLAT)
        for key, loci in self.dict.items():
            sub_array = np.empty(len(loci), ComparativeBins._DTYPE_FLAT)
            sub_array.fill((*key, 0, 0))
            sub_array['start'] = loci['start']
            sub_array['stop'] = loci['stop']
            array = np.append(array, sub_array)
        return array
