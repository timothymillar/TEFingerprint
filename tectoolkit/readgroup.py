#! /usr/bin/env python

import pysam
import numpy as np
from itertools import product

locus = np.dtype([('start', np.int64),
                  ('stop', np.int64)])

_read = np.dtype([('tip', np.int64),
                  ('tail', np.int64),
                  ('reverse', np.bool),
                  ('name', np.str_, 254)])


def _flag_attributes(flag):
    attributes = {"read paired": False,
                  "read mapped in proper pair": False,
                  "read unmapped mate unmapped": False,
                  "read reverse strand": False,
                  "mate reverse strand": False,
                  "first in pair": False,
                  "second in pair": False,
                  "not primary alignment": False,
                  "read fails platform / vendor quality checks": False,
                  "read is PCR or optical duplicate": False,
                  "supplementary alignment": False}
    return list(map(lambda x: bool(int(x)), tuple(bin(int(flag)))[2:]))



class ReadGroup(object):
    """"""
    def __init__(self, reads):
        self.reads = reads

    def _parse_sam_string(self, sam_string):
        """

        :param sam_strings:
        :param strand:
        :return:
        """
        attr = sam_string.split("\t")
        start = int(attr[3])
        length = len(attr[9])
        end = start + length
        if strand == '+':
            reverse = True
            tip = end
            tail = start
            return tip, tail, length, reverse
        elif strand == '-':
            reverse = False
            tip = start
            tail = end
            return tip, tail, length, reverse



    @staticmethod
    def from_sam_strings(sam_strings):
        reads = (parse_sam_string(string, strand) for string in sam_strings)
        reads = np.fromiter(reads, dtype=sam_read)
        reads.sort(order=('tip', 'tail'))
        return reads