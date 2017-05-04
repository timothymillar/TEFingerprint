#! /usr/bin/env python

import pytest
import numpy as np
import numpy.testing as npt
from tectoolkit import loci


@pytest.mark.parametrize("query,answer",
                         # empty values unhandled in this function
                         # single locus spanning single base
                         [([(13, 13)],
                           [(13, 13)]),
                          # nested loci
                          ([(15, 25), (16, 17), (19, 20)],
                           [(15, 25)]),
                          # adjacent loci
                          ([(7, 9), (10, 12)],
                           [(7, 9), (10, 12)]),
                          # combined
                          ([(3, 6), (6, 8), (7, 9), (10, 12), (13, 13), (15, 25), (16, 17), (19, 20)],
                           [(3, 9), (10, 12), (13, 13), (15, 25)])],
                         ids=['single', 'nested', 'adjacent', 'combined'])
def test_loci_melter(query, answer):
    """
    Test for function _loci_melter.
    Test includes following edge cases:
     * Long locus completely overlaps short loci: (15, 25) & (16, 17) & (19, 20) --> (15, 25)
     * Adjacent loci do not get merged: (7, 9) & (10, 12) -->  (*, 9) & (10, *)
     * Locus may span a single base: (13, 13) --> (13, 13)
    """
    query = np.array(query, dtype=loci.GenomeLoci._DTYPE_LOCI)
    query = np.fromiter(loci._loci_melter(query), dtype=loci.GenomeLoci._DTYPE_LOCI)
    answer = np.array(answer, dtype=loci.GenomeLoci._DTYPE_LOCI)
    npt.assert_array_equal(query, answer)


class TestGenomeLoci:
    """Tests for class GenomeLoci"""
    def test_keys(self):
        """Test that keys are created from dict correctly and returned correctly"""
        query = {('chr1:0-100', '+', 'Copia'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                        dtype=loci.GenomeLoci._DTYPE_LOCI),
                 ('chr1:0-100', '-', 'Copia'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                        dtype=loci.GenomeLoci._DTYPE_LOCI)}
        query = loci.GenomeLoci.from_dict(query)
        answer = {loci.LociKey('chr1:0-100', '+', 'Copia'), loci.LociKey('chr1:0-100', '-', 'Copia')}
        assert set(query.keys()) == answer

    def test_split(self):
        """Split into separate GenomeLoci objects based on keys"""
        dictionary = {('chr1:0-100', '+', 'Copia'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                             dtype=loci.GenomeLoci._DTYPE_LOCI),
                      ('chr1:0-100', '-', 'Copia'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                             dtype=loci.GenomeLoci._DTYPE_LOCI),
                      ('chr1:0-100', '+', 'Gypsy'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                             dtype=loci.GenomeLoci._DTYPE_LOCI),
                      ('chr1:0-100', '-', 'Gypsy'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                             dtype=loci.GenomeLoci._DTYPE_LOCI),
                      ('chr2:0-100', '+', 'Copia'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                             dtype=loci.GenomeLoci._DTYPE_LOCI),
                      ('chr2:0-100', '-', 'Copia'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                             dtype=loci.GenomeLoci._DTYPE_LOCI),
                      ('chr2:0-100', '+', 'Gypsy'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                             dtype=loci.GenomeLoci._DTYPE_LOCI),
                      ('chr2:0-100', '-', 'Gypsy'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                             dtype=loci.GenomeLoci._DTYPE_LOCI)}
        query = loci.GenomeLoci.from_dict(dictionary)
        assert len(list(query.split())) == 8

    def test_append(self):
        """Tests that loci arrays are added or appended in the result of a key clash"""
        query1 = {('chr1:0-100', '+', 'Copia'): np.array([(1, 2), (2, 3), (3, 4)],  # clashing key
                                                         dtype=loci.GenomeLoci._DTYPE_LOCI),
                  ('chr1:0-100', '-', 'Gypsy'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                         dtype=loci.GenomeLoci._DTYPE_LOCI)}
        query1 = loci.GenomeLoci.from_dict(query1)

        query2 = {('chr1:0-100', '+', 'Copia'): np.array([(4, 5), (5, 6)],  # clashing key
                                                         dtype=loci.GenomeLoci._DTYPE_LOCI),
                  ('chr1:0-100', '-', 'vLINE'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                         dtype=loci.GenomeLoci._DTYPE_LOCI)}
        query2 = loci.GenomeLoci.from_dict(query2)

        query = query1.append(query2)

        answer = {('chr1:0-100', '+', 'Copia'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                         dtype=loci.GenomeLoci._DTYPE_LOCI),
                  ('chr1:0-100', '-', 'Gypsy'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                         dtype=loci.GenomeLoci._DTYPE_LOCI),
                  ('chr1:0-100', '-', 'vLINE'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                         dtype=loci.GenomeLoci._DTYPE_LOCI)}
        answer = loci.GenomeLoci.from_dict(answer)

        assert set(query.keys()) == set(answer.keys())
        for key in query.keys():
            npt.assert_array_equal(query[key], answer[key])

    def test_append_inplace(self):
        """Tests that loci arrays are added or appended inplace"""
        query = {('chr1:0-100', '+', 'Copia'): np.array([(1, 2), (2, 3), (3, 4)],  # clashing key
                                                        dtype=loci.GenomeLoci._DTYPE_LOCI),
                 ('chr1:0-100', '-', 'Gypsy'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                        dtype=loci.GenomeLoci._DTYPE_LOCI)}
        query = loci.GenomeLoci.from_dict(query)

        query2 = {('chr1:0-100', '+', 'Copia'): np.array([(4, 5), (5, 6)],  # clashing key
                                                         dtype=loci.GenomeLoci._DTYPE_LOCI),
                  ('chr1:0-100', '-', 'vLINE'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                         dtype=loci.GenomeLoci._DTYPE_LOCI)}
        query2 = loci.GenomeLoci.from_dict(query2)

        query.append(query2, inplace=True)

        answer = {('chr1:0-100', '+', 'Copia'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                         dtype=loci.GenomeLoci._DTYPE_LOCI),
                  ('chr1:0-100', '-', 'Gypsy'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                         dtype=loci.GenomeLoci._DTYPE_LOCI),
                  ('chr1:0-100', '-', 'vLINE'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                         dtype=loci.GenomeLoci._DTYPE_LOCI)}
        answer = loci.GenomeLoci.from_dict(answer)

        assert set(query.keys()) == set(answer.keys())
        for key in query.keys():
            npt.assert_array_equal(query[key], answer[key])

    def test_buffer(self):
        """Test that loci are buffered correctly"""
        query = loci.GenomeLoci.from_dict({('chr1:0-300',
                                            '+',
                                            'Gypsy'): np.array([(10, 20),
                                                                (51, 60),
                                                                (160, 180),
                                                                (200, 300)],
                                                               dtype=loci.GenomeLoci._DTYPE_LOCI),
                                           ('chr2:0-300',
                                            '+',
                                            'Gypsy'): np.array([(10, 20)],
                                                               dtype=loci.GenomeLoci._DTYPE_LOCI),
                                           ('chr3:0-300',
                                            '+',
                                            'Gypsy'): np.array([],
                                                               dtype=loci.GenomeLoci._DTYPE_LOCI)})
        query.buffer(20)
        answer = loci.GenomeLoci.from_dict({('chr1:0-300',
                                             '+',
                                             'Gypsy'): np.array([(0, 35),
                                                                 (36, 80),
                                                                 (140, 190),
                                                                 (191, 300)],
                                                                dtype=loci.GenomeLoci._DTYPE_LOCI),
                                            ('chr2:0-300',
                                             '+',
                                             'Gypsy'): np.array([(0, 40)],
                                                                dtype=loci.GenomeLoci._DTYPE_LOCI),
                                            ('chr3:0-300',
                                             '+',
                                             'Gypsy'): np.array([],
                                                                dtype=loci.GenomeLoci._DTYPE_LOCI)})
        assert set(query.keys()) == set(answer.keys())
        for key in query.keys():
            npt.assert_array_equal(query[key], answer[key])

    def test_as_array(self):
        """"""
        dictionary = {('chr1:0-100', '+', 'Copia'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.GenomeLoci._DTYPE_LOCI),
                      ('chr1:0-100', '-', 'Copia'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.GenomeLoci._DTYPE_LOCI),
                      ('chr1:0-100', '+', 'Gypsy'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.GenomeLoci._DTYPE_LOCI),
                      ('chr1:0-100', '-', 'Gypsy'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.GenomeLoci._DTYPE_LOCI),
                      ('chr2:0-100', '+', 'Copia'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.GenomeLoci._DTYPE_LOCI),
                      ('chr2:0-100', '-', 'Copia'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.GenomeLoci._DTYPE_LOCI),
                      ('chr2:0-100', '+', 'Gypsy'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.GenomeLoci._DTYPE_LOCI),
                      ('chr2:0-100', '-', 'Gypsy'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.GenomeLoci._DTYPE_LOCI)}
        query = loci.GenomeLoci.from_dict(dictionary)

        answer = np.array([('chr1:0-100', '+', 'Copia', 1, 2),
                           ('chr1:0-100', '+', 'Copia', 2, 3),
                           ('chr1:0-100', '-', 'Copia', 1, 2),
                           ('chr1:0-100', '-', 'Copia', 2, 3),
                           ('chr1:0-100', '+', 'Gypsy', 1, 2),
                           ('chr1:0-100', '+', 'Gypsy', 2, 3),
                           ('chr1:0-100', '-', 'Gypsy', 1, 2),
                           ('chr1:0-100', '-', 'Gypsy', 2, 3),
                           ('chr2:0-100', '+', 'Copia', 1, 2),
                           ('chr2:0-100', '+', 'Copia', 2, 3),
                           ('chr2:0-100', '-', 'Copia', 1, 2),
                           ('chr2:0-100', '-', 'Copia', 2, 3),
                           ('chr2:0-100', '+', 'Gypsy', 1, 2),
                           ('chr2:0-100', '+', 'Gypsy', 2, 3),
                           ('chr2:0-100', '-', 'Gypsy', 1, 2),
                           ('chr2:0-100', '-', 'Gypsy', 2, 3)], dtype=loci.GenomeLoci._DTYPE_ARRAY)
        npt.assert_array_equal(np.sort(query.as_array()), np.sort(answer))

    def test_as_gff(self):
        """"""
        dictionary = {('chr1:0-100', '+', 'Copia'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.GenomeLoci._DTYPE_LOCI),
                      ('chr1:0-100', '-', 'Copia'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.GenomeLoci._DTYPE_LOCI),
                      ('chr1:0-100', '+', 'Gypsy'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.GenomeLoci._DTYPE_LOCI),
                      ('chr1:0-100', '-', 'Gypsy'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.GenomeLoci._DTYPE_LOCI),
                      ('chr2:0-100', '+', 'Copia'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.GenomeLoci._DTYPE_LOCI)}
        query = loci.GenomeLoci.from_dict(dictionary)

        answer = '\n'.join(['chr1\t.\t.\t1\t2\t.\t+\t.\tID=chr1_+_Copia_1;category=Copia',
                            'chr1\t.\t.\t1\t2\t.\t+\t.\tID=chr1_+_Gypsy_1;category=Gypsy',
                            'chr1\t.\t.\t1\t2\t.\t-\t.\tID=chr1_-_Copia_1;category=Copia',
                            'chr1\t.\t.\t1\t2\t.\t-\t.\tID=chr1_-_Gypsy_1;category=Gypsy',
                            'chr1\t.\t.\t2\t3\t.\t+\t.\tID=chr1_+_Copia_2;category=Copia',
                            'chr1\t.\t.\t2\t3\t.\t+\t.\tID=chr1_+_Gypsy_2;category=Gypsy',
                            'chr1\t.\t.\t2\t3\t.\t-\t.\tID=chr1_-_Copia_2;category=Copia',
                            'chr1\t.\t.\t2\t3\t.\t-\t.\tID=chr1_-_Gypsy_2;category=Gypsy',
                            'chr2\t.\t.\t1\t2\t.\t+\t.\tID=chr2_+_Copia_1;category=Copia',
                            'chr2\t.\t.\t2\t3\t.\t+\t.\tID=chr2_+_Copia_2;category=Copia'])
        assert query.as_gff() == answer


class TestReadLoci:
    """Tests for class ReadLoci"""
    def test_from_bam(self):
        """"""
        pass

    def test_tips(self):
        """Test for method to extract the tips of loci as a dict"""
        dictionary = {('chr1:0-100', '+', 'Copia', 'bam1'): np.array([(1, 2, 'a'), (2, 3, 'b')],
                                                                     dtype=loci.ReadLoci._DTYPE_LOCI),
                      ('chr1:0-100', '-', 'Copia', 'bam1'): np.array([(1, 2, 'c'), (2, 3, 'd')],
                                                                     dtype=loci.ReadLoci._DTYPE_LOCI),
                      ('chr1:0-100', '+', 'Gypsy', 'bam1'): np.array([(1, 2, 'e'), (2, 3, 'f')],
                                                                     dtype=loci.ReadLoci._DTYPE_LOCI),
                      ('chr1:0-100', '-', 'Gypsy', 'bam1'): np.array([(1, 2, 'g'), (2, 3, 'h')],
                                                                     dtype=loci.ReadLoci._DTYPE_LOCI),
                      ('chr2:0-100', '+', 'Copia', 'bam1'): np.array([(1, 2, 'i'), (2, 3, 'j')],
                                                                     dtype=loci.ReadLoci._DTYPE_LOCI)}
        query = loci.ReadLoci.from_dict(dictionary)
        query = query.tips()

        answer = {loci.LociKey('chr1:0-100', '+', 'Copia', 'bam1'): np.array([2, 3]),
                  loci.LociKey('chr1:0-100', '-', 'Copia', 'bam1'): np.array([1, 2]),
                  loci.LociKey('chr1:0-100', '+', 'Gypsy', 'bam1'): np.array([2, 3]),
                  loci.LociKey('chr1:0-100', '-', 'Gypsy', 'bam1'): np.array([1, 2]),
                  loci.LociKey('chr2:0-100', '+', 'Copia', 'bam1'): np.array([2, 3])}
        assert set(query.keys()) == set(answer.keys())
        for key in query.keys():
            npt.assert_array_equal(query[key], answer[key])

    def test_fingerprint(self):
        """Test for fingerprinting of read loci when one strand is empty (no reads)"""
        query = {('chr1:0-3000', '-', 'Gypsy', 'bam1'): np.array([], dtype=loci.ReadLoci._DTYPE_LOCI),
                 ('chr1:0-3000', '+', 'Gypsy', 'bam1'): np.array([(0,    0, 'Gypsy27_uniqueID'),
                                                                  (0,    0, 'Gypsy27_uniqueID'),
                                                                  (0,   60, 'Gypsy27_uniqueID'),
                                                                  (0,   61, 'Gypsy27_uniqueID'),
                                                                  (0,   61, 'Gypsy27_uniqueID'),
                                                                  (0,   61, 'Gypsy27_uniqueID'),
                                                                  (0,   76, 'Gypsy27_uniqueID'),
                                                                  (0,   78, 'Gypsy27_uniqueID'),
                                                                  (0,  122, 'Gypsy27_uniqueID'),
                                                                  (0,  122, 'Gypsy27_uniqueID'),
                                                                  (0,  141, 'Gypsy27_uniqueID'),
                                                                  (0,  183, 'Gypsy27_uniqueID'),
                                                                  (0,  251, 'Gypsy27_uniqueID'),
                                                                  (0,  260, 'Gypsy27_uniqueID'),
                                                                  (0,  260, 'Gypsy27_uniqueID'),
                                                                  (0,  263, 'Gypsy27_uniqueID'),
                                                                  (0,  263, 'Gypsy27_uniqueID'),
                                                                  (0,  267, 'Gypsy27_uniqueID'),
                                                                  (0,  267, 'Gypsy27_uniqueID'),
                                                                  (0,  288, 'Gypsy27_uniqueID'),
                                                                  (0,  288, 'Gypsy27_uniqueID'),
                                                                  (0,  295, 'Gypsy27_uniqueID'),
                                                                  (0,  300, 'Gypsy27_uniqueID'),
                                                                  (0,  310, 'Gypsy27_uniqueID'),
                                                                  (0,  310, 'Gypsy27_uniqueID'),
                                                                  (0,  317, 'Gypsy27_uniqueID'),
                                                                  (0,  317, 'Gypsy27_uniqueID'),
                                                                  (0,  334, 'Gypsy27_uniqueID'),
                                                                  (0,  334, 'Gypsy27_uniqueID'),
                                                                  (0,  335, 'Gypsy27_uniqueID'),
                                                                  (0,  338, 'Gypsy27_uniqueID'),
                                                                  (0,  338, 'Gypsy27_uniqueID'),
                                                                  (0,  338, 'Gypsy27_uniqueID'),
                                                                  (0,  338, 'Gypsy27_uniqueID'),
                                                                  (0,  340, 'Gypsy27_uniqueID'),
                                                                  (0,  342, 'Gypsy27_uniqueID'),
                                                                  (0,  342, 'Gypsy27_uniqueID'),
                                                                  (0,  344, 'Gypsy27_uniqueID'),
                                                                  (0,  344, 'Gypsy27_uniqueID'),
                                                                  (0,  358, 'Gypsy27_uniqueID'),
                                                                  (0,  367, 'Gypsy27_uniqueID'),
                                                                  (0,  370, 'Gypsy27_uniqueID'),
                                                                  (0,  370, 'Gypsy27_uniqueID'),
                                                                  (0,  377, 'Gypsy27_uniqueID'),
                                                                  (0,  387, 'Gypsy27_uniqueID'),
                                                                  (0,  402, 'Gypsy27_uniqueID'),
                                                                  (0,  403, 'Gypsy27_uniqueID'),
                                                                  (0,  410, 'Gypsy27_uniqueID'),
                                                                  (0,  410, 'Gypsy27_uniqueID'),
                                                                  (0,  410, 'Gypsy27_uniqueID'),
                                                                  (0,  418, 'Gypsy27_uniqueID'),
                                                                  (0,  418, 'Gypsy27_uniqueID'),
                                                                  (0,  424, 'Gypsy27_uniqueID'),
                                                                  (0,  424, 'Gypsy27_uniqueID'),
                                                                  (0,  577, 'Gypsy27_uniqueID'),
                                                                  (0,  857, 'Gypsy27_uniqueID'),
                                                                  (0,  879, 'Gypsy27_uniqueID'),
                                                                  (0,  921, 'Gypsy27_uniqueID'),
                                                                  (0,  921, 'Gypsy27_uniqueID'),
                                                                  (0, 1007, 'Gypsy27_uniqueID'),
                                                                  (0, 1031, 'Gypsy27_uniqueID'),
                                                                  (0, 1051, 'Gypsy27_uniqueID'),
                                                                  (0, 1051, 'Gypsy27_uniqueID'),
                                                                  (0, 1059, 'Gypsy27_uniqueID'),
                                                                  (0, 1071, 'Gypsy27_uniqueID'),
                                                                  (0, 1071, 'Gypsy27_uniqueID'),
                                                                  (0, 1080, 'Gypsy27_uniqueID'),
                                                                  (0, 1094, 'Gypsy27_uniqueID'),
                                                                  (0, 1094, 'Gypsy27_uniqueID'),
                                                                  (0, 1110, 'Gypsy27_uniqueID'),
                                                                  (0, 1110, 'Gypsy27_uniqueID'),
                                                                  (0, 1113, 'Gypsy27_uniqueID'),
                                                                  (0, 1113, 'Gypsy27_uniqueID'),
                                                                  (0, 1183, 'Gypsy27_uniqueID'),
                                                                  (0, 1189, 'Gypsy27_uniqueID'),
                                                                  (0, 1200, 'Gypsy27_uniqueID'),
                                                                  (0, 1200, 'Gypsy27_uniqueID'),
                                                                  (0, 1217, 'Gypsy27_uniqueID'),
                                                                  (0, 1234, 'Gypsy27_uniqueID'),
                                                                  (0, 1234, 'Gypsy27_uniqueID'),
                                                                  (0, 1591, 'Gypsy27_uniqueID'),
                                                                  (0, 1620, 'Gypsy27_uniqueID'),
                                                                  (0, 1620, 'Gypsy27_uniqueID'),
                                                                  (0, 1662, 'Gypsy27_uniqueID'),
                                                                  (0, 1686, 'Gypsy27_uniqueID'),
                                                                  (0, 1707, 'Gypsy27_uniqueID'),
                                                                  (0, 1755, 'Gypsy27_uniqueID'),
                                                                  (0, 1828, 'Gypsy27_uniqueID'),
                                                                  (0, 1828, 'Gypsy27_uniqueID'),
                                                                  (0, 1848, 'Gypsy27_uniqueID'),
                                                                  (0, 1848, 'Gypsy27_uniqueID'),
                                                                  (0, 1848, 'Gypsy27_uniqueID'),
                                                                  (0, 1848, 'Gypsy27_uniqueID'),
                                                                  (0, 1851, 'Gypsy27_uniqueID'),
                                                                  (0, 1851, 'Gypsy27_uniqueID'),
                                                                  (0, 1852, 'Gypsy27_uniqueID'),
                                                                  (0, 1917, 'Gypsy27_uniqueID')],
                                                                 dtype=loci.ReadLoci._DTYPE_LOCI)}
        query = loci.ReadLoci.from_dict(query)
        query = query.fingerprint(10, eps=200, min_eps=10, hierarchical=True)
        answer = loci.FingerPrint.from_dict({('chr1:0-3000',
                                              '-',
                                              'Gypsy',
                                              'bam1'): np.array([], dtype=loci.FingerPrint._DTYPE_LOCI),
                                             ('chr1:0-3000',
                                              '+',
                                              'Gypsy',
                                              'bam1'): np.array([(0, 577),
                                                                 (879, 1234),
                                                                 (1662, 1917)], dtype=loci.FingerPrint._DTYPE_LOCI)})
        assert set(query.keys()) == set(answer.keys())
        for key in query.keys():
            npt.assert_array_equal(query[key], answer[key])

    def test_as_array(self):
        query = loci.ReadLoci.from_dict({('chr1:0-3000',
                                          '+',
                                          'Gypsy',
                                          'bam1'): np.array([(0, 20, 'name1'),
                                                             (21, 30, 'name2'),
                                                             (42, 100, 'name3')], dtype=loci.ReadLoci._DTYPE_LOCI),
                                         ('chr1:0-3000',
                                          '-',
                                          'Gypsy',
                                          'bam1'): np.array([(0, 20, 'name4'),
                                                             (21, 30, 'name5'),
                                                             (42, 100, 'name6')], dtype=loci.ReadLoci._DTYPE_LOCI)})
        answer = np.array([('chr1:0-3000', '+', 'Gypsy', 'bam1', 0, 20, 'name1'),
                           ('chr1:0-3000', '-', 'Gypsy', 'bam1', 0, 20, 'name4'),
                           ('chr1:0-3000', '+', 'Gypsy', 'bam1', 21, 30, 'name2'),
                           ('chr1:0-3000', '-', 'Gypsy', 'bam1', 21, 30, 'name5'),
                           ('chr1:0-3000', '+', 'Gypsy', 'bam1', 42, 100, 'name3'),
                           ('chr1:0-3000', '-', 'Gypsy', 'bam1', 42, 100, 'name6')],
                          dtype=loci.ReadLoci._DTYPE_ARRAY)

        npt.assert_array_equal(np.sort(query.as_array()), np.sort(answer))

    def test_as_gff(self):
        """"""
        query = loci.ReadLoci.from_dict({('chr1:0-3000',
                                          '+',
                                          'Gypsy',
                                          'bam1'): np.array([(0, 20, 'name1'),
                                                             (21, 30, 'name2'),
                                                             (42, 100, 'name3')], dtype=loci.ReadLoci._DTYPE_LOCI),
                                         ('chr1:0-3000',
                                          '-',
                                          'Gypsy',
                                          'bam1'): np.array([(0, 20, 'name4'),
                                                             (21, 30, 'name5'),
                                                             (42, 100, 'name6')], dtype=loci.ReadLoci._DTYPE_LOCI)})
        answer = '\n'.join(['chr1\t.\t.\t0\t20\t.\t+\t.\tID=name1;category=Gypsy;source=bam1',
                            'chr1\t.\t.\t0\t20\t.\t-\t.\tID=name4;category=Gypsy;source=bam1',
                            'chr1\t.\t.\t21\t30\t.\t+\t.\tID=name2;category=Gypsy;source=bam1',
                            'chr1\t.\t.\t21\t30\t.\t-\t.\tID=name5;category=Gypsy;source=bam1',
                            'chr1\t.\t.\t42\t100\t.\t+\t.\tID=name3;category=Gypsy;source=bam1',
                            'chr1\t.\t.\t42\t100\t.\t-\t.\tID=name6;category=Gypsy;source=bam1'])
        assert query.as_gff() == answer


class TestFingerPrint:
    """Tests for class FingerPrint"""
    def test_as_array(self):
        """"""
        query = loci.FingerPrint.from_dict({('chr1:0-3000',
                                             '+',
                                             'Gypsy',
                                             'bam1'): np.array([(0, 577),
                                                                (879, 1234),
                                                                (1662, 1917)], dtype=loci.FingerPrint._DTYPE_LOCI),
                                            ('chr1:0-3000',
                                             '-',
                                             'Gypsy',
                                             'bam1'): np.array([(20, 570),
                                                                (870, 1230),
                                                                (1662, 1917)], dtype=loci.FingerPrint._DTYPE_LOCI),
                                            ('chr2:0-4000',
                                             '+',
                                             'Gypsy',
                                             'bam1'): np.array([(0, 577),
                                                                (879, 1234),
                                                                (1662, 1917)], dtype=loci.FingerPrint._DTYPE_LOCI),
                                            ('chr2:0-4000',
                                             '-',
                                             'Gypsy',
                                             'bam1'): np.array([], dtype=loci.FingerPrint._DTYPE_LOCI)
                                            })
        answer = np.array([('chr1:0-3000', '+', 'Gypsy', 'bam1', 0, 577),
                           ('chr1:0-3000', '-', 'Gypsy', 'bam1', 20, 570),
                           ('chr1:0-3000', '-', 'Gypsy', 'bam1', 870, 1230),
                           ('chr1:0-3000', '+', 'Gypsy', 'bam1', 879, 1234),
                           ('chr1:0-3000', '+', 'Gypsy', 'bam1', 1662, 1917),
                           ('chr1:0-3000', '-', 'Gypsy', 'bam1', 1662, 1917),
                           ('chr2:0-4000', '+', 'Gypsy', 'bam1', 0, 577),
                           ('chr2:0-4000', '+', 'Gypsy', 'bam1', 879, 1234),
                           ('chr2:0-4000', '+', 'Gypsy', 'bam1', 1662, 1917)],
                          dtype=loci.FingerPrint._DTYPE_ARRAY)
        npt.assert_array_equal(np.sort(query.as_array()), np.sort(answer))

    def test_as_gff(self):
        """"""
        query = loci.FingerPrint.from_dict({('chr1:0-3000',
                                             '+',
                                             'Gypsy',
                                             'bam1'): np.array([(0, 577),
                                                                (879, 1234),
                                                                (1662, 1917)], dtype=loci.FingerPrint._DTYPE_LOCI),
                                            ('chr1:0-3000',
                                             '-',
                                             'Gypsy',
                                             'bam1'): np.array([(20, 570),
                                                                (870, 1230),
                                                                (1662, 1917)], dtype=loci.FingerPrint._DTYPE_LOCI),
                                            ('chr2:0-4000',
                                             '+',
                                             'Gypsy',
                                             'bam1'): np.array([(0, 577),
                                                                (879, 1234),
                                                                (1662, 1917)], dtype=loci.FingerPrint._DTYPE_LOCI),
                                            ('chr2:0-4000',
                                             '-',
                                             'Gypsy',
                                             'bam1'): np.array([], dtype=loci.FingerPrint._DTYPE_LOCI)
                                            })
        answer = '\n'.join(['chr1\t.\t.\t0\t577\t.\t+\t.\tID=chr1_+_Gypsy_0;category=Gypsy;source=bam1',
                            'chr1\t.\t.\t20\t570\t.\t-\t.\tID=chr1_-_Gypsy_20;category=Gypsy;source=bam1',
                            'chr1\t.\t.\t870\t1230\t.\t-\t.\tID=chr1_-_Gypsy_870;category=Gypsy;source=bam1',
                            'chr1\t.\t.\t879\t1234\t.\t+\t.\tID=chr1_+_Gypsy_879;category=Gypsy;source=bam1',
                            'chr1\t.\t.\t1662\t1917\t.\t+\t.\tID=chr1_+_Gypsy_1662;category=Gypsy;source=bam1',
                            'chr1\t.\t.\t1662\t1917\t.\t-\t.\tID=chr1_-_Gypsy_1662;category=Gypsy;source=bam1',
                            'chr2\t.\t.\t0\t577\t.\t+\t.\tID=chr2_+_Gypsy_0;category=Gypsy;source=bam1',
                            'chr2\t.\t.\t879\t1234\t.\t+\t.\tID=chr2_+_Gypsy_879;category=Gypsy;source=bam1',
                            'chr2\t.\t.\t1662\t1917\t.\t+\t.\tID=chr2_+_Gypsy_1662;category=Gypsy;source=bam1'])
        assert query.as_gff() == answer


class TestComparativeBins:
    """Tests for class ComparativeBins"""
    def test_from_fingerprints(self):
        """"""
        fprint1 = loci.FingerPrint.from_dict({('chr1:0-300',
                                               '+',
                                               'Gypsy',
                                               'bam1'): np.array([(0, 10),
                                                                  (21, 30)],
                                                                 dtype=loci.FingerPrint._DTYPE_LOCI)})
        fprint2 = loci.FingerPrint.from_dict({('chr1:0-300',
                                               '+',
                                               'Gypsy',
                                               'bam2'): np.array([(10, 20),
                                                                  (35, 40)],
                                                                 dtype=loci.FingerPrint._DTYPE_LOCI)})
        query = loci.ComparativeBins.from_fingerprints(fprint1, fprint2)
        answer = loci.ComparativeBins.from_dict({('chr1:0-300',
                                                  '+',
                                                  'Gypsy'): np.array([(0, 20),
                                                                      (21, 30),
                                                                      (35, 40)],
                                                                     dtype=loci.ComparativeBins._DTYPE_LOCI)})
        assert set(query.keys()) == set(answer.keys())
        for key in query.keys():
            npt.assert_array_equal(query[key], answer[key])

    def test_compare(self):
        """Test the comparison of read counts based on comparative bins"""
        bins = loci.ComparativeBins.from_dict({('chr1:0-3000',
                                                '+',
                                                'Gypsy'): np.array([(0, 20),
                                                                    (21, 30),
                                                                    (35, 40)],
                                                                   dtype=loci.ComparativeBins._DTYPE_LOCI),
                                               ('chr1:0-3000',
                                                '+',
                                                'Copia'): np.array([(0, 20),
                                                                    (21, 30),
                                                                    (35, 40)],
                                                                   dtype=loci.ComparativeBins._DTYPE_LOCI)
                                               })
        reads = {('chr1:0-3000', '+', 'Gypsy', 'bam1'): np.array([(0,  0, 'Gypsy_uniqueID'),
                                                                  (0,  2, 'Gypsy_uniqueID'),
                                                                  (0, 20, 'Gypsy_uniqueID'),
                                                                  (0, 21, 'Gypsy_uniqueID'),
                                                                  (0, 27, 'Gypsy_uniqueID'),
                                                                  (0, 28, 'Gypsy_uniqueID'),
                                                                  (0, 28, 'Gypsy_uniqueID'),
                                                                  (0, 30, 'Gypsy_uniqueID'),
                                                                  (0, 33, 'Gypsy_uniqueID'),
                                                                  (0, 35, 'Gypsy_uniqueID'),
                                                                  (0, 36, 'Gypsy_uniqueID'),
                                                                  (0, 40, 'Gypsy_uniqueID'),
                                                                  (0, 70, 'Gypsy_uniqueID')],
                                                                 dtype=loci.ReadLoci._DTYPE_LOCI),
                 ('chr1:0-3000', '+', 'Gypsy', 'bam2'): np.array([(0,  2, 'Gypsy_uniqueID'),
                                                                  (0, 21, 'Gypsy_uniqueID'),
                                                                  (0, 27, 'Gypsy_uniqueID'),
                                                                  (0, 28, 'Gypsy_uniqueID'),
                                                                  (0, 28, 'Gypsy_uniqueID'),
                                                                  (0, 30, 'Gypsy_uniqueID'),
                                                                  (0, 70, 'Gypsy_uniqueID')],
                                                                 dtype=loci.ReadLoci._DTYPE_LOCI),
                 ('chr1:0-3000', '+', 'Copia', 'bam1'): np.array([(0,  0, 'Gypsy_uniqueID'),
                                                                  (0,  2, 'Gypsy_uniqueID'),
                                                                  (0, 20, 'Gypsy_uniqueID'),
                                                                  (0, 21, 'Gypsy_uniqueID'),
                                                                  (0, 27, 'Gypsy_uniqueID'),
                                                                  (0, 28, 'Gypsy_uniqueID'),
                                                                  (0, 28, 'Gypsy_uniqueID'),
                                                                  (0, 30, 'Gypsy_uniqueID'),
                                                                  (0, 33, 'Gypsy_uniqueID'),
                                                                  (0, 35, 'Gypsy_uniqueID'),
                                                                  (0, 36, 'Gypsy_uniqueID'),
                                                                  (0, 40, 'Gypsy_uniqueID'),
                                                                  (0, 70, 'Gypsy_uniqueID')],
                                                                 dtype=loci.ReadLoci._DTYPE_LOCI),
                 ('chr1:0-3000', '+', 'Copia', 'bam2'): np.array([],
                                                                 dtype=loci.ReadLoci._DTYPE_LOCI)}
        reads = loci.ReadLoci.from_dict(reads)
        query = bins.compare(reads)

        answer = {('chr1:0-3000', '+', 'Copia'): np.array([(0, 20, ['bam1', 'bam2'], [3, 0], [1.0, 0.0]),
                                                           (21, 30, ['bam1', 'bam2'], [5, 0], [1.0, 0.0]),
                                                           (35, 40, ['bam1', 'bam2'], [3, 0], [1.0, 0.0])],
                                                          dtype=loci.Comparison._DTYPE_LOCI),
                  ('chr1:0-3000', '+', 'Gypsy'): np.array([(0, 20, ['bam1', 'bam2'], [3, 1], [0.75, 0.25]),
                                                           (21, 30, ['bam1', 'bam2'], [5, 5], [0.5, 0.5]),
                                                           (35, 40, ['bam1', 'bam2'], [3, 0], [1.0, 0.0])],
                                                          dtype=loci.Comparison._DTYPE_LOCI)}
        answer = loci.Comparison.from_dict(answer)

        assert set(query.keys()) == set(answer.keys())
        for key in query.keys():
            npt.assert_array_equal(query[key], answer[key])

    def test_as_array(self):
        """"""
        dictionary = {('chr1:0-100', '+', 'Copia'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.ComparativeBins._DTYPE_LOCI),
                      ('chr1:0-100', '-', 'Copia'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.ComparativeBins._DTYPE_LOCI),
                      ('chr1:0-100', '+', 'Gypsy'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.ComparativeBins._DTYPE_LOCI),
                      ('chr1:0-100', '-', 'Gypsy'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.ComparativeBins._DTYPE_LOCI),
                      ('chr2:0-100', '+', 'Copia'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.ComparativeBins._DTYPE_LOCI),
                      ('chr2:0-100', '-', 'Copia'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.ComparativeBins._DTYPE_LOCI),
                      ('chr2:0-100', '+', 'Gypsy'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.ComparativeBins._DTYPE_LOCI),
                      ('chr2:0-100', '-', 'Gypsy'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.ComparativeBins._DTYPE_LOCI)}
        query = loci.ComparativeBins.from_dict(dictionary)

        answer = np.array([('chr1:0-100', '+', 'Copia', 1, 2),
                           ('chr1:0-100', '+', 'Copia', 2, 3),
                           ('chr1:0-100', '-', 'Copia', 1, 2),
                           ('chr1:0-100', '-', 'Copia', 2, 3),
                           ('chr1:0-100', '+', 'Gypsy', 1, 2),
                           ('chr1:0-100', '+', 'Gypsy', 2, 3),
                           ('chr1:0-100', '-', 'Gypsy', 1, 2),
                           ('chr1:0-100', '-', 'Gypsy', 2, 3),
                           ('chr2:0-100', '+', 'Copia', 1, 2),
                           ('chr2:0-100', '+', 'Copia', 2, 3),
                           ('chr2:0-100', '-', 'Copia', 1, 2),
                           ('chr2:0-100', '-', 'Copia', 2, 3),
                           ('chr2:0-100', '+', 'Gypsy', 1, 2),
                           ('chr2:0-100', '+', 'Gypsy', 2, 3),
                           ('chr2:0-100', '-', 'Gypsy', 1, 2),
                           ('chr2:0-100', '-', 'Gypsy', 2, 3)], dtype=loci.ComparativeBins._DTYPE_ARRAY)
        npt.assert_array_equal(np.sort(query.as_array()), np.sort(answer))

    def test_as_gff(self):
        """"""
        dictionary = {('chr1:0-100', '+', 'Copia'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.ComparativeBins._DTYPE_LOCI),
                      ('chr1:0-100', '-', 'Copia'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.ComparativeBins._DTYPE_LOCI),
                      ('chr1:0-100', '+', 'Gypsy'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.ComparativeBins._DTYPE_LOCI),
                      ('chr1:0-100', '-', 'Gypsy'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.ComparativeBins._DTYPE_LOCI),
                      ('chr2:0-100', '+', 'Copia'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci.ComparativeBins._DTYPE_LOCI)}
        query = loci.ComparativeBins.from_dict(dictionary)

        answer = '\n'.join(['chr1\t.\t.\t1\t2\t.\t+\t.\tID=chr1_+_Copia_1;category=Copia',
                            'chr1\t.\t.\t1\t2\t.\t+\t.\tID=chr1_+_Gypsy_1;category=Gypsy',
                            'chr1\t.\t.\t1\t2\t.\t-\t.\tID=chr1_-_Copia_1;category=Copia',
                            'chr1\t.\t.\t1\t2\t.\t-\t.\tID=chr1_-_Gypsy_1;category=Gypsy',
                            'chr1\t.\t.\t2\t3\t.\t+\t.\tID=chr1_+_Copia_2;category=Copia',
                            'chr1\t.\t.\t2\t3\t.\t+\t.\tID=chr1_+_Gypsy_2;category=Gypsy',
                            'chr1\t.\t.\t2\t3\t.\t-\t.\tID=chr1_-_Copia_2;category=Copia',
                            'chr1\t.\t.\t2\t3\t.\t-\t.\tID=chr1_-_Gypsy_2;category=Gypsy',
                            'chr2\t.\t.\t1\t2\t.\t+\t.\tID=chr2_+_Copia_1;category=Copia',
                            'chr2\t.\t.\t2\t3\t.\t+\t.\tID=chr2_+_Copia_2;category=Copia'])
        assert query.as_gff() == answer


class TestComparison:
    """Tests for class Comparison"""
    def test_as_array(self):
        """"""
        query = {('chr1:0-3000', '+', 'Copia'): np.array([(0, 20, ['bam1', 'bam2'], [3, 0], [1.0, 0.0]),
                                                          (21, 30, ['bam1', 'bam2'], [5, 0], [1.0, 0.0]),
                                                          (35, 40, ['bam1', 'bam2'], [3, 0], [1.0, 0.0])],
                                                         dtype=loci.Comparison._DTYPE_LOCI),
                 ('chr1:0-3000', '+', 'Gypsy'): np.array([(0, 20, ['bam1', 'bam2'], [3, 1], [0.75, 0.25]),
                                                          (21, 30, ['bam1', 'bam2'], [5, 5], [0.5, 0.5]),
                                                          (35, 40, ['bam1', 'bam2'], [3, 0], [1.0, 0.0])],
                                                         dtype=loci.Comparison._DTYPE_LOCI)}
        query = loci.Comparison.from_dict(query)

        answer = np.array([('chr1:0-3000', '+', 'Copia', 0, 20, ['bam1', 'bam2'], [3, 0], [1.0, 0.0]),
                           ('chr1:0-3000', '+', 'Gypsy', 0, 20, ['bam1', 'bam2'], [3, 1], [0.75, 0.25]),
                           ('chr1:0-3000', '+', 'Copia', 21, 30, ['bam1', 'bam2'], [5, 0], [1.0, 0.0]),
                           ('chr1:0-3000', '+', 'Gypsy', 21, 30, ['bam1', 'bam2'], [5, 5], [0.5, 0.5]),
                           ('chr1:0-3000', '+', 'Copia', 35, 40, ['bam1', 'bam2'], [3, 0], [1.0, 0.0]),
                           ('chr1:0-3000', '+', 'Gypsy', 35, 40, ['bam1', 'bam2'], [3, 0], [1.0, 0.0])],
                          dtype=loci.Comparison._DTYPE_ARRAY)

        npt.assert_array_equal(query.as_array(), answer)

    def test_as_gff(self):
        """"""
        query = {('chr1:0-3000', '+', 'Copia'): np.array([(0, 20, ['bam1', 'bam2'], [3, 0], [1.0, 0.0]),
                                                          (21, 30, ['bam1', 'bam2'], [5, 0], [1.0, 0.0]),
                                                          (35, 40, ['bam1', 'bam2'], [3, 0], [1.0, 0.0])],
                                                         dtype=loci.Comparison._DTYPE_LOCI),
                 ('chr1:0-3000', '+', 'Gypsy'): np.array([(0, 20, ['bam1', 'bam2'], [3, 1], [0.75, 0.25]),
                                                          (21, 30, ['bam1', 'bam2'], [5, 5], [0.5, 0.5]),
                                                          (35, 40, ['bam1', 'bam2'], [3, 0], [1.0, 0.0])],
                                                         dtype=loci.Comparison._DTYPE_LOCI)}
        query = loci.Comparison.from_dict(query)

        answer = '\n'.join(['chr1\t.\t.\t0\t20\t.\t+\t.\tID=chr1_+_Copia_0;category=Copia;sources=bam1,bam2;counts=3,0;proportions=1.0,0.0;color=#08306B',
                            'chr1\t.\t.\t0\t20\t.\t+\t.\tID=chr1_+_Gypsy_0;category=Gypsy;sources=bam1,bam2;counts=3,1;proportions=0.75,0.25;color=#6BAED6',
                            'chr1\t.\t.\t21\t30\t.\t+\t.\tID=chr1_+_Copia_21;category=Copia;sources=bam1,bam2;counts=5,0;proportions=1.0,0.0;color=#08306B',
                            'chr1\t.\t.\t21\t30\t.\t+\t.\tID=chr1_+_Gypsy_21;category=Gypsy;sources=bam1,bam2;counts=5,5;proportions=0.5,0.5;color=#C6DBEF',
                            'chr1\t.\t.\t35\t40\t.\t+\t.\tID=chr1_+_Copia_35;category=Copia;sources=bam1,bam2;counts=3,0;proportions=1.0,0.0;color=#08306B',
                            'chr1\t.\t.\t35\t40\t.\t+\t.\tID=chr1_+_Gypsy_35;category=Gypsy;sources=bam1,bam2;counts=3,0;proportions=1.0,0.0;color=#08306B'])
        assert query.as_gff() == answer

    def test_flat_array(self):
        """"""
        query = {('chr1:0-3000', '+', 'Copia'): np.array([(0, 20, ['bam1', 'bam2'], [3, 0], [1.0, 0.0]),
                                                          (21, 30, ['bam1', 'bam2'], [5, 0], [1.0, 0.0]),
                                                          (35, 40, ['bam1', 'bam2'], [3, 0], [1.0, 0.0])],
                                                         dtype=loci.Comparison._DTYPE_LOCI),
                 ('chr1:0-3000', '+', 'Gypsy'): np.array([(0, 20, ['bam1', 'bam2'], [3, 1], [0.75, 0.25]),
                                                          (21, 30, ['bam1', 'bam2'], [5, 5], [0.5, 0.5]),
                                                          (35, 40, ['bam1', 'bam2'], [3, 0], [1.0, 0.0])],
                                                         dtype=loci.Comparison._DTYPE_LOCI)}
        query = loci.Comparison.from_dict(query)

        answer = np.array([('chr1:0-3000', '+', 'Copia', 0, 20, 'bam1', 3, 1.0),
                           ('chr1:0-3000', '+', 'Copia', 0, 20, 'bam2', 0, 0.0),
                           ('chr1:0-3000', '+', 'Gypsy', 0, 20, 'bam1', 3, 0.75),
                           ('chr1:0-3000', '+', 'Gypsy', 0, 20, 'bam2', 1, 0.25),
                           ('chr1:0-3000', '+', 'Copia', 21, 30, 'bam1', 5, 1.0),
                           ('chr1:0-3000', '+', 'Copia', 21, 30, 'bam2', 0, 0.0),
                           ('chr1:0-3000', '+', 'Gypsy', 21, 30, 'bam1', 5, 0.5),
                           ('chr1:0-3000', '+', 'Gypsy', 21, 30, 'bam2', 5, 0.5),
                           ('chr1:0-3000', '+', 'Copia', 35, 40, 'bam1', 3, 1.0),
                           ('chr1:0-3000', '+', 'Copia', 35, 40, 'bam2', 0, 0.0),
                           ('chr1:0-3000', '+', 'Gypsy', 35, 40, 'bam1', 3, 1.0),
                           ('chr1:0-3000', '+', 'Gypsy', 35, 40, 'bam2', 0, 0.0)],
                          dtype=loci.Comparison._DTYPE_FLAT_ARRAY)

        npt.assert_array_equal(query.as_flat_array(), answer)

    def test_as_flat_gff(self):
        """"""
        query = {('chr1:0-3000', '+', 'Copia'): np.array([(0, 20, ['bam1', 'bam2'], [3, 0], [1.0, 0.0]),
                                                          (21, 30, ['bam1', 'bam2'], [5, 0], [1.0, 0.0]),
                                                          (35, 40, ['bam1', 'bam2'], [3, 0], [1.0, 0.0])],
                                                         dtype=loci.Comparison._DTYPE_LOCI),
                 ('chr1:0-3000', '+', 'Gypsy'): np.array([(0, 20, ['bam1', 'bam2'], [3, 1], [0.75, 0.25]),
                                                          (21, 30, ['bam1', 'bam2'], [5, 5], [0.5, 0.5]),
                                                          (35, 40, ['bam1', 'bam2'], [3, 0], [1.0, 0.0])],
                                                         dtype=loci.Comparison._DTYPE_LOCI)}
        query = loci.Comparison.from_dict(query)

        answer = '\n'.join(['chr1\t.\t.\t0\t20\t.\t+\t.\tID=chr1_+_Copia_0_0;category=Copia;source=bam1;count=3;proportion=1.0;color=#08306B',
                            'chr1\t.\t.\t0\t20\t.\t+\t.\tID=chr1_+_Copia_0_1;category=Copia;source=bam2;count=0;proportion=0.0;color=#C6DBEF',
                            'chr1\t.\t.\t0\t20\t.\t+\t.\tID=chr1_+_Gypsy_0_0;category=Gypsy;source=bam1;count=3;proportion=0.75;color=#6BAED6',
                            'chr1\t.\t.\t0\t20\t.\t+\t.\tID=chr1_+_Gypsy_0_1;category=Gypsy;source=bam2;count=1;proportion=0.25;color=#C6DBEF',
                            'chr1\t.\t.\t21\t30\t.\t+\t.\tID=chr1_+_Copia_21_0;category=Copia;source=bam1;count=5;proportion=1.0;color=#08306B',
                            'chr1\t.\t.\t21\t30\t.\t+\t.\tID=chr1_+_Copia_21_1;category=Copia;source=bam2;count=0;proportion=0.0;color=#C6DBEF',
                            'chr1\t.\t.\t21\t30\t.\t+\t.\tID=chr1_+_Gypsy_21_0;category=Gypsy;source=bam1;count=5;proportion=0.5;color=#C6DBEF',
                            'chr1\t.\t.\t21\t30\t.\t+\t.\tID=chr1_+_Gypsy_21_1;category=Gypsy;source=bam2;count=5;proportion=0.5;color=#C6DBEF',
                            'chr1\t.\t.\t35\t40\t.\t+\t.\tID=chr1_+_Copia_35_0;category=Copia;source=bam1;count=3;proportion=1.0;color=#08306B',
                            'chr1\t.\t.\t35\t40\t.\t+\t.\tID=chr1_+_Copia_35_1;category=Copia;source=bam2;count=0;proportion=0.0;color=#C6DBEF',
                            'chr1\t.\t.\t35\t40\t.\t+\t.\tID=chr1_+_Gypsy_35_0;category=Gypsy;source=bam1;count=3;proportion=1.0;color=#08306B',
                            'chr1\t.\t.\t35\t40\t.\t+\t.\tID=chr1_+_Gypsy_35_1;category=Gypsy;source=bam2;count=0;proportion=0.0;color=#C6DBEF'])
        assert query.as_flat_gff() == answer
