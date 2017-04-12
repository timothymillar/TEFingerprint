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
    query = np.array(query, dtype=loci._Loci._DTYPE_LOCI)
    query = np.fromiter(loci._loci_melter(query), dtype=loci._Loci._DTYPE_LOCI)
    answer = np.array(answer, dtype=loci._Loci._DTYPE_LOCI)
    npt.assert_array_equal(query, answer)


class TestLoci:
    """Tests for class _Loci"""
    def test_split(self):
        """Split into separate _Loci objects based on keys"""
        dictionary = {('chr1:0-100', '+', 'Copia'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                             dtype=loci._Loci._DTYPE_LOCI),
                      ('chr1:0-100', '-', 'Copia'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                             dtype=loci._Loci._DTYPE_LOCI),
                      ('chr1:0-100', '+', 'Gypsy'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                             dtype=loci._Loci._DTYPE_LOCI),
                      ('chr1:0-100', '-', 'Gypsy'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                             dtype=loci._Loci._DTYPE_LOCI),
                      ('chr2:0-100', '+', 'Copia'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                             dtype=loci._Loci._DTYPE_LOCI),
                      ('chr2:0-100', '-', 'Copia'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                             dtype=loci._Loci._DTYPE_LOCI),
                      ('chr2:0-100', '+', 'Gypsy'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                             dtype=loci._Loci._DTYPE_LOCI),
                      ('chr2:0-100', '-', 'Gypsy'): np.array([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
                                                             dtype=loci._Loci._DTYPE_LOCI)}
        query = loci._Loci.from_dict(dictionary)
        assert len(list(query.split())) == 8

    def test_as_array(self):
        """"""
        dictionary = {('chr1:0-100', '+', 'Copia'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci._Loci._DTYPE_LOCI),
                      ('chr1:0-100', '-', 'Copia'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci._Loci._DTYPE_LOCI),
                      ('chr1:0-100', '+', 'Gypsy'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci._Loci._DTYPE_LOCI),
                      ('chr1:0-100', '-', 'Gypsy'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci._Loci._DTYPE_LOCI),
                      ('chr2:0-100', '+', 'Copia'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci._Loci._DTYPE_LOCI),
                      ('chr2:0-100', '-', 'Copia'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci._Loci._DTYPE_LOCI),
                      ('chr2:0-100', '+', 'Gypsy'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci._Loci._DTYPE_LOCI),
                      ('chr2:0-100', '-', 'Gypsy'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci._Loci._DTYPE_LOCI)}
        query = loci._Loci.from_dict(dictionary)

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
                           ('chr2:0-100', '-', 'Gypsy', 2, 3)], dtype=loci._Loci._DTYPE_ARRAY)
        npt.assert_array_equal(np.sort(query.as_array()), np.sort(answer))

    def test_as_gff(self):
        """"""
        dictionary = {('chr1:0-100', '+', 'Copia'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci._Loci._DTYPE_LOCI),
                      ('chr1:0-100', '-', 'Copia'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci._Loci._DTYPE_LOCI),
                      ('chr1:0-100', '+', 'Gypsy'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci._Loci._DTYPE_LOCI),
                      ('chr1:0-100', '-', 'Gypsy'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci._Loci._DTYPE_LOCI),
                      ('chr2:0-100', '+', 'Copia'): np.array([(1, 2), (2, 3)],
                                                             dtype=loci._Loci._DTYPE_LOCI)}
        query = loci._Loci.from_dict(dictionary)

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
