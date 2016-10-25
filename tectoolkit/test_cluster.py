#! /usr/bin/env python

import numpy as np
import numpy.testing as npt
from tectoolkit.cluster import _UnivariateLoci, FlatUnivariateDensityCluster, HierarchicalUnivariateDensityCluster



class TestUL:
    """"""
    def test_melt_uloci(self):
        """"""
        input = np.array([(3, 6),
                          (6, 8),
                          (7, 9),
                          (10, 12),
                          (13, 13),
                          (15, 25),
                          (16, 17),
                          (19, 20)], dtype=_UnivariateLoci._ulocus)
        answer = np.array([(3, 9),
                           (10, 12),
                           (13, 13),
                           (15, 25)], dtype=_UnivariateLoci._ulocus)
        output = _UnivariateLoci._melt_uloci(input)
        npt.assert_array_equal(output, answer)

    def test_sort_uloci(self):
        """"""
        query = _UnivariateLoci()
        query.loci = np.array([(2, 4),
                               (3, 4),
                               (3, 3),
                               (4, 4),
                               (3, 99),
                               (1, 1)], dtype=_UnivariateLoci._ulocus)
        answer = _UnivariateLoci()
        answer.loci = np.array([(1, 1),
                                (2, 4),
                                (3, 3),
                                (3, 4),
                                (3, 99),
                                (4, 4)], dtype=_UnivariateLoci._ulocus)
        query._sort_uloci()
        npt.assert_array_equal(query.loci, answer.loci)
