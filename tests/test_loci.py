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

