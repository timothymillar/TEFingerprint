#! /usr/bin/env python3
import os
import pytest
import subprocess
from operator import itemgetter

DATA = os.path.dirname(os.path.realpath(__file__)) + '/data/testFingerprint-2017-12-06.gff'


with open(DATA) as f:
    SAMPLES = [i.strip() for i in f]


@pytest.mark.parametrize('query,answer',
                         [(["tef-filter-gff", DATA,
                            "--all"],
                           SAMPLES),
                          (["tef-filter-gff", DATA,
                            "--all", "end>40000"],
                           list(itemgetter(6, 7, 8, 9)(SAMPLES))),
                          (["tef-filter-gff", DATA,
                            "--all", "end>40000", "strand==-"],
                           list(itemgetter(7, 9)(SAMPLES))),
                          (["tef-filter-gff", DATA,
                            "--any", "end>40000", "strand==-"],
                           list(itemgetter(1, 3, 5, 6, 7, 8, 9)(SAMPLES))),
                          (["tef-filter-gff", DATA,
                            "--all", "sample_0_element_*_count>=8"],
                           list(itemgetter(4, 7)(SAMPLES))),
                          (["tef-filter-gff", DATA,
                            "--any", "sample_0_element_*_count>=8"],
                           list(itemgetter(0, 2, 3, 4, 5, 6, 7, 8, 9)(SAMPLES))),
                          (["tef-filter-gff", DATA,
                            "--all", "pair!=.", "strand=+",
                            "--any", "sample_?_count<50"],
                           list(itemgetter(0, 2, 8)(SAMPLES)))
                          ],
                         ids=['no-filter',
                              'end>40000',
                              '--all end>40000 strand==-',
                              '--any end>40000 strand==-',
                              '--all sample_0_element_*_count>=8',
                              '--any sample_0_element_*_count>=8',
                              '--all and --any'])
def test_filter_gff(query, answer):
    output, error = subprocess.Popen(query, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE).communicate()

    assert output.decode().splitlines() == answer



