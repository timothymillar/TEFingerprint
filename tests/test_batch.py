#! /usr/bin/env python

import os
import pytest
import numpy as np
import numpy.testing as npt
from tefingerprint import loci
from tefingerprint import batch

DATA_PATH = os.path.dirname(os.path.realpath(__file__)) + '/data/'


class TestFingerprint:
    """Integration tests for batch fingerprinting"""
    CORES = 1

    @pytest.mark.parametrize("parameters, answer",
                             [({'bams': [DATA_PATH + 'testB-2017-06-08.bam'],
                                'categories':['Gypsy', 'PIF'],
                                'references':['chr1'],
                                'quality':30,
                                'transposon_tag':'ME',
                                'min_reads':5,
                                'eps':250,
                                'min_eps':0,
                                'hierarchical':True},
                               {('chr1:1-25000', '+', 'PIF', 'testB-2017-06-08.bam'):
                                    np.array([],
                                             dtype=loci.FingerPrint._DTYPE_LOCI),
                                ('chr1:1-25000', '+', 'Gypsy', 'testB-2017-06-08.bam'):
                                    np.array([(2452, 2577), (2841, 3138), (24065, 24217)],
                                             dtype=loci.FingerPrint._DTYPE_LOCI),
                                ('chr1:1-25000', '-', 'Gypsy', 'testB-2017-06-08.bam'):
                                    np.array([(24850, 24919)],
                                             dtype=loci.FingerPrint._DTYPE_LOCI),
                                ('chr1:1-25000', '-', 'PIF', 'testB-2017-06-08.bam'):
                                    np.array([(21834, 21982)],
                                                          dtype=loci.FingerPrint._DTYPE_LOCI)}),
                              ({'bams': [DATA_PATH + 'testB-2017-06-08.bam'],
                                'categories': ['Gypsy', 'PIF'],
                                'references': ['chr1'],
                                'quality': 60,
                                'transposon_tag': 'ME',
                                'min_reads': 5,
                                'eps': 250,
                                'min_eps': 0,
                                'hierarchical': True},
                               {('chr1:1-25000', '+', 'PIF', 'testB-2017-06-08.bam'):
                                    np.array([],
                                             dtype=loci.FingerPrint._DTYPE_LOCI),
                                ('chr1:1-25000', '+', 'Gypsy', 'testB-2017-06-08.bam'):
                                    np.array([(2452, 2577), (2841, 3138), (24065, 24217)],
                                             dtype=loci.FingerPrint._DTYPE_LOCI),
                                ('chr1:1-25000', '-', 'Gypsy', 'testB-2017-06-08.bam'):
                                    np.array([],
                                             dtype=loci.FingerPrint._DTYPE_LOCI),
                                ('chr1:1-25000', '-', 'PIF', 'testB-2017-06-08.bam'):
                                    np.array([(21834, 21982)],
                                             dtype=loci.FingerPrint._DTYPE_LOCI)}),
                              ({'bams': [DATA_PATH + 'testB-2017-06-08.bam'],
                                'categories': ['INVALID'],
                                'references': ['chr1'],
                                'quality': 30,
                                'transposon_tag': 'ME',
                                'min_reads': 5,
                                'eps': 250,
                                'min_eps': 0,
                                'hierarchical': True},
                               {('chr1:1-25000', '+', 'INVALID', 'testB-2017-06-08.bam'):
                                    np.array([],
                                             dtype=loci.FingerPrint._DTYPE_LOCI),
                                ('chr1:1-25000', '-', 'INVALID', 'testB-2017-06-08.bam'):
                                    np.array([],
                                             dtype=loci.FingerPrint._DTYPE_LOCI)})
                              ])
    def test_fingerprint(self, parameters, answer):
        parameters['cores'] = self.CORES
        query = batch.fingerprint(**parameters)
        answer = loci.FingerPrint.from_dict(answer)

        assert set(query.keys()) == set(answer.keys())
        for key in query.keys():
            npt.assert_array_equal(query[key], answer[key])


class TestFingerprintMulticore(TestFingerprint):
    """Multi-core versions of integration tests for batch fingerprinting"""
    CORES = 2


class TestComparison:
    """Integration tests for batch fingerprinting"""
    CORES = 1

    @pytest.mark.parametrize("parameters, answer",
                             [({'bams': [DATA_PATH + 'testA-2017-06-08.bam',
                                         DATA_PATH + 'testB-2017-06-08.bam'],
                                'categories': ['Gypsy', 'PIF'],
                                'references': ['chr1'],
                                'quality': 30,
                                'transposon_tag': 'ME',
                                'min_reads': 5,
                                'eps': 250,
                                'min_eps': 0,
                                'hierarchical': True,
                                'fingerprint_buffer': 50,
                                'bin_buffer': 0},
                               {('chr1:1-25000', '-', 'PIF'): np.array([(21784, 22032,
                                                                         ['testA-2017-06-08.bam',
                                                                          'testB-2017-06-08.bam'],
                                                                         [4, 6],
                                                                         [0.4, 0.6])],
                                                                       dtype=loci.Comparison._DTYPE_LOCI),
                                ('chr1:1-25000', '-', 'Gypsy'): np.array([(24737, 24969,
                                                                           ['testA-2017-06-08.bam',
                                                                            'testB-2017-06-08.bam'],
                                                                           [12, 10],
                                                                           [0.545, 0.455])],
                                                                         dtype=loci.Comparison._DTYPE_LOCI),
                                ('chr1:1-25000', '+', 'Gypsy'): np.array([(2402, 2627,
                                                                           ['testA-2017-06-08.bam',
                                                                            'testB-2017-06-08.bam'],
                                                                           [6, 6],
                                                                           [0.5, 0.5]),
                                                                          (2791, 3188,
                                                                           ['testA-2017-06-08.bam',
                                                                            'testB-2017-06-08.bam'],
                                                                           [8, 8],
                                                                           [0.5, 0.5]),
                                                                          (24015, 24267,
                                                                           ['testA-2017-06-08.bam',
                                                                            'testB-2017-06-08.bam'],
                                                                           [4, 6],
                                                                           [0.4, 0.6])],
                                                                         dtype=loci.Comparison._DTYPE_LOCI),
                                ('chr1:1-25000', '+', 'PIF'): np.array([],
                                                                       dtype=loci.Comparison._DTYPE_LOCI)}),
                              ({'bams': [DATA_PATH + 'testA-2017-06-08.bam',
                                         DATA_PATH + 'testB-2017-06-08.bam'],
                                'categories': ['Gypsy', 'PIF'],
                                'references': ['chr1'],
                                'quality': 60,
                                'transposon_tag': 'ME',
                                'min_reads': 5,
                                'eps': 250,
                                'min_eps': 0,
                                'hierarchical': True,
                                'fingerprint_buffer': 50,
                                'bin_buffer': 0},
                               {('chr1:1-25000', '-', 'PIF'): np.array([(21784, 22032,
                                                                         ['testA-2017-06-08.bam',
                                                                          'testB-2017-06-08.bam'],
                                                                         [4, 6],
                                                                         [0.4, 0.6])],
                                                                       dtype=loci.Comparison._DTYPE_LOCI),
                                ('chr1:1-25000', '-', 'Gypsy'): np.array([(24737, 24944,
                                                                           ['testA-2017-06-08.bam',
                                                                            'testB-2017-06-08.bam'],
                                                                           [6, 0],
                                                                           [1.0, 0.0])],
                                                                         dtype=loci.Comparison._DTYPE_LOCI),
                                ('chr1:1-25000', '+', 'Gypsy'): np.array([(2402, 2627,
                                                                           ['testA-2017-06-08.bam',
                                                                            'testB-2017-06-08.bam'],
                                                                           [6, 6], [0.5, 0.5]),
                                                                          (2791, 3188,
                                                                           ['testA-2017-06-08.bam',
                                                                            'testB-2017-06-08.bam'],
                                                                           [8, 8],
                                                                           [0.5, 0.5]),
                                                                          (24015, 24267,
                                                                           ['testA-2017-06-08.bam',
                                                                            'testB-2017-06-08.bam'],
                                                                           [4, 6],
                                                                           [0.4, 0.6])],
                                                                         dtype=loci.Comparison._DTYPE_LOCI),
                                ('chr1:1-25000', '+', 'PIF'): np.array([],
                                                                       dtype=loci.Comparison._DTYPE_LOCI)})
                              ])
    def test_comparison(self, parameters, answer):
        parameters['cores'] = self.CORES
        query = batch.comparison(**parameters)
        answer = loci.Comparison.from_dict(answer)

        assert set(query.keys()) == set(answer.keys())
        for key in query.keys():
            npt.assert_array_equal(query[key], answer[key])


class TestComparisonMulticore(TestComparison):
    """Multi-core versions of integration tests for batch fingerprinting"""
    CORES = 2
