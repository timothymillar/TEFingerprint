#! /usr/bin/env python

import os
import numpy as np
from tefingerprint import loci
from tefingerprint import bamio

DATA_PATH = os.path.dirname(os.path.realpath(__file__)) + '/data/'


def test_extract_informative_read_tips():
    """
    Test extraction of informative reads.
    Not all families of reads extracted.
    Family with no reads ('NOT-A-FAMILY') extracted.
    """
    bam = DATA_PATH + 'testA-2017-06-08.bam'

    query = bamio.extract_informative_read_tips(bam,
                                                'chr1',
                                                ['Gypsy',
                                                 'PIF-Harbinger',
                                                 'NOT-A-FAMILY'],
                                                quality=0,
                                                tag='ME')

    dtype_loci = np.dtype([('tip', np.int64),
                           ('element', 'O')])

    answer = loci.ContigSet(loci.Contig(loci.Header(reference='chr1',
                                                    strand='+',
                                                    category='Gypsy',
                                                    source='testA-2017-06-08.bam'),
                                        np.array([( 2452, 'Gypsy_Gypsy26_chr15_18793972'),
                                                  ( 2506, 'Gypsy_Gypsy26_chr15_18793972'),
                                                  ( 2553, 'Gypsy_Gypsy26_chr15_18793972'),
                                                  ( 2566, 'Gypsy_Gypsy26_chr15_18793972'),
                                                  ( 2577, 'Gypsy_Gypsy26_chr15_18793972'),
                                                  ( 2577, 'Gypsy_Gypsy26_chr15_18793972'),
                                                  ( 2841, 'Gypsy_Gypsy26_chr15_18793972'),
                                                  ( 2841, 'Gypsy_Gypsy26_chr15_18793972'),
                                                  ( 2841, 'Gypsy_Gypsy26_chr8_2502854'),
                                                  ( 2973, 'Gypsy_Gypsy26_chr18_27801424'),
                                                  ( 3024, 'Gypsy_Gypsy26_chr8_5114633'),
                                                  ( 3062, 'Gypsy_Gypsy26_chr8_5114633'),
                                                  ( 3039, 'Gypsy_Gypsy26_chr2_1987286'),
                                                  ( 3138, 'Gypsy_Gypsy26_chr18_27801424'),
                                                  (24065, 'Gypsy_Gypsy12_chr1_12715223'),
                                                  (24184, 'Gypsy_Gypsy7_chr4_10302390'),
                                                  (24195, 'Gypsy_Gypsy12_chr1_12715223'),
                                                  (24217, 'Gypsy_Gypsy12_chr1_12715223')],
                                                 dtype=dtype_loci)),
                            loci.Contig(loci.Header(reference='chr1',
                                                    strand='-',
                                                    category='Gypsy',
                                                    source='testA-2017-06-08.bam'),
                                        np.array([( 3217, 'Gypsy_Gypsy26_chr15_18793972'),
                                                  ( 3226, 'Gypsy_Gypsy26_chr15_18793972'),
                                                  ( 3246, 'Gypsy_Gypsy26_chr15_18793972'),
                                                  ( 3405, 'Gypsy_Gypsy26_chr2_1987286'),
                                                  ( 3646, 'Gypsy_Gypsy26_chr15_18793972'),
                                                  ( 3776, 'Gypsy_Gypsy26_chr18_27801424'),
                                                  ( 3779, 'Gypsy_Gypsy26_chr8_5114633'),
                                                  ( 3800, 'Gypsy_Gypsy26_chr8_5114633'),
                                                  (24787, 'Gypsy_Gypsy7_chr4_10302390'),
                                                  (24799, 'Gypsy_Gypsy29_chr11_13193899'),
                                                  (24850, 'Gypsy_Gypsy7_chr4_10302390'),
                                                  (24854, 'Gypsy_Gypsy12_chr1_12715223'),
                                                  (24857, 'Gypsy_Gypsy23_chr15_8310356'),
                                                  (24860, 'Gypsy_Gypsy23_chrUn_38723460'),
                                                  (24872, 'Gypsy_Gypsy23_chrUn_38723460'),
                                                  (24877, 'Gypsy_GYVIT1_chr6_13115950'),
                                                  (24894, 'Gypsy_Gypsy23_chrUn_38723460'),
                                                  (24895, 'Gypsy_Gypsy12_chr1_12715223'),
                                                  (24910, 'Gypsy_Gypsy23_chr14_11656393'),
                                                  (24919, 'Gypsy_Gypsy23_chrUn_38723460')],
                                                 dtype=dtype_loci)),
                            loci.Contig(loci.Header(reference='chr1',
                                                    strand='+',
                                                    category='PIF-Harbinger',
                                                    source='testA-2017-06-08.bam'),
                                        np.array([(21282, 'PIF-Harbinger_Harbinger-3N3_chr16_20723579'),
                                                  (21308, 'PIF-Harbinger_Harbinger-3_chr2_4407914'),
                                                  (21435, 'PIF-Harbinger_Harbinger-3N3_chr16_20723579'),
                                                  (21448, 'PIF-Harbinger_Harbinger-3N3_chr16_20723579')],
                                                 dtype=dtype_loci)),
                            loci.Contig(loci.Header(reference='chr1',
                                                    strand='-',
                                                    category='PIF-Harbinger',
                                                    source='testA-2017-06-08.bam'),
                                        np.array([(21834, 'PIF-Harbinger_Harbinger-3N3_chr16_20723579'),
                                                  (21945, 'PIF-Harbinger_Harbinger-3N3_chr16_20723579'),
                                                  (21968, 'PIF-Harbinger_Harbinger-3N3_chr16_20723579'),
                                                  (21982, 'PIF-Harbinger_Harbinger-3N3_chr16_20723579')],
                                                 dtype=dtype_loci)),
                            loci.Contig(loci.Header(reference='chr1',
                                                    strand='+',
                                                    category='NOT-A-FAMILY',
                                                    source='testA-2017-06-08.bam'),
                                        np.array([],
                                                 dtype=dtype_loci)),
                            loci.Contig(loci.Header(reference='chr1',
                                                    strand='-',
                                                    category='NOT-A-FAMILY',
                                                    source='testA-2017-06-08.bam'),
                                        np.array([],
                                                 dtype=dtype_loci)))

    assert query == answer


def test_extract_gff_intervals():
    gff = DATA_PATH + 'testAnnotation-2017-11-27.gff'

    query = bamio.extract_gff_intervals(gff, 'chr1', ['Gypsy', 'Copia'])

    dtype_loci = np.dtype([('start', np.int64),
                           ('stop', np.int64),
                           ('element', '<O')])

    answer = loci.ContigSet(loci.Contig(loci.Header(reference='chr1',
                                                    category='Gypsy',
                                                    source='testAnnotation-2017-11-27.gff'),
                                        np.array([(3150, 3200, 'Gypsy-21_ClassI;chr1:3150-3200'),
                                                  (24250, 24700, 'Gypsy-21_ClassI;chr1:24250-24700')],
                                                 dtype=dtype_loci)),
                            loci.Contig(loci.Header(reference='chr1',
                                                    category='Copia',
                                                    source='testAnnotation-2017-11-27.gff'),
                                        np.array([(98260, 98322, 'Copia-10_ClassI;chr1:98260-98322')],
                                                 dtype=dtype_loci)))

    assert query == answer


def test_extract_anchor_intervals():
    pass  # TODO: this test
