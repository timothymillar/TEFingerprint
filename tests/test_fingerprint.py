#! /usr/bin/env python

import numpy as np
from tefingerprint import loci
from tefingerprint import fingerprint


def test_create_contig_ids():
    dtype_loci_query = np.dtype([('start', np.int64),
                                 ('stop', np.int64),
                                 ('element', 'O')])

    query = loci.ContigSet(loci.Contig(loci.Header(reference='chr1',
                                                   strand='+',
                                                   category='gypsy'),
                                       np.array([(1, 5, 'gypsy1'),
                                                 (7, 9, 'gypsy4')],
                                                dtype=dtype_loci_query)),
                           loci.Contig(loci.Header(reference='chr1',
                                                   strand='-',
                                                   category='gypsy'),
                                       np.array([(3, 8, 'gypsy7'),
                                                 (9, 12, 'gypsy1')],
                                                dtype=dtype_loci_query)))

    dtype_loci_answer = np.dtype([('start', np.int64),
                                  ('stop', np.int64),
                                  ('element', 'O'),
                                  ('ID', 'O')])

    answer = loci.ContigSet(loci.Contig(loci.Header(reference='chr1',
                                                    strand='+',
                                                    category='gypsy'),
                                        np.array([(1, 5, 'gypsy1', 'gypsy_chr1_+_5'),
                                                    (7, 9, 'gypsy4', 'gypsy_chr1_+_9')],
                                                 dtype=dtype_loci_answer)),
                            loci.Contig(loci.Header(reference='chr1',
                                                    strand='-',
                                                    category='gypsy'),
                                        np.array([(3, 8, 'gypsy7', 'gypsy_chr1_-_3'),
                                                    (9, 12, 'gypsy1', 'gypsy_chr1_-_9')],
                                                 dtype=dtype_loci_answer)))

    assert query.map(fingerprint.create_contig_ids) == answer


def test_count_reads_n0():
    dtype_loci_reads = np.dtype([('tip', np.int64),
                                 ('element', 'O')])

    reads = loci.ContigSet(loci.Contig(loci.Header(reference='chr1',
                                                   strand='+',
                                                   category='gypsy',
                                                   source='bam1'),
                                       np.array([(2, 'gypsy1'),
                                                 (4, 'gypsy1'),
                                                 (5, 'gypsy4'),
                                                 (7, 'gypsy4'),
                                                 (7, 'gypsy7'),
                                                 (7, 'gypsy1'),
                                                 (8, 'gypsy1'),
                                                 (8, 'gypsy1')],
                                                dtype=dtype_loci_reads)),
                           loci.Contig(loci.Header(reference='chr1',
                                                   strand='+',
                                                   category='gypsy',
                                                   source='bam2'),
                                       np.array([(3, 'gypsy1'),
                                                 (4, 'gypsy1'),
                                                 (6, 'gypsy1'),
                                                 (7, 'gypsy1'),
                                                 (7, 'gypsy1'),
                                                 (7, 'gypsy1'),
                                                 (7, 'gypsy1'),
                                                 (50, 'gypsy7')],
                                                dtype=dtype_loci_reads)))

    dtype_loci_query = np.dtype([('start', np.int64),
                                 ('stop', np.int64)])

    query = loci.ContigSet(loci.Contig(loci.Header(reference='chr1',
                                                   strand='+',
                                                   category='gypsy'),
                                       np.array([(1, 15),
                                                 (30, 60)],
                                                dtype=dtype_loci_query)))

    dtype_loci_answer = np.dtype([('start', np.int64),
                                  ('stop', np.int64),
                                  ('median', np.int64),
                                  ('sample', [('0', [('name', 'O'),
                                                     ('count', np.int64)]),
                                              ('1', [('name', 'O'),
                                                     ('count', np.int64)])])])

    answer = loci.ContigSet(loci.Contig(loci.Header(reference='chr1',
                                                    strand='+',
                                                    category='gypsy'),
                                        np.array([(2,
                                                   8,
                                                   7,
                                                   (('bam1', 8),
                                                    ('bam2', 7))),
                                                  (50,
                                                   50,
                                                   50,
                                                   (('bam1', 0),
                                                    ('bam2', 1)))],
                                                 dtype=dtype_loci_answer)))

    assert fingerprint.count_reads(query,
                                   reads,
                                   trim=True,
                                   n_common_elements=0) == answer


def test_count_reads_n2():
    dtype_loci_reads = np.dtype([('tip', np.int64),
                                 ('element', 'O')])

    reads = loci.ContigSet(loci.Contig(loci.Header(reference='chr1',
                                                   strand='+',
                                                   category='gypsy',
                                                   source='bam1'),
                                       np.array([(2, 'gypsy1'),
                                                 (4, 'gypsy1'),
                                                 (5, 'gypsy4'),
                                                 (7, 'gypsy4'),
                                                 (7, 'gypsy7'),
                                                 (7, 'gypsy1'),
                                                 (8, 'gypsy1'),
                                                 (8, 'gypsy1')],
                                                dtype=dtype_loci_reads)),
                           loci.Contig(loci.Header(reference='chr1',
                                                   strand='+',
                                                   category='gypsy',
                                                   source='bam2'),
                                       np.array([(3, 'gypsy1'),
                                                 (4, 'gypsy1'),
                                                 (6, 'gypsy1'),
                                                 (7, 'gypsy1'),
                                                 (7, 'gypsy1'),
                                                 (7, 'gypsy1'),
                                                 (7, 'gypsy1'),
                                                 (50, 'gypsy7')],
                                                dtype=dtype_loci_reads)))

    dtype_loci_query = np.dtype([('start', np.int64),
                                 ('stop', np.int64)])

    query = loci.ContigSet(loci.Contig(loci.Header(reference='chr1',
                                                   strand='+',
                                                   category='gypsy'),
                                       np.array([(1, 15),
                                                 (30, 60)],
                                                dtype=dtype_loci_query)))

    dtype_loci_answer = np.dtype([('start', np.int64),
                                  ('stop', np.int64),
                                  ('median', np.int64),
                                  ('sample', [('0', [('name', 'O'),
                                                     ('count', np.int64),
                                                     ('element', [('0', [('name', 'O'),
                                                                         ('count', np.int64)]),
                                                                  ('1', [('name', 'O'),
                                                                         ('count', np.int64)])])]),
                                              ('1', [('name', 'O'),
                                                     ('count', np.int64),
                                                     ('element', [('0', [('name', 'O'),
                                                                         ('count', np.int64)]),
                                                                  ('1', [('name', 'O'),
                                                                         ('count', np.int64)])])])])])

    answer = loci.ContigSet(loci.Contig(loci.Header(reference='chr1',
                                                    strand='+',
                                                    category='gypsy'),
                                        np.array([(2,
                                                   8,
                                                   7,
                                                   (('bam1', 8, (('gypsy1', 5),
                                                                 ('gypsy4', 2))),
                                                    ('bam2', 7, (('gypsy1', 7),
                                                                 ('.', 0))))),
                                                  (50,
                                                   50,
                                                   50,
                                                   (('bam1', 0, (('.', 0),
                                                                 ('.', 0))),
                                                    ('bam2', 1, (('gypsy7', 1),
                                                                 ('.', 0)))))],
                                                 dtype=dtype_loci_answer)))

    assert fingerprint.count_reads(query,
                                   reads,
                                   trim=True,
                                   n_common_elements=2) == answer


def test_match_known_insertions():
    pass


def test_pair_clusters():
    pass


if __name__ == "__main__":
    pass
