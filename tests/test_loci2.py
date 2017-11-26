#! /usr/bin/env python

import pytest
import numpy as np
import numpy.testing as npt
from tefingerprint import loci2


class TestHeader:
    """Tests for class Header"""

    @pytest.mark.parametrize("header, tuple_",
                             [(loci2.Header(reference='chr1',
                                            strand='-',
                                            category='family',
                                            source='bam'),
                               ('chr1', '-', 'family', 'bam')),
                              (loci2.Header(reference='chr1',
                                            strand='-'),
                               ('chr1', '-'))])
    def test_tuple(self, header, tuple_):
        assert header.tuple == tuple_

    @pytest.mark.parametrize("header, dtype",
                             [(loci2.Header(reference='chr1',
                                            strand='-',
                                            category='family',
                                            source='bam'),
                               (np.dtype([('reference', 'O'),
                                          ('strand', '<U1'),
                                          ('category', 'O'),
                                          ('source', 'O')]))),
                              (loci2.Header(reference='chr1',
                                            strand='-'),
                               (np.dtype([('reference', 'O'),
                                          ('strand', '<U1')])))])
    def test_dtype(self, header, dtype):
        assert header.dtype == dtype

    def test_mutate(self):
        origional = loci2.Header(reference='chr1',
                                 strand='-',
                                 category=None,
                                 source='bam')

        answer = loci2.Header(reference='chr1',
                              strand='.',
                              category='family',
                              source=None)

        assert origional.mutate(strand='.',
                                category='family',
                                source=None) == answer


def test_sort():
    loci_dtype = np.dtype([('tip', np.int64), ('element', 'O')])

    header = loci2.Header(reference='chr1',
                          strand='-',
                          category='family')

    query_loci = np.array([(5, 'element1'),
                           (1, 'element2'),
                           (7, 'element3'),
                           (9, 'element4'),
                           (6, 'element5'),
                           (2, 'element6')], dtype=loci_dtype)
    query = loci2.Contig(header, query_loci)

    answer_loci = np.array([(1, 'element2'),
                            (2, 'element6'),
                            (5, 'element1'),
                            (6, 'element5'),
                            (7, 'element3'),
                            (9, 'element4')], dtype=loci_dtype)
    answer = loci2.Contig(header, answer_loci)

    assert loci2.sort(query, order='tip') == answer
    assert loci2.sort(query, order='tip') != query


def test_iter_values():
    loci_dtype = np.dtype([('tip', np.int64),
                           ('element', 'O')])
    loci = np.array([(5, 'element1'),
                     (1, 'element2'),
                     (7, 'element3')], dtype=loci_dtype)

    query_header = loci2.Header(reference='chr1',
                                strand='-',
                                category=None,
                                source='bam')
    query = loci2.Contig(query_header, loci)

    answer = [('chr1', '-', 'bam', 5, 'element1'),
              ('chr1', '-', 'bam', 1, 'element2'),
              ('chr1', '-', 'bam', 7, 'element3')]

    assert list(loci2.iter_values(query)) == answer


def test_as_array():
    """Test conversion of nested loci data to flat array"""
    dtype_element_count = np.dtype([('name', 'O'),
                                    ('count', np.int64)])
    dtype_elements = np.dtype([(str(i), dtype_element_count)
                               for i in range(2)])
    dtype_sample_count = np.dtype([('name', 'O'),
                                   ('count', np.int64),
                                   ('element', dtype_elements)])
    dtype_samples = np.dtype([(str(i), dtype_sample_count)
                              for i in range(2)])
    dtype_loci = np.dtype([('start', np.int64),
                           ('stop', np.int64),
                           ('sample', dtype_samples)])

    # 3 element array with nested structured data
    loci = np.array([(10, 15, (('bam1', 9, (('gypsy7', 5), ('gypsy3', 3))),
                               ('bam2', 8, (('gypsy7', 7), ('gypsy3', 1))))),
                     (21, 32, (('bam1', 7, (('gypsy3', 5), ('gypsy1', 2))),
                               ('bam2', 7, (('gypsy3', 7), (None, 0))))),
                     (43, 61, (('bam1', 5, (('gypsy9', 3), ('gypsy3', 2))),
                               ('bam2', 6, (('gypsy3', 5), ('gypsy9', 1)))))
                     ],
                    dtype=dtype_loci)

    header = loci2.Header(reference='chr1',
                          strand='-',
                          category='gypsy')

    query = loci2.Contig(header, loci)

    # array dtype includes header fields
    answer_dtype = np.dtype([('reference', 'O'),
                             ('strand', '<U1'),
                             ('category', 'O'),
                             ('start', np.int64),
                             ('stop', np.int64),
                             ('sample_0_name', 'O'),
                             ('sample_0_count', np.int64),
                             ('sample_0_element_0_name', 'O'),
                             ('sample_0_element_0_count', np.int64),
                             ('sample_0_element_1_name', 'O'),
                             ('sample_0_element_1_count', np.int64),
                             ('sample_1_name', 'O'),
                             ('sample_1_count', np.int64),
                             ('sample_1_element_0_name', 'O'),
                             ('sample_1_element_0_count', np.int64),
                             ('sample_1_element_1_name', 'O'),
                             ('sample_1_element_1_count', np.int64)])

    # 3 element array with flat structured data
    answer = np.array([('chr1', '-', 'gypsy', 10, 15, 'bam1', 9, 'gypsy7', 5,
                        'gypsy3', 3, 'bam2', 8, 'gypsy7', 7, 'gypsy3', 1),
                       ('chr1', '-', 'gypsy', 21, 32, 'bam1', 7, 'gypsy3', 5,
                        'gypsy1', 2, 'bam2', 7, 'gypsy3', 7, None, 0),
                       ('chr1', '-', 'gypsy', 43, 61, 'bam1', 5, 'gypsy9', 3,
                        'gypsy3', 2, 'bam2', 6, 'gypsy3', 5, 'gypsy9', 1)],
                      dtype=answer_dtype)

    npt.assert_array_equal(loci2.as_array(query), answer)


def test_mutate_header():
    loci_dtype = np.dtype([('tip', np.int64), ('element', 'O')])
    loci = np.array([(5, 'element1'),
                     (1, 'element2'),
                     (7, 'element3')], dtype=loci_dtype)

    query_header = loci2.Header(reference='chr1',
                                strand='-',
                                category=None,
                                source='bam')
    query = loci2.Contig(query_header, loci)

    answer_header = loci2.Header(reference='chr1',
                                 strand='.',
                                 category='family',
                                 source=None)
    answer = loci2.Contig(answer_header, loci)

    assert loci2.mutate_header(query,
                               strand='.',
                               category='family',
                               source=None) == answer
    assert loci2.mutate_header(query,
                               strand='.',
                               category='family',
                               source=None) != query


def test_append():
    loci_dtype = np.dtype([('tip', np.int64), ('element', 'O')])

    query_1 = loci2.Contig(loci2.Header(reference='chr1',
                                        strand='-',
                                        category='family',
                                        source='bam'),
                           np.array([(5, 'element1'),
                                     (1, 'element2'),
                                     (7, 'element3')], dtype=loci_dtype))

    query_2 = loci2.Contig(loci2.Header(reference='chr1',
                                        strand='-',
                                        category='family',
                                        source='bam'),
                           np.array([(9, 'element4'),
                                     (6, 'element5'),
                                     (2, 'element6')], dtype=loci_dtype))

    answer = loci2.Contig(loci2.Header(reference='chr1',
                                       strand='-',
                                       category='family',
                                       source='bam'),
                          np.array([(5, 'element1'),
                                    (1, 'element2'),
                                    (7, 'element3'),
                                    (9, 'element4'),
                                    (6, 'element5'),
                                    (2, 'element6')], dtype=loci_dtype))

    assert loci2.append(query_1, query_2) == answer


def test_append_header_miss_match():
    loci_dtype = np.dtype([('tip', np.int64), ('element', 'O')])

    query_1 = loci2.Contig(loci2.Header(reference='chr1',
                                        strand='-',
                                        category='family',
                                        source='bam'),
                           np.array([(5, 'element1'),
                                     (1, 'element2'),
                                     (7, 'element3')], dtype=loci_dtype))

    query_2 = loci2.Contig(loci2.Header(reference='chr2', # miss-matched
                                        strand='-',
                                        category='family',
                                        source='bam'),
                           np.array([(9, 'element4'),
                                     (6, 'element5'),
                                     (2, 'element6')], dtype=loci_dtype))

    try:
        loci2.append(query_1, query_2)
    except ValueError:
        pass
    else:
        assert False


def test_drop_field():
    query = loci2.Contig(loci2.Header(reference='chr1',
                                      strand='-',
                                      category='family',
                                      source='bam'),
                         np.array([(5, 'element1'),
                                   (1, 'element2'),
                                   (7, 'element3')],
                                  dtype=np.dtype([('tip', np.int64),
                                                  ('element', 'O')])))

    answer = loci2.Contig(loci2.Header(reference='chr1',
                                       strand='-',
                                       category='family',
                                       source='bam'),
                          np.array([5, 1, 7]))
    answer.loci = np.array(answer.loci, np.dtype([('tip', np.int64)]))

    assert loci2.drop_field(query, 'element') == answer
    assert loci2.drop_field(query, 'element') != query


def test_add_field():
    query = loci2.Contig(loci2.Header(reference='chr1',
                                      strand='-',
                                      category='family',
                                      source='bam'),
                         np.array([(5, 'element1'),
                                   (1, 'element2'),
                                   (7, 'element3')],
                                  dtype=np.dtype([('tip', np.int64),
                                                  ('element', 'O')])))

    answer = loci2.Contig(loci2.Header(reference='chr1',
                                       strand='-',
                                       category='family',
                                       source='bam'),
                          np.array([(5, 'element1', 0),
                                    (1, 'element2', 0),
                                    (7, 'element3', 0)],
                                   dtype=np.dtype([('tip', np.int64),
                                                   ('element', 'O'),
                                                   ('field', np.int64)])))

    assert loci2.add_field(query, np.dtype([('field', np.int64)])) == answer
    assert loci2.add_field(query, np.dtype([('field', np.int64)])) != query


def test_cluster():
    header = loci2.Header(reference='chr1',
                          strand='-',
                          category='Gypsy',
                          source='bam')
    query_loci = np.array([(   0, 'Gypsy'),
                           (   0, 'Gypsy'),
                           (  60, 'Gypsy'),
                           (  61, 'Gypsy'),
                           (  61, 'Gypsy'),
                           (  61, 'Gypsy'),
                           (  76, 'Gypsy'),
                           (  78, 'Gypsy'),
                           ( 122, 'Gypsy'),
                           ( 122, 'Gypsy'),
                           ( 141, 'Gypsy'),
                           ( 183, 'Gypsy'),
                           ( 251, 'Gypsy'),
                           ( 260, 'Gypsy'),
                           ( 260, 'Gypsy'),
                           ( 263, 'Gypsy'),
                           ( 263, 'Gypsy'),
                           ( 267, 'Gypsy'),
                           ( 267, 'Gypsy'),
                           ( 288, 'Gypsy'),
                           ( 288, 'Gypsy'),
                           ( 295, 'Gypsy'),
                           ( 300, 'Gypsy'),
                           ( 310, 'Gypsy'),
                           ( 310, 'Gypsy'),
                           ( 317, 'Gypsy'),
                           ( 317, 'Gypsy'),
                           ( 334, 'Gypsy'),
                           ( 334, 'Gypsy'),
                           ( 335, 'Gypsy'),
                           ( 338, 'Gypsy'),
                           ( 338, 'Gypsy'),
                           ( 338, 'Gypsy'),
                           ( 338, 'Gypsy'),
                           ( 340, 'Gypsy'),
                           ( 342, 'Gypsy'),
                           ( 342, 'Gypsy'),
                           ( 344, 'Gypsy'),
                           ( 344, 'Gypsy'),
                           ( 358, 'Gypsy'),
                           ( 367, 'Gypsy'),
                           ( 370, 'Gypsy'),
                           ( 370, 'Gypsy'),
                           ( 377, 'Gypsy'),
                           ( 387, 'Gypsy'),
                           ( 402, 'Gypsy'),
                           ( 403, 'Gypsy'),
                           ( 410, 'Gypsy'),
                           ( 410, 'Gypsy'),
                           ( 410, 'Gypsy'),
                           ( 418, 'Gypsy'),
                           ( 418, 'Gypsy'),
                           ( 424, 'Gypsy'),
                           ( 424, 'Gypsy'),
                           ( 577, 'Gypsy'),
                           ( 857, 'Gypsy'),
                           ( 879, 'Gypsy'),
                           ( 921, 'Gypsy'),
                           ( 921, 'Gypsy'),
                           (1007, 'Gypsy'),
                           (1031, 'Gypsy'),
                           (1051, 'Gypsy'),
                           (1051, 'Gypsy'),
                           (1059, 'Gypsy'),
                           (1071, 'Gypsy'),
                           (1071, 'Gypsy'),
                           (1080, 'Gypsy'),
                           (1094, 'Gypsy'),
                           (1094, 'Gypsy'),
                           (1110, 'Gypsy'),
                           (1110, 'Gypsy'),
                           (1113, 'Gypsy'),
                           (1113, 'Gypsy'),
                           (1183, 'Gypsy'),
                           (1189, 'Gypsy'),
                           (1200, 'Gypsy'),
                           (1200, 'Gypsy'),
                           (1217, 'Gypsy'),
                           (1234, 'Gypsy'),
                           (1234, 'Gypsy'),
                           (1591, 'Gypsy'),
                           (1620, 'Gypsy'),
                           (1620, 'Gypsy'),
                           (1662, 'Gypsy'),
                           (1686, 'Gypsy'),
                           (1707, 'Gypsy'),
                           (1755, 'Gypsy'),
                           (1828, 'Gypsy'),
                           (1828, 'Gypsy'),
                           (1848, 'Gypsy'),
                           (1848, 'Gypsy'),
                           (1848, 'Gypsy'),
                           (1848, 'Gypsy'),
                           (1851, 'Gypsy'),
                           (1851, 'Gypsy'),
                           (1852, 'Gypsy'),
                           (1917, 'Gypsy')],
                          dtype=np.dtype([('tip', np.int64),
                                          ('element', 'O')]))
    query = loci2.Contig(header, query_loci)

    answer_loci = np.array([(0, 577),
                            (879, 1234),
                            (1662, 1917)],
                           dtype=np.dtype([('start', np.int64),
                                           ('stop', np.int64)]))
    answer = loci2.Contig(header, answer_loci)

    assert loci2.cluster(query,
                         'tip',
                         10,
                         epsilon=200,
                         minimum_epsilon=10,
                         hierarchical=True) == answer


def test_cluster_empty():
    header = loci2.Header(reference='chr1',
                          strand='-',
                          category='Gypsy',
                          source='bam')
    query_loci = np.array([],
                          dtype=np.dtype([('tip', np.int64),
                                          ('element', 'O')]))
    query = loci2.Contig(header, query_loci)

    answer_loci = np.array([],
                           dtype=np.dtype([('start', np.int64),
                                           ('stop', np.int64)]))
    answer = loci2.Contig(header, answer_loci)

    assert loci2.cluster(query,
                         'tip',
                         10,
                         epsilon=200,
                         minimum_epsilon=10,
                         hierarchical=True) == answer


@pytest.mark.parametrize("query, answer",
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
                          ([(3, 6), (6, 8), (7, 9), (10, 12), (13, 13),
                            (15, 25), (16, 17), (19, 20)],
                           [(3, 9), (10, 12), (13, 13), (15, 25)])],
                         ids=['single', 'nested', 'adjacent', 'combined'])
def test_unions(query, answer):
    """
    Test includes following edge cases:
     * Long locus completely overlaps short loci:
        (15, 25) & (16, 17) & (19, 20) --> (15, 25)
     * Adjacent loci do not get merged:
        (7, 9) & (10, 12) -->  (*, 9) & (10, *)
     * Locus may span a single base:
        (13, 13) --> (13, 13)
    """
    header = loci2.Header(reference='chr1',
                          strand='-',
                          source='bam')

    dtype = np.dtype([('start', np.int64),
                      ('stop', np.int64)])

    query = loci2.Contig(header,
                         np.array(query, dtype=dtype))

    answer = loci2.Contig(header,
                          np.array(answer, dtype=dtype))

    assert loci2.unions(query) == answer


def test_unions_buffered():
    header = loci2.Header(reference='chr1',
                          strand='-',
                          source='bam')

    dtype = np.dtype([('start', np.int64),
                      ('stop', np.int64)])

    query = loci2.Contig(header,
                         np.array([(3, 6),
                                   (6, 8),
                                   (7, 9),
                                   (10, 12),
                                   (13, 13),
                                   (15, 25),
                                   (16, 17),
                                   (19, 20)], dtype=dtype))

    answer = loci2.Contig(header,
                          np.array([(-2, 9),
                                    (10, 12),
                                    (13, 14),
                                    (15, 30)], dtype=dtype))

    assert loci2.unions_buffered(query, 5) == answer


class TestContigSet:
    """Tests for class ContigSet"""

    def test_init_empty(self):
        """"""
        answer = loci2.ContigSet()
        assert type(answer) == loci2.ContigSet

    def test_init_different_headers(self):
        """"""
        dtype_loci = np.dtype([('tip', np.int64), ('element', 'O')])

        header_1 = loci2.Header(reference='chr1',
                                strand='+',
                                category='gypsy')
        contig_1=loci2.Contig(header_1,
                              np.array([(1, 'gypsy1'),
                                        (7, 'gypsy4')],
                                       dtype=dtype_loci))

        header_2 = loci2.Header(reference='chr2',
                                strand='+',
                                category='gypsy')
        contig_2=loci2.Contig(header_2,
                              np.array([(3, 'gypsy7'),
                                        (9, 'gypsy1')],
                                       dtype=dtype_loci))

        query = loci2.ContigSet(contig_1, contig_2)

        assert len(query) == 4
        assert len(list(query.contigs())) == 2
        assert len(query.headers()) == 2

    def test_init_clashing_headers(self):
        """Contigs with same header should cause ValueError"""
        dtype_loci = np.dtype([('tip', np.int64), ('element', 'O')])

        header = loci2.Header(reference='chr1',
                              strand='+',
                              category='gypsy')

        contig_1 = loci2.Contig(header,
                                np.array([(1, 'gypsy1'),
                                          (7, 'gypsy4')],
                                         dtype=dtype_loci))

        contig_2 = loci2.Contig(header,
                                np.array([(3, 'gypsy7'),
                                          (9, 'gypsy1')],
                                         dtype=dtype_loci))

        try:
            loci2.ContigSet(contig_1, contig_2)
        except ValueError:
            assert True
        else:
            assert False

    def test_init_append_headers(self):
        """Contigs with same header should be appended"""
        dtype_loci = np.dtype([('tip', np.int64), ('element', 'O')])

        header = loci2.Header(reference='chr1',
                              strand='+',
                              category='gypsy')

        contig_1 = loci2.Contig(header,
                                np.array([(1, 'gypsy1'),
                                          (7, 'gypsy4')],
                                         dtype=dtype_loci))

        contig_2 = loci2.Contig(header,
                                np.array([(3, 'gypsy7'),
                                          (9, 'gypsy1')],
                                         dtype=dtype_loci))

        query = loci2.ContigSet(contig_1,
                                contig_2,
                                append_duplicate_headers=True)

        assert len(query) == 4
        assert len(list(query.contigs())) == 1
        assert len(query.headers()) == 1

        query_loci = list(query.contigs())[0].loci

        answer_loci = np.array([(1, 'gypsy1'),
                                (7, 'gypsy4'),
                                (3, 'gypsy7'),
                                (9, 'gypsy1')],
                               dtype=dtype_loci)

        npt.assert_array_equal(query_loci, answer_loci)

    def test_dtype_headers(self):
        dtype_loci = np.dtype([('tip', np.int64), ('element', 'O')])

        header_1 = loci2.Header(reference='chr1',
                                strand='+',
                                category='gypsy')
        contig_1=loci2.Contig(header_1,
                              np.array([(1, 'gypsy1'),
                                        (7, 'gypsy4')],
                                       dtype=dtype_loci))

        header_2 = loci2.Header(reference='chr2',
                                strand='+',
                                category='gypsy')
        contig_2=loci2.Contig(header_2,
                              np.array([(3, 'gypsy7'),
                                        (9, 'gypsy1')],
                                       dtype=dtype_loci))

        query = loci2.ContigSet(contig_1, contig_2)

        assert query.dtype_headers() == contig_1.header.dtype
        assert query.dtype_headers() == contig_1.header.dtype

    def test_dtype_loci(self):
        dtype_loci = np.dtype([('tip', np.int64), ('element', 'O')])

        header_1 = loci2.Header(reference='chr1',
                                strand='+',
                                category='gypsy')
        contig_1 = loci2.Contig(header_1,
                                np.array([(1, 'gypsy1'),
                                          (7, 'gypsy4')],
                                         dtype=dtype_loci))

        header_2 = loci2.Header(reference='chr2',
                                strand='+',
                                category='gypsy')
        contig_2 = loci2.Contig(header_2,
                                np.array([(3, 'gypsy7'),
                                          (9, 'gypsy1')],
                                         dtype=dtype_loci))

        query = loci2.ContigSet(contig_1, contig_2)

        assert query.dtype_loci() == contig_1.loci.dtype
        assert query.dtype_loci() == contig_1.loci.dtype

    def test_headers(self):
        dtype_loci = np.dtype([('tip', np.int64), ('element', 'O')])

        header_1 = loci2.Header(reference='chr1',
                                strand='+',
                                category='gypsy')
        contig_1=loci2.Contig(header_1,
                              np.array([(1, 'gypsy1'),
                                        (7, 'gypsy4')],
                                       dtype=dtype_loci))

        header_2 = loci2.Header(reference='chr2',
                                strand='+',
                                category='gypsy')
        contig_2=loci2.Contig(header_2,
                              np.array([(3, 'gypsy7'),
                                        (9, 'gypsy1')],
                                       dtype=dtype_loci))

        query = set(loci2.ContigSet(contig_1, contig_2).headers())
        answer = {header_1, header_2}

        assert query == answer

    def test_add_different_headers(self):
        """"""
        dtype_loci = np.dtype([('tip', np.int64), ('element', 'O')])

        header_1 = loci2.Header(reference='chr1',
                                strand='+',
                                category='gypsy')
        contig_1=loci2.Contig(header_1,
                              np.array([(1, 'gypsy1'),
                                        (7, 'gypsy4')],
                                       dtype=dtype_loci))

        header_2 = loci2.Header(reference='chr2',
                                strand='+',
                                category='gypsy')
        contig_2=loci2.Contig(header_2,
                              np.array([(3, 'gypsy7'),
                                        (9, 'gypsy1')],
                                       dtype=dtype_loci))

        query = loci2.ContigSet(contig_1)
        query.add(contig_2)

        assert len(query) == 4
        assert len(list(query.contigs())) == 2
        assert len(query.headers()) == 2

    def test_add_clashing_headers(self):
        """Contigs with same header should cause ValueError"""
        dtype_loci = np.dtype([('tip', np.int64), ('element', 'O')])

        header = loci2.Header(reference='chr1',
                              strand='+',
                              category='gypsy')

        contig_1 = loci2.Contig(header,
                                np.array([(1, 'gypsy1'),
                                          (7, 'gypsy4')],
                                         dtype=dtype_loci))

        contig_2 = loci2.Contig(header,
                                np.array([(3, 'gypsy7'),
                                          (9, 'gypsy1')],
                                         dtype=dtype_loci))

        query = loci2.ContigSet(contig_1)
        try:
            query.add(contig_2)
        except ValueError:
            assert True
        else:
            assert False

    def test_add_append_headers(self):
        """Contigs with same header should be appended"""
        dtype_loci = np.dtype([('tip', np.int64), ('element', 'O')])

        header = loci2.Header(reference='chr1',
                              strand='+',
                              category='gypsy')

        contig_1 = loci2.Contig(header,
                                np.array([(1, 'gypsy1'),
                                          (7, 'gypsy4')],
                                         dtype=dtype_loci))

        contig_2 = loci2.Contig(header,
                                np.array([(3, 'gypsy7'),
                                          (9, 'gypsy1')],
                                         dtype=dtype_loci))

        query = loci2.ContigSet(contig_1)
        query.add(contig_2, append_duplicate_headers=True)

        assert len(query) == 4
        assert len(list(query.contigs())) == 1
        assert len(query.headers()) == 1

        query_loci = list(query.contigs())[0].loci

        answer_loci = np.array([(1, 'gypsy1'),
                                (7, 'gypsy4'),
                                (3, 'gypsy7'),
                                (9, 'gypsy1')],
                               dtype=dtype_loci)

        npt.assert_array_equal(query_loci, answer_loci)

    def test_update(self):
        dtype_loci = np.dtype([('tip', np.int64), ('element', 'O')])

        header_1 = loci2.Header(reference='chr1',
                                strand='+',
                                category='gypsy')
        contig_1 = loci2.Contig(header_1,
                              np.array([(1, 'gypsy1'),
                                        (7, 'gypsy4')],
                                       dtype=dtype_loci))

        header_2 = loci2.Header(reference='chr2',
                                strand='+',
                                category='gypsy')
        contig_2 = loci2.Contig(header_2,
                              np.array([(3, 'gypsy7'),
                                        (9, 'gypsy1')],
                                       dtype=dtype_loci))

        query = loci2.ContigSet(contig_1)
        query_2 = loci2.ContigSet(contig_2)
        query.update(query_2.contigs())

        assert len(query) == 4
        assert len(list(query.contigs())) == 2
        assert len(query.headers()) == 2

    def test_map(self):
        dtype_loci = np.dtype([('tip', np.int64), ('element', 'O')])

        header_1 = loci2.Header(reference='chr1',
                                strand='+',
                                category='gypsy')
        contig_1 = loci2.Contig(header_1,
                              np.array([(1, 'gypsy1'),
                                        (7, 'gypsy4')],
                                       dtype=dtype_loci))
        contig_alt_1 = loci2.Contig(header_1,
                              np.array([(101, 'gypsy1'),
                                        (107, 'gypsy4')],
                                       dtype=dtype_loci))

        header_2 = loci2.Header(reference='chr2',
                                strand='+',
                                category='gypsy')
        contig_2 = loci2.Contig(header_2,
                              np.array([(3, 'gypsy7'),
                                        (9, 'gypsy1')],
                                       dtype=dtype_loci))
        contig_alt_2 = loci2.Contig(header_2,
                              np.array([(103, 'gypsy7'),
                                        (109, 'gypsy1')],
                                       dtype=dtype_loci))

        query = loci2.ContigSet(contig_1, contig_2)

        def func(contig):
            """dummy function that adds 100 to contig loci 'tip's"""
            loci = np.copy(contig.loci)
            loci['tip'] += 100
            return loci2.Contig(contig.header, loci)

        query = query.map(func)

        answer = loci2.ContigSet(contig_alt_1, contig_alt_2)

        assert query == answer

    def test_iter_values(self):
        dtype_loci = np.dtype([('tip', np.int64), ('element', 'O')])

        header_1 = loci2.Header(reference='chr1',
                                strand='+',
                                category='gypsy')
        contig_1 = loci2.Contig(header_1,
                              np.array([(1, 'gypsy1'),
                                        (7, 'gypsy4')],
                                       dtype=dtype_loci))

        header_2 = loci2.Header(reference='chr2',
                                strand='+',
                                category='gypsy')
        contig_2 = loci2.Contig(header_2,
                              np.array([(3, 'gypsy7'),
                                        (9, 'gypsy1')],
                                       dtype=dtype_loci))

        query = loci2.ContigSet(contig_1, contig_2)

        answer = [('chr1', '+', 'gypsy', 1, 'gypsy1'),
                  ('chr1', '+', 'gypsy', 7, 'gypsy4'),
                  ('chr2', '+', 'gypsy', 3, 'gypsy7'),
                  ('chr2', '+', 'gypsy', 9, 'gypsy1')]

        assert list(query.iter_values()) == answer

    def test_as_array(self):
        pass

    def test_as_tabular_lines(self):
        pass

    def test_as_gff_lines(self):
        pass

