#! /usr/bin/env python

import pytest
import gffutils
import tectoolkit.filter_gff

GFF_UNSORTED = """chr11\t.\t.\t126799\t127584\t.\t+\t.\tName=Copia;ID=bin_Copia_chr11_+_126799;read_count_max=98;read_count_min=57;cluster_absence=0;read_presence=3;cluster_presence=3;read_absence=0
chr11\t.\t.\t126915\t127484\t.\t+\t.\tsample=Regen_49_danglers.vitis.bwa_mem.bam;Name=Copia;ID=0_Copia_chr11_+_126915;Parent=bin_Copia_chr11_+_126799
chr11\t.\t.\t126969\t127466\t.\t+\t.\tsample=Regen_84_danglers.vitis.bwa_mem.bam;Name=Copia;ID=1_Copia_chr11_+_126969;Parent=bin_Copia_chr11_+_126799
chr11\t.\t.\t126899\t127476\t.\t+\t.\tsample=Regen_105_danglers.vitis.bwa_mem.bam;Name=Copia;ID=2_Copia_chr11_+_126899;Parent=bin_Copia_chr11_+_126799
chr11\t.\t.\t10275156\t10276194\t.\t-\t.\tName=Gypsy;ID=bin_Gypsy_chr11_-_10275156;read_count_max=31;read_count_min=22;cluster_absence=0;read_presence=3;cluster_presence=3;read_absence=0
chr11\t.\t.\t10275256\t10276045\t.\t-\t.\tsample=Regen_49_danglers.vitis.bwa_mem.bam;Name=Gypsy;ID=0_Gypsy_chr11_-_10275256;Parent=bin_Gypsy_chr11_-_10275156
chr11\t.\t.\t10276179\t10276939\t.\t-\t.\tsample=Regen_49_danglers.vitis.bwa_mem.bam;Name=Gypsy;ID=0_Gypsy_chr11_-_10276179;Parent=bin_Gypsy_chr11_-_10275156
chr11\t.\t.\t10275392\t10276024\t.\t-\t.\tsample=Regen_84_danglers.vitis.bwa_mem.bam;Name=Gypsy;ID=1_Gypsy_chr11_-_10275392;Parent=bin_Gypsy_chr11_-_10275156
chr11\t.\t.\t10275513\t10276094\t.\t-\t.\tsample=Regen_105_danglers.vitis.bwa_mem.bam;Name=Gypsy;ID=2_Gypsy_chr11_-_10275513;Parent=bin_Gypsy_chr11_-_10275156
chr11\t.\t.\t671652\t672285\t.\t-\t.\tName=hAT;ID=bin_hAT_chr11_-_671652;read_count_max=11;read_count_min=0;cluster_absence=2;read_presence=1;cluster_presence=1;read_absence=2
chr11\t.\t.\t671752\t672185\t.\t-\t.\tName=hAT;ID=0_hAT_chr11_-_671752;sample=Regen_49_danglers.vitis.bwa_mem.bam;Parent=bin_hAT_chr11_-_671652
chr11\t.\t.\t1177277\t1177878\t.\t-\t.\tName=hAT;ID=bin_hAT_chr11_-_1177277;read_count_max=45;read_count_min=0;cluster_absence=2;read_presence=1;cluster_presence=1;read_absence=2
chr11\t.\t.\t1177377\t1177778\t.\t-\t.\tName=hAT;ID=0_hAT_chr11_-_1177377;sample=Regen_49_danglers.vitis.bwa_mem.bam;Parent=bin_hAT_chr11_-_1177277
chr11\t.\t.\t1177894\t1178429\t.\t-\t.\tName=hAT;ID=bin_hAT_chr11_-_1177894;read_count_max=20;read_count_min=0;cluster_absence=2;read_presence=1;cluster_presence=1;read_absence=2
chr11\t.\t.\t1177994\t1178329\t.\t-\t.\tName=hAT;ID=0_hAT_chr11_-_1177994;sample=Regen_49_danglers.vitis.bwa_mem.bam;Parent=bin_hAT_chr11_-_1177894"""

GFF_SORTED = """chr11\t.\t.\t126799\t127584\t.\t+\t.\tName=Copia;ID=bin_Copia_chr11_+_126799;read_count_max=98;read_count_min=57;cluster_absence=0;read_presence=3;cluster_presence=3;read_absence=0
chr11\t.\t.\t126915\t127484\t.\t+\t.\tName=Copia;ID=0_Copia_chr11_+_126915;sample=Regen_49_danglers.vitis.bwa_mem.bam;Parent=bin_Copia_chr11_+_126799
chr11\t.\t.\t126969\t127466\t.\t+\t.\tName=Copia;ID=1_Copia_chr11_+_126969;sample=Regen_84_danglers.vitis.bwa_mem.bam;Parent=bin_Copia_chr11_+_126799
chr11\t.\t.\t126899\t127476\t.\t+\t.\tName=Copia;ID=2_Copia_chr11_+_126899;sample=Regen_105_danglers.vitis.bwa_mem.bam;Parent=bin_Copia_chr11_+_126799
chr11\t.\t.\t10275156\t10276194\t.\t-\t.\tName=Gypsy;ID=bin_Gypsy_chr11_-_10275156;read_count_max=31;read_count_min=22;cluster_absence=0;read_presence=3;cluster_presence=3;read_absence=0
chr11\t.\t.\t10275256\t10276045\t.\t-\t.\tName=Gypsy;ID=0_Gypsy_chr11_-_10275256;sample=Regen_49_danglers.vitis.bwa_mem.bam;Parent=bin_Gypsy_chr11_-_10275156
chr11\t.\t.\t10276179\t10276939\t.\t-\t.\tName=Gypsy;ID=0_Gypsy_chr11_-_10276179;sample=Regen_49_danglers.vitis.bwa_mem.bam;Parent=bin_Gypsy_chr11_-_10275156
chr11\t.\t.\t10275392\t10276024\t.\t-\t.\tName=Gypsy;ID=1_Gypsy_chr11_-_10275392;sample=Regen_84_danglers.vitis.bwa_mem.bam;Parent=bin_Gypsy_chr11_-_10275156
chr11\t.\t.\t10275513\t10276094\t.\t-\t.\tName=Gypsy;ID=2_Gypsy_chr11_-_10275513;sample=Regen_105_danglers.vitis.bwa_mem.bam;Parent=bin_Gypsy_chr11_-_10275156
chr11\t.\t.\t671652\t672285\t.\t-\t.\tName=hAT;ID=bin_hAT_chr11_-_671652;read_count_max=11;read_count_min=0;cluster_absence=2;read_presence=1;cluster_presence=1;read_absence=2
chr11\t.\t.\t671752\t672185\t.\t-\t.\tName=hAT;ID=0_hAT_chr11_-_671752;sample=Regen_49_danglers.vitis.bwa_mem.bam;Parent=bin_hAT_chr11_-_671652
chr11\t.\t.\t1177277\t1177878\t.\t-\t.\tName=hAT;ID=bin_hAT_chr11_-_1177277;read_count_max=45;read_count_min=0;cluster_absence=2;read_presence=1;cluster_presence=1;read_absence=2
chr11\t.\t.\t1177377\t1177778\t.\t-\t.\tName=hAT;ID=0_hAT_chr11_-_1177377;sample=Regen_49_danglers.vitis.bwa_mem.bam;Parent=bin_hAT_chr11_-_1177277
chr11\t.\t.\t1177894\t1178429\t.\t-\t.\tName=hAT;ID=bin_hAT_chr11_-_1177894;read_count_max=20;read_count_min=0;cluster_absence=2;read_presence=1;cluster_presence=1;read_absence=2
chr11\t.\t.\t1177994\t1178329\t.\t-\t.\tName=hAT;ID=0_hAT_chr11_-_1177994;sample=Regen_49_danglers.vitis.bwa_mem.bam;Parent=bin_hAT_chr11_-_1177894"""

GFF_SINGLE_CLUSTERS_WITH_20_READS = """chr11\t.\t.\t1177277\t1177878\t.\t-\t.\tName=hAT;ID=bin_hAT_chr11_-_1177277;read_count_max=45;read_count_min=0;cluster_absence=2;read_presence=1;cluster_presence=1;read_absence=2
chr11\t.\t.\t1177377\t1177778\t.\t-\t.\tName=hAT;ID=0_hAT_chr11_-_1177377;sample=Regen_49_danglers.vitis.bwa_mem.bam;Parent=bin_hAT_chr11_-_1177277
chr11\t.\t.\t1177894\t1178429\t.\t-\t.\tName=hAT;ID=bin_hAT_chr11_-_1177894;read_count_max=20;read_count_min=0;cluster_absence=2;read_presence=1;cluster_presence=1;read_absence=2
chr11\t.\t.\t1177994\t1178329\t.\t-\t.\tName=hAT;ID=0_hAT_chr11_-_1177994;sample=Regen_49_danglers.vitis.bwa_mem.bam;Parent=bin_hAT_chr11_-_1177894"""


def test_descendants():
    """
    Test for function descendants.
    This function returns a generator of :class:`gffutils.Feature` which should be coerced to a
    list of feature ids for comparison.
    """
    db = gffutils.create_db(GFF_UNSORTED,
                            dbfn=':memory:',
                            keep_order=True,
                            from_string=True,
                            merge_strategy='merge',
                            sort_attribute_values=True)
    feature = db['bin_Gypsy_chr11_-_10275156']  # Get feature object by ID
    answer = ['0_Gypsy_chr11_-_10275256',
              '0_Gypsy_chr11_-_10276179',
              '1_Gypsy_chr11_-_10275392',
              '2_Gypsy_chr11_-_10275513']
    assert [descendant.id for descendant in tectoolkit.filter_gff.descendants(feature, db)] == answer


def test_ancestors():
    """
    Test for function ancestors.
    This function returns a generator of :class:`gffutils.Feature` which should be coerced to a
    list of feature ids for comparison.
    """
    db = gffutils.create_db(GFF_UNSORTED,
                            dbfn=':memory:',
                            keep_order=True,
                            from_string=True,
                            merge_strategy='merge',
                            sort_attribute_values=True)
    feature = db['1_Gypsy_chr11_-_10275392']  # Get feature object by ID
    answer = ['bin_Gypsy_chr11_-_10275156']
    assert [ancestor.id for ancestor in tectoolkit.filter_gff.ancestors(feature, db)] == answer


@pytest.mark.parametrize('feature,filt,answer',
                         # Check if string values are equivalent they are ('Gypsy'=='Gypsy')
                         [('bin_Gypsy_chr11_-_10275156',
                           {'attribute': 'Name', 'operator': '==', 'value': 'Gypsy'},
                           True),
                          # Check if string values are not equivalent they are not not equivalent ('Gypsy'!='Gypsy')
                          ('bin_Gypsy_chr11_-_10275156',
                           {'attribute': 'Name', 'operator': '!=', 'value': 'Gypsy'},
                           False),
                          # Check if string values are not equivalent they are not equivalent ('Gypsy'!='Copia')
                          ('bin_Gypsy_chr11_-_10275156',
                           {'attribute': 'Name', 'operator': '!=', 'value': 'Copia'},
                           True),
                          # Check if numerical values are equal, they are (31==31.0)
                          ('bin_Gypsy_chr11_-_10275156',
                           {'attribute': 'read_count_max', 'operator': '==', 'value': '31.0'},
                           True),
                          # Check if numerical values are equal, they are not (31==32)
                          ('bin_Gypsy_chr11_-_10275156',
                           {'attribute': 'read_count_max', 'operator': '==', 'value': '32'},
                           False),
                          # Check if numerical value is greater, it is (31>30)
                          ('bin_Gypsy_chr11_-_10275156',
                           {'attribute': 'read_count_max', 'operator': '>', 'value': '30'},
                           True),
                          # Check if numerical value is greater, it is not (31>31)
                          ('bin_Gypsy_chr11_-_10275156',
                           {'attribute': 'read_count_max', 'operator': '>', 'value': '31'},
                           False),
                          # Check if numerical value is less, it is (31<32)
                          ('bin_Gypsy_chr11_-_10275156',
                           {'attribute': 'read_count_max', 'operator': '<', 'value': '32'},
                           True),
                          # Check if numerical value is less, it is not (31<30)
                          ('bin_Gypsy_chr11_-_10275156',
                           {'attribute': 'read_count_max', 'operator': '<', 'value': '30'},
                           False)],
                         ids=['', '', '', '', '', '', '', '', ''])
def test_matches_filter(feature, filt, answer):
    """
    Test factory for function matches_filter.
    Tests many cases with different filters (operators and attributes).

    :param feature: The name of a feature in the test db
    :type feature: str
    :param filt: A filter to test the feature with
    :type filt: dict[str, str]
    :param answer: Boolean specifying whether the feature matches all filters
    :type answer: bool
    """
    db = gffutils.create_db(GFF_UNSORTED,
                            dbfn=':memory:',
                            keep_order=True,
                            from_string=True,
                            merge_strategy='merge',
                            sort_attribute_values=True)
    feature = db[feature]
    assert tectoolkit.filter_gff.matches_filter(feature, filt) is answer


@pytest.mark.parametrize('feature,filters,answer',
                         # feature matches all filters
                         [('bin_Gypsy_chr11_-_10275156',
                           [{'attribute': 'Name', 'operator': '==', 'value': 'Gypsy'},
                            {'attribute': 'Name', 'operator': '!=', 'value': 'Copia'},
                            {'attribute': 'read_count_max', 'operator': '>=', 'value': '31'},
                            {'attribute': 'read_count_max', 'operator': '<', 'value': '32'}],
                           True),
                          # feature does not matches all filters
                          ('bin_Gypsy_chr11_-_10275156',
                           [{'attribute': 'Name', 'operator': '==', 'value': 'Gypsy'},
                            {'attribute': 'Name', 'operator': '!=', 'value': 'Copia'},
                            {'attribute': 'read_count_max', 'operator': '>=', 'value': '32'},  # 32!=31
                            {'attribute': 'read_count_max', 'operator': '<', 'value': '32'}],
                           False)
                          ],
                         ids=["", ""])
def test_matches_filters(feature, filters, answer):
    """
    Test for function matches_filters.
    Tests a case where feature matches all filters and a case where feature does not match all filters.

    :param feature: The name of a feature in the test db
    :type feature: str
    :param filters: A list of filters to test the feature with
    :type filters: list[dict[str, str]]
    :param answer: Boolean specifying whether the feature matches all filters
    :type answer: bool
    """
    db = gffutils.create_db(GFF_UNSORTED,
                            dbfn=':memory:',
                            keep_order=True,
                            from_string=True,
                            merge_strategy='merge',
                            sort_attribute_values=True)
    feature = db[feature]
    assert tectoolkit.filter_gff.matches_filters(feature, filters) is answer


@pytest.mark.parametrize('feature,filters,answer',
                         # '0_Gypsy_chr11_-_10276179' is a child of 'bin_Gypsy_chr11_-_10275156'
                         # 'bin_Gypsy_chr11_-_10275156' matches all filters
                         [('0_Gypsy_chr11_-_10276179',
                           [{'attribute': 'Name', 'operator': '==', 'value': 'Gypsy'},
                            {'attribute': 'Name', 'operator': '!=', 'value': 'Copia'},
                            {'attribute': 'read_count_max', 'operator': '>=', 'value': '31'},
                            {'attribute': 'read_count_max', 'operator': '<', 'value': '32'}],
                           True),
                          # '0_Copia_chr11_+_126915' is a child of 'bin_Copia_chr11_+_126799'
                          # neither of which matches all filters
                          ('0_Copia_chr11_+_126915',
                           [{'attribute': 'Name', 'operator': '==', 'value': 'Gypsy'},
                            {'attribute': 'Name', 'operator': '!=', 'value': 'Copia'},
                            {'attribute': 'read_count_max', 'operator': '>=', 'value': '31'},
                            {'attribute': 'read_count_max', 'operator': '<', 'value': '32'}],
                           False),
                          # 'bin_Gypsy_chr11_-_10275156' is a parent of '0_Gypsy_chr11_-_10276179'
                          # '0_Gypsy_chr11_-_10276179' matches all filters
                          ('bin_Gypsy_chr11_-_10275156',
                           [{'attribute': 'Name', 'operator': '==', 'value': 'Gypsy'},
                            {'attribute': 'sample', 'operator': '!=',
                             'value': 'Regen_49_danglers.vitis.bwa_mem.bam'}],
                           True)],
                         ids=['', '', ''])
def test_relative_matches_filters(feature, filters, answer):
    """
    Test for function relative_matches_filters.
    Tests whether a feature or at least one of its relatives matches a all of list of filters.

    :param feature: The name of a feature in the test db
    :type feature: str
    :param filters: A list of filters to test the feature with
    :type filters: list[dict[str, str]]
    :param answer: Boolean specifying whether the feature matches all filters
    :type answer: bool
    """
    db = gffutils.create_db(GFF_UNSORTED,
                            dbfn=':memory:',
                            keep_order=True,
                            from_string=True,
                            merge_strategy='merge',
                            sort_attribute_values=True)
    feature = db[feature]
    assert tectoolkit.filter_gff.relative_matches_filters(feature, filters, db) is answer


def test_filter_by_attributes():
    """
    Test for function filter_by_attributes.
    This function mutates the :class:`gffutils.FeatureDB` in place.
    This test selects features that meet the filters 'cluster_presence=1' and 'read_count_max>=20'
    as well as relatives (descendants or ancestors) of those features.
    """
    db = gffutils.create_db(GFF_UNSORTED,
                            dbfn=':memory:',
                            keep_order=True,
                            from_string=True,
                            merge_strategy='merge',
                            sort_attribute_values=True)
    filters = [{'attribute': 'cluster_presence', 'operator': '==', 'value': '1'},
               {'attribute': 'read_count_max', 'operator': '>=', 'value': '20'}]
    tectoolkit.filter_gff.filter_by_attributes(db, filters)
    assert '\n'.join([str(feature) for feature in db.all_features()]) == GFF_SINGLE_CLUSTERS_WITH_20_READS
