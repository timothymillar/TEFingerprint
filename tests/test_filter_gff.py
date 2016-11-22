#! /usr/bin/env python

import gffutils
from tectoolkit.filter_gff import GffFilterDB

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


class TestGffFilterDB:
    """Tests for class GffFilterDB"""
    def test__str__(self):
        """
        Test for __str__ method of :class:`GffFilterDB`.
        This method should return an identical GFF formatted string to the input to the database
        with the exception that attributes should be ordered consistently.
        """
        db = GffFilterDB(gffutils.create_db(GFF_UNSORTED,
                                            dbfn=':memory:',
                                            keep_order=True,
                                            from_string=True,
                                            merge_strategy='merge',
                                            sort_attribute_values=True))
        assert str(db) == GFF_SORTED

    def test_descendants(self):
        """
        Test for method descendants.
        This method returns a generator of :class:`gffutils.Feature` which should be coerced to a
        list of feature ids for comparison.
        """
        filter_db = GffFilterDB(gffutils.create_db(GFF_UNSORTED,
                                                   dbfn=':memory:',
                                                   keep_order=True,
                                                   from_string=True,
                                                   merge_strategy='merge',
                                                   sort_attribute_values=True))
        query = filter_db.db['bin_Gypsy_chr11_-_10275156']  # Get feature object by ID
        answer = ['0_Gypsy_chr11_-_10275256',
                  '0_Gypsy_chr11_-_10276179',
                  '1_Gypsy_chr11_-_10275392',
                  '2_Gypsy_chr11_-_10275513']
        assert [descendant.id for descendant in filter_db.descendants(query)] == answer

    def test_ancestors(self):
        """
        Test for method ancestors.
        This method returns a generator of :class:`gffutils.Feature` which should be coerced to a
        list of feature ids for comparison.
        """
        filter_db = GffFilterDB(gffutils.create_db(GFF_UNSORTED,
                                                   dbfn=':memory:',
                                                   keep_order=True,
                                                   from_string=True,
                                                   merge_strategy='merge',
                                                   sort_attribute_values=True))
        query = filter_db.db['1_Gypsy_chr11_-_10275392']  # Get feature object by ID
        answer = ['bin_Gypsy_chr11_-_10275156']
        assert [ancestor.id for ancestor in filter_db.ancestors(query)] == answer

    def test_matches_filter(self):
        """
        Test for method matches_filter.
        Tests many cases with different filters (operators and attributes).
        """
        filter_db = GffFilterDB(gffutils.create_db(GFF_UNSORTED,
                                                   dbfn=':memory:',
                                                   keep_order=True,
                                                   from_string=True,
                                                   merge_strategy='merge',
                                                   sort_attribute_values=True))

        # Check if string values are equivalent they are ('Gypsy'=='Gypsy')
        query = filter_db.db['bin_Gypsy_chr11_-_10275156']  # Get feature object by ID
        filt = {'attribute': 'Name', 'operator': '=', 'value': 'Gypsy'}  # Name=Gypsy
        assert filter_db.matches_filter(query, filt) is True

        # Check if string values are equal, they are not('31'=='31.0')
        query = filter_db.db['bin_Gypsy_chr11_-_10275156']
        filt = {'attribute': 'read_count_max', 'operator': '=', 'value': '31.0'}  # Note floating point
        assert filter_db.matches_filter(query, filt) is False

        # Check if string values are not equivalent they are not not equivalent ('Gypsy'!='Gypsy')
        query = filter_db.db['bin_Gypsy_chr11_-_10275156']
        filt = {'attribute': 'Name', 'operator': '!=', 'value': 'Gypsy'}  # Name=Gypsy
        assert filter_db.matches_filter(query, filt) is False

        # Check if string values are not equivalent they are are not equivalent ('Gypsy'!='Copia')
        query = filter_db.db['bin_Gypsy_chr11_-_10275156']
        filt = {'attribute': 'Name', 'operator': '!=', 'value': 'Copia'}  # Name=Copia
        assert filter_db.matches_filter(query, filt) is True

        # Check if numerical values are equal, they are (31==31.0)
        query = filter_db.db['bin_Gypsy_chr11_-_10275156']
        filt = {'attribute': 'read_count_max', 'operator': '==', 'value': '31.0'}  # Note floating point
        assert filter_db.matches_filter(query, filt) is True

        # Check if numerical values are equal, they are not (31==32)
        query = filter_db.db['bin_Gypsy_chr11_-_10275156']
        filt = {'attribute': 'read_count_max', 'operator': '==', 'value': '32'}
        assert filter_db.matches_filter(query, filt) is False

        # Check if numerical value is greater, it is (31>30)
        query = filter_db.db['bin_Gypsy_chr11_-_10275156']
        filt = {'attribute': 'read_count_max', 'operator': '>', 'value': '30'}
        assert filter_db.matches_filter(query, filt) is True

        # Check if numerical value is greater, it is not (31>31)
        query = filter_db.db['bin_Gypsy_chr11_-_10275156']
        filt = {'attribute': 'read_count_max', 'operator': '>', 'value': '31'}
        assert filter_db.matches_filter(query, filt) is False

        # Check if numerical value is less, it is (31<32)
        query = filter_db.db['bin_Gypsy_chr11_-_10275156']
        filt = {'attribute': 'read_count_max', 'operator': '<', 'value': '32'}
        assert filter_db.matches_filter(query, filt) is True

        # Check if numerical value is less, it is not (31<30)
        query = filter_db.db['bin_Gypsy_chr11_-_10275156']
        filt = {'attribute': 'read_count_max', 'operator': '<', 'value': '30'}
        assert filter_db.matches_filter(query, filt) is False

    def test_matches_filters(self):
        """
        Test for method matches_filters.
        Tests a case where feature matches all filters and a case where feature does not match all filters.
        """
        filter_db = GffFilterDB(gffutils.create_db(GFF_UNSORTED,
                                                   dbfn=':memory:',
                                                   keep_order=True,
                                                   from_string=True,
                                                   merge_strategy='merge',
                                                   sort_attribute_values=True))

        # feature matches all filters
        feature = filter_db.db['bin_Gypsy_chr11_-_10275156']
        filters = [{'attribute': 'Name', 'operator': '=', 'value': 'Gypsy'},
                   {'attribute': 'Name', 'operator': '!=', 'value': 'Copia'},
                   {'attribute': 'read_count_max', 'operator': '>=', 'value': '31'},
                   {'attribute': 'read_count_max', 'operator': '<', 'value': '32'}]
        assert filter_db.matches_filters(feature, filters) is True

        # feature does not matches all filters
        feature = filter_db.db['bin_Gypsy_chr11_-_10275156']
        filters = [{'attribute': 'Name', 'operator': '=', 'value': 'Gypsy'},
                   {'attribute': 'Name', 'operator': '!=', 'value': 'Copia'},
                   {'attribute': 'read_count_max', 'operator': '>=', 'value': '32'},  # feature does not matches filter
                   {'attribute': 'read_count_max', 'operator': '<', 'value': '32'}]
        assert filter_db.matches_filters(feature, filters) is False

    def test_relative_matches_filters(self):
        """relative_matches_filters"""
        pass

    def test_filter_by_attributes(self):
        """filter_by_attributes"""
        pass