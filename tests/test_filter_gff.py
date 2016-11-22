#! /usr/bin/env python

import gffutils
from tectoolkit.filter_gff import GffFilterDB


class TestGffFilterDB:
    """Tests for class GffFilterDB"""
    def test__str__(self):
        query = ("chr11\t.\t.\t671652\t672285\t.\t-\t.\tName=hAT;ID=bin_hAT_chr11_-_671652;read_count_max=11;read_count_min=0;cluster_absence=2;read_presence=1;cluster_presence=1;read_absence=2\n"
                 "chr11\t.\t.\t671752\t672185\t.\t-\t.\tName=hAT;ID=0_hAT_chr11_-_671752;sample=Regen_49_danglers.vitis.bwa_mem.bam;Parent=bin_hAT_chr11_-_671652\n"
                 "chr11\t.\t.\t1177277\t1177878\t.\t-\t.\tName=hAT;ID=bin_hAT_chr11_-_1177277;read_count_max=45;read_count_min=0;cluster_absence=2;read_presence=1;cluster_presence=1;read_absence=2\n"
                 "chr11\t.\t.\t1177377\t1177778\t.\t-\t.\tName=hAT;ID=0_hAT_chr11_-_1177377;sample=Regen_49_danglers.vitis.bwa_mem.bam;Parent=bin_hAT_chr11_-_1177277\n"
                 "chr11\t.\t.\t1177894\t1178429\t.\t-\t.\tName=hAT;ID=bin_hAT_chr11_-_1177894;read_count_max=20;read_count_min=0;cluster_absence=2;read_presence=1;cluster_presence=1;read_absence=2\n"
                 "chr11\t.\t.\t1177994\t1178329\t.\t-\t.\tName=hAT;ID=0_hAT_chr11_-_1177994;sample=Regen_49_danglers.vitis.bwa_mem.bam;Parent=bin_hAT_chr11_-_1177894")
        db = GffFilterDB(gffutils.create_db(query,
                                            dbfn=':memory:',
                                            keep_order=True,
                                            from_string=True,
                                            merge_strategy='merge',
                                            sort_attribute_values=True))
        print(db)
        assert str(db) == query

    def test_descendants(self):
        """descendants"""
        pass

    def test_ancestors(self):
        """ancestors"""
        pass

    def test_matches_filter(self):
        """matches_filter"""
        pass

    def test_matches_filters(self):
        """matches_filters"""
        pass

    def test_relative_matches_filters(self):
        """relative_matches_filters"""
        pass

    def test_filter_by_attributes(self):
        """filter_by_attributes"""
        pass