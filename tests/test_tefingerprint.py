#! /usr/bin/env python3
import os
import pytest
import subprocess

DATA_PATH = os.path.dirname(os.path.realpath(__file__)) + '/data/'


class TestTEFingerprint:
    """Integration tests for batch fingerprinting"""
    CORES = 1

    @pytest.mark.parametrize('query,answer',
                             [(['tefingerprint', DATA_PATH + 'testA-2017-06-08.bam', '-f', 'Gypsy', 'PIF', '-e', '30000', '-m', '5', '-n', '0', '--gff', '-'],
                               ('chr1\tTEFingerprint\t.\t2452\t24217\t.\t+\t.\tcategory=Gypsy;median=2907;sample_0_name=testA-2017-06-08.bam;sample_0_count=18;ID=Gypsy_chr1_+_24217;pair=Gypsy_chr1_-_3217\n'
                                'chr1\tTEFingerprint\t.\t3217\t3800\t.\t-\t.\tcategory=Gypsy;median=3525;sample_0_name=testA-2017-06-08.bam;sample_0_count=8;ID=Gypsy_chr1_-_3217;pair=Gypsy_chr1_+_24217\n'
                                'chr1\tTEFingerprint\t.\t24787\t24919\t.\t-\t.\tcategory=Gypsy;median=24866;sample_0_name=testA-2017-06-08.bam;sample_0_count=12;ID=Gypsy_chr1_-_24787;pair=.\n')
                               ),
                              (['tefingerprint', DATA_PATH + 'testA-2017-06-08.bam', '-f', 'Gypsy', 'PIF', '-e', '30000', '-m', '5', '-n', '0', '-s', 'IDBCAN', '--gff', '-'],
                               ('chr1\tTEFingerprint\t.\t2452\t24217\t.\t+\t.\tcategory=Gypsy;median=2907;sample_0_name=testA-2017-06-08.bam;sample_0_count=18;ID=Gypsy_chr1_+_24217;pair=Gypsy_chr1_-_3217\n'
                                'chr1\tTEFingerprint\t.\t3217\t24919\t.\t-\t.\tcategory=Gypsy;median=24824;sample_0_name=testA-2017-06-08.bam;sample_0_count=20;ID=Gypsy_chr1_-_3217;pair=Gypsy_chr1_+_24217\n')
                               )
                              ],
                             ids=['hierarchical', 'non-hierarchical'])
    def test_hierarchical(self, query, answer):
        output, error = subprocess.Popen(query, stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE).communicate()
        assert output.decode() == answer

    @pytest.mark.parametrize('query,answer',
                             [(['tefingerprint', DATA_PATH + '/testA-2017-06-08.bam',  DATA_PATH + 'testB-2017-06-08.bam', '-f', 'Gypsy', 'PIF', '-e', '30000', '-m', '5', '-b', '10000', '-n', '0', '--no-colour', '--gff', '-'],
                               ('chr1\tTEFingerprint\t.\t2452\t24217\t.\t+\t.\tcategory=Gypsy;median=2973;sample_0_name=testA-2017-06-08.bam;sample_0_count=18;sample_1_name=testB-2017-06-08.bam;sample_1_count=20;ID=Gypsy_chr1_+_24217;pair=Gypsy_chr1_-_3217\n'
                                'chr1\tTEFingerprint\t.\t3217\t3800\t.\t-\t.\tcategory=Gypsy;median=3525;sample_0_name=testA-2017-06-08.bam;sample_0_count=8;sample_1_name=testB-2017-06-08.bam;sample_1_count=0;ID=Gypsy_chr1_-_3217;pair=Gypsy_chr1_+_24217\n'
                                'chr1\tTEFingerprint\t.\t21834\t21982\t.\t-\t.\tcategory=PIF;median=21956;sample_0_name=testA-2017-06-08.bam;sample_0_count=4;sample_1_name=testB-2017-06-08.bam;sample_1_count=6;ID=PIF_chr1_-_21834;pair=.\n'
                                'chr1\tTEFingerprint\t.\t24787\t24919\t.\t-\t.\tcategory=Gypsy;median=24872;sample_0_name=testA-2017-06-08.bam;sample_0_count=12;sample_1_name=testB-2017-06-08.bam;sample_1_count=10;ID=Gypsy_chr1_-_24787;pair=.\n')
                               ),
                              (['tefingerprint', DATA_PATH + '/testA-2017-06-08.bam', DATA_PATH + 'testB-2017-06-08.bam', '-f', 'Gypsy', 'PIF', '-e', '30000', '-m', '5', '-b', '10000', '-n', '0', '--gff', '-'],
                              ('chr1\tTEFingerprint\t.\t2452\t24217\t.\t+\t.\tcategory=Gypsy;median=2973;sample_0_name=testA-2017-06-08.bam;sample_0_count=18;sample_1_name=testB-2017-06-08.bam;sample_1_count=20;color=#C6DBEF;ID=Gypsy_chr1_+_24217;pair=Gypsy_chr1_-_3217\n'
                               'chr1\tTEFingerprint\t.\t3217\t3800\t.\t-\t.\tcategory=Gypsy;median=3525;sample_0_name=testA-2017-06-08.bam;sample_0_count=8;sample_1_name=testB-2017-06-08.bam;sample_1_count=0;color=#08306B;ID=Gypsy_chr1_-_3217;pair=Gypsy_chr1_+_24217\n'
                               'chr1\tTEFingerprint\t.\t21834\t21982\t.\t-\t.\tcategory=PIF;median=21956;sample_0_name=testA-2017-06-08.bam;sample_0_count=4;sample_1_name=testB-2017-06-08.bam;sample_1_count=6;color=#9ECAE1;ID=PIF_chr1_-_21834;pair=.\n'
                               'chr1\tTEFingerprint\t.\t24787\t24919\t.\t-\t.\tcategory=Gypsy;median=24872;sample_0_name=testA-2017-06-08.bam;sample_0_count=12;sample_1_name=testB-2017-06-08.bam;sample_1_count=10;color=#C6DBEF;ID=Gypsy_chr1_-_24787;pair=.\n')
                               ),
                              (['tefingerprint', DATA_PATH + '/testA-2017-06-08.bam',  DATA_PATH + 'testB-2017-06-08.bam', '-f', 'Gypsy', 'PIF', '-e', '30000', '-m', '5', '-b', '100000', '-n', '0', '--no-colour', '--gff', '-'],
                               ('chr1\tTEFingerprint\t.\t2452\t24217\t.\t+\t.\tcategory=Gypsy;median=2973;sample_0_name=testA-2017-06-08.bam;sample_0_count=18;sample_1_name=testB-2017-06-08.bam;sample_1_count=20;ID=Gypsy_chr1_+_24217;pair=Gypsy_chr1_-_3217\n'
                                'chr1\tTEFingerprint\t.\t3217\t24919\t.\t-\t.\tcategory=Gypsy;median=24857;sample_0_name=testA-2017-06-08.bam;sample_0_count=20;sample_1_name=testB-2017-06-08.bam;sample_1_count=10;ID=Gypsy_chr1_-_3217;pair=Gypsy_chr1_+_24217\n'
                                'chr1\tTEFingerprint\t.\t21834\t21982\t.\t-\t.\tcategory=PIF;median=21956;sample_0_name=testA-2017-06-08.bam;sample_0_count=4;sample_1_name=testB-2017-06-08.bam;sample_1_count=6;ID=PIF_chr1_-_21834;pair=.\n')
                               )
                              ],
                             ids=['no-overlap', 'colour', 'overlap'])
    def test_buffer(self, query, answer):
        output, error = subprocess.Popen(query, stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE).communicate()
        assert output.decode() == answer


class TestTEFingerprintMulticore(TestTEFingerprint):
    """Multi-core versions of integration tests for batch fingerprinting"""
    CORES = 2
