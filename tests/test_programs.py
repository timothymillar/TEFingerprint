#! /usr/bin/env python
import os
import pytest
import subprocess

DATA_PATH = os.path.dirname(os.path.realpath(__file__)) + '/data/'


@pytest.mark.parametrize('query,answer',
                         [(['tef', 'fingerprint', DATA_PATH + 'testA-2017-06-08.bam', '-f', 'Gypsy', 'PIF', '-e', '30000', '-m', '5'],
                           ('chr1\t.\t.\t2452\t24217\t.\t+\t.\tID=chr1_+_Gypsy_2452;category=Gypsy;source=testA-2017-06-08.bam\n'
                            'chr1\t.\t.\t3217\t3800\t.\t-\t.\tID=chr1_-_Gypsy_3217;category=Gypsy;source=testA-2017-06-08.bam\n'
                            'chr1\t.\t.\t24787\t24919\t.\t-\t.\tID=chr1_-_Gypsy_24787;category=Gypsy;source=testA-2017-06-08.bam\n')
                           ),
                          (['tef', 'fingerprint', DATA_PATH + 'testA-2017-06-08.bam', '-f', 'Gypsy', 'PIF', '-e', '30000', '-m', '5', '--non-hierarchical'],
                           ('chr1\t.\t.\t2452\t24217\t.\t+\t.\tID=chr1_+_Gypsy_2452;category=Gypsy;source=testA-2017-06-08.bam\n'
                            'chr1\t.\t.\t3217\t24919\t.\t-\t.\tID=chr1_-_Gypsy_3217;category=Gypsy;source=testA-2017-06-08.bam\n')
                           )
                          ],
                         ids=['hierarchical', 'non-hierarchical'])
def test_fingerprint(query, answer):
    """
    Tests the tef fingerprint command-line tool.

    :param query: list of arguments to run the tool
    :type query: list[str]
    :param answer: resulting gff string
    :type answer: str
    """
    output, error = subprocess.Popen(query, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    assert output.decode() == answer


@pytest.mark.parametrize('query,answer',
                         [(['tef', 'compare', DATA_PATH + 'testA-2017-06-08.bam', DATA_PATH + 'testB-2017-06-08.bam', '-f', 'Gypsy', 'PIF', '-e', '30000', '-m', '5', '-b 50'],
                           ('chr1\t.\t.\t2402\t24267\t.\t+\t.\tID=chr1_+_Gypsy_2402;category=Gypsy;sources=testA-2017-06-08.bam,testB-2017-06-08.bam;counts=18,20;proportions=0.474,0.526;color=#C6DBEF\n'
                            'chr1\t.\t.\t3167\t3850\t.\t-\t.\tID=chr1_-_Gypsy_3167;category=Gypsy;sources=testA-2017-06-08.bam,testB-2017-06-08.bam;counts=8,0;proportions=1.0,0.0;color=#08306B\n'
                            'chr1\t.\t.\t21784\t22032\t.\t-\t.\tID=chr1_-_PIF_21784;category=PIF;sources=testA-2017-06-08.bam,testB-2017-06-08.bam;counts=4,6;proportions=0.4,0.6;color=#9ECAE1\n'
                            'chr1\t.\t.\t24737\t24969\t.\t-\t.\tID=chr1_-_Gypsy_24737;category=Gypsy;sources=testA-2017-06-08.bam,testB-2017-06-08.bam;counts=12,10;proportions=0.545,0.455;color=#C6DBEF\n')
                           ),
                          (['tef', 'compare', DATA_PATH + 'testA-2017-06-08.bam', DATA_PATH + 'testB-2017-06-08.bam', '-f', 'Gypsy', 'PIF', '-e', '30000', '-m', '5', '-b 50', '--non-hierarchical'],
                           ('chr1\t.\t.\t2402\t24267\t.\t+\t.\tID=chr1_+_Gypsy_2402;category=Gypsy;sources=testA-2017-06-08.bam,testB-2017-06-08.bam;counts=18,20;proportions=0.474,0.526;color=#C6DBEF\n'
                            'chr1\t.\t.\t3167\t24969\t.\t-\t.\tID=chr1_-_Gypsy_3167;category=Gypsy;sources=testA-2017-06-08.bam,testB-2017-06-08.bam;counts=20,10;proportions=0.667,0.333;color=#9ECAE1\n'
                            'chr1\t.\t.\t21784\t22032\t.\t-\t.\tID=chr1_-_PIF_21784;category=PIF;sources=testA-2017-06-08.bam,testB-2017-06-08.bam;counts=4,6;proportions=0.4,0.6;color=#9ECAE1\n')
                           ),
                          (['tef', 'compare', DATA_PATH + 'testA-2017-06-08.bam', DATA_PATH + 'testB-2017-06-08.bam', '-f', 'Gypsy', 'PIF', '-e', '30000', '-m', '5', '-b 50', '--long-form-gff'],
                           ('chr1\t.\t.\t2402\t24267\t.\t+\t.\tID=chr1_+_Gypsy_2402_0;category=Gypsy;source=testA-2017-06-08.bam;count=18;proportion=0.474;color=#C6DBEF\n'
                            'chr1\t.\t.\t2402\t24267\t.\t+\t.\tID=chr1_+_Gypsy_2402_1;category=Gypsy;source=testB-2017-06-08.bam;count=20;proportion=0.526;color=#C6DBEF\n'
                            'chr1\t.\t.\t3167\t3850\t.\t-\t.\tID=chr1_-_Gypsy_3167_0;category=Gypsy;source=testA-2017-06-08.bam;count=8;proportion=1.0;color=#08306B\n'
                            'chr1\t.\t.\t3167\t3850\t.\t-\t.\tID=chr1_-_Gypsy_3167_1;category=Gypsy;source=testB-2017-06-08.bam;count=0;proportion=0.0;color=#C6DBEF\n'
                            'chr1\t.\t.\t21784\t22032\t.\t-\t.\tID=chr1_-_PIF_21784_0;category=PIF;source=testA-2017-06-08.bam;count=4;proportion=0.4;color=#C6DBEF\n'
                            'chr1\t.\t.\t21784\t22032\t.\t-\t.\tID=chr1_-_PIF_21784_1;category=PIF;source=testB-2017-06-08.bam;count=6;proportion=0.6;color=#9ECAE1\n'
                            'chr1\t.\t.\t24737\t24969\t.\t-\t.\tID=chr1_-_Gypsy_24737_0;category=Gypsy;source=testA-2017-06-08.bam;count=12;proportion=0.545;color=#C6DBEF\n'
                            'chr1\t.\t.\t24737\t24969\t.\t-\t.\tID=chr1_-_Gypsy_24737_1;category=Gypsy;source=testB-2017-06-08.bam;count=10;proportion=0.455;color=#C6DBEF\n'))
                          ],
                         ids=['hierarchical', 'non-hierarchical', 'long-form-gff'])
def test_compare(query, answer):
    """
    Tests the tef compare command-line tool.

    :param query: list of arguments to run the tool
    :type query: list[str]
    :param answer: resulting gff string
    :type answer: str
    """
    output, error = subprocess.Popen(query, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    assert output.decode() == answer
