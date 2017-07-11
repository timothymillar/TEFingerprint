#! /usr/bin/env python
import os
import pytest
import subprocess

DATA_PATH = os.path.dirname(os.path.realpath(__file__)) + '/data/'

@pytest.mark.parametrize('query,answer',
                         [(['tef', 'filter-gff', DATA_PATH + 'testC-2017-06-29.gff', '-a', "proportions>=0.6"],
                           ('chr1\t.\t.\t3167\t3850\t.\t-\t.\tID=chr1_-_Gypsy_3167;category=Gypsy;sources=testA-2017-06-08.bam,testB-2017-06-08.bam;counts=8,0;proportions=1.0,0.0;color=#08306B\n'
                            'chr1\t.\t.\t21784\t22032\t.\t-\t.\tID=chr1_-_PIF_21784;category=PIF;sources=testA-2017-06-08.bam,testB-2017-06-08.bam;counts=4,6;proportions=0.4,0.6;color=#9ECAE1\n')
                           ),
                          (['tef', 'filter-gff', DATA_PATH + 'testC-2017-06-29.gff', '-a', "proportions>=0.6", 'category=Gypsy'],
                           ('chr1\t.\t.\t3167\t3850\t.\t-\t.\tID=chr1_-_Gypsy_3167;category=Gypsy;sources=testA-2017-06-08.bam,testB-2017-06-08.bam;counts=8,0;proportions=1.0,0.0;color=#08306B\n')
                           ),
                          (['tef', 'filter-gff', DATA_PATH + 'testC-2017-06-29.gff', '-c', 'start<20000', 'strand=+'],
                           ('chr1\t.\t.\t2402\t24267\t.\t+\t.\tID=chr1_+_Gypsy_2402;category=Gypsy;sources=testA-2017-06-08.bam,testB-2017-06-08.bam;counts=18,20;proportions=0.474,0.526;color=#C6DBEF\n')
                           ),
                          (['tef', 'filter-gff', DATA_PATH + 'testC-2017-06-29.gff', '-c', 'start>20000', '-a', 'category=Gypsy'],
                           ('chr1\t.\t.\t24737\t24969\t.\t-\t.\tID=chr1_-_Gypsy_24737;category=Gypsy;sources=testA-2017-06-08.bam,testB-2017-06-08.bam;counts=12,10;proportions=0.545,0.455;color=#C6DBEF\n')
                           ),
                          (['tef', 'filter-gff', DATA_PATH + 'testC-2017-06-29.gff', '-a', 'non=nan'],
                           ('\n')
                           )
                          ],
                         ids=['attribute*1', 'attribute*2', 'column*2', 'attribute+column', 'attribute=nan'])
def test_filter_gff(query, answer):
    """
    Tests the tef filter-gff command-line tool.

    :param query: list of arguments to run the tool
    :type query: list[str]
    :param answer: resulting gff string
    :type answer: str
    """
    output, error = subprocess.Popen(query, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    assert output.decode() == answer