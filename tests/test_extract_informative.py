#! /usr/bin/env python3
import tempfile
import os
import shutil
import gzip
import pytest
from tefingerprint._applications import extract_informative


@pytest.mark.parametrize('query,answer',
                         [('AGATC', 'GATCT'),
                          ('aGtcn', 'NGACT')])
def test_reverse_complement(query, answer):
    assert extract_informative.reverse_complement(query) == answer


FULL_LENGTH_READS = [
    {'element': 'MULE_1',
     'name': 'read01_forward_dangler',
     'quality': 'HHHHHGHHHHHHHHHHHHHFCHHHHGIEGGHHHHHHEHHHHHHHHHF@HEGGGGGFHHAEHHDFGHHHHHHHFHHHHHHE?D?GFHHEHHEEFFHGF;FE',
     'sequence': 'TCTGAGCTTAATATCGCCGGTCAACGGTCAAAATGGAGCTTTTTTCTTCATGCTGTTGGGGGGATTCACCCAACAAAAGATTTCCACTTCAGGCCCATTT'},
    {'element': 'Gypsy_2',
     'name': 'read02_forward_dangler',
     'quality': 'HHHHHHHHHHHHHHHHHHHHHHHHHHH?HHHHEFHGEGHGFGGFFHHHHHFHHHFGFHHHHHGEHHHFHGHHHFECFHFHHHHHH',
     'sequence': 'TATATAGCACGGATATATCGCCTGGTCAACGGTTGGTTAACACGGGCAATCAACAGTCAAAGCTCAAAATTGAGCCTTTTTTTTT'},
    {'element': 'MULE_3',
     'name': 'read03_forward_dangler',
     'quality': 'CBCFFFFFHHHHHJJJJJJJHIGGHIJJHIIJJGIJJIIIIJIJJJIIJJJJJJHHH==?CDF@BBCDACDDCCDCDDDCDEDDDDDDDCDB<ADDDDDD',
     'sequence': 'GAAATTAGGCGAAATATCGGTGTGTCAGTACCCACTGACACGACAAATAGCATGGATATATCGTGTGGTCAAAGCTTGATCAACATAGGTGGTCAACGGT'},
    {'element': 'Gypsy_1',
     'name': 'read04_forward_dangler',
     'quality': '@@@FFFDAHFFHH?FHFI=GHGGGI<FFGGGHEEE<FGE@<7;A5@FEFGBH:<ACECAC@CC>CCACCDDD:3>5@CCBBCDE34@B@CB<@B',
     'sequence': 'ATCGGTGTGTCGGTGGCTGTTATATAGCACGGATATATCGGTGGTCAACGCGGTCAAAGCTCAAAACTGAGCTATTTTTTGCTACTGTTGGGGG'},
    {'element': 'MULE_2',
     'name': 'read05_reverse_dangler',
     'quality': 'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHFHHHHHHHHGHHHHHHHHHFHEHHFHFFHHHHFEHEEHHHHHHBFHFHGHEHHHHFFGHHHE',
     'sequence': 'TTATATTGTGATATTTGATTATAATATATTTCGGTGGATATTTTGATACAAAATATCGGTAGACCTAAAATTGATCAAAACTTATAAAAATATAAAATAA'},
    {'element': 'Gypsy_3',
     'name': 'read06_reverse_dangler',
     'quality': '@@@FFFFFDHHHHJJJIJIIIGEDHEHIGIIIIJIJIJIHHIJJFIHGFGGEHJIJJJJGJIIIIHIIJIIIJJHFGHAHFFFFFFFEEEDEEDDDDDEE',
     'sequence': 'TATTTGATTATAATATATTTCGGTGGATATTTTGATACAAAATATCGGTAGACCTAAAATTGATCAAAACTTATAAAAATATAAAATAAACTCTTATAAA'},
    {'element': 'MULE_1',
     'name': 'read07_reverse_dangler',
     'quality': 'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHGHGHHGHHGHHHHHHHHHHHHHFHHGHEHHHGHGHHHHHGHHGFFHHHHGG=',
     'sequence': 'CCTTTAATTTGTTATTATTTTTATAAAATTTGTCTCAAAATTTTATTAGGAGATAATAATTACTATAATCAACTGAAAATTGATTTATATAACAAAATTC'},
    {'element': 'Gypsy_2',
     'name': 'read08_reverse_dangler',
     'quality': 'CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJIJJJJJHIGJJGJJJJJJJGIIJJJJJJJJJJJJJJJJJJJJJJJJJJJIJJIJJJJJJHHHHHHHHFF',
     'sequence': 'ACATAACTTATCATATTTGATAATAATATCCTATACGTCAATAAAAATATAAATTTTATAAATATATATTTATTATTAAGTTGCATTATATATTAATTTA'}
]


def test_extract_informative_full_reads():
    """
    Tests read extraction parsing and filtering.
    """
    answer = FULL_LENGTH_READS
    data_path = os.path.dirname(os.path.realpath(__file__)) + '/data/testPreprocessInput-2017-07-28.bam'
    query = list(extract_informative.extract_informative_full_reads(data_path))
    assert query == answer


# soft-tip results for combining

# tail clips of length >= 38 (default length)
TAIL_CLIPS = [
    {'element': 'MULE_3',
     'name': 'read09_forward_tail_F5',
     'quality': 'CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJIJJJJJJ',
     'sequence': 'TTATGATATATTTTTATTAATTGAAAATTTTTCATTCAATT'},
    {'element': 'Gypsy_3',
     'name': 'read12_reverse_tail_R5',
     'quality': '@CCFFFDDHFFDFIHIJJ=FGIJIIIIJGIJIIIFGCHIJIJIJJ',
     'sequence': 'TCCACTGTTTGGGGATTCGAACCCCCAACAAAAGGTTCAACATAT'},
    {'element': 'Gypsy_3',
     'name': 'read15_reverse_not_proper_tail_R5',  # included from v0.4
     'quality': '@CCFFFDDHFFDFIHIJJ=FGIJIIIIJGIJIIIFGCHIJIJIJJ',
     'sequence': 'TCCACTGTTTGGGGATTCGAACCCCCAACAAAAGGTTCAACATAT'}
]

# tip clips of length >= 38
TIP_CLIPS = [
    {'name': 'read17_mate_unmapped_F3',
     'element': 'MULE_2',
     'sequence': 'ATATAAATATAACAAAATATTTAAAGAATATTTTAATAAAC',
     'quality': 'DBBEFHDFHHFHIGG@IHF@IGHCFDGIHGDIGGHHGFBD@'},
    {'name': 'read18_mate_unmapped_F3',
     'element': 'MULE_2',
     'sequence': 'ATATAAATATAACAAAATATTTAAAGAATATTTTAATAAAC',
     'quality': 'DBBEFHDFHHFHIGG@IHF@IGHCFDGIHGDIGGHHGFBD@'}
]

# tail clips of length 38 > len >= 10
SHORT_TAIL_CLIPS = [
    {'element': 'MULE_2',
     'name': 'read11_forward_short_tail_F5',
     'quality': 'HHHHHHHHHHH',
     'sequence': 'AATTAAATGAA'},
    {'element': 'MULE_2',
     'name': 'read14_reverse_short_tail_R5',
     'quality': '@CCFFFDDHF',
     'sequence': 'TCCACTGTTT'}
]


SHORT_TIP_CLIPS = []


@pytest.mark.parametrize('include_tails,include_tips,soft_clip_minimum_length,answer',
                         [(False, False, 38, []),
                          (True, False, 38, TAIL_CLIPS),
                          (True, False, 10, TAIL_CLIPS + SHORT_TAIL_CLIPS),
                          (False, True, 38, TIP_CLIPS),
                          (False, True, 10, TIP_CLIPS + SHORT_TIP_CLIPS),
                          (True, True, 38, TAIL_CLIPS + TIP_CLIPS)
                          ])
def test_extract_informative_clips(include_tails,
                                   include_tips,
                                   soft_clip_minimum_length,
                                   answer):
    """
    Tests read extraction parsing and filtering.
    """
    data_path = os.path.dirname(os.path.realpath(__file__)) + '/data/testPreprocessInput-2017-07-28.bam'
    query = list(extract_informative.extract_informative_clips(
        data_path,
        include_soft_tails=include_tails,
        include_soft_tips=include_tips,
        minimum_soft_clip_len=soft_clip_minimum_length))

    # order not important
    # use sets of tuples to allow out of order comparision
    query = {(read['element'],
              read['name'],
              read['quality'],
              read['sequence']) for read in query}
    answer = {(read['element'],
               read['name'],
               read['quality'],
               read['sequence']) for read in answer}
    assert query == answer


def test_capture_elements_and_write_fastq():
    query = [{'element': 'MULE_1',
              'name': 'read01_forward_dangler',
              'quality': 'HHHHHGHHHHHHHHHHHHHFCHHHHGIEGGHHHHHHEHHHHHHHHHF@HEGGGGGFHHAEHHDFGHHHHHHHFHHHHHHE?D?GFHHEHHEEFFHGF;FE',
              'sequence': 'TCTGAGCTTAATATCGCCGGTCAACGGTCAAAATGGAGCTTTTTTCTTCATGCTGTTGGGGGGATTCACCCAACAAAAGATTTCCACTTCAGGCCCATTT'},
             {'element': 'Gypsy_2',
              'name': 'read02_forward_dangler',
              'quality': 'HHHHHHHHHHHHHHHHHHHHHHHHHHH?HHHHEFHGEGHGFGGFFHHHHHFHHHFGFHHHHHGEHHHFHGHHHFECFHFHHHHHH',
              'sequence': 'TATATAGCACGGATATATCGCCTGGTCAACGGTTGGTTAACACGGGCAATCAACAGTCAAAGCTCAAAATTGAGCCTTTTTTTTT'},
             {'element': 'MULE_3',
              'name': 'read03_forward_dangler',
              'quality': 'CBCFFFFFHHHHHJJJJJJJHIGGHIJJHIIJJGIJJIIIIJIJJJIIJJJJJJHHH==?CDF@BBCDACDDCCDCDDDCDEDDDDDDDCDB<ADDDDDD',
              'sequence': 'GAAATTAGGCGAAATATCGGTGTGTCAGTACCCACTGACACGACAAATAGCATGGATATATCGTGTGGTCAAAGCTTGATCAACATAGGTGGTCAACGGT'},
             {'element': 'Gypsy_1',
              'name': 'read04_forward_dangler',
              'quality': '@@@FFFDAHFFHH?FHFI=GHGGGI<FFGGGHEEE<FGE@<7;A5@FEFGBH:<ACECAC@CC>CCACCDDD:3>5@CCBBCDE34@B@CB<@B',
              'sequence': 'ATCGGTGTGTCGGTGGCTGTTATATAGCACGGATATATCGGTGGTCAACGCGGTCAAAGCTCAAAACTGAGCTATTTTTTGCTACTGTTGGGGG'},
             {'element': 'MULE_2',
              'name': 'read05_reverse_dangler',
              'quality': 'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHFHHHHHHHHGHHHHHHHHHFHEHHFHFFHHHHFEHEEHHHHHHBFHFHGHEHHHHFFGHHHE',
              'sequence': 'TTATATTGTGATATTTGATTATAATATATTTCGGTGGATATTTTGATACAAAATATCGGTAGACCTAAAATTGATCAAAACTTATAAAAATATAAAATAA'},
             {'element': 'Gypsy_3',
              'name': 'read06_reverse_dangler',
              'quality': '@@@FFFFFDHHHHJJJIJIIIGEDHEHIGIIIIJIJIJIHHIJJFIHGFGGEHJIJJJJGJIIIIHIIJIIIJJHFGHAHFFFFFFFEEEDEEDDDDDEE',
              'sequence': 'TATTTGATTATAATATATTTCGGTGGATATTTTGATACAAAATATCGGTAGACCTAAAATTGATCAAAACTTATAAAAATATAAAATAAACTCTTATAAA'},
             {'element': 'MULE_1',
              'name': 'read07_reverse_dangler',
              'quality': 'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHGHGHHGHHGHHHHHHHHHHHHHFHHGHEHHHGHGHHHHHGHHGFFHHHHGG=',
              'sequence': 'CCTTTAATTTGTTATTATTTTTATAAAATTTGTCTCAAAATTTTATTAGGAGATAATAATTACTATAATCAACTGAAAATTGATTTATATAACAAAATTC'},
             {'element': 'Gypsy_2',
              'name': 'read08_reverse_dangler',
              'quality': 'CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJIJJJJJHIGJJGJJJJJJJGIIJJJJJJJJJJJJJJJJJJJJJJJJJJJIJJIJJJJJJHHHHHHHHFF',
              'sequence': 'ACATAACTTATCATATTTGATAATAATATCCTATACGTCAATAAAAATATAAATTTTATAAATATATATTTATTATTAAGTTGCATTATATATTAATTTA'},
             {'element': 'MULE_3',
              'name': 'read09_forward_tail',
              'quality': 'CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJIJJJJJJ',
              'sequence': 'TTATGATATATTTTTATTAATTGAAAATTTTTCATTCAATT'},
             {'element': 'Gypsy_3',
              'name': 'read12_reverse_tail',
              'quality': '@CCFFFDDHFFDFIHIJJ=FGIJIIIIJGIJIIIFGCHIJIJIJJ',
              'sequence': 'TCCACTGTTTGGGGATTCGAACCCCCAACAAAAGGTTCAACATAT'}]

    answer_elements = {'read01_forward_dangler': 'MULE_1',
                       'read02_forward_dangler': 'Gypsy_2',
                       'read03_forward_dangler': 'MULE_3',
                       'read04_forward_dangler': 'Gypsy_1',
                       'read05_reverse_dangler': 'MULE_2',
                       'read06_reverse_dangler': 'Gypsy_3',
                       'read07_reverse_dangler': 'MULE_1',
                       'read08_reverse_dangler': 'Gypsy_2',
                       'read09_forward_tail': 'MULE_3',
                       'read12_reverse_tail': 'Gypsy_3'}

    answer_fastq = """@read01_forward_dangler
TCTGAGCTTAATATCGCCGGTCAACGGTCAAAATGGAGCTTTTTTCTTCATGCTGTTGGGGGGATTCACCCAACAAAAGATTTCCACTTCAGGCCCATTT
+
HHHHHGHHHHHHHHHHHHHFCHHHHGIEGGHHHHHHEHHHHHHHHHF@HEGGGGGFHHAEHHDFGHHHHHHHFHHHHHHE?D?GFHHEHHEEFFHGF;FE
@read02_forward_dangler
TATATAGCACGGATATATCGCCTGGTCAACGGTTGGTTAACACGGGCAATCAACAGTCAAAGCTCAAAATTGAGCCTTTTTTTTT
+
HHHHHHHHHHHHHHHHHHHHHHHHHHH?HHHHEFHGEGHGFGGFFHHHHHFHHHFGFHHHHHGEHHHFHGHHHFECFHFHHHHHH
@read03_forward_dangler
GAAATTAGGCGAAATATCGGTGTGTCAGTACCCACTGACACGACAAATAGCATGGATATATCGTGTGGTCAAAGCTTGATCAACATAGGTGGTCAACGGT
+
CBCFFFFFHHHHHJJJJJJJHIGGHIJJHIIJJGIJJIIIIJIJJJIIJJJJJJHHH==?CDF@BBCDACDDCCDCDDDCDEDDDDDDDCDB<ADDDDDD
@read04_forward_dangler
ATCGGTGTGTCGGTGGCTGTTATATAGCACGGATATATCGGTGGTCAACGCGGTCAAAGCTCAAAACTGAGCTATTTTTTGCTACTGTTGGGGG
+
@@@FFFDAHFFHH?FHFI=GHGGGI<FFGGGHEEE<FGE@<7;A5@FEFGBH:<ACECAC@CC>CCACCDDD:3>5@CCBBCDE34@B@CB<@B
@read05_reverse_dangler
TTATATTGTGATATTTGATTATAATATATTTCGGTGGATATTTTGATACAAAATATCGGTAGACCTAAAATTGATCAAAACTTATAAAAATATAAAATAA
+
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHFHHHHHHHHGHHHHHHHHHFHEHHFHFFHHHHFEHEEHHHHHHBFHFHGHEHHHHFFGHHHE
@read06_reverse_dangler
TATTTGATTATAATATATTTCGGTGGATATTTTGATACAAAATATCGGTAGACCTAAAATTGATCAAAACTTATAAAAATATAAAATAAACTCTTATAAA
+
@@@FFFFFDHHHHJJJIJIIIGEDHEHIGIIIIJIJIJIHHIJJFIHGFGGEHJIJJJJGJIIIIHIIJIIIJJHFGHAHFFFFFFFEEEDEEDDDDDEE
@read07_reverse_dangler
CCTTTAATTTGTTATTATTTTTATAAAATTTGTCTCAAAATTTTATTAGGAGATAATAATTACTATAATCAACTGAAAATTGATTTATATAACAAAATTC
+
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHGHGHHGHHGHHHHHHHHHHHHHFHHGHEHHHGHGHHHHHGHHGFFHHHHGG=
@read08_reverse_dangler
ACATAACTTATCATATTTGATAATAATATCCTATACGTCAATAAAAATATAAATTTTATAAATATATATTTATTATTAAGTTGCATTATATATTAATTTA
+
CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJIJJJJJHIGJJGJJJJJJJGIIJJJJJJJJJJJJJJJJJJJJJJJJJJJIJJIJJJJJJHHHHHHHHFF
@read09_forward_tail
TTATGATATATTTTTATTAATTGAAAATTTTTCATTCAATT
+
CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJIJJJJJJ
@read12_reverse_tail
TCCACTGTTTGGGGATTCGAACCCCCAACAAAAGGTTCAACATAT
+
@CCFFFDDHFFDFIHIJJ=FGIJIIIIJGIJIIIFGCHIJIJIJJ
"""
    tempdir = tempfile.mkdtemp()
    fastq = tempdir + '/temp.fastq.gz'
    elements = extract_informative.capture_elements_and_write_fastq(query, fastq)
    assert elements == answer_elements

    with gzip.open(fastq, 'r') as handle:
        assert handle.read().decode() == answer_fastq
    shutil.rmtree(tempdir)


def test_tag_danglers():
    pass
