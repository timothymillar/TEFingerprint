import numpy as np
from tectoolkit.reads import StrandReads
from tectoolkit.fingerprint import Fingerprint


class TestFingerprint:
    """"""
    def test_fingerprint_empty(self):
        """
        Test for fingerprint of an empty read group.
        This should run but not produce any gff lines.
        """
        reads = StrandReads(np.array([], dtype=StrandReads.DTYPE_READ),
                            reference='chr1',
                            grouping='Gypsy',
                            source='file.bam')
        fingerprint = Fingerprint(reads, [200, 10], 10)
        gff = format(fingerprint, 'gff')

        answer = ''
        assert gff == answer

    def test_fingerprint(self):
        """
        Test for automatic calculation of :class:`Fingerprint`.
        Effectively integration test for :class:`ReadGroup`,
        :class:`Fingerprint`, :class:`HUDC` and :class:`NestedFeature`.
        """
        # note that the second value for each read is it's tail position and is not used in any calculation
        reads = StrandReads(np.array([(0, 0, '+', 'Gypsy27_uniqueID'),
                                      (   0, 0, '+', 'Gypsy27_uniqueID'),
                                      (  60, 0, '+', 'Gypsy27_uniqueID'),
                                      (  61, 0, '+', 'Gypsy27_uniqueID'),
                                      (  61, 0, '+', 'Gypsy27_uniqueID'),
                                      (  61, 0, '+', 'Gypsy27_uniqueID'),
                                      (  76, 0, '+', 'Gypsy27_uniqueID'),
                                      (  78, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 122, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 122, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 141, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 183, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 251, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 260, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 260, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 263, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 263, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 267, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 267, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 288, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 288, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 295, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 300, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 310, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 310, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 317, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 317, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 334, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 334, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 335, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 338, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 338, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 338, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 338, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 340, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 342, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 342, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 344, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 344, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 358, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 367, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 370, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 370, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 377, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 387, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 402, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 403, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 410, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 410, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 410, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 418, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 418, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 424, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 424, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 577, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 857, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 879, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 921, 0, '+', 'Gypsy27_uniqueID'),
                                      ( 921, 0, '+', 'Gypsy27_uniqueID'),
                                      (1007, 0, '+', 'Gypsy27_uniqueID'),
                                      (1031, 0, '+', 'Gypsy27_uniqueID'),
                                      (1051, 0, '+', 'Gypsy27_uniqueID'),
                                      (1051, 0, '+', 'Gypsy27_uniqueID'),
                                      (1059, 0, '+', 'Gypsy27_uniqueID'),
                                      (1071, 0, '+', 'Gypsy27_uniqueID'),
                                      (1071, 0, '+', 'Gypsy27_uniqueID'),
                                      (1080, 0, '+', 'Gypsy27_uniqueID'),
                                      (1094, 0, '+', 'Gypsy27_uniqueID'),
                                      (1094, 0, '+', 'Gypsy27_uniqueID'),
                                      (1110, 0, '+', 'Gypsy27_uniqueID'),
                                      (1110, 0, '+', 'Gypsy27_uniqueID'),
                                      (1113, 0, '+', 'Gypsy27_uniqueID'),
                                      (1113, 0, '+', 'Gypsy27_uniqueID'),
                                      (1183, 0, '+', 'Gypsy27_uniqueID'),
                                      (1189, 0, '+', 'Gypsy27_uniqueID'),
                                      (1200, 0, '+', 'Gypsy27_uniqueID'),
                                      (1200, 0, '+', 'Gypsy27_uniqueID'),
                                      (1217, 0, '+', 'Gypsy27_uniqueID'),
                                      (1234, 0, '+', 'Gypsy27_uniqueID'),
                                      (1234, 0, '+', 'Gypsy27_uniqueID'),
                                      (1591, 0, '+', 'Gypsy27_uniqueID'),
                                      (1620, 0, '+', 'Gypsy27_uniqueID'),
                                      (1620, 0, '+', 'Gypsy27_uniqueID'),
                                      (1662, 0, '+', 'Gypsy27_uniqueID'),
                                      (1686, 0, '+', 'Gypsy27_uniqueID'),
                                      (1707, 0, '+', 'Gypsy27_uniqueID'),
                                      (1755, 0, '+', 'Gypsy27_uniqueID'),
                                      (1828, 0, '+', 'Gypsy27_uniqueID'),
                                      (1828, 0, '+', 'Gypsy27_uniqueID'),
                                      (1848, 0, '+', 'Gypsy27_uniqueID'),
                                      (1848, 0, '+', 'Gypsy27_uniqueID'),
                                      (1848, 0, '+', 'Gypsy27_uniqueID'),
                                      (1848, 0, '+', 'Gypsy27_uniqueID'),
                                      (1851, 0, '+', 'Gypsy27_uniqueID'),
                                      (1851, 0, '+', 'Gypsy27_uniqueID'),
                                      (1852, 0, '+', 'Gypsy27_uniqueID'),
                                      (1917, 0, '+', 'Gypsy27_uniqueID')], dtype=StrandReads.DTYPE_READ),
                            reference='chr1', grouping='Gypsy', source='file.bam')
        fingerprint = Fingerprint(reads, [200, 10], 10)
        gff = format(fingerprint, 'gff')

        answer = '\n'.join(['chr1\t.\t.\t0\t577\t.\t+\t.\tID=Gypsy_chr1_+_0;Name=Gypsy;sample=file.bam',
                            'chr1\t.\t.\t879\t1234\t.\t+\t.\tID=Gypsy_chr1_+_879;Name=Gypsy;sample=file.bam',
                            'chr1\t.\t.\t1662\t1917\t.\t+\t.\tID=Gypsy_chr1_+_1662;Name=Gypsy;sample=file.bam'])
        assert gff == answer

