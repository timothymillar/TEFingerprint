#! /usr/bin/env python

import sys
import argparse
from itertools import product
from multiprocessing import Pool
from tectoolkit import bam_io
from tectoolkit.reads import ReadGroup, UnivariateLoci
from tectoolkit.gff_io import NestedFeature
from tectoolkit.cluster import UDC, HUDC


class FingerprintProgram(object):
    """Main class for the fingerprint program"""
    def __init__(self, arguments):
        """
        Init method for :class:`FingerprintProgram`.

        :param arguments: A list of commandline arguments to be parsed for the fingerprint program
        """
        self.args = self.parse_args(arguments)
        if self.args.references == ['']:
            self.references = bam_io.read_bam_references(self.args.input_bam[0])
        else:
            self.references = self.args.references

    def parse_args(self, args):
        """
        Defines an argument parser to handle commandline inputs for the fingerprint program.

        :param args: A list of commandline arguments for the fingerprint program

        :return: A dictionary like object of arguments and values for the fingerprint program
        """
        parser = argparse.ArgumentParser('Identify potential TE flanking regions')
        parser.add_argument('input_bam',
                            nargs=1,
                            help='A single bam file to be fingerprinted')
        parser.add_argument('-r', '--references',
                            type=str,
                            nargs='*',
                            default=[''],
                            help='The reference sequence(s) (e.g. chromosome) to be fingerprinted. ' + \
                                 'If left blank all references sequences in the input file will be used.')
        parser.add_argument('-f', '--families',
                            type=str,
                            nargs='*',
                            default=[''],
                            help='TE grouping(s) to be used. '
                                 'These must be exact string match\'s to start of read name and are used to split '
                                 'reads into categories for analysis')
        parser.add_argument('--mate_element_tag',
                            type=str,
                            default='ME',
                            help='Tag used in bam file to indicate the element type mate read')
        parser.add_argument('-s', '--strands',
                            type=str,
                            nargs='+',
                            choices=set("+-"),
                            default=['+', '-'],
                            help='Strand(s) to be analysed. Use + for forward or - for reverse. Default is to analyse '
                                 'both strands (separately).')
        parser.add_argument('-m', '--min_reads',
                            type=int,
                            default=[5],
                            nargs=1,
                            help='Minimum number of read tips required to be considered a cluster. '
                                 'This values is used in combination with epsilon to describe the density of '
                                 'read tips that is required for identification of a clusters. '
                                 'For every set of <min_reads> reads tips, if those reads are within epsilon range of '
                                 'one another, they are classified as a subcluster. '
                                 'Overlapping sets of subclusters are then merged to form clusters.')
        parser.add_argument('-e', '--eps',
                            type=int,
                            default=[100],
                            nargs='+',
                            help='Epsilon is the maximum allowable distance among a set of read tips to be '
                                 'considered a (sub)cluster. '
                                 'If a single value is given, the UDC algorithm will be used to identify all '
                                 'clusters at the specified density (defined by epsilon and min_points). '
                                 'If two values are given, they will be interpreted as maximum and minimum epsilon '
                                 'values using the Hierarchical HUDC algorithm.'
                                 'The maximum (or only) epsilon value given should be larger than the insert size, and '
                                 'the minimum epsilon (if used) should be much smaller (often zero) in order to find '
                                 'adequate support for child clusters. '
                                 'HUDC identifies all clusters at the '
                                 'maximum specified density and then attempts to split them into logical '
                                 'child clusters at all values of epsilon between maximum and minimum. '
                                 'The robustness of each parent cluster is compared to it\'s children. '
                                 'If the parent is more robust it is selected, otherwise the process is repeated for '
                                 'child cluster recursively until a parent or terminal (cluster with no children) '
                                 'is selected. ')
        parser.add_argument('-t', '--threads',
                            type=int,
                            default=1,
                            help='Maximum number of cpu threads to be used')
        try:
            arguments = parser.parse_args(args)
        except:
            parser.print_help()
            sys.exit(0)
        else:
            return arguments

    def _build_jobs(self):
        """
        Builds a collection of tuples with parameters for fingerprint jobs from class attributes.

        :return: A generator like object of tuples of job parameters
        :rtype: :class:`itertools.product`
        """
        return product(self.args.input_bam,
                       self.references,
                       self.args.strands,
                       [self.args.families],
                       [self.args.mate_element_tag],
                       [self.args.eps],
                       self.args.min_reads)

    @staticmethod
    def _fingerprint(input_bam, reference, strand, families, mate_element_tag, eps, min_reads):
        """
        Creates an instance of  :class:`Fingerprint` and prints the GFF3 formatted results to stdout.

        :param input_bam: A bam file path
        :type input_bam: str
        :param reference: The target reference name to be fingerprinted
        :type reference: str
        :param families: The target families/categories of TE's to be fingerprinted
        :type families: list[str]
        :param mate_element_tag: The sam tag that contains the mates element name in the bam file
        :type mate_element_tag: str
        :param strand: The target strand ('+' or '-') to fingerprinted
        :type strand: str
        :param eps: The eps value(s) to be used in the cluster analysis (:class:`FUDC` for one values
        or :class:`HUDC` for two values)
        :type eps: list[int]
        :param min_reads: Minimum number of reads required to form cluster in cluster analysis
        :type min_reads: int
        """
        read_groups = bam_io.read_bam_into_groups(input_bam, reference, strand, mate_element_tag, families)
        fingerprints = (Fingerprint(reads, eps, min_reads) for reads in read_groups)
        for fingerprint in fingerprints:
            print(format(fingerprint, 'gff'))

    def run(self):
        """

        :return:
        """
        jobs = self._build_jobs()
        if self.args.threads == 1:
            for job in jobs:
                self._fingerprint(*job)
        else:
            with Pool(self.args.threads) as pool:
                pool.starmap(self._fingerprint, jobs)


class Fingerprint(object):
    """Fingerprint a bam file"""
    def __init__(self, reads, eps, min_reads):
        """
        Init method for :class:`Fingerprint`.

        :param reads: A collection of reads
        :type reads: :class:`ReadGroup`
        :param eps: The eps value(s) to be used in the cluster analysis (:class:`FUDC` for one values
        or :class:`HUDC` for two values)
        :type eps: list[int]
        :param min_reads: Minimum number of reads required to form cluster in cluster analysis
        :type min_reads: int
        """
        self.eps = eps
        self.min_reads = min_reads
        self.reads = reads
        self.reference = self.reads.reference
        self.family = self.reads.grouping
        self.strand = self.reads.strand()
        self.source = self.reads.source
        self.loci = self._fit()

    def _fit(self):
        """
        Run a clustering algorithm on the targeted reads specified in an instance of :class:`Fingerprint`.
        If a single eps value was passed to the instance of :class:`Fingerprint` the :class:`FUDC` algorithm is used.
        If two eps values were passed to the instance of :class:`Fingerprint` the :class:`HUDC` algorithm is used.

        :return: Loci identified by the clustering algorithm
        :rtype: :class:`UnivariateLoci`
        """
        if len(self.eps) == 2:
            # use hierarchical clustering method
            max_eps, min_eps = max(self.eps), min(self.eps)
            hudc = HUDC(self.min_reads, max_eps, min_eps)
            hudc.fit(self.reads['tip'])
            return UnivariateLoci.from_iter(hudc.cluster_extremities())
        elif len(self.eps) == 1:
            # use flat clustering method
            eps = max(self.eps)
            fudc = UDC(self.min_reads, eps)
            fudc.fit(self.reads['tip'])
            return UnivariateLoci.from_iter(fudc.cluster_extremities())
        else:
            pass

    def __format__(self, code):
        assert code in {'gff'}
        return '\n'.join([str(l) for l in list(self._to_gff())])

    def _to_gff(self):
        """
        Creates :class:`GffFeature` object for each loci in :class:`Fingerprint`.

        :return: A generator of :class:`GffFeature` objects
        :rtype: generator[:class:`GffFeature`]
        """
        for start, end in self.loci:
            yield NestedFeature(seqid=self.reference,
                                start=start,
                                end=end,
                                strand=self.strand,
                                ID="{0}_{1}_{2}_{3}".format(self.family, self.reference, self.strand, start),
                                Name=self.family,
                                sample=self.source)

if __name__ == '__main__':
    pass