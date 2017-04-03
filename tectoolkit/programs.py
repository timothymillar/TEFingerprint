#! /usr/bin/env python

import argparse
from tectoolkit import batch
from tectoolkit import bamio


class FingerprintProgram(object):
    """Class for the fingerprint program"""
    def __init__(self, arguments):
        self.args = self.parse_args(arguments)
        result = batch.fingerprint(bams=self.args.input_bams,
                                   references=self.args.references,
                                   categories=self.args.families,
                                   transposon_tag=self.args.mate_element_tag[0],
                                   min_reads=self.args.min_reads[0],
                                   eps=self.args.epsilon[0],
                                   min_eps=self.args.min_eps[0],
                                   hierarchical=self.args.hierarchical_clustering[0],
                                   cores=self.args.threads[0])
        print(result.as_gff())

    @staticmethod
    def parse_args(args):
        """
        Defines an argument parser to handle commandline inputs for the fingerprint program.

        :param args: A list of commandline arguments for the fingerprint program

        :return: A dictionary like object of arguments and values for the fingerprint program
        """
        parser = argparse.ArgumentParser('Identify potential TE flanking regions')
        parser.add_argument('input_bams',
                            nargs='+',
                            help='One or more bam files to be fingerprinted')
        parser.add_argument('-r', '--references',
                            type=str,
                            nargs='+',
                            default=[None],
                            help='The reference sequence(s) (e.g. chromosome) to be fingerprinted. ' + \
                                 'If left blank all references sequences in the input file will be used.')
        parser.add_argument('-f', '--families',
                            type=str,
                            nargs='*',
                            default=[''],
                            help='TE categories to be used. '
                                 'These must be exact string match\'s to start of read name and are used to split '
                                 'reads into categories for analysis')
        parser.add_argument('--mate_element_tag',
                            type=str,
                            default=['ME'],
                            nargs=1,
                            help='Tag used in bam file to indicate the element type mate read')
        parser.add_argument('-m', '--min_reads',
                            type=int,
                            default=[10],
                            nargs=1,
                            help='Minimum number of read tips required to be considered a cluster. '
                                 'This values is used in combination with epsilon to describe the density of '
                                 'read tips that is required for identification of a clusters. '
                                 'For every set of <min_reads> reads tips, if those reads are within epsilon range of '
                                 'one another, they are classified as a subcluster. '
                                 'Overlapping sets of subclusters are then merged to form clusters.')
        parser.add_argument('-e', '--epsilon',
                            type=int,
                            default=[250],
                            nargs=1,
                            help='Epsilon is the maximum allowable distance among a set of read tips to be '
                                 'considered a (sub)cluster. '
                                 'The epsilon value given should be larger than the insert size. '
                                 'HUDC identifies all clusters at the '
                                 'maximum specified density and then attempts to split them into logical '
                                 'child clusters at all values of epsilon between maximum and minimum. '
                                 'The robustness of each parent cluster is compared to it\'s children. '
                                 'If the parent is more robust it is selected, otherwise the process is repeated for '
                                 'child cluster recursively until a parent or terminal (cluster with no children) '
                                 'is selected. ')
        parser.add_argument('--min_eps',
                            type=int,
                            default=[0],
                            nargs=1,
                            help='Minimum eps values used by the HUDC algorithm when calculating support for clusters. '
                                 'This should usually be left as the default value of 0.')
        parser.add_argument('--hierarchical_clustering',
                            type=bool,
                            default=[True],
                            nargs=1,
                            help='By default hierarchical HUDC algorithm is used. If this is set to False, '
                                 'the non-hierarchical UDC algorithm is used and min_eps is inored.')
        parser.add_argument('-j', '--join_distance',
                            type=int,
                            default=[0],
                            nargs=1,
                            help='The maximum distance allowable between neighbouring clusters (of the same family '
                                 'and opposite strands) to be associated with one another as a pair')
        parser.add_argument('-t', '--threads',
                            type=int,
                            default=[1],
                            nargs=1,
                            help='Maximum number of cpu threads to be used')
        args = parser.parse_args(args)
        if args.references == [None]:
            args.references = bamio.extract_bam_references(*args.input_bams)
        return args


class ComparisonProgram(object):
    """Class for the comparison program"""
    def __init__(self, arguments):
        self.args = self.parse_args(arguments)
        result = batch.comparison(bams=self.args.input_bams,
                                  references=self.args.references,
                                  categories=self.args.families,
                                  transposon_tag=self.args.mate_element_tag[0],
                                  min_reads=self.args.min_reads[0],
                                  eps=self.args.epsilon[0],
                                  min_eps=self.args.min_eps[0],
                                  hierarchical=self.args.hierarchical_clustering[0],
                                  bin_buffer=self.args.bin_buffer[0],
                                  cores=self.args.threads[0])
        if self.args.long_form[0] is True:
            print(result.as_flat_gff())
        else:
            print(result.as_gff())

    @staticmethod
    def parse_args(args):
        """
        Defines an argument parser to handle commandline inputs for the comparison program.

        :param args: A list of commandline arguments for the comparison program

        :return: A dictionary like object of arguments and values for the comparison program
        """
        parser = argparse.ArgumentParser('Compare potential TE flanking regions from multiple samples')
        parser.add_argument('input_bams',
                            nargs='+',
                            help='A list of two or more bam files to be compared')
        parser.add_argument('-r', '--references',
                            type=str,
                            nargs='*',
                            default=[None],
                            help='The reference sequence(s) (e.g. chromosome) to be fingerprinted. ' + \
                                 'If left blank all references sequences in the input file will be used.')
        parser.add_argument('-f', '--families',
                            type=str,
                            nargs='*',
                            default=[''],
                            help='TE categories to be used. '
                                 'These must be exact string match\'s to start of read name and are used to split '
                                 'reads into categories for analysis')
        parser.add_argument('--mate_element_tag',
                            type=str,
                            default=['ME'],
                            nargs=1,
                            help='Tag used in bam file to indicate the element type mate read')
        parser.add_argument('-m', '--min_reads',
                            type=int,
                            default=[10],
                            nargs=1,
                            help='Minimum number of read tips required to be considered a cluster. '
                                 'This values is used in combination with epsilon to describe the density of '
                                 'read tips that is required for identification of a clusters. '
                                 'For every set of <min_reads> reads tips, if those reads are within epsilon range of '
                                 'one another, they are classified as a subcluster. '
                                 'Overlapping sets of subclusters are then merged to form clusters.')
        parser.add_argument('-e', '--epsilon',
                            type=int,
                            default=[250],
                            nargs=1,
                            help='Epsilon is the maximum allowable distance among a set of read tips to be '
                                 'considered a (sub)cluster. '
                                 'The epsilon value given should be larger than the insert size. '
                                 'HUDC identifies all clusters at the '
                                 'maximum specified density and then attempts to split them into logical '
                                 'child clusters at all values of epsilon between maximum and minimum. '
                                 'The robustness of each parent cluster is compared to it\'s children. '
                                 'If the parent is more robust it is selected, otherwise the process is repeated for '
                                 'child cluster recursively until a parent or terminal (cluster with no children) '
                                 'is selected. ')
        parser.add_argument('--min_eps',
                            type=int,
                            default=[0],
                            nargs=1,
                            help='Minimum eps values used by the HUDC algorithm when calculating support for clusters. '
                                 'This should usually be left as the default value of 0.')
        parser.add_argument('--hierarchical_clustering',
                            type=bool,
                            default=[True],
                            nargs=1,
                            help='By default hierarchical HUDC algorithm is used. If this is set to False, '
                                 'the non-hierarchical UDC algorithm is used and min_eps is inored.')
        parser.add_argument('-j', '--join_distance',
                            type=int,
                            default=[0],
                            nargs=1,
                            help='The maximum distance allowable between neighbouring clusters (of the same family '
                                 'and opposite strands) to be associated with one another as a pair')
        parser.add_argument('-b', '--bin_buffer',
                            type=int,
                            default=[None],
                            nargs=1,
                            help='Additional buffer to be added to margins of comparative bins. '
                                 'This is used avoid identifying small clusters as unique, when these is only '
                                 'slight miss-match in read positions across samples (i.e. false positives). '
                                 'By default this will use the same value as epsilon.')
        parser.add_argument('--long_form',
                            type=bool,
                            default=[False],
                            nargs=1,
                            help='If True, the resulting gff file will contain one feature per sample per bin. '
                                 'This avoids nested lists in the feature attributes but results in many more features')
        parser.add_argument('-t', '--threads',
                            type=int,
                            default=[1],
                            nargs=1,
                            help='Maximum number of cpu threads to be used')
        args = parser.parse_args(args)
        if args.references == [None]:
            args.references = bamio.extract_bam_references(*args.input_bams)
        if args.bin_buffer == [None]:
            args.bin_buffer = args.epsilon
        return args


if __name__ == '__main__':
    pass
