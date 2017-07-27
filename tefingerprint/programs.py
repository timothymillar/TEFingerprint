#! /usr/bin/env python

import argparse
from tefingerprint import batch
from tefingerprint import bamio


class FingerprintProgram(object):
    """Class for the fingerprint program"""
    def __init__(self, arguments):
        self.args = self.parse_args(arguments)
        result = batch.fingerprint(bams=self.args.bams,
                                   references=self.args.references,
                                   categories=self.args.families,
                                   quality=self.args.mapping_quality[0],
                                   transposon_tag=self.args.mate_element_tag[0],
                                   min_reads=self.args.minimum_reads[0],
                                   eps=self.args.epsilon[0],
                                   min_eps=self.args.minimum_epsilon[0],
                                   hierarchical=self.args.hierarchical_clustering,
                                   cores=self.args.threads[0])
        print(result.as_gff())
        if self.args.feature_csv[0]:
            with open(self.args.feature_csv[0], 'w') as csv:
                csv.write(result.as_csv())

    @staticmethod
    def parse_args(args):
        """
        Defines an argument parser to handle commandline inputs for the fingerprint program.

        :param args: A list of commandline arguments for the fingerprint program

        :return: A dictionary like object of arguments and values for the fingerprint program
        """
        parser = argparse.ArgumentParser('Identify potential TE flanking regions.')
        parser.add_argument('bams',
                            nargs='+',
                            help='One or more bam files to be fingerprinted.')
        parser.add_argument('-r', '--references',
                            type=str,
                            nargs='+',
                            default=[None],
                            help='The reference sequence(s) (e.g. chromosome) to be fingerprinted. '
                                 'If left blank all references sequences in the input file will be used.')
        parser.add_argument('-f', '--families',
                            type=str,
                            nargs='*',
                            default=[''],
                            help='TE categories to be used. '
                                 'These must be exact string match\'s to start of read name and are used to split '
                                 'reads into categories for analysis.')
        parser.add_argument('-q', '--mapping-quality',
                            type=int,
                            nargs=1,
                            default=[30],
                            help='Minimum mapping quality of reads.')
        parser.add_argument('--mate-element-tag',
                            type=str,
                            default=['ME'],
                            nargs=1,
                            help='Tag used in bam file to indicate the element type of the mate read.')
        parser.add_argument('-m', '--minimum-reads',
                            type=int,
                            default=[10],
                            nargs=1,
                            help='Minimum number of read tips required to be considered a cluster.')
        parser.add_argument('-e', '--epsilon',
                            type=int,
                            default=[250],
                            nargs=1,
                            help='The maximum allowable distance among a set of read tips to be considered a cluster.')
        parser.add_argument('--minimum-epsilon',
                            type=int,
                            default=[0],
                            nargs=1,
                            help='Minimum epsilon values used when calculating support for clusters. '
                                 'This is only used in hierarchical clustering and should usually be left '
                                 'as the default value of 0.')
        parser.set_defaults(hierarchical_clustering=True)
        parser.add_argument('--non-hierarchical',
                            dest='hierarchical_clustering',
                            action='store_false',
                            help='Use non-hierarchical clustering.')
        parser.add_argument('--feature-csv',
                            type=str,
                            default=[False],
                            nargs=1,
                            help='Optionally write a csv file of features.')
        parser.add_argument('-t', '--threads',
                            type=int,
                            default=[1],
                            nargs=1,
                            help='Maximum number of cpu threads to be used.')
        args = parser.parse_args(args)
        if args.references == [None]:
            args.references = bamio.extract_bam_reference_strings(*args.bams)
        return args


class ComparisonProgram(object):
    """Class for the comparison program"""
    def __init__(self, arguments):
        self.args = self.parse_args(arguments)
        result = batch.comparison(bams=self.args.bams,
                                  references=self.args.references,
                                  categories=self.args.families,
                                  quality=self.args.mapping_quality[0],
                                  transposon_tag=self.args.mate_element_tag[0],
                                  min_reads=self.args.minimum_reads[0],
                                  eps=self.args.epsilon[0],
                                  min_eps=self.args.minimum_epsilon[0],
                                  hierarchical=self.args.hierarchical_clustering,
                                  fingerprint_buffer=self.args.buffer_fingerprints[0],
                                  bin_buffer=self.args.buffer_comparative_bins[0],
                                  cores=self.args.threads[0])
        if self.args.long_form_gff is True:
            print(result.as_flat_gff())
        else:
            print(result.as_gff())
        if self.args.character_csv[0]:
            with open(self.args.character_csv[0], 'w') as csv:
                csv.write(result.as_character_csv())
        if self.args.feature_csv[0]:
            with open(self.args.feature_csv[0], 'w') as csv:
                csv.write(result.as_flat_csv())

    @staticmethod
    def parse_args(args):
        """
        Defines an argument parser to handle commandline inputs for the comparison program.

        :param args: A list of commandline arguments for the comparison program

        :return: A dictionary like object of arguments and values for the comparison program
        """
        parser = argparse.ArgumentParser('Compare potential TE flanking regions from multiple samples')
        parser.add_argument('bams',
                            nargs='+',
                            help='A list of two or more bam files to be compared.')
        parser.add_argument('-r', '--references',
                            type=str,
                            nargs='*',
                            default=[None],
                            help='The reference sequence(s) (e.g. chromosome) to be fingerprinted. '
                                 'If left blank all references sequences in the input file will be used.')
        parser.add_argument('-f', '--families',
                            type=str,
                            nargs='*',
                            default=[''],
                            help='TE categories to be used. '
                                 'These must be exact string match\'s to start of read name and are used to split '
                                 'reads into categories for analysis.')
        parser.add_argument('-q', '--mapping-quality',
                            type=int,
                            nargs=1,
                            default=[30],
                            help='Minimum mapping quality of reads.')
        parser.add_argument('--mate-element-tag',
                            type=str,
                            default=['ME'],
                            nargs=1,
                            help='Tag used in bam file to indicate the element type of the mate read.')
        parser.add_argument('-m', '--minimum-reads',
                            type=int,
                            default=[10],
                            nargs=1,
                            help='Minimum number of read tips required to be considered a cluster.')
        parser.add_argument('-e', '--epsilon',
                            type=int,
                            default=[250],
                            nargs=1,
                            help='The maximum allowable distance among a set of read tips to be considered a cluster.')
        parser.add_argument('--minimum-epsilon',
                            type=int,
                            default=[0],
                            nargs=1,
                            help='Minimum epsilon values used when calculating support for clusters. '
                                 'This is only used in hierarchical clustering and should usually be left '
                                 'as the default value of 0.')
        parser.set_defaults(hierarchical_clustering=True)
        parser.add_argument('--non-hierarchical',
                            dest='hierarchical_clustering',
                            action='store_false',
                            help='Use non-hierarchical clustering.')
        parser.add_argument('-b', '--buffer-fingerprints',
                            type=int,
                            default=[0],
                            nargs=1,
                            help='Additional buffer to be added to margins of fingerprints. '
                                 'This is used avoid identifying small clusters as unique, when these is only '
                                 'slight miss-match in read positions across samples (i.e. false positives).')
        parser.add_argument('--buffer-comparative-bins',
                            type=int,
                            default=[0],
                            nargs=1,
                            help='Additional buffer to be added to margins of comparative bins.')
        parser.set_defaults(long_form_gff=False)
        parser.add_argument('--long-form-gff',
                            dest='long_form_gff',
                            action='store_true',
                            help='The resulting gff output will contain one feature per sample per bin. '
                                 'This avoids nested lists in the feature attributes.')
        parser.add_argument('--feature-csv',
                            type=str,
                            default=[False],
                            nargs=1,
                            help='Optionally write a csv file of features (one row per sample per comparative bin).')
        parser.add_argument('--character-csv',
                            type=str,
                            default=[False],
                            nargs=1,
                            help='Optionally write a csv file of data as character states '
                                 '(rows of samples * columns of comparative bins).')
        parser.add_argument('-t', '--threads',
                            type=int,
                            default=[1],
                            nargs=1,
                            help='Maximum number of cpu threads to be used.')
        args = parser.parse_args(args)
        if args.references == [None]:
            args.references = bamio.extract_bam_reference_strings(*args.bams)
        return args


if __name__ == '__main__':
    pass
