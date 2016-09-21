#! /usr/bin/env python

import sys
import argparse
from itertools import product
from multiprocessing import Pool
from tectoolkit import io
from tectoolkit.classes import ReadGroup, ReferenceLoci, GffFeature



class Fingerprint:
    """"""
    def __init__(self, arguments):
        self.args = self.parse_args(arguments)

    def parse_args(self, args):
        """

        :param args:
        :return:
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
                            help='TE grouping(s) to be used. Must be exact string match(s) to start of read name')
        parser.add_argument('-s', '--strands',
                            type=str,
                            nargs='+',
                            choices=set("+-."),
                            default=['+', '-'],
                            help='Strand(s) to be analysed. Use + for forward, - for reverse and . for both')
        parser.add_argument('-e', '--eps',
                            type=int,
                            default=[100],
                            nargs=1,
                            help='Maximum allowable distance among read tips to be considered a cluster')
        parser.add_argument('-m', '--min_reads',
                            type=int,
                            default=[5],
                            nargs=1,
                            help='Minimum allowable number of read tips to be considered a cluster')
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

    def _build_jobs(self, input_bam, references, families, strands, eps, min_reads):
        """

        :param input_bam:
        :param references:
        :param families:
        :param strands:
        :param eps:
        :param min_reads:
        :return:
        """
        if references == ['']:
            references = io.read_bam_references(input_bam)
        else:
            pass
        return product(input_bam,
                       references,
                       families,
                       strands,
                       eps,
                       min_reads)

    def _fingerprint(self, input_bam, reference, family, strand, eps, min_reads):
        """

        :param input_bam:
        :param reference:
        :param family:
        :param strand:
        :param eps:
        :param min_reads:
        :return:
        """
        strings = io.read_bam_strings(input_bam, reference=reference, family=family, strand=strand)
        reads = ReadGroup.from_sam_strings(strings, strand=strand)
        reads.sort()
        clusters = ReferenceLoci.from_simple_cluster(reads['tip'], min_reads, eps)
        for start, end in clusters:
            gff = GffFeature(reference,
                             start=start,
                             end=end,
                             strand=strand,
                             ID="{0}_{1}_{2}_{3}".format(family, reference, strand, start),
                             Name=family)
            print(str(gff))

    def run(self):
        """

        :return:
        """
        jobs = self._build_jobs(self.args.input_bam,
                                self.args.references,
                                self.args.families,
                                self.args.strands,
                                self.args.eps,
                                self.args.min_reads)
        with Pool(self.args.threads) as pool:
            results = pool.starmap(self._fingerprint, jobs)
        for result in results:
            if result:
                print(result)
            else:
                pass


if __name__ == '__main__':
    pass