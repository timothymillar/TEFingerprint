#! /usr/bin/env python

import sys
import argparse
from itertools import product
from tectoolkit import libtec
from multiprocessing import Pool


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
            references = libtec.read_bam_references(input_bam)
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
        reads = libtec.read_sam_reads(input_bam, reference, family, strand)
        clusters = libtec.simple_cluster(reads, min_reads, eps)
        for cluster in clusters:
            gff = libtec.GffFeature(reference,
                                    start=cluster['start'],
                                    end=cluster['stop'],
                                    strand=strand,
                                    ID="{0}_{1}_{2}_{3}".format(family, reference, strand, cluster['start']),
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
        with Pool(self.args.threads) as p:
            results = p.starmap(self._fingerprint, jobs)
        for result in results:
            if result:
                print(result)
            else:
                pass


if __name__ == '__main__':
    pass