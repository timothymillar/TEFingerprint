#! /usr/bin/env python

import sys
import argparse
from functools import reduce
from itertools import product
from multiprocessing import Pool
from tectoolkit import io
from tectoolkit.classes import ReadLoci
from tectoolkit.gff import GffFeature
from tectoolkit.fingerprint import Fingerprint


class CompareProgram(object):
    """"""
    def __init__(self, arguments):
        self.args = self.parse_args(arguments)

    def parse_args(self, args):
        """

        :param args:
        :return:
        """
        parser = argparse.ArgumentParser('Compare potential TE flanking regions')
        parser.add_argument('input_bams',
                            nargs='+',
                            help='A list of two or more bam files to be compared')
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
                            nargs='+',
                            help='Maximum allowable distance among read tips to be considered a cluster')
        parser.add_argument('-m', '--min_reads',
                            type=int,
                            default=[5],
                            nargs=1,
                            help='Minimum allowable number of reads (tips) to be considered a cluster')
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

    def _build_jobs(self, input_bams, references, families, strands, eps, min_reads):
        """

        :param input_bams:
        :param references:
        :param families:
        :param strands:
        :param eps:
        :param min_reads:
        :return:
        """
        if references == ['']:
            bam_refs = [set(io.read_bam_references(bam)) for bam in input_bams]
            references = bam_refs[0]
            for refs in bam_refs[1:]:
                references = references.intersection(refs)
            references = list(references)
        else:
            pass
        return product([input_bams],
                       references,
                       families,
                       strands,
                       [eps],
                       min_reads)

    def _run_comparison(self, input_bams, reference, family, strand, eps, min_reads):
        fingerprints = (Fingerprint(bam, reference, family, strand, eps, min_reads) for bam in input_bams)
        comparison = FingerprintComparison(tuple(fingerprints))
        for feature in comparison.to_gff():
            print(feature)

    def run(self):
        """

        :return:
        """
        jobs = self._build_jobs(self.args.input_bams,
                                self.args.references,
                                self.args.families,
                                self.args.strands,
                                self.args.eps,
                                self.args.min_reads)
        if self.args.threads == 1:
            for job in jobs:
                self._run_comparison(*job)
        else:
            with Pool(self.args.threads) as pool:
                pool.starmap(self._run_comparison, jobs)


class FingerprintComparison(object):

    def __init__(self, fingerprints):
        self.fingerprints = fingerprints
        for f in fingerprints:
            assert type(f) == Fingerprint
        assert len({(f.reference,
                     f.family,
                     f.strand,
                     f.eps,
                     f.min_reads) for f in fingerprints}) == 1
        self.reference = self.fingerprints[0].reference
        self.family = self.fingerprints[0].family
        self.strand = self.fingerprints[0].strand
        self.eps = self.fingerprints[0].eps
        self.min_reads = self.fingerprints[0].min_reads
        self.bin_loci = self._identify_bins()

    def _identify_bins(self):
        loci = reduce(ReadLoci.append, [f.loci for f in self.fingerprints])
        loci.melt()
        return loci

    def to_gff(self):
        for start, end in self.bin_loci:
            yield GffFeature(self.reference,
                             start=start,
                             end=end,
                             strand=self.strand,
                             ID="{0}_{1}_{2}_{3}".format(self.family, self.reference, self.strand, start),
                             Name=self.family)

if __name__ == '__main__':
    pass
