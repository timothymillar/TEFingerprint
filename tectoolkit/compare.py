#! /usr/bin/env python

import sys
import argparse
import numpy as np
from functools import reduce
from scipy.stats import kruskal
from itertools import product
from multiprocessing import Pool
from tectoolkit import io
from tectoolkit.classes import ReadGroup, GffFeature
from tectoolkit.cluster import UnivariateLoci
from tectoolkit.cluster import FlatUnivariateDensityCluster as FUDC
from tectoolkit.cluster import HierarchicalUnivariateDensityCluster as HUDC


class Compare:
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

    def _robust_kruskal(self, *args):
        try:
            kruskal_h, kruskal_p = kruskal(*args)
        except ValueError:
            kruskal_h, kruskal_p = 0, 1
        finally:
            return kruskal_h, kruskal_p

    def _fcluster(self, reads, min_reads, eps):
        fudc = FUDC(min_reads, eps)
        fudc.fit(reads['tip'])
        return fudc.loci

    def _hcluster(self, reads, min_reads, max_eps, min_eps):
        hudc = HUDC(min_reads, max_eps, min_eps)
        hudc.fit(reads['tip'])
        return hudc.loci

    def _compare(self, input_bams, reference, family, strand, eps, min_reads):
        """

        :param input_bams:
        :param reference:
        :param family:
        :param strand:
        :param eps:
        :param min_reads:
        :return:
        """
        sams = (io.read_bam_strings(bam, reference=reference, family=family, strand=strand) for bam in input_bams)
        read_groups = (ReadGroup.from_sam_strings(sam, strand=strand) for sam in sams)
        if len(eps) == 1:
            eps = max(eps)
            loci_groups = (self._fcluster(reads, min_reads, eps) for reads in read_groups)
        if len(eps) == 2:
            max_eps, min_eps = max(eps), min(eps)
            loci_groups = (self._hcluster(reads, min_reads, max_eps, min_eps) for reads in read_groups)
        else:
            pass  # throw error
        loci = UnivariateLoci(reduce(np.append, loci_groups))
        loci.melt()
        for start, end in loci:
            #locus_read_groups = tuple(group.sub_group_by_locus(start, end, margin=eps) for group in read_groups)
            #kruskal_h, kruskal_p = self._robust_kruskal(*[group['tip'] for group in locus_read_groups])
            #counts = np.fromiter(map(len, locus_read_groups), dtype=int)
            #count_stdev = (counts/np.max(counts)).std()
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
        jobs = self._build_jobs(self.args.input_bams,
                                self.args.references,
                                self.args.families,
                                self.args.strands,
                                self.args.eps,
                                self.args.min_reads)
        if self.args.threads == 1:
            for job in jobs:
                self._compare(*job)
        else:
            with Pool(self.args.threads) as pool:
                pool.starmap(self._compare, jobs)


if __name__ == '__main__':
    pass
