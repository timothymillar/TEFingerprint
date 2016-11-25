#! /usr/bin/env python

import sys
import argparse
import numpy as np
from functools import reduce
from itertools import product
from multiprocessing import Pool
from tectoolkit import io
from tectoolkit.cluster import UnivariateLoci
from tectoolkit.gff import GffFeature
from tectoolkit.fingerprint import Fingerprint


class CompareProgram(object):
    """Main class for the compare program"""
    def __init__(self, arguments):
        """
        Init method for :class:`CompareProgram`.

        :param arguments: A list of commandline arguments to be parsed for the compare program
        """
        self.args = self.parse_args(arguments)
        if self.args.references == ['']:
            self.references = self._return_references(self.args.input_bams)
        else:
            self.references = self.args.references
        self.reference_lengths = self._return_reference_lengths(self.references, self.args.input_bams)

    def parse_args(self, args):
        """
        Defines an argument parser to handle commandline inputs for the compare program.

        :param args: A list of commandline arguments for the compare program

        :return: A dictionary like object of arguments and values for the compare program
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
        parser.add_argument('-b', '--bin_buffer',
                            type=int,
                            default=[0],
                            nargs=1,
                            help='Additional buffer to be added to margins of comparative bins')
        try:
            arguments = parser.parse_args(args)
        except:
            parser.print_help()
            sys.exit(0)
        else:
            return arguments

    def _return_references(self, input_bams):
        """
        Returns the intersection of references names in sam header for a collection of bam files.

        :param input_bams: A list of bam file paths
        :type input_bams: list[str]

        :return: A list of reference names
        :rtype: list[str]
        """
        bam_refs = [set(io.read_bam_references(bam)) for bam in input_bams]
        references = bam_refs[0]
        for refs in bam_refs[1:]:
            references = references.intersection(refs)
        references = list(references)
        return references

    def _return_reference_lengths(self, references, input_bams):
        """
        Returns a dictionary of reference lengths for a collection of bam files.
        Reference lengths must be identical in all bam files.

        :param references: A collection of reference names
        :type references: list[str]
        :param input_bams: A list of bam file paths
        :type input_bams: list[str]

        :return: A dictionary of reference lengths with names as keys
        :rtype: dict[str, int]
        """
        reference_dicts = [io.read_bam_reference_lengths(bam) for bam in input_bams]
        reference_lengths = {tuple(d[ref] for ref in references) for d in reference_dicts}
        assert len(reference_lengths) == 1
        return {ref: length for ref, length in zip(references, reference_lengths.pop())}

    def _build_jobs(self):
        """
        Builds a collection of tuples with parameters for compare jobs from class attributes.

        :return: A generator like object of tuples of job parameters
        :rtype: :class:`itertools.product`
        """
        return product([self.args.input_bams],
                       self.references,
                       self.args.families,
                       self.args.strands,
                       [self.args.eps],
                       self.args.min_reads,
                       self.args.bin_buffer)

    def _run_comparison(self, input_bams, reference, family, strand, eps, min_reads, bin_buffer):
        """
        Creates a :class:`Fingerprint` for each of a set of bam files and then combines these in an
        instance of :class:`FingerprintComparison` and prints the GFF3 formatted results to stdout.

        :param input_bams: A list of bam file paths
        :type input_bams: list[str]
        :param reference: The target reference name common to all input bam files
        :type reference: str
        :param family: The target family/category of TE's to be fingerprinted and compared
        :type family: str
        :param strand: The target strand to compared ('+' or '-')
        :type strand: str
        :param eps: The eps value(s) to be used in the cluster analysis (:class:`FUDC` for one values
        or :class:`HUDC` for two values)
        :type eps: list[int]
        :param min_reads: Minimum number of reads required to form cluster in cluster analysis
        :type min_reads: int
        :param bin_buffer: A value to extend the comparative bins by in both directions
        :type bin_buffer: int
        """
        fingerprints = (Fingerprint(bam, reference, family, strand, eps, min_reads) for bam in input_bams)
        reference_length = self.reference_lengths[reference]
        comparison = FingerprintComparison(tuple(fingerprints), bin_buffer, reference_length)
        for feature in comparison.to_gff():
            print(format(feature, 'nested'))

    def run(self):
        """
        Method to run the compare program.
        """
        jobs = self._build_jobs()
        if self.args.threads == 1:
            for job in jobs:
                self._run_comparison(*job)
        else:
            with Pool(self.args.threads) as pool:
                pool.starmap(self._run_comparison, jobs)


class FingerprintComparison(object):
    """Compare a collection of :class:`fingerprint` objects."""
    def __init__(self, fingerprints, buffer, reference_length):
        """
        Init method for :class:`FingerprintComparison`.

        :param fingerprints: a collection of :class:`fingerprint` objects
        :type fingerprints: tuple[:class:`FingerprintComparison`]
        :param buffer: A value to extend the comparative bins by in both directions
        :type buffer: int
        :param reference_length: The length of the reference used in fingerprints
        :type reference_length: int
        """
        self.fingerprints = fingerprints
        for f in fingerprints:
            assert type(f) == Fingerprint
        assert len({(f.reference,
                     f.family,
                     f.strand,
                     f.eps[0],
                     f.eps[-1],
                     f.min_reads) for f in fingerprints}) == 1
        self.reference = self.fingerprints[0].reference
        self.family = self.fingerprints[0].family
        self.strand = self.fingerprints[0].strand
        self.eps = self.fingerprints[0].eps
        self.min_reads = self.fingerprints[0].min_reads
        self.bin_loci = self._identify_bins()
        self._buffer_bins(buffer, reference_length)

    def _identify_bins(self):
        """
        Identifies bins (loci) in which to compare fingerprints from different samples (bam files).
        Bins are calculated by merging the overlapping fingerprints from all samples.

        :return: The union of loci among compared fingerprints
        :rtype: :class:`UnivariateLoci`
        """
        loci = reduce(UnivariateLoci.append, [f.loci for f in self.fingerprints])
        loci.melt()
        return loci

    def _buffer_bins(self, buffer, reference_length):
        """
        Expands comparative bins by buffer zone.
        """
        if buffer == 0:
            pass
        else:
            self.bin_loci.loci['start'] -= buffer
            self.bin_loci.loci['stop'] += buffer
            self.bin_loci.loci['start'][self.bin_loci.loci['start'] <= 0] = 0
            self.bin_loci.loci['stop'][self.bin_loci.loci['stop'] >= reference_length] = reference_length

    def _compare_bin(self, start, end):
        """
        Calculates basic statistics for a given bin (loci) to be compared across samples.
        Identifies (nested) fingerprint loci that fall within the bin bounds.
        Collects bin and nested fingerprint loci attributes into dictionaries of parameters for :class:`GffFeature`

        :param start: Start location of comparative bin
        :type start: int
        :param end: End location of comparative bin
        :type end: int

        :return: Bin and nested fingerprint loci attributes
        :rtype: (dict[str, any], list[dict[str, any]])
        """
        local_reads = tuple(f.reads.subset_by_locus(start, end) for f in self.fingerprints)
        local_clusters = tuple(f.loci.subset_by_locus(start, end) for f in self.fingerprints)
        sources = tuple(f.source for f in self.fingerprints)
        local_read_counts = np.array([len(reads) for reads in local_reads])
        read_count_max = max(local_read_counts)
        read_count_min = min(local_read_counts)
        read_presence = sum(local_read_counts != 0)
        read_absence = sum(local_read_counts == 0)
        local_cluster_counts = np.array([len(loci) for loci in local_clusters])
        cluster_presence = sum(local_cluster_counts != 0)
        cluster_absence = sum(local_cluster_counts == 0)
        bin_dict = {'seqid': self.reference,
                    'start': start,
                    'end': end,
                    'strand': self.strand,
                    'ID': "bin_{0}_{1}_{2}_{3}".format(self.family, self.reference, self.strand, start),
                    'Name': self.family,
                    'read_count_min': read_count_min,
                    'read_count_max': read_count_max,
                    'read_presence': read_presence,
                    'read_absence': read_absence,
                    'cluster_presence': cluster_presence,
                    'cluster_absence': cluster_absence}
        sample_dicts = []
        for number, sample in enumerate(zip(sources, local_clusters)):
            source, loci = sample
            if len(loci) == 0:
                pass
            else:
                sample_dicts += [{'seqid': self.reference,
                                  'start': start,
                                  'end': end,
                                  'strand': self.strand,
                                  'ID': "{0}_{1}_{2}_{3}_{4}".format(number,
                                                                     self.family,
                                                                     self.reference,
                                                                     self.strand,
                                                                     start),
                                  'Name': self.family,
                                  'sample': source} for start, end in loci]
        return bin_dict, sample_dicts

    def to_gff(self):
        """
        Creates :class:`GffFeature` object for each comparative bins in :class:`FingerprintComparison`.
        Component :class:`Fingerprint` loci found within the comparative bin are included as nested features.

        :return: A generator of nested :class:`GffFeature` objects
        :rtype: generator[:class:`GffFeature`]
        """
        for start, end in self.bin_loci:
            bin_dict, sample_dicts = self._compare_bin(start, end)
            feature = GffFeature(**bin_dict)
            feature.add_children(*[GffFeature(**d) for d in sample_dicts])
            yield feature

if __name__ == '__main__':
    pass
