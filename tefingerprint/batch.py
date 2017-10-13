#! /usr/bin/env python

from itertools import product
from multiprocessing import Pool
from tefingerprint import loci


def _fingerprint_worker(bams,
                        categories,
                        reference,
                        quality,
                        transposon_tag,
                        min_reads,
                        eps,
                        min_eps,
                        hierarchical):
    """
    Worker function for batch fingerprinting.
    Runs a single job dispatched from fingerprint.
    """
    reads = loci.append(*[loci.InformativeReadLoci.from_bam(bam,
                                                            categories,
                                                            references=reference,
                                                            quality=quality,
                                                            tag=transposon_tag) for bam in bams])
    return reads.fingerprint(min_reads,
                             eps,
                             min_eps=min_eps,
                             hierarchical=hierarchical)


def fingerprint(bams=None,
                categories=None,
                references=None,
                quality=30,
                transposon_tag='ME',
                min_reads=None,
                eps=None,
                min_eps=0,
                hierarchical=True,
                cores=1):
    """
    Multi-processes pipeline for batch fingerprinting one or more bam file(s) by transposon categories.

    :param bams: a list of bam file(s) to be fingerprinted
    :type bams: list[str]
    :param categories: targeted transposon categories
    :type categories: list[str]
    :param references: targeted (slices of) references with format 'name' or 'name:minimum-maximum'
    :type references: list[str]
    :param quality: minimum mapping quality
    :type quality: int
    :param transposon_tag: the two letter sam tag containing the transposon associated with each read
    :type transposon_tag: str
    :param min_reads: minimum number of read tips required to form a fingerprint locus
    :type min_reads: int
    :param eps: maximum distance allowed among each set of read tips to form a fingerprint locus
    :type eps: int
    :param min_eps: minimum epsilon used when calculating support for child clusters in the HUDC algorithm
    :type min_eps: int
    :param hierarchical: if false, the non-hierarchical UDC algorithm is used and min_eps ignored
    :type hierarchical: bool
    :param cores: number of processes to use for the batch job
    :type cores: int

    :return: Fingerprints of bam files
    :rtype: :class:`loci.GenomicBins`
    """
    jobs = product([bams],
                   [categories],
                   references,  # job per reference
                   [quality],
                   [transposon_tag],
                   [min_reads],
                   [eps],
                   [min_eps],
                   [hierarchical])

    if cores == 1:
        return loci.append(*[_fingerprint_worker(*job) for job in jobs])
    else:
        with Pool(cores) as pool:
            return loci.append(*pool.starmap(_fingerprint_worker, jobs))


def _comparison_worker(bams,
                       categories,
                       reference,
                       quality,
                       transposon_tag,
                       min_reads,
                       eps,
                       min_eps,
                       number_of_elements,
                       hierarchical,
                       fingerprint_buffer):
    """
    Worker function for batch fingerprint comparisons.
    Runs a single job dispatched from comparison.
    """
    reads = loci.append(*[loci.InformativeReadLoci.from_bam(bam,
                                                            categories,
                                                            references=reference,
                                                            quality=quality,
                                                            tag=transposon_tag)
                          for bam in bams])
    bins = reads.cluster(min_reads, eps, min_eps=min_eps, hierarchical=hierarchical)
    return bins.merge_sources().buffered_melt(fingerprint_buffer).count_reads(reads,
                                                                              n_common_elements=number_of_elements)


def comparison(bams=None,
               categories=None,
               references=None,
               quality=30,
               transposon_tag='ME',
               min_reads=None,
               eps=None,
               min_eps=0,
               number_of_elements=3,
               hierarchical=True,
               fingerprint_buffer=0,
               cores=1):
    """
    Multi-processes pipeline for batch comparison of multiple bam file(s) by transposon fingerprints.

    :param bams: a list of bam file(s) to be fingerprinted
    :type bams: list[str]
    :param categories: targeted transposon categories
    :type categories: list[str]
    :param references: targeted (slices of) references with format 'name' or 'name:minimum-maximum'
    :type references: list[str]
    :param quality: minimum mapping quality
    :type quality: int
    :param transposon_tag: the two letter sam tag containing the transposon associated with each read
    :type transposon_tag: str
    :param min_reads: minimum number of read tips required to form a fingerprint locus
    :type min_reads: int
    :param eps: maximum distance allowed among each set of read tips to form a fingerprint locus
    :type eps: int
    :param min_eps: minimum epsilon used when calculating support for child clusters in the HUDC algorithm
    :type min_eps: int
    :param number_of_elements: number of most common elements to count for ech cluster.
    :type number_of_elements: int
    :param hierarchical: if false, the non-hierarchical UDC algorithm is used and min_eps ignored
    :type hierarchical: bool
    :param fingerprint_buffer: value to buffer fingerprints by, defaults to 0
    :type fingerprint_buffer: int
    :param cores: number of processes to use for the batch job
    :type cores: int

    :return: Fingerprint comparisons of bam files
    :rtype: :class:`loci.BinCounts`
    """
    jobs = product([bams],
                   [categories],
                   references,  # job per reference
                   [quality],
                   [transposon_tag],
                   [min_reads],
                   [eps],
                   [min_eps],
                   [number_of_elements],
                   [hierarchical],
                   [fingerprint_buffer])

    if cores == 1:
        return loci.append(*[_comparison_worker(*job) for job in jobs])
    else:
        with Pool(cores) as pool:
            return loci.append(*pool.starmap(_comparison_worker, jobs))


if __name__ == '__main__':
    pass
