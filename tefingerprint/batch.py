#! /usr/bin/env python

from itertools import product
from multiprocessing import Pool
from tefingerprint import loci


def _fingerprint_worker(bams, categories, reference, transposon_tag, min_reads, eps, min_eps, hierarchical):
    """Worker function for batch fingerprinting"""
    reads = loci.append(*[loci.ReadLoci.from_bam(bam, categories, references=reference, tag=transposon_tag)
                          for bam in bams])
    return reads.fingerprint(min_reads, eps, min_eps=min_eps, hierarchical=hierarchical)


def fingerprint(bams=None,
                categories=None,
                references=None,
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
    :rtype: :class:`loci.FingerPrint`
    """
    jobs = product([bams], [categories], references, [transposon_tag], [min_reads], [eps], [min_eps], [hierarchical])

    if cores == 1:
        return loci.append(*[_fingerprint_worker(*job) for job in jobs])
    else:
        with Pool(cores) as pool:
            return loci.append(*pool.starmap(_fingerprint_worker, jobs))


def _comparison_worker(bams, categories, reference, transposon_tag, min_reads, eps, min_eps, hierarchical, fingerprint_buffer, bin_buffer):
    """Worker function for batch fingerprint comparisons"""
    reads = loci.append(*[loci.ReadLoci.from_bam(bam, categories, references=reference, tag=transposon_tag)
                          for bam in bams])
    fprint = reads.fingerprint(min_reads, eps, min_eps=min_eps, hierarchical=hierarchical)
    fprint.buffer(fingerprint_buffer)
    bins = loci.ComparativeBins.from_fingerprints(fprint)
    bins.buffer(bin_buffer)
    return bins.compare(reads)


def comparison(bams=None,
               categories=None,
               references=None,
               transposon_tag='ME',
               min_reads=None,
               eps=None,
               min_eps=0,
               hierarchical=True,
               fingerprint_buffer=0,
               bin_buffer=0,
               cores=1):
    """
    Multi-processes pipeline for batch comparison of multiple bam file(s) by transposon fingerprints.

    :param bams: a list of bam file(s) to be fingerprinted
    :type bams: list[str]
    :param categories: targeted transposon categories
    :type categories: list[str]
    :param references: targeted (slices of) references with format 'name' or 'name:minimum-maximum'
    :type references: list[str]
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
    :param fingerprint_buffer: value to buffer fingerprints by, defaults to 0
    :type fingerprint_buffer: int
    :param bin_buffer: value to buffer comparative-bins by, defaults to 0
    :type bin_buffer: int
    :param cores: number of processes to use for the batch job
    :type cores: int

    :return: Fingerprint comparisons of bam files
    :rtype: :class:`loci.Comparison`
    """
    jobs = product([bams], [categories], references, [transposon_tag], [min_reads], [eps], [min_eps], [hierarchical], [fingerprint_buffer],  [bin_buffer])

    if cores == 1:
        return loci.append(*[_comparison_worker(*job) for job in jobs])
    else:
        with Pool(cores) as pool:
            return loci.append(*pool.starmap(_comparison_worker, jobs))


if __name__ == '__main__':
    pass
