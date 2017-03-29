#! /usr/bin/env python

from itertools import product
from multiprocessing import Pool
from tectoolkit.loci import *


def _fingerprint_worker(bams, reference, categories, mate_element_tag, min_reads, eps, min_eps, hierarchical):
    reads = merge(*[ReadLoci.from_bam(bam, reference, categories, tag=mate_element_tag) for bam in bams])
    return reads.fingerprint(min_reads, eps, min_eps=min_eps, hierarchical=hierarchical)


def fingerprint(bams=None,
                references=None,
                categories=None,
                mate_element_tag='ME',
                min_reads=None,
                eps=None,
                min_eps=0,
                hierarchical=True,
                cores=1):
    jobs = product([bams], references, [categories], [mate_element_tag], [min_reads], [eps], [min_eps], [hierarchical])

    if cores == 1:
        return merge(*[_fingerprint_worker(*job) for job in jobs])
    else:
        with Pool(cores) as pool:
            return merge(*pool.starmap(_fingerprint_worker, jobs))


def _comparison_worker(bams, reference, categories, mate_element_tag, min_reads, eps, min_eps, hierarchical):
    reads = merge(*[ReadLoci.from_bam(bam, reference, categories, tag=mate_element_tag) for bam in bams])
    fprint = reads.fingerprint(min_reads, eps, min_eps=min_eps, hierarchical=hierarchical)
    return ComparativeBins.from_union(fprint).compare(reads)


def comparison(bams=None,
               references=None,
               categories=None,
               mate_element_tag='ME',
               min_reads=None,
               eps=None,
               min_eps=0,
               hierarchical=True,
               cores=1):
    jobs = product([bams], references, [categories], [mate_element_tag], [min_reads], [eps], [min_eps], [hierarchical])

    if cores == 1:
        return merge(*[_comparison_worker(*job) for job in jobs])
    else:
        with Pool(cores) as pool:
            return merge(*pool.starmap(_comparison_worker, jobs))