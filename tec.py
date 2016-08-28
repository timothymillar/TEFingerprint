import sys
import pysam
import argparse
from multiprocessing import Pool
from itertools import product


def parse_program_arg(arg):
    parser = argparse.ArgumentParser('Identify program to run')
    parser.add_argument('program',
                        type=str,
                        choices=("fingerprint", "compare"))


def parse_fingerprint_args(args):
    parser = argparse.ArgumentParser('Identify transposon flanking regions')
    parser.add_argument('input_bam',
                        nargs='1')
    parser.add_argument('-r', '--references',
                        type=str,
                        nargs='*',
                        default=[''])
    parser.add_argument('-f', '--families',
                        type=str,
                        nargs='*',
                        default=[''])
    parser.add_argument('-s', '--strands',
                        type=str,
                        nargs='+',
                        choices=set("+-."),
                        default=['+', '-'])
    parser.add_argument('-e', '--eps',
                        type=int,
                        default=100,
                        nargs='1')
    parser.add_argument('-m', '--min_reads',
                        type=int,
                        default=5,
                        nargs='1')
    parser.add_argument('-c', '--cores',
                        type=int,
                        default=1)
    return parser.parse_args(args)


def parse_compare_args(args):
    parser = argparse.ArgumentParser('Identify transposon flanking regions')
    parser.add_argument('input_bams',
                        nargs='+')
    parser.add_argument('-r', '--references',
                        type=str,
                        nargs='*',
                        default=[''])
    parser.add_argument('-f', '--families',
                        type=str,
                        nargs='*',
                        default=[''])
    parser.add_argument('-s', '--strands',
                        type=str,
                        nargs='+',
                        choices=set("+-."),
                        default=['+', '-'])
    parser.add_argument('-e', '--eps',
                        type=int,
                        default=100,
                        nargs='1')
    parser.add_argument('-m', '--min_reads',
                        type=int,
                        default=5,
                        nargs='1')
    parser.add_argument('-c', '--cores',
                        type=int,
                        default=1)
    return parser.parse_args(args)


def bam_references(input_bam):
    bam = pysam.AlignmentFile(input_bam, 'rb')
    references = bam.references
    bam.close()
    return references


def build_fingerprint_jobs(input_bam, references, families, strands, eps, min_reads):
    if references == ['']:
        references = bam_references(input_bam)
    else:
        pass
    return product(input_bam,
                   references,
                   families,
                   strands,
                   eps,
                   min_reads)


def build_compare_jobs(input_bams, references, families, strands, eps, min_reads):
    if references == ['']:
        bam_refs = [set(bam_references(bam)) for bam in input_bams]
        references = bam_refs[0]
        for refs in bam_refs[1:]:
            references = references.intersection(refs)
        references = list(references)
    else:
        pass
    return product(input_bams,
                   references,
                   families,
                   strands,
                   eps,
                   min_reads)


def main():
    program = parse_program_arg(sys.argv[1]).program
    if program == "fingerprint":
        args = parse_fingerprint_args(sys.argv[2:])
        jobs = build_fingerprint_jobs(args.input_bam,
                                      args.references,
                                      args.families,
                                      args.strands,
                                      args.eps,
                                      args.min_reads)
        with Pool(args.cores) as p:
            results = p.starmap(fingerprint, jobs)
    if program == "compare":
        args = parse_fingerprint_args(sys.argv[2:])
        jobs = build_compare_jobs(args.input_bam,
                                  args.references,
                                  args.families,
                                  args.strands,
                                  args.eps,
                                  args.min_reads)
        with Pool(args.cores) as p:
            results = p.starmap(compare, jobs)
    for result in results:
        print(result)








