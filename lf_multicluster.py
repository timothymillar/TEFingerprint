import sys
import pysam
import argparse
from multiprocessing import Pool
from lf_splitcluster import split_cluster
from itertools import product


def parse_args(args):
    parser = argparse.ArgumentParser('Identify transposon flanking regions')
    parser.add_argument('input_bam')
    parser.add_argument('--references', type=str, nargs='*', default=[''])
    parser.add_argument('--families', type=str, nargs='*', default=[''])
    parser.add_argument('--strands', type=str, nargs='+',
                        choices=set("+-."), default=['+', '-'])
    parser.add_argument('--eps', type=int, default=100)
    parser.add_argument('--min_reads', type=int, default=5)
    parser.add_argument('--cores', type=int, default=1)
    return parser.parse_args(args)


def sam_references(input_bam):
    sam = pysam.AlignmentFile(input_bam, 'rb')
    references = sam.references
    sam.close()
    return references


def build_jobs(input_bam, references, families, strands, eps, min_reads):
    if references == ['']:
        references = sam_references(input_bam)
    else:
        pass
    return product([input_bam],
                   references,
                   families,
                   strands,
                   [eps],
                   [min_reads])


def main():
    args = parse_args(sys.argv[1:])
    jobs = build_jobs(args.input_bam,
                      args.references,
                      args.families,
                      args.strands,
                      args.eps,
                      args.min_reads)
    with Pool(args.cores) as p:
        results = p.starmap(split_cluster, jobs)
    for result in results:
        print(result)

if __name__ == '__main__':
    main()
