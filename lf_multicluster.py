import sys
import pysam
import argparse
from multiprocessing import Pool
from lf_splitcluster import split_cluster
from itertools import product


def parse_args(args):
    parser = argparse.ArgumentParser('Identify transposon flanking regions')
    parser.add_argument('input_bam')
    parser.add_argument('--references', nargs='?', default=False)
    parser.add_argument('--read_groups', nargs='?', default=False)
    parser.add_argument('--strands', nargs='?', default=True)
    parser.add_argument('--eps', type=int, default=100)
    parser.add_argument('--min_tips', type=int, default=5)
    parser.add_argument('--cores', type=int, default=1)
    return parser.parse_args(args)


def sam_references(input_bam):
    sam = pysam.AlignmentFile(input_bam, 'rb')
    references = sam.references
    sam.close()
    return references


def build_jobs(input_bam, references, read_groups, strands, eps, min_tips):
    if references:
        if isinstance(references, str):
            references = [references]
        else:
            pass
    else:
        references = sam_references(input_bam)
    if read_groups:
        if isinstance(read_groups, str):
            read_groups = [read_groups]
        else:
            pass
    else:
        read_groups = [False]
    if strands is True:
        strands = ['+', '-']
    elif strands:
        if isinstance(strands, str):
            strands = [strands]
        else:
            pass
    else:
        strands = [False]
    return product([input_bam],
                   references,
                   read_groups,
                   strands,
                   [eps],
                   [min_tips])


def main():
    args = parse_args(sys.argv[1:])
    jobs = build_jobs(args.input_bam,
                      args.references,
                      args.read_groups,
                      args.strands,
                      args.eps,
                      args.min_tips)
    with Pool(args.cores) as p:
        results = p.starmap(split_cluster, jobs)
    for result in results:
        print(result)

if __name__ == '__main__':
    main()
