import sys
import argparse
import pysam
from subprocess import Popen, PIPE, STDOUT


def parse_args(args):
    parser = argparse.ArgumentParser('Identify transposon flanking regions')
    parser.add_argument('input_bam')
    parser.add_argument('--reference', type=str)
    parser.add_argument('--read_group', type=str, nargs='?', default='')
    parser.add_argument('--strand', type=str, choices=set("+-."))
    parser.add_argument('--eps', type=int, default=100)
    parser.add_argument('--min_tips', type=int, default=5)
    return parser.parse_args(args)


def strand2flag(strand):
    if strand == '+':
        return ('-F', '20')
    elif strand == '-':
        return ('-f', '16')
    elif strand == '.':
        return ('-F', '4')
    else:
        pass  # throw error


def split_cluster(bamfile, reference, read_group, strand, eps, min_tips):
    flag = strand2flag(strand)
    sam = pysam.view(flag[0], flag[1], bamfile, reference)
    if strand:
        sam = [line for line in sam if line.startswith(read_group)]
    head = pysam.view('-H', flag[0], flag[1], bamfile, reference)
    sam = head + sam
    bam = bytearray(''.join(sam), 'utf-8')
    proc = Popen(['python',
                  'lf_cluster.py',
                  '--reference', reference,
                  '--read_group', read_group,
                  '--strand', strand,
                  '--eps', str(eps),
                  '--min_tips', str(min_tips)],
                 stdout=PIPE,
                 stdin=PIPE,
                 stderr=STDOUT)
    result = proc.communicate(input=bam)[0]
    return result.decode()


def main():
    args = parse_args(sys.argv[1:])
    print(split_cluster(args.input_bam,
                        args.reference,
                        args.read_group,
                        args.strand,
                        args.eps,
                        args.min_tips))

if __name__ == '__main__':
    main()
