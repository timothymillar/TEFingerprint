import io
import os
import sys
import argparse
import pysam
from subprocess import Popen, PIPE, STDOUT


def parse_args(args):
    parser = argparse.ArgumentParser('Identify transposon flanking regions')
    parser.add_argument('input_bam')
    parser.add_argument('-r', '--reference',
                        type=str)
    parser.add_argument('-f', '--family',
                        type=str,
                        nargs='?',
                        default='')
    parser.add_argument('-s', '--strand',
                        type=str,
                        choices=set("+-."))
    parser.add_argument('-e', '--eps',
                        type=int,
                        default=100)
    parser.add_argument('-m', '--min_reads',
                        type=int,
                        default=5)
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


def split_cluster(bamfile, reference, family, strand, eps, min_reads):
    flag = strand2flag(strand)
    sam_buffer = io.StringIO(pysam.view(flag[0], flag[1], bamfile, reference))
    sam = [line for line in sam_buffer if line.startswith(family)]
    head = [pysam.view('-H', flag[0], flag[1], bamfile, reference)]
    sam = head + sam
    bam = bytearray(''.join(sam), 'utf-8')
    cluster_script = os.path.dirname(os.path.realpath(__file__)) + '/lf_cluster.py'
    proc = Popen(['python',
                  cluster_script,
                  '--reference', reference,
                  '--family', family,
                  '--strand', strand,
                  '--eps', str(eps),
                  '--min_reads', str(min_reads)],
                 stdout=PIPE,
                 stdin=PIPE,
                 stderr=STDOUT)
    result = proc.communicate(input=bam)[0]
    return result.decode()


def main():
    args = parse_args(sys.argv[1:])
    print(split_cluster(args.input_bam,
                        args.reference,
                        args.family,
                        args.strand,
                        args.eps,
                        args.min_reads))

if __name__ == '__main__':
    main()
