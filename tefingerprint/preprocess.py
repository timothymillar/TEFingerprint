#! /usr/bin/env python

import argparse
import pysam
import subprocess
import os
import shutil
import sys
import re
from tempfile import mkdtemp


class PreProcessProgram(object):
    """"""
    def __init__(self, fastq_1, fastq_2,
                 reference_fasta,
                 repeats_fasta,
                 output_bam,
                 mate_element_tag='ME',
                 temp_dir=False, threads=1):
        self.fastq_1 = fastq_1
        self.fastq_2 = fastq_2
        self.reference_fasta = reference_fasta
        self.repeats_fasta = repeats_fasta
        self.output_bam = output_bam
        self.mate_element_tag = mate_element_tag
        self.temp_dir = temp_dir
        self.threads = str(threads)

    @staticmethod
    def _cli(args):
        """
        Argument parser to handle commandline inputs for the preprocessing program.

        :param args: List of commandline arguments for the preprocess program
        :type args: []

        :return: Dictionary like object of arguments and values
        """
        parser = argparse.ArgumentParser('Identify potential TE flanking regions')
        parser.add_argument('reads',
                            type=str,
                            nargs=2,
                            help='A pair of fastq files containing paired end reads')
        parser.add_argument('--reference',
                            type=str,
                            nargs=1,
                            help="Reference genome in fasta format with bwa index")
        parser.add_argument('--repeats',
                            type=str,
                            nargs=1,
                            help="Library of repeat elements in fasta format with bwa index")
        parser.add_argument('-o', '--output',
                            type=str,
                            nargs=1,
                            help="Name of output bam file")
        parser.add_argument('--mate_element_tag',
                            type=str,
                            nargs=1,
                            default=['ME'],
                            help='Tag used in bam file to indicate the element type mate read')
        parser.add_argument('--tempdir',
                            type=str,
                            nargs=1,
                            default=[False],
                            help="Optional directory to store temp files in")
        parser.add_argument('-t', '--threads',
                            type=int,
                            nargs=1,
                            default=[1],
                            help='Maximum number of cpu threads to be used')
        return parser.parse_args(args)

    @staticmethod
    def from_cli(args):
        arguments = PreProcessProgram._cli(args)
        return PreProcessProgram(arguments.reads[0],
                                 arguments.reads[1],
                                 arguments.reference[0],
                                 arguments.repeats[0],
                                 arguments.output[0],
                                 mate_element_tag=arguments.mate_element_tag[0],
                                 temp_dir=arguments.tempdir[0],
                                 threads=arguments.threads[0])

    def _run_pipeline(self, scratch_dir):
        """"""
        # map reads to repeats and store as temp file
        temp_bam_1 = os.path.join(scratch_dir, '1_pairedReadsMappedToRepeats.bam')
        print('>>> Mapping paired reads to repeat library at: {0}'.format(temp_bam_1))
        map_pairs_to_repeat_elements(self.fastq_1,
                                     self.fastq_2,
                                     self.repeats_fasta,
                                     temp_bam_1,
                                     self.threads)

        # index bam
        print('>>> Indexing: {0}'.format(temp_bam_1))
        index_bam(temp_bam_1)

        # extract danglers from bam
        print('>>> Extracting dangler reads from: {0}'.format(temp_bam_1))
        dangler_strings = extract_forward_danglers(temp_bam_1)

        # convert dangler strings to fastq
        fastq_lines = '\n'.join((sam_strings_as_fastq(dangler_strings)))

        # extract soft clipped tails from bam
        print('>>> Extracting soft clipped tails from: {0}'.format(temp_bam_1))
        soft_clipped_forward_tails = extract_forward_soft_clipped_tails(temp_bam_1)
        soft_clipped_reverse_tails = extract_reverse_soft_clipped_tails(temp_bam_1)

        # add soft clipped tail strings to fastq
        fastq_lines += '\n' + '\n'.join(forward_soft_clipped_tails_as_fastq(soft_clipped_forward_tails))
        fastq_lines += '\n' + '\n'.join(reverse_soft_clipped_tails_as_fastq(soft_clipped_reverse_tails))

        # write temp fastq
        temp_fastq_danglers = os.path.join(scratch_dir, '2_danglerReads.fastq')
        print('>>> Writing dangler reads to: {0}'.format(temp_fastq_danglers))
        with open(temp_fastq_danglers, 'w') as f:
            for line in fastq_lines:
                f.write(line)

        # map danglers to reference
        temp_bam_2 = os.path.join(scratch_dir, '3_danglerReadsMappedToReference.bam')
        print('>>> Mapping dangler reads to reference at: {0}'.format(temp_bam_2))
        map_danglers_to_reference(temp_fastq_danglers,
                                  self.reference_fasta,
                                  temp_bam_2,
                                  self.threads)

        # index bam
        print('>>> Indexing: {0}'.format(temp_bam_2))
        index_bam(temp_bam_2)

        # tag danglers and write output file
        sam_strings = dangler_strings + soft_clipped_forward_tails + soft_clipped_reverse_tails
        print('>>> Adding tags to mapped danglers at: {0}'.format(self.output_bam))
        tag_danglers(temp_bam_2, sam_strings, self.output_bam, self.mate_element_tag)

        # index bam
        print('>>> Indexing: {0}'.format(self.output_bam))
        index_bam(self.output_bam)

    def run(self):
        """"""
        # check if samtools and bwa are available
        _check_programs_installed('bwa', 'samtools')

        # create temp dir for intermediate files unless one is suplied by user
        if self.temp_dir:
            temp_dir = self.temp_dir
        else:
            temp_dir = mkdtemp()
        print('>>> Creating temporary directory at: {0}'.format(temp_dir))

        # attempt running pipeline
        try:
            self._run_pipeline(temp_dir)
        except:
            # remove temp dir unless it was supplied by user
            if self.temp_dir:
                pass
            else:
                print('>>> Error Encountered... Removing temporary directory at: {0}'.format(temp_dir))
                shutil.rmtree(temp_dir)
            # re-raise the error for debugging
            raise

        # remove temp dir unless it was supplied by user
        if self.temp_dir:
            pass
        else:
            print('>>> Removing temporary directory at: {0}'.format(temp_dir))
            shutil.rmtree(temp_dir)


def _check_programs_installed(*args):
    """
    Check if a program exists on the users path.

    :param args: names of programs
    :type args: str

    :return: true if all of the program are on the users path
    :rtype: bool
    """
    for program in args:
        if shutil.which(program) is None:
            raise EnvironmentError("could not find '{0}' in $PATH: {1}".format(program, os.environ['PATH']))
    return True


def index_fasta(fasta):
    """
    Create index files for bwa mem

    :param fasta: fasta file to index
    :type fasta: str

    :return: subprocess code
    :rtype: int
    """
    _check_programs_installed('bwa')
    return subprocess.call(['bwa index -a is ' + fasta], shell=True)


def index_bam(bam):
    """
    Create index files for bwa mem

    :param bam: bam file to index
    :type bam: str

    :return: subprocess code
    :rtype: int
    """
    _check_programs_installed('samtools')
    return subprocess.call(['samtools index ' + bam], shell=True)


def map_pairs_to_repeat_elements(fastq_1, fastq_2, repeats_fasta, output_bam, threads=1):
    """
    Maps paired end reads to repeat-element reference and writes bam file.

    :param fastq_1: paired end reads file 1
    :type fastq_1: str
    :param fastq_2: paired end reads file 2
    :type fastq_2: str
    :param repeats_fasta: repeat-element reference
    :type repeats_fasta: str
    :param output_bam: bam file of paired end reads aligned to repeat-element reference
    :type output_bam: str
    :param threads: number of threads to use in bwa
    :type threads: str

    :return: subprocess code
    :rtype: int
    """
    _check_programs_installed('bwa', 'samtools')
    command = ' '.join(['bwa', 'mem', '-t', threads, repeats_fasta, fastq_1, fastq_2,
                        '| samtools view -Su - | samtools sort - -o', output_bam])
    return subprocess.call([command], shell=True)


def reverse_complement(sequence):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement.get(base, base) for base in reversed(sequence))


def parse_sam_string(string):
    attributes = string.split("\t")
    return {'QNAME': attributes[0],
            'FLAG': attributes[1],
            'RNAME': attributes[2],
            'POS': attributes[3],
            'MAPQ': attributes[4],
            'CIGAR': attributes[5],
            'RNEXT': attributes[6],
            'PNEXT': attributes[7],
            'TLEN': attributes[8],
            'SEQ': attributes[9],
            'QUAL': attributes[10]}


def extract_forward_soft_clipped_tails(bam):
    """

    :param bam:
    :return:
    """
    return pysam.view('-f', '2', '-F', '16', '-F', '4', '-F', '8', bam).splitlines()


def forward_soft_clipped_tails_as_fastq(sam_strings, min_length=20):
    reads = [parse_sam_string(read) for read in sam_strings]
    # extract soft clipped tails
    clips = [re.findall(r"^[0-9]+[S]", read['CIGAR']) for read in reads]
    # subset and pair reads with soft clipped tails
    reads = [case for case in zip(reads, clips) if case[1]]
    # clean up clip lengths
    reads = [(read, int(clip[0].strip('S'))) for read, clip in reads]
    # filter out short clips
    reads = [(read, clip) for read, clip in reads if clip >= min_length]

    def format_fastq(read, clip):
        return '@{0}\n{1}\n+\n{2}'.format(read['QNAME'],
                                          read['SEQ'][0: clip + 1],
                                          read['QUAL'][0: clip + 1])

    return (format_fastq(read, clip) for read, clip in reads)


def extract_reverse_soft_clipped_tails(bam):
    """

    :param bam:
    :return:
    """
    return pysam.view('-f', '2', '-f', '16', '-F', '4', '-F', '8', bam).splitlines()


def reverse_soft_clipped_tails_as_fastq(sam_strings, min_length=20):
    reads = [parse_sam_string(read) for read in sam_strings]
    # extract soft clipped tails
    clips = [re.findall(r"[0-9]+[S]$", read['CIGAR']) for read in reads]
    # subset and pair reads with soft clipped tails
    reads = [case for case in zip(reads, clips) if case[1]]
    # clean up clip lengths
    reads = [(read, int(clip[0].strip('S'))) for read, clip in reads]
    # filter out short clips
    reads = [(read, clip) for read, clip in reads if clip >= min_length]

    def format_fastq(read, clip):
        return '@{0}\n{1}\n+\n{2}'.format(read['QNAME'],
                                          reverse_complement(read['SEQ'][-clip:]),
                                          read['QUAL'][-clip:])

    return (format_fastq(read, clip) for read, clip in reads)


def extract_forward_danglers(bam):
    """
    Extracts unmapped reads with a mapped pair.

    :param bam: bam file of paired end reads aligned to repeat-element reference
    :type bam: str

    :return: sam formatted strings
    :rtype: str
    """
    return pysam.samtools.view('-F', '0x800', '-F', '0x100', '-F', '8', '-f', '4', bam).splitlines()


def sam_strings_as_fastq(sam_strings):
    """
    Extracts unmapped reads with a mapped pair and writes them to a fastq file.

    :param sam_strings: sam formatted strings
    :type sam_strings: str
    :param output_fastq: fastq file of non-mapped reads with pair mapping to repeat-element
    :type output_fastq: str
    """
    reads = iter(parse_sam_string(string) for string in sam_strings)
    fastq_lines = ("@{0}\n{1}\n+\n{2}".format(read['QNAME'],
                                              read['SEQ'],
                                              read['QUAL']) for read in reads)
    return fastq_lines


def map_danglers_to_reference(fastq, reference_fasta, output_bam, threads):
    """
    Maps dangler reads (non-mapped reads with pair mapping to repeat-element) to reference genome and writes a bam file.

    :param fastq: fastq file of non-mapped reads with pair mapping to repeat-element
    :type fastq: str
    :param reference_fasta: fasta file containing reference genome
    :type reference_fasta: str
    :param output_bam: bam file of dangler reads aligned to reference genome
    :type output_bam: str
    :param threads: number of threads to use in bwa
    :type threads: str

    :return: subprocess code
    :rtype: int
    """
    _check_programs_installed('bwa', 'samtools')
    command = ' '.join(['bwa', 'mem', '-t', threads, reference_fasta, fastq,
                        '| samtools view -Su - | samtools sort - -o', output_bam])
    return subprocess.call([command], shell=True)


def tag_danglers(dangler_bam, sam_strings, output_bam, tag):
    """
    Tags mapped danglers with the id of the element which their mate mapped to and writes to new bam.

    :param dangler_bam: bam file of dangler reads aligned to reference genome
    :type dangler_bam: str
    :param sam_strings: sam formatted strings of non-mapped reads with pair mapping to repeat-element
    :type sam_strings: str
    :param output_bam: bam file of dangler reads aligned to reference genome and tagged with mate element
    :type output_bam: str
    """
    sam_attributes = iter(line.split('\t') for line in sam_strings)
    read_mate_elements = {attrs[0]: attrs[2] for attrs in sam_attributes}
    danglers_untagged = pysam.Samfile(dangler_bam, "rb")
    danglers_tagged = pysam.Samfile(output_bam, "wb", template=danglers_untagged)
    for read in danglers_untagged.fetch():
        read.tags += [(tag, read_mate_elements[read.qname])]
        danglers_tagged.write(read)
    danglers_tagged.close()


if __name__ == '__main__':
    PreProcessProgram.from_cli(sys.argv)
