#! /usr/bin/env python

import argparse
import pysam
import subprocess
import os
import shutil
import sys
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
        try:
            arguments = parser.parse_args(args)
        except:
            parser.print_help()
            sys.exit(0)
        else:
            return arguments

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
        dangler_strings = extract_danglers(temp_bam_1)

        # convert dangler strings to temp fastq
        temp_fastq_danglers = os.path.join(scratch_dir, '2_danglerReads.fastq')
        print('>>> Writing dangler reads to: {0}'.format(temp_fastq_danglers))
        sam_strings_to_fastq(dangler_strings,
                             temp_fastq_danglers)

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
        print('>>> Adding tags to mapped danglers at: {0}'.format(self.output_bam))
        tag_danglers(temp_bam_2, dangler_strings, self.output_bam, self.mate_element_tag)

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


def extract_danglers(bam):
    """
    Extracts unmapped reads with a mapped pair.

    :param bam: bam file of paired end reads aligned to repeat-element reference
    :type bam: str

    :return: sam formatted strings
    :rtype: str
    """
    return pysam.samtools.view('-F', '0x800', '-F', '0x100', '-F', '8', '-f', '4', bam)


def sam_strings_to_fastq(sam_strings, output_fastq):
    """
    Extracts unmapped reads with a mapped pair and writes them to a fastq file.

    :param sam_strings: sam formatted strings
    :type sam_strings: str
    :param output_fastq: fastq file of non-mapped reads with pair mapping to repeat-element
    :type output_fastq: str
    """
    sam_attributes = iter(line.split('\t') for line in sam_strings.splitlines())
    fastq_lines = ("@{0}\n{1}\n+\n{2}\n".format(items[0], items[9], items[10]) for items in sam_attributes)
    with open(output_fastq, 'w') as f:
        for line in fastq_lines:
            f.write(line)


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


def tag_danglers(dangler_bam, dangler_strings, output_bam, tag):
    """
    Tags mapped danglers with the id of the element which their mate mapped to and writes to new bam.

    :param dangler_bam: bam file of dangler reads aligned to reference genome
    :type dangler_bam: str
    :param dangler_strings: sam formatted strings of non-mapped reads with pair mapping to repeat-element
    :type dangler_strings: str
    :param output_bam: bam file of dangler reads aligned to reference genome and tagged with mate element
    :type output_bam: str
    """
    sam_attributes = iter(line.split('\t') for line in dangler_strings.splitlines())
    read_mate_elements = {attrs[0]: attrs[2] for attrs in sam_attributes}
    danglers_untagged = pysam.Samfile(dangler_bam, "rb")
    danglers_tagged = pysam.Samfile(output_bam, "wb", template=danglers_untagged)
    for read in danglers_untagged.fetch():
        read.tags += [(tag, read_mate_elements[read.qname])]
        danglers_tagged.write(read)
    danglers_tagged.close()


if __name__ == '__main__':
    pass
