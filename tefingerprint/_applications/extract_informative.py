#! /usr/bin/env python3

import argparse
import pysam
import subprocess
import os
import glob
import shutil
import sys
import gzip
from tempfile import mkdtemp
import numpy


class Program(object):
    """"""
    def __init__(self, input_bam,
                 reference_fasta,
                 output_bam,
                 include_tips=False,
                 include_tails=True,
                 include_full_length_reads=True,
                 soft_clip_minimum_length=38,
                 mate_element_tag='ME',
                 use_temp_dir=False,
                 keep_temp_files=False,
                 threads=1):
        self.input_bam = input_bam
        self.reference_fasta = reference_fasta
        self.output_bam = output_bam
        self.include_tips = include_tips
        self.include_tails = include_tails
        self.include_full_length_reads = include_full_length_reads
        self.soft_clip_minimum_length = soft_clip_minimum_length
        self.mate_element_tag = mate_element_tag
        self.threads = str(threads)
        self.use_temp_dir = use_temp_dir
        self.temp_dir = None
        if self.use_temp_dir:
            self.keep_temp_files = False
        else:
            self.keep_temp_files = keep_temp_files
        self.temp_file_prefix = None

    @staticmethod
    def cli(args):
        """
        Argument parser to handle commandline inputs for the
        extract-informative program.

        :param args: List of commandline arguments for the
        extract-informative program
        :type args: []

        :return: Dictionary like object of arguments and values
        """
        parser = argparse.ArgumentParser('Extract reads whose pair was '
                                         'mapped to a repeat element, '
                                         'map them to a reference genome '
                                         'and tag them with their mates '
                                         'element name.')
        parser.add_argument('bam',
                            type=str,
                            nargs=1,
                            help='A bam file containing paired end reads '
                                 'mapped to a library of repeat elements.')
        parser.add_argument('-r', '--reference',
                            type=str,
                            nargs=1,
                            help="Reference genome in fasta format indexed "
                                 "with bwa.")
        parser.add_argument('-o', '--output',
                            type=str,
                            nargs=1,
                            help="Name of output bam file.")
        parser.set_defaults(include_tips=True)
        parser.add_argument('--exclude-tips',
                            dest='include_tips',
                            action='store_false',
                            help="Exclude soft-clipped tips and instead "
                                 "include their associated full length mate "
                                 "read. This may reduce precision but can "
                                 "increase depth of informative clusters.")
        parser.set_defaults(include_tails=True)
        parser.add_argument('--exclude-tails',
                            dest='include_tails',
                            action='store_false',
                            help="Don't include soft-clipped tails.")
        parser.set_defaults(include_full_length_reads=True)
        parser.add_argument('--exclude-full-length-reads',
                            dest='include_full_length_reads',
                            action='store_false',
                            help="Don't include full-length informative "
                                 "reads. With high quality data and "
                                 "annotations the exclusion of full-length "
                                 "'hanging' informative reads can improve "
                                 "precision of insertion sites at the "
                                 "expense of much lower cluster depth.")
        parser.add_argument('--soft-clip-minimum-length',
                            nargs=1,
                            type=int,
                            default=[38],
                            help="Minimum required length for inclusion of "
                                 "soft-clipped tips and tails.")
        parser.add_argument('--mate-element-tag',
                            type=str,
                            nargs=1,
                            default=['ME'],
                            help='Tag used in bam file to indicate the '
                                 'element type mate read.')
        parser.set_defaults(keep_temp_files=False)
        parser.add_argument('--keep-temp-files',
                            dest='keep_temp_files',
                            action='store_true',
                            help='Temporary intermediate files will '
                                 'not be deleted.')
        parser.set_defaults(use_os_temp=False)
        parser.add_argument('--use-os-temp',
                            dest='use_os_temp',
                            action='store_true',
                            help='Optional writes temp data to a temp '
                                 'file provided by the OS.')
        parser.add_argument('-t', '--threads',
                            type=int,
                            nargs=1,
                            default=[1],
                            help='Maximum number of cpu threads to be used.')
        return parser.parse_args(args)

    @staticmethod
    def from_cli(args):
        arguments = Program.cli(args)
        return Program(arguments.bam[0],
                       arguments.reference[0],
                       arguments.output[0],
                       include_tips=arguments.include_tips,
                       include_tails=arguments.include_tails,
                       include_full_length_reads=arguments.include_full_length_reads,
                       soft_clip_minimum_length=arguments.soft_clip_minimum_length[0],
                       mate_element_tag=arguments.mate_element_tag[0],
                       use_temp_dir=arguments.use_os_temp,
                       keep_temp_files=arguments.keep_temp_files,
                       threads=arguments.threads[0])

    def _run_pipeline(self):
        """"""
        if self.include_tips or self.include_tails:
            print('>>> Extracting soft clips from : {0}'.format(self.input_bam))
            informative_clips = extract_informative_clips(
                self.input_bam,
                include_soft_tips=self.include_tips,
                include_soft_tails=self.include_tails,
                minimum_soft_clip_len=self.soft_clip_minimum_length)

            print('>>> Writing soft clips to : {0}'.format(self.temp_fastq))
            # parameters for next segment
            clip_tags = capture_elements_and_write_fastq(informative_clips,
                                                         self.temp_fastq)
            # get names of original reads to skip in danglers
            # clip names have 3 characters appended so remove these
            skip_names = {name[:-3] for name in clip_tags.keys()}
            fastq_mode = 'a'  # append to file just written
        else:
            # defaults if no soft clips used
            clip_tags = dict()
            skip_names = set()
            fastq_mode = 'w'  # write new file

        if self.include_full_length_reads:
            print('>>> Extracting reads from : {0}'.format(self.input_bam))
            informative_reads = extract_informative_full_reads(self.input_bam)

            print('>>> Writing reads to : {0}'.format(self.temp_fastq))
            element_tags = capture_elements_and_write_fastq(informative_reads,
                                                            self.temp_fastq,
                                                            skip=skip_names,
                                                            mode=fastq_mode)
        else:
            element_tags = {}

        # update tags, soft-clip tags overwrite regular read tags
        element_tags.update(clip_tags)

        print('>>> Mapping reads to {0} in {1}'.format(self.reference_fasta,
                                                       self.temp_bam))
        map_danglers_to_reference(self.temp_fastq,
                                  self.reference_fasta,
                                  self.temp_bam,
                                  temp_file_prefix=self.temp_sam_sort,
                                  threads=self.threads)

        print('>>> Indexing: {0}'.format(self.temp_bam))
        pysam.index(self.temp_bam)

        print('>>> Tagging reads in: {0}'.format(self.output_bam))
        tag_danglers(self.temp_bam,
                     element_tags,
                     self.output_bam,
                     self.mate_element_tag)

        print('>>> Indexing: {0}'.format(self.output_bam))
        pysam.index(self.output_bam)

    def _cleanup(self):
        # if temp directory used then always remove files and dir
        if self.use_temp_dir:
            print('>>> Removing temporary directory: {0}'.format(self.temp_dir))
            shutil.rmtree(self.temp_dir)
        # else if keeping files requested then don't remove them
        elif self.keep_temp_files:
            temp_files = glob.glob(self.temp_file_prefix + '*')
            print('>>> Keeping temporary files: {0}'.format(temp_files))
        # else remove temp files
        else:
            temp_files = glob.glob(self.temp_file_prefix + '*')
            print('>>> Removing temporary files: {0}'.format(temp_files))
            for f in temp_files:
                os.remove(f)

    def _generate_temp_file_path(self, prefix):
        """Generate unique prefix for temp files.
        Primarily defined by output bam name.
        In addition a random number used (checked to be unique).
        This is not bullet-proof but should be good enough.
        """
        searching = True
        while searching:
            suffix = '.TEFTMP' + str(numpy.random.randint(1000000)).zfill(7)
            prefix = prefix + suffix
            if glob.glob(prefix + '*'):
                # this prefix is not unique
                pass
            else:
                # prefix is unique
                return prefix

    def run(self):
        """"""
        # check if samtools and bwa are available
        check_programs_installed('bwa', 'samtools')

        # check reference exists and is indexed
        check_reference_fasta(self.reference_fasta)

        # create temp dir for intermediate files if requested
        if self.use_temp_dir:
            self.temp_dir = mkdtemp()
            print('>>> Creating temporary directory: {0}'.format(self.temp_dir))
            # generate tempfile prefix
            self.temp_file_prefix = self._generate_temp_file_path(self.temp_dir)
            self.temp_fastq = (self.temp_file_prefix + '.0.fastq.gz')
            self.temp_bam = (self.temp_file_prefix + '.1.bam')
            self.temp_sam_sort = (self.temp_file_prefix + '.sort')

        else:
            # generate tempfile prefix
            self.temp_file_prefix = self._generate_temp_file_path(self.output_bam)
            self.temp_fastq = (self.temp_file_prefix + '.0.fastq.gz')
            self.temp_bam = (self.temp_file_prefix + '.1.bam')
            self.temp_sam_sort = (self.temp_file_prefix + '.sort')

        # attempt running pipeline
        print('>>> Running Pipeline')
        try:
            self._run_pipeline()
        except:
            # clean up and re-raise the error for debugging
            print('>>> Error Encountered, cleaning up')
            self._cleanup()
            raise

        # remove temp dir unless it was supplied by user
        self._cleanup()


def extract_informative_full_reads(bam):
    """Extracts informative reads from a bam file of paired end reads
    aligned to a library of transposons.

    :param bam: path to a bam file of paired end reads aligned to a library
        of transposons
    :type bam: str

    :return: parsed and filtered reads relating information about location
        of transposon insertions
    :rtype: generator[dict[str: str]]
    """
    alignment = pysam.AlignmentFile(bam, 'r')
    for read in alignment:
        if read.is_secondary or read.is_supplementary or read.mate_is_unmapped:
            # read is bad
            pass

        elif read.is_unmapped:
            # read is informative
            # ensure unique name
            qname = read.qname + (':R2' if read.is_read2 else ':R1')
            if read.is_reverse:
                # read is a reverse informative
                yield {'name': qname,
                       'element': read.next_reference_name,
                       'sequence': reverse_complement(read.seq),
                       'quality': read.qual[::-1]}

            else:
                # read is a forward informative
                yield {'name': qname,
                       'element': read.next_reference_name,
                       'sequence': read.seq,
                       'quality': read.qual}


def extract_informative_clips(bam,
                              include_soft_tips=True,
                              include_soft_tails=True,
                              minimum_soft_clip_len=38):
    """Extracts informative soft clipped ends of reads from a bam file
    of paired end reads aligned to a library of transposons.

    :param bam: path to a bam file of paired end reads aligned to a library
        of transposons
    :type bam: str
    :param include_soft_tips: whether or not to include soft clipped tips
    :type include_soft_tips: bool
    :param include_soft_tails: whether or not to include soft clipped tails
    :type include_soft_tails: bool
    :param minimum_soft_clip_len: minimum length of tail soft clips to include
    :type minimum_soft_clip_len: int

    :return: parsed and filtered reads relating information about location
        of transposon insertions
    :rtype: generator[dict[str: str]]
    """
    alignment = pysam.AlignmentFile(bam, 'r')
    for read in alignment:

        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            # read is bad, soft clipped reads are always mapped
            pass

        else:

            tip_lower = read.cigar[0]
            if tip_lower[0] == 4 and tip_lower[1] >= minimum_soft_clip_len:
                # soft clip at lower end of read
                clip = tip_lower[1]
                # ensure unique name
                qname = read.qname + (':R2' if read.is_read2 else ':R1')
                if read.is_reverse:
                    # reverse tip clip
                    if include_soft_tips:
                        yield {'name': qname + ':-3',
                               'element': read.reference_name,
                               'sequence': read.seq[0: clip + 1],
                               'quality': read.qual[0: clip + 1]}
                else:
                    # forward tail clip
                    if include_soft_tails:
                        yield {'name': qname + ':+5',
                               'element': read.reference_name,
                               'sequence': read.seq[0: clip + 1],
                               'quality': read.qual[0: clip + 1]}

            tip_upper = read.cigar[-1]
            if tip_upper[0] == 4 and tip_upper[1] >= minimum_soft_clip_len:
                # soft clip at upper end of read
                # ensure unique name
                qname = read.qname + (':R2' if read.is_read2 else ':R1')
                clip = tip_upper[1]
                if read.is_reverse:
                    # reverse tail clip
                    if include_soft_tails:
                        yield {'name': qname + ':-5',
                               'element': read.reference_name,
                               'sequence': reverse_complement(read.seq[-clip:]),
                               'quality': read.qual[-clip:][::-1]}
                else:
                    # forward tip clip
                    if include_soft_tips:
                        yield {'name': qname + ':+3',
                               'element': read.reference_name,
                               'sequence': reverse_complement(read.seq[-clip:]),
                               'quality': read.qual[-clip:][::-1]}


def capture_elements_and_write_fastq(reads, fastq, skip=(), mode='w'):
    """
    Captures writes parsed reads to a fastq file and returns a dictionary of
    read-element pairs.

    :param reads: sequence of dictionaries of reads with fields 'name',
        'element', sequence' and 'quality'
    :type reads: iterator[dict[str: str]]
    :param fastq: path to fastq to write
    :type fastq: str
    :param skip: set of element names to not include
    :type skip: set(str)
    :param mode: mode for opening fastq file 'w' or 'a'
    :type mode: str

    :return: dictionary of read-element pairs
    :rtype: dict[str, str]
    """
    element_tags = {}
    template = "@{0}\n{1}\n+\n{2}\n"

    with gzip.open(fastq, mode) as fq:
        for read in reads:
            # filter
            if read['name'] in skip:
                pass
            else:

                # add read name and element name to dict
                element_tags[read['name']] = read['element']

                # write read to fastq
                fq.write(template.format(read['name'],
                                         read['sequence'],
                                         read['quality']).encode('utf-8'))
    return element_tags


def reverse_complement(sequence):
    """
    Computes the reverse compliment of a nucleotide string.

    :param sequence: capitalised string of nucleotides in the set
        {'A', 'T', 'G', 'C', 'N'}
    :type sequence: str

    :return: capitalised reverse compliment string of nucleotides
    :rtype: str
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base)
                   for base in reversed(sequence.upper()))


def check_reference_fasta(path):
    """Check reference fasta path is a valid file path and that the
    `.pac` and `.sa` index files are valid file paths.

    This is a simple check to save ambiguous errors after the initial fastq
    is generated which can take some time.
    This does not check that the reference and index files are correct,
    just that they exist.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError("Reference file not found: '{0}'".format(path))
    for suffix in ['.pac', '.sa']:
        idx_path = path + suffix
        if not os.path.isfile(idx_path):
            message = "Reference index file not found: '{0}'".format(idx_path)
            raise FileNotFoundError(message)


def check_programs_installed(*args):
    """
    Check if a program exists on the users path.

    :param args: names of programs
    :type args: str

    :return: true if all of the program are on the users path
    :rtype: bool
    """
    for program in args:
        if shutil.which(program) is None:
            message="could not find '{0}' in $PATH: {1}"
            raise EnvironmentError(message.format(program,
                                                  os.environ['PATH']))
    return True


def map_danglers_to_reference(fastq,
                              reference_fasta,
                              output_bam,
                              temp_file_prefix,
                              threads=1):
    """
    Maps dangler reads (non-mapped reads with pair mapping to repeat-element)
    to reference genome and writes a bam file.

    :param fastq: fastq file of non-mapped reads with pair mapping to
        repeat-element
    :type fastq: str
    :param reference_fasta: fasta file containing reference genome
    :type reference_fasta: str
    :param output_bam: bam file of dangler reads aligned to reference
        genome
    :type output_bam: str
    :param threads: number of threads to use in bwa
    :type threads: str

    :return: subprocess code
    :rtype: int
    """
    check_programs_installed('bwa', 'samtools')
    command = ('bwa mem -t {threads} {reference} {fastq} ' +
               '| samtools view -Su - | samtools sort - -T {tmp} -o {bam}')
    command = command.format(threads=threads,
                             reference=reference_fasta,
                             fastq=fastq,
                             tmp=temp_file_prefix,
                             bam=output_bam)
    return subprocess.call([command], shell=True)


def tag_danglers(dangler_bam, mate_element_dict, output_bam, tag):
    """
    Tags mapped danglers with the id of the element which their mate mapped
    to and writes to new bam.

    :param dangler_bam: bam file of dangler reads aligned to reference genome
    :type dangler_bam: str
    :param mate_element_dict: dict of read-name, mate-element pairs
    :type mate_element_dict: dict[str, str]
    :param output_bam: bam file of dangler reads aligned to reference genome
        and tagged with mate element
    :type output_bam: str
    """
    danglers_untagged = pysam.Samfile(dangler_bam, "rb")
    danglers_tagged = pysam.Samfile(output_bam, "wb",
                                    template=danglers_untagged)
    for read in danglers_untagged.fetch():
        read.tags += [(tag, mate_element_dict[read.qname])]
        danglers_tagged.write(read)
    danglers_tagged.close()


if __name__ == '__main__':
    Program.from_cli(sys.argv)
