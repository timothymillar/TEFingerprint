#! /usr/bin/env python

import pysam
import subprocess


def map_pairs_to_repeat_elements(fastq_1, fastq_2, reference_fasta, output_bam):
    """
    Maps paired end reads to repeat-element reference and writes bam file.

    :param fastq_1: paired end reads file 1
    :type fastq_1: str
    :param fastq_2: paired end reads file 2
    :type fastq_2: str
    :param reference_fasta: repeat-element reference
    :type reference_fasta: str
    :param output_bam: bam file of paired end reads aligned to repeat-element reference
    :type output_bam: str
    """
    pass


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


def map_danglers_to_reference(fastq, output_bam):
    """
    Maps dangler reads (non-mapped reads with pair mapping to repeat-element) to reference genome and writes a bam file.

    :param fastq: fastq file of non-mapped reads with pair mapping to repeat-element
    :type fastq: str
    :param output_bam: bam file of dangler reads aligned to reference genome
    :type output_bam: str
    """
    pass


def tag_danglers(dangler_bam, dangler_strings, output_bam):
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
        read.tags += [('ME', read_mate_elements[read.qname])]
        danglers_tagged.write(read)
    danglers_tagged.close()


if __name__ == '__main__':
    pass
