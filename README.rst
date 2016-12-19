TEFingerprint
=============

A Python library with CLI for fingerprinting transposable elements based
on paired-end reads. This project is currently in a pe-alpha state and
rapidly changing.

Development Environment
-----------------------

A development environment can be built using the conda tool available in
the `Anaconda Python
distribution <https://www.continuum.io/downloads>`__. An `environment
file <http://conda.pydata.org/docs/using/envs.html#share-an-environment>`__
"environment.yml" is available in the root directory of this project.

::

    conda env create -f environment.yml

This will create a `virtual
environment <http://conda.pydata.org/docs/using/envs.html>`__ called
'tectoolkit'.

Dependencies
------------

The current development version supports Python 3.5 or higher and pysam
0.9 or higher.

Building and Installation
-------------------------

-  Activate the development environment:

::

    source activate tectoolkit

-  Uninstall any previous version:

::

    pip uninstall tectoolkit

-  Navigate to the projects root directory build and install:

   -  (note the ``*`` wildcard should be replaced with the current
      version number to install)

::

    python setup.py sdist
    pip install dist/tectoolkit-*.tar.gz

-  Alternatively, The project can be installed in editable mode from the
   projects root directory:

::

    pip install -e ./

Testing
-------

Unit testing using `pytest <http://doc.pytest.org/en/latest/>`__. To run
tests from the projects root directory:

::

    python -m pytest -v ./

Usage
-----

Currently the required pre-processing of paired end reads is not
supported natively by the tool (see the `Leapfrog
Readme <https://github.com/mfiers/leapfrog/blob/master/README>`__ for
details). The basic 'tec' tool is used as a wrapper for included
programs like ``fingerprint``, ``compare`` and ``filter_gff``. If no
aruments are given the help text will be displayed:

::

     tec
    usage: Identify program to run [-h] {fingerprint,compare,filter_gff}
    Identify program to run: error: the following arguments are required: program
    usage: Identify program to run [-h] {fingerprint,compare,filter_gff}

    positional arguments:
      {fingerprint,compare,filter_gff}

    optional arguments:
      -h, --help            show this help message and exit

Fingerprint
~~~~~~~~~~~

::

    $ tec fingerprint
    usage: Identify potential TE flanking regions [-h]
                                                  [-r [REFERENCES [REFERENCES ...]]]
                                                  [-f [FAMILIES [FAMILIES ...]]]
                                                  [-s {-,+} [{-,+} ...]]
                                                  [-m MIN_READS]
                                                  [-e EPSILON [EPSILON ...]]
                                                  [-t THREADS]
                                                  input_bam
    Identify potential TE flanking regions: error: the following arguments are required: input_bam
    usage: Identify potential TE flanking regions [-h]
                                                  [-r [REFERENCES [REFERENCES ...]]]
                                                  [-f [FAMILIES [FAMILIES ...]]]
                                                  [-s {-,+} [{-,+} ...]]
                                                  [-m MIN_READS]
                                                  [-e EPSILON [EPSILON ...]]
                                                  [-t THREADS]
                                                  input_bam

    positional arguments:
      input_bam             A single bam file to be fingerprinted

    optional arguments:
      -h, --help            show this help message and exit
      -r [REFERENCES [REFERENCES ...]], --references [REFERENCES [REFERENCES ...]]
                            The reference sequence(s) (e.g. chromosome) to be
                            fingerprinted. If left blank all references sequences
                            in the input file will be used.
      -f [FAMILIES [FAMILIES ...]], --families [FAMILIES [FAMILIES ...]]
                            TE grouping(s) to be used. These must be exact string
                            match's to start of read name and are used to split
                            reads into categories for analysis
      -s {-,+} [{-,+} ...], --strands {-,+} [{-,+} ...]
                            Strand(s) to be analysed. Use + for forward or - for
                            reverse. Default is to analyse both strands
                            (separately).
      -m MIN_READS, --min_reads MIN_READS
                            Minimum number of read tips required to be considered
                            a cluster. This values is used in combination with
                            epsilon to describe the density of read tips that is
                            required for identification of a clusters. For every
                            set of <min_reads> reads tips, if those reads are
                            within epsilon range of one another, they are
                            classified as a subcluster. Overlapping sets of
                            subclusters are then merged to form clusters.
      -e EPSILON [EPSILON ...], --epsilon EPSILON [EPSILON ...]
                            Epsilon is the maximum allowable distance among a set
                            of read tips to be considered a (sub)cluster. If a
                            single value is given, the UDC algorithm will be used
                            to identify all clusters at the specified density
                            (defined by epsilon and min_points). If two values are
                            given, they will be interpreted as maximum and minimum
                            epsilon values using the Hierarchical HUDC
                            algorithm.The maximum (or only) epsilon value given
                            should be larger than the insert size, and the minimum
                            epsilon (if used) should be much smaller (often zero)
                            in order to find adequate support for child clusters.
                            HUDC identifies all clusters at the maximum specified
                            density and then attempts to split them into logical
                            child clusters at all values of epsilon between
                            maximum and minimum. The robustness of each parent
                            cluster is compared to it's children. If the parent is
                            more robust it is selected, otherwise the process is
                            repeated for child cluster recursively until a parent
                            or terminal (cluster with no children) is selected.
      -t THREADS, --threads THREADS
                            Maximum number of cpu threads to be used

Compare
~~~~~~~

::

    $ tec compare
    usage: Compare potential TE flanking regions [-h]
                                                 [-r [REFERENCES [REFERENCES ...]]]
                                                 [-f [FAMILIES [FAMILIES ...]]]
                                                 [-s {-,+} [{-,+} ...]]
                                                 [-m MIN_READS]
                                                 [-e EPSILON [EPSILON ...]]
                                                 [-b BIN_BUFFER] [-t THREADS]
                                                 input_bams [input_bams ...]
    Compare potential TE flanking regions: error: the following arguments are required: input_bams
    usage: Compare potential TE flanking regions [-h]
                                                 [-r [REFERENCES [REFERENCES ...]]]
                                                 [-f [FAMILIES [FAMILIES ...]]]
                                                 [-s {-,+} [{-,+} ...]]
                                                 [-m MIN_READS]
                                                 [-e EPSILON [EPSILON ...]]
                                                 [-b BIN_BUFFER] [-t THREADS]
                                                 input_bams [input_bams ...]

    positional arguments:
      input_bams            A list of two or more bam files to be compared

    optional arguments:
      -h, --help            show this help message and exit
      -r [REFERENCES [REFERENCES ...]], --references [REFERENCES [REFERENCES ...]]
                            The reference sequence(s) (e.g. chromosome) to be
                            fingerprinted. If left blank all references sequences
                            in the input file will be used.
      -f [FAMILIES [FAMILIES ...]], --families [FAMILIES [FAMILIES ...]]
                            TE grouping(s) to be used. These must be exact string
                            match's to start of read name and are used to split
                            reads into categories for analysis
      -s {-,+} [{-,+} ...], --strands {-,+} [{-,+} ...]
                            Strand(s) to be analysed. Use + for forward or - for
                            reverse. Default is to analyse both strands
                            (separately).
      -m MIN_READS, --min_reads MIN_READS
                            Minimum number of read tips required to be considered
                            a cluster. This values is used in combination with
                            epsilon to describe the density of read tips that is
                            required for identification of a clusters. For every
                            set of <min_reads> reads tips, if those reads are
                            within epsilon range of one another, they are
                            classified as a subcluster. Overlapping sets of
                            subclusters are then merged to form clusters.
      -e EPSILON [EPSILON ...], --epsilon EPSILON [EPSILON ...]
                            Epsilon is the maximum allowable distance among a set
                            of read tips to be considered a (sub)cluster. If a
                            single value is given, the UDC algorithm will be used
                            to identify all clusters at the specified density
                            (defined by epsilon and min_points). If two values are
                            given, they will be interpreted as maximum and minimum
                            epsilon values using the Hierarchical HUDC
                            algorithm.The maximum (or only) epsilon value given
                            should be larger than the insert size, and the minimum
                            epsilon (if used) should be much smaller (often zero)
                            in order to find adequate support for child clusters.
                            HUDC identifies all clusters at the maximum specified
                            density and then attempts to split them into logical
                            child clusters at all values of epsilon between
                            maximum and minimum. The robustness of each parent
                            cluster is compared to it's children. If the parent is
                            more robust it is selected, otherwise the process is
                            repeated for child cluster recursively until a parent
                            or terminal (cluster with no children) is selected.
      -b BIN_BUFFER, --bin_buffer BIN_BUFFER
                            Additional buffer to be added to margins of
                            comparative bins. This is used avoid identifying small
                            clusters as unique, when these is only slight miss-
                            match in read positions across samples (i.e. false
                            positives). A value of 20-50 should be sufficient in
                            most cases
      -t THREADS, --threads THREADS
                            Maximum number of cpu threads to be used

Filter GFF
~~~~~~~~~~

::

    $ tec filter_gff
    usage: Identify potential TE flanking regions [-h] [-f FILTERS [FILTERS ...]]
                                                  input_gff
    Identify potential TE flanking regions: error: the following arguments are required: input_gff
    usage: Identify potential TE flanking regions [-h] [-f FILTERS [FILTERS ...]]
                                                  input_gff

    positional arguments:
      input_gff             A single gff file to be filtered

    optional arguments:
      -h, --help            show this help message and exit
      -f FILTERS [FILTERS ...], --filters FILTERS [FILTERS ...]
                            List of filters to apply. A valid filter takes the
                            form '<attribute><operator><value>'where <attribute>
                            is the name of a GFF attribute, <operator> is one of
                            '=', '==', '!=', '>=', '<=', '>' or '<' and the value
                            of the GFF attribute is compared to <value> using the
                            operator The list of filters is applied additively
                            (i.e. a feature must meet all filters) and, if a
                            feature is selected, all of it's ancestors and
                            descendants will also be included in the output.
                            Operators '=', '==' and '!=' will attempt to compare
                            values as floating point numbers if possible and
                            otherwise compare values as strings. Operators '>=',
                            '<=', '>' and '<' will coerce values to floating point
                            numbers before comparison.
