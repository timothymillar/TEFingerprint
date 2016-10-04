# TEFingerprint
A toolkit for fingerprinting transposable elements based on paired-end reads.
This project is currently in a pe-alpha state and rapidly changing.


## Development Environment
A development environment can be built using the conda tool available in the [Anaconda Python distribution](https://www.continuum.io/downloads).
An [environment file](http://conda.pydata.org/docs/using/envs.html#share-an-environment) "environment.yml" is available in the root directory of this project. 

```
conda env create -f environment.yml
```

This will create a [virtual environment](http://conda.pydata.org/docs/using/envs.html) called 'tectoolkit'.


## Dependencies

The current development version supports Python 3.5 or higher and pysam 0.9 or higher.


## Building and Installation

Activate the development environment:

```
source activate tectoolkit
```

Uninstall any previous version:

```
pip uninstall tectoolkit
```

Navigate to the projects root directory build and install:

```
python setup.py sdist
pip install dist/tectoolkit-*.tar.gz
```
(note the `*` wildcard should be replaced with the current version number to install)


## Usage
 
Currently the required pre-processing of paired end reads is not supported natively by the tool (see the [Leapfrog Readme](https://github.com/mfiers/leapfrog/blob/master/README) for details).
The basic 'tec' tool is used as a wrapper for included programs like 'fingerprint' and 'compare'.
If no aruments are given the help text will be displayed:

```
$ tec
usage: Identify program to run [-h] {fingerprint,compare}
Identify program to run: error: the following arguments are required: program
usage: Identify program to run [-h] {fingerprint,compare}

positional arguments:
  {fingerprint,compare}

optional arguments:
  -h, --help            show this help message and exit
```

or:

```
$ tec fingerprint
usage: Identify potential TE flanking regions [-h]
                                              [-r [REFERENCES [REFERENCES ...]]]
                                              [-f [FAMILIES [FAMILIES ...]]]
                                              [-s {+,.,-} [{+,.,-} ...]]
                                              [-e EPS [EPS ...]]
                                              [-m MIN_READS] [-t THREADS]
                                              input_bam
Identify potential TE flanking regions: error: the following arguments are required: input_bam
usage: Identify potential TE flanking regions [-h]
                                              [-r [REFERENCES [REFERENCES ...]]]
                                              [-f [FAMILIES [FAMILIES ...]]]
                                              [-s {+,.,-} [{+,.,-} ...]]
                                              [-e EPS [EPS ...]]
                                              [-m MIN_READS] [-t THREADS]
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
                        TE grouping(s) to be used. Must be exact string
                        match(s) to start of read name
  -s {+,.,-} [{+,.,-} ...], --strands {+,.,-} [{+,.,-} ...]
                        Strand(s) to be analysed. Use + for forward, - for
                        reverse and . for both
  -e EPS [EPS ...], --eps EPS [EPS ...]
                        Maximum allowable distance among read tips to be
                        considered a cluster
  -m MIN_READS, --min_reads MIN_READS
                        Minimum allowable number of read tips to be considered
                        a cluster
  -t THREADS, --threads THREADS
                        Maximum number of cpu threads to be used
```
