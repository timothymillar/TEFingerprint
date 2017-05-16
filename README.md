# TEFingerprint
A Python library with CLI for fingerprinting transposable elements based on paired-end reads.
This project is currently in a pe-alpha state and rapidly changing.


## Development Environment
A development environment can be built using the conda tool available in the [Anaconda Python distribution](https://www.continuum.io/downloads).
An [environment file](http://conda.pydata.org/docs/using/envs.html#share-an-environment) "environment.yml" is available in the root directory of this project. 

```
conda env create -f environment.yml
```

This will create a [virtual environment](http://conda.pydata.org/docs/using/envs.html) called `tefingerprint`.


## Dependencies

The current development version supports Python 3.5 or higher and pysam 0.9 or higher.


## Building and Installation

* Activate the development environment:

```
source activate tefingerprint
```

* Uninstall any previous version:

```
pip uninstall tefingerprint
```

* Navigate to the projects root directory then build and install:
	* (note the `*` wildcard should be replaced with the current version number to install)

```
python setup.py sdist
pip install dist/tefingerprint-*.tar.gz
```

* Alternatively, The project can be installed in editable mode from the projects root directory:

```
pip install -e ./
```

## Testing

Unit testing using [pytest](http://doc.pytest.org/en/latest/).
To run tests from the projects root directory:

```
python -m pytest -v ./
```

## Usage
 
The basic 'tef' tool is used as a wrapper for included programs: `preprocess`, `fingerprint`, `compare` and `filter_gff`.
The help text can be displayed for `tef` or any of the wrapped tols with the `-h` flag e.g:

```
tef -h 
tef preprocess -h 
tef fingerprint -h 
tef compare -h 
tef filter_gff -h
```

### Preprocess

The preprocessing pipeline requires a paired-end reads in fastq format, a library of tranposable elements in fasta format
and a reference genome in fasta format. Fasta files should be pre indexed using BWA. This pipeline also requires that 
BWA and Samtools are installed and on the users `$PATH`.

Paired end reads are aligned to the library of transposable elements.
Unmapped reads with mapped mates (dangler reads) are identified and extracted.
Dangler reads are mapped against the reference genome and tagged with the ID of their mates transposable element.

You should have the following four files:

* `reads1.fastq` and `reads2.fastq`: paired-end reads in fastq format.
* `reference.fasta` a reference genome in fasta format.
* `repeats.fasta ` a library of repeat elements in fasta format.


Index fasta files with bwa:

```
bwa index -a bwtsw reference.fasta
bwa index -a is repeats.fasta
```

Map paired end reads to repeat elements using `BWA MEM` and convert to a sorted bam file with `samtools`:

```
bwa mem -t 8 \
	repeats.fasta \
	reads1.fastq \
	reads2.fastq \
	| samtools view -Su - \
	| samtools sort - \
	-o reads_mapped_to_repeats.bam
```

Index the resulting bam file with `samtools`:

```
samtools index reads_mapped_to_repeats.bam
```

Runing the preprocess program on the indexed bam file:

```
tef preprocess \
    reads_mapped_to_repeats.bam
    --reference reference.fasta \
    --output danglers.bam \
    --threads 8
```

The output file `danglers.bam` contains reads mapped to the reference genome. Each of these reads is tagged with the repeat element that their pair was mapped to.

Arguments:

* The `threads` argument specifies how many threads to use for alignment in bwa (python components of the pipeline are currently single threaded).
* By default, the intermediate files are written to a temporary directory that is automatically removed when the pipeline is completed. These files can be saved by manually specifying a directory with the `--tempdir` option.
* By default, the same tag used to store repeat element names associated with each read is `ME` (Mate Element). This can be changed with the `--mate_element_tag` option.

### Fingerprint

Example usage:

```
tef fingerprint danglers.bam \
	-f family1 family2 ... \
	-m 20 \
	-e 500 \
	-t 4 \
	> fingerprint.gff
```

Where `danglers.bam` is the bam file being fingerprinted and `fingerprint.gff` is the output gff file.

Arguments:

* `-r/references` May optionally be used to specify a subset of chromosomes to fingerprint. By default all reference chromosomes are fingerprinted (based on the bam header).
* `-f/--families` Specifies the (super) families or grouping of repeated elements to fingerprint. These names are matched against the start of the mate element name i.e. the name `Gypsy` would treat reads with tagged with a mate element called `Gypsy3`, `Gypsy27` or `GypsyX` as the same.
* `-m/--minreads` Specifies the minimum number of read (tips) required to form a cluster. It is used in combination with `-e/epsilon`.
* `-e/epsilon` Specifies the maximum allowable distance among a set of read tips to be considered a (sub) cluster. Sub-clusters are calculated based on `-m/--minreads` and `-e/epsilon` and then overlapping sub-clusters are combined to create cluster.
* `-t/--threads` Specifies the number of CPU threads to use. The maximum number of threads that may be used is the same as the number of references specified.
 
Additional arguments:

 * `--min_eps` The minimum value of epsilon to be used in hierarchical clustering. Defaults to `0`.
 * `--hierarchical_clustering` Specifies wether or not to use the hierarchical clustering algorithm in order to differentiate between proximate clusters. Defaults to `True`.
 *  `--mate_element_tag` The sam tag used to specify the name of each reads mate element. Defaults to `ME`.


### Compare

Example usage:

```
tef compare danglers1.bam danglers2.bam ... \
	-f family1 family2 ... \
	-m 20 \
	-e 500 \
	-b 50 \
	-t 4 \
	> comparison.gff
```

Where `danglers1.bam ...` are the bam files being compared and `comparison.gff` is the output gff file.

Arguments:

* `-r/references` May optionally be used to specify a subset of chromosomes to fingerprint. By default all reference chromosomes are fingerprinted (based on the bam header).
* `-f/--families` Specifies the (super) families or grouping of repeated elements to fingerprint. These names are matched against the start of the mate element name i.e. the name `Gypsy` would treat reads with tagged with a mate element called `Gypsy3`, `Gypsy27` or `GypsyX` as the same.
* `-m/--minreads` Specifies the minimum number of read (tips) required to form a cluster. It is used in combination with `-e/epsilon`.
* `-e/epsilon` Specifies the maximum allowable distance among a set of read tips to be considered a (sub) cluster. Sub-clusters are calculated based on `-m/--minreads` and `-e/epsilon` and then overlapping sub-clusters are combined to create cluster.
* `-b/--fingerprint_buffer` Specifies a distance (in base pairs) to buffer fingerprints by before combining them into comparative bins. This is used to ensure that small clusters, that are slightly offset in different samples, are treated as a single comparative bin. It also improves the robustness of comparisons by allowing more reads to be included in each bin. Defaults to `0`
* `-t/--threads` Specifies the number of CPU threads to use. The maximum number of threads that may be used is the same as the number of references specified.
 
Additional arguments:
 
 * `--long_form` Option to produce a GFF file in which each comparative bin is duplicated for each input bam file. This produces a gff file that does not contatin nested lists of counts or source names. Defaults to `False`
 * `--min_eps` The minimum value of epsilon to be used in hierarchical clustering. Defaults to `0`.
 * `--hierarchical_clustering` Specifies wether or not to use the hierarchical clustering algorithm in order to differentiate between proximate clusters. Defaults to `True`.
 * `--bin_buffer` The same as `--fingerprint_buffer` but buffering is performed after fingerprints are combined, therefore less likely to combine slightly offset clusters. Defaults to `0`
 *  `--mate_element_tag` The sam tag used to specify the name of each reads mate element. Defaults to `ME`.


### Filter GFF

This script can be used to filter down the results of `fingerprint` or `compare`. Filters can be applied to attributes in the attribute column or to the first 8 standard gff3 columns. 

Multiple filters may be combined, in which case a feature must pass all of them to be kept. 

If an attribute contains a comma separated list of values e.g. `proportions=0.9,0.1,0.0` only one of the values must pass the filter for the feature to be retained.

Filters take the form `'<column/attribute><operator><value>'` where:

* `<column/attribute>` is the name of the column or attribute that the filter is applied to.
* `<operator>` is one of the following operators `=`, `==`, `!=`, `<` `>`, `>=`, `<=` that describes the comparason being performed.
* `<value>` is the value the each feature is compared to.

Filters should be contained within quotes `''` so that the operator is not interpreted as a shell command.

The following operators are only used for numerical comparisons: `<` `>`, `>=`, `<=`. 

The operators `=`, `==` and `!=` will try to compare values as numerical (floating points) but will also check for equivalence or non-equivalence of string values. Note that `=`, `==` are identical.


Example usage with two attribute filters:

```
tef filter_gff comparison.gff \
	-a 'counts>=10' 'proportions>0.95' \
	> comparison_filtered.gff
```
Where `comparison.gff ` is a gff file and `comparison_filtered.gff` is a filtered version of that file.

Arguments:

* `-c/--column_filters`: filters to apply to the first 8 standard gff3 columns. These should take the form `'<column><operator><value>'`
* `-a/--attribute_filters`: filters to apply to the attributes column. These should take the form `'<attribute><operator><value>'`