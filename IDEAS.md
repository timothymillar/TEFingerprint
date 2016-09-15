# TODO

## Split libtec into components

* tectoolkit.gff
    * Gff3 class for output
* tectoolkit.readgroup
    * ReadGroup class
    * Loci class
    * read/parse sam/bam formats
* tectoolkit.cluster
    * clustering algorithms etc


## Read Collection Class
A class for a colection of reads that could replace/wrap the sam read dtype and the locus dtype.

The class would be based on a dtype with slots for start stop length and reversed. It would reperesent a collection of reads which could include a cluster  defined by its density

### Components:

* A numpy array of reads using a custom dtype
	* tip
	* tail
	* name

* Attributes
	* reference
	* family
	* strand

* Built in methods
	* Len
	* contains
	* iter
	* getitem
	* add?
	
* Methods
	* start
	* stop
	* Constructor from Bam
	* sort (calls np.sort())
	* simple density cluster

### Functions for a collection

* Sort
	* collection.sort(key = start)

* Merge 
	* Iterating through sorted clusters
	* remove duplicating reads?
	* use start stop
	* can be done with np.union1d if name is used
	* pseudo code:
```
def merge(collection):
	collection.sort(key=(start, stop))
	current = collection[0]
	for next in collection[1:]:
		if current.stop() > next.start():
			current = union(current, next)
	else: 
		yield current
		current = next
	yield current
```






