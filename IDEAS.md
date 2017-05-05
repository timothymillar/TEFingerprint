# TEFingerprint Ideas

This is record of feature/roadmap ideas 

1. Preprocessing:
	- [x] A simple python wrapper for the established process to identifies the “dangler” reads and  map them against the reference. 
	- [x] Store the name of the TE each dangler represents as a SAM tag rather than appended to the name.
	- 	[ ] Use soft clipped reads to get more accurate location of insertion

2. Fingerprinting:
	- [x] Basic “flat” clustering
	- [x] Hierarchical clustering to split close/nested clusters
	- [x] Output to GFF3 format (recorded statistics could be improved)
	- [ ] Associate clusters on opposite strands representing ends of the same TE insertion (removed until required)
	- [ ] Use of anchor reads to assess homozygosity of insertions (for re-sequence data)

3. Fingerprint comparisons:
	- [x] Identify “comparative bins” where clusters are found in at least one sample
	- [x] summary statistics for comparison across each bin (basic counts)
4. Output filetypes:
	- [x] GFF3
	- [x] GFF3 long form (no nested counts)
	- [x] Tabular data (available in python using pandas)
4. Filtering output:
	- [x] Script to filter the GFF output (could be improved)
