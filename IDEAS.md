# TEFingerprint Ideas

This is record of feature/roadmap ideas 

1. Preprocessing:
	- [ ] A simple python wrapper for the established process to identifies the “dangler” reads and  map them against the reference. Ideally the name of the TE each dangler represents is stored as a SAM tag rather than appended to the name (partially implemented)

2. Fingerprinting:
	- [x] Basic “flat” clustering
	- [x] Hierarchical clustering to split close/nested clusters
	- [x] Output to GFF3 format (recorded statistics could be improved)
	- [ ] Associate clusters on opposite strands representing ends of the same TE insertion (Not implemented)
	- [ ] Use soft clipped reads to get more accurate location of insertion (Not implemented)
	- [ ] Use of anchor reads to assess homozygosity of insertions

3. Fingerprint comparisons:
	- [x] Identify “comparative bins” where clusters are found in at least one sample
	- [x] summary statistics for comparison across each bin (very basic although still on par with competitors)
4. Output filetypes:
	- [x] GFF3 (not sorted properly)
	- [ ] SQLite using GffUtils
4. Filtering output:
	- [x] Script to filter the GFF output (slowish and crude)
	- [ ] Script to filter SQLite output