TEFingerprint Ideas
===================

This is record of feature/roadmap ideas

1. Preprocessing:

   -  [x] A simple python wrapper for the established process to
      identifies the “dangler” reads and map them against the reference.
   -  [x] Store the name of the TE each dangler represents as a SAM tag
      rather than appended to the name.
   -  [x] Use soft clipped reads to get more accurate location of
      insertion

2. Fingerprinting:

   -  [x] Basic “flat” clustering
   -  [x] Hierarchical clustering to split close/nested clusters
   -  [X] Associate clusters on opposite strands representing ends of
      the same TE insertion (removed until required)
   -  [x] Tag known transposons from annotation GFF
   -  [ ] Use of anchor reads to assess homozygosity of insertions (for
      re-sequence data)

3. Fingerprint comparisons:

   -  [x] Identify “comparative bins” where clusters are found in at
      least one sample
   -  [x] summary statistics for comparison across each bin (basic
      counts)

4. Output filetypes:

   -  [x] GFF3
   -  [x] Tabular data (available in python using pandas)

5. Filtering output:

   -  [x] Script to filter the GFF output (could be improved)
