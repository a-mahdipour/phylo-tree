# phylo-tree
Making phylogenetics tree starting with a transcription data.



so here is my algorithm to produce a phylogeny for a differential expression matrix:
1. Applying t-test (unpaired case) between two subtypes (test and control)
2. Finding the p-value and significance for each probe by using gene
filtering FDR of level 0.01 (or 0.05)
3. Finding the average expression matrix
4. Finding the distance matrix among subtypes 
5. Constructing most probable trees based by considering the large number of
bootstrap (1000 or more).
to do all these steps I used "ape" R package in which I used the following functions:
- maphylo_bootstrapmaphylo_reconstruct functions
- and then maphylo_consensus_phylip function
-and finally maphylo_consensus_phylip from Phylip package to derive the consensus phylogeny
