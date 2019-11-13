# phylo-tree

_Making phylogenetics tree starting with an expression data_

---

A suggesting algorithm to produce a phylogeny for a transcription data is as follows:
1. Applying (unpaired) t-test between two subtypes (test and control)
2. Finding the p-value and significance for each probe by using gene
filtering FDR of level 0.01 (or 0.05)
3. Finding the average expression matrix
4. Finding the distance matrix among subtypes 
5. Constructing most probable trees based by considering the large number of
bootstrap (100, 1000 or even more for smaller trees).

Note: To apply all these steps you may try _ape_ R package in which the following functions are used:
* maphylo_bootstrapmaphylo_reconstruct function,
* then maphylo_consensus_phylip function
* and finally maphylo_consensus_phylip from _Phylip_ package to derive the consensus phylogeny

