########################################
## Constructing most probable trees   ##
## from expression data               ##
########################################

library(maphylogeny)
library(Phylogenetics)
library(ape)

# Reading the transcription data  
mat <- read.csv(file.choose(), header=TRUE)

# 
boots.mat <- maphylo_bootstrap(mat, 
                               group=c(),  
                               bootstrap = 100, 
                               bootstrap_groups=FALSE,
                               r=1, 
                               dm="pearson",  # "euclidean" 
                               snow=FALSE,
                               na.rm=FALSE)

# Defining the distance matrix
require(phangorn)
dist.mat <- dist.logDet(boots.mat)

# Constructing potential phylogenies
recons.mat <- maphylo_reconstruct(boots.mat, 
                                  snow=FALSE, 
                                  pm=nj  # fast
                                 )

# Checking the method by sketching one of the trees
plot(recons.mat [[1]])
4)-----------------------------------
phy <-read.csv(file.choose(), header=T)

maphylo_consensus_phylip(recons.mat, 
                         outfile="N-euc-nj", 
                         outgroup=1)

________________________________________________


# Using Neighbor-Joining to reconstruct trees
dist.mat <-  maphylo_bootstrap(mat, 
                               r=c(1), 
                               bootstrap=1000,
                               snow=TRUE)
phylogenies <-  maphylo_reconstruct(dist.mat)

# Reconstructing trees for slower methods and large trees
phylogenies <-  maphylo_reconstruct(boots.mat,
                                    snow=TRUE,
                                    pm=fastme.bal)
phylogenies <-  maphylo_reconstruct(dist.mat,
                                    snow=TRUE,
                                    pm=fastme.bal)

# plot the first tree
plot(phylogenies[[1]])

# Write the phylogenies in NEXUS format. Now analyze the tree (e.g. Dendroscope)
write.nexus(phylogenies, file="N-pear-fast.nxs")


# PHYLIP implementation in the EMBOSS extension EMBASSY.
maphylo_consensus_phylip(recons.mat,outfile="PhylipRun-euc-nj")
maphylo_consensus_phylip(phylogenies, outfile="hylipRun-euc-fast")


################ PVCLUST ###########################################
# Hierarchical clustering with p-values using pvclust R package
library("pvclust")
pvcl <- pvclust(M, method.hclust="complete", method.dist="correlation", use.cor="pairwise.complete.obs", nboot=1000, r=1, store=FALSE)

pdf(file.path(file, "hclust_spca_ge_all.pdf"), width=35)
plot(pvcl)
dev.off()




