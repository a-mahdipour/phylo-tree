

install.views("Phylogenetics")




%--------------------------------------
library(devtools)
install.packages("devtools")


install_github("maphylogeny","lima1")
?


install.packages("biocGenerics")


install.packages("bioDist")
library("bioDist")

biocLite("affy")


source("http://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")

source("http://bioconductor.org/biocLite.R")
biocLite()

install.packages("ape");require(ape)


plot(cos, -10*pi,10*pi)

source("http://bioconductor.org/biocLite.R")
biocLite("bioDist")


source("http://bioconductor.org/biocLite.R")
biocLite()

source("http://bioconductor.org/biocLite.R")
biocLite("sva")
------------------------------------
library(maphylogeny)

M <- read.csv(file.choose(), header=T)


maphylo_bootstrap(M, group=c(),  bootstrap = 100,bootstrap_groups=FALSE,r=seq(.5,1.4,by=.1),dm="euclidean", snow=FALSE,na.rm=FALSE)

maphylo_reconstruct(dists,snow=FALSE,pm=nj, ...)


maphylo_reconstruct(dists,snow=FALSE,pm=nj)


maphylo_consensus_phylip(trees, outfile="phy", outgroup=1)


signature <- signatureFinder(which(colnames(geData) == "HIF1a"),  cpuCluster = aMakeCluster)

install.packages("ctv")

library("ctv")

library(genefilter)
install.views("Phylogenetics")

library("ape")
1)--------------------------------------------------
M <- read.csv(file.choose(), header=TRUE)

2)-----------------------
maphylo_bootstrap(M, group=c(),  bootstrap = 2,bootstrap_groups=FALSE,r=1, dm="pearson", snow=FALSE,na.rm=FALSE)
A <- maphylo_bootstrap(M, group=c(),  bootstrap = 100,bootstrap_groups=FALSE,r=1, dm="pearson", snow=FALSE,na.rm=FALSE)

A <- maphylo_bootstrap(M, group=c(),  bootstrap = 100,bootstrap_groups=FALSE,r=1, dm="euclidean", snow=FALSE,na.rm=FALSE)

install.packages("phangorn"); require(phangorn)
dm<- dist.logDet(M)
3)-----------------------------------
maphylo_reconstruct(A,snow=FALSE,pm=nj)?

B <- maphylo_reconstruct(A,snow=TRUE,pm=nj)
B <- maphylo_reconstruct(A,snow=TRUE)

plot(B[[1]])
4)-----------------------------------
phy <-read.csv(file.choose(), header=T)

maphylo_consensus_phylip(B, outfile="N-euc-nj", outgroup=1)

________________________________________________

dists = maphylo_bootstrap(M,r=c(1),bootstrap=1000,snow=TRUE)
stopCluster(cl)

## End(Not run)

# use Neighbor-Joining to reconstruct trees
trees = maphylo_reconstruct(dists)

## Not run: 
# it is also possible to reconstruct trees in parallel, especially useful 
# for slower methods and large trees

trees = maphylo_reconstruct(A,snow=TRUE,pm=fastme.bal)

trees = maphylo_reconstruct(dists,snow=TRUE,pm=fastme.bal)

## End(Not run)

# plot the first tree
plot(trees[[1]])

# Write the trees in NEXUS format. Now analyze the tree (in Dendroscope for
# example). 
write.nexus(trees, file="N-pear-fast.nxs")

## Not run: 
# PHYLIP is supported as a quick and easy way to analyze the trees. Requires
# the PHYLIP implementation in the EMBOSS extension EMBASSY.

maphylo_consensus_phylip(B, outfile="NwPMLMYHRUN-euc-nj")

maphylo_consensus_phylip(trees, outfile="NwPMLMYHRUN-euc-fast")
## End(Not run)

################ PVCLUST ###########################################

install.packages("pvclust")

library("pvclust")

pvcl <- pvclust(M, method.hclust="complete", method.dist="correlation", use.cor="pairwise.complete.obs", nboot=1000, r=1, store=FALSE)

pdf(file.path(saveres, "hclust_spca_ge_all.pdf"), width=35)
plot(pvcl)
dev.off()




