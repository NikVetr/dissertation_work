setwd("/Volumes/1TB/Harvati")

library(ape)
library(coda)
library(rwty)
library(phangorn)
library(lattice)
library(latticeExtra)
library(MCMCpack)
library(mvMORPH)
library(phytools)

traitNumIter <- (15*(2^(0:4)))*3+1
# traitNumIter <- c(46)
nreps <- 100

rfDists <- matrix(0, nreps, 7); 
colnames(rfDists) <- c("mcc", "upgma", "nj", "pars", "pars-pca", "mk-mvc", "mk-mvc-pca")
rfDists <- lapply(1:length(traitNumIter), function(x) rfDists)

for (j in 1:length(traitNumIter)) { #iterates over trait number
  
  for(i in 1:nreps){ #iterates over replicates
    
    print(i)
    
    replicateNum <- i
    ntraits <- traitNumIter[j]
    fileName <- paste0("simulations_", ntraits, "traits_", "replicate_", replicateNum)
    
    trueTree <- read.nexus(file = paste0("trees/", fileName, "_tree.nex"))

    trees1 <- read.tree(paste0("output/", fileName, "_trees_run_1.trees") )
    trees2 <- read.tree(paste0("output/", fileName, "_trees_run_2.trees") )
    trees <- c(trees1, trees2)

    mccTree <- maxCladeCred(trees)
    njTree <- read.tree(file = paste0("output_distance/", fileName, "_neighborJoining.txt"))
    upgmaTree <- read.tree(file = paste0("output_distance/", fileName, "_upgma.txt"))
    parsTree <- read.tree(file = paste0("output_parsimony/", fileName, "_divergenceCodingParsimony.txt"))
    parsTree_pca <- read.tree(file = paste0("output_parsimony_PCA/", fileName, "_divergenceCodingParsimony.txt"))
    
    mktrees1 <- read.tree(paste0("output_discrete/", fileName, "_MVC_raw_trees_run_1.trees") )
    mktrees2 <- read.tree(paste0("output_discrete/", fileName, "_MVC_raw_trees_run_2.trees") )
    mktrees <- c(mktrees1, mktrees2)
    mccTree_mkRAW <- maxCladeCred(mktrees)
    mkPCAtrees1 <- read.tree(paste0("output_discrete/", fileName, "_MVC_PCs_trees_run_1.trees") )
    mkPCAtrees2 <- read.tree(paste0("output_discrete/", fileName, "_MVC_PCs_trees_run_2.trees") )
    mkPCAtrees <- c(mkPCAtrees1, mkPCAtrees2)
    mccTree_mkPCA <- maxCladeCred(mkPCAtrees)
    
    rfDists[[j]][i,1] <- RF.dist(trueTree, mccTree)
    rfDists[[j]][i,2] <- RF.dist(trueTree, upgmaTree)
    rfDists[[j]][i,3] <- RF.dist(trueTree, njTree)
    if(class(parsTree) == "multiPhylo"){
      rfDists[[j]][i,4] <- mean(sapply(1:length(parsTree), function(x) RF.dist(trueTree, parsTree[[x]])))
    } else if(class(parsTree) == "phylo"){
      rfDists[[j]][i,4] <- RF.dist(trueTree, parsTree)
    }    
    if(class(parsTree_pca) == "multiPhylo"){
      rfDists[[j]][i,5] <- mean(sapply(1:length(parsTree_pca), function(x) RF.dist(trueTree, parsTree_pca[[x]])))
    } else if(class(parsTree_pca) == "phylo"){
      rfDists[[j]][i,5] <- RF.dist(trueTree, parsTree_pca)
    }
    rfDists[[j]][i,6] <- RF.dist(trueTree, mccTree_mkRAW)
    rfDists[[j]][i,7] <- RF.dist(trueTree, mccTree_mkPCA)
  }
}

save(rfDists, file =  "rfDists.txt")
load(file =  "rfDists.txt")
par(mfrow = c(5,1))
hist(rfDists[,1], breaks = 2*0:25, main = "bayesian mcc"); abline(v = mean(rfDists[,1]), col = 2, lwd = 2)
hist(rfDists[,3], breaks = 2*0:25, main = "neighbor joining"); abline(v = mean(rfDists[,3]), col = 2, lwd = 2)
hist(rfDists[,2], breaks = 2*0:25, main = "upgma"); abline(v = mean(rfDists[,2]), col = 2, lwd = 2)
hist(rfDists[,4], breaks = 2*0:25, main = "parsimony"); abline(v = mean(rfDists[,4]), col = 2, lwd = 2)
hist(rfDists[,5], breaks = 2*0:25, main = "parsimony-PCA"); abline(v = mean(rfDists[,5]), col = 2, lwd = 2)

par(mfrow = c(1,1))
library(RColorBrewer)
cols <- brewer.pal(5, "Dark2")
bw <- 1.5
alph <- 0.25
plot(density(rfDists[,1], bw = bw), xlim = c(0,50), ylim = c(0,0.1), col = "white", main = "rf-dists to true, data-generating tree", xlab = "Robinson-Foulds Distance")
for(i in 1:5){
  polygon(density(rfDists[,i], bw = bw), col = add.alpha(cols[i], alpha = alph))
  lines(density(rfDists[,i], bw = bw), col = cols[i]) 
}
legend(legend = colnames(rfDists), fill = cols, x = "topright")


par(mfrow = c(1,1))
plot(rfDists[,3] - rfDists[,1], type = "l", ylim = c(-10, 48))
lines(rfDists[,2] - rfDists[,1], type = "l", col = 2)
lines(rfDists[,4] - rfDists[,1], type = "l", col = "blue")
lines(rfDists[,5] - rfDists[,1], type = "l", col = "green")
abline(h = 0, lty = 2)
legend(legend = c("neighbor-joining", "upgma", "parsimony", "parsimony-pca"), fill = c(1,2,"blue","green"), x = "topright"); title("difference in RF-dists from MCC Tree RF-dist")
