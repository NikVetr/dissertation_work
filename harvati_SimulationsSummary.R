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


traitNumIter <- c(46)
nreps <- 500

rfDists <- matrix(0, nreps, 4); colnames(rfDists) <- c("mcc", "upgma", "nj", "pars")

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

    rfDists[i,1] <- RF.dist(trueTree, mccTree)
    rfDists[i,2] <- RF.dist(trueTree, upgmaTree)
    rfDists[i,3] <- RF.dist(trueTree, njTree)
    if(class(parsTree) == "multiPhylo"){
      rfDists[i,4] <- mean(sapply(1:length(parsTree), function(x) RF.dist(trueTree, parsTree[[x]])))
    } else if(class(parsTree) == "phylo"){
      rfDists[i,4] <- RF.dist(trueTree, parsTree)
    }
  }
}

save(rfDists, file =  "rfDists.txt")
load(file =  "rfDists.txt")
par(mfrow = c(4,1))
hist(rfDists[,1], breaks = 2*0:25, main = "bayesian mcc"); abline(v = mean(rfDists[,1]), col = 2, lwd = 2)
hist(rfDists[,3], breaks = 2*0:25, main = "neighbor joining"); abline(v = mean(rfDists[,3]), col = 2, lwd = 2)
hist(rfDists[,2], breaks = 2*0:25, main = "upgma"); abline(v = mean(rfDists[,2]), col = 2, lwd = 2)
hist(rfDists[,4], breaks = 2*0:25, main = "parsimony"); abline(v = mean(rfDists[,4]), col = 2, lwd = 2)

par(mfrow = c(1,1))
plot(rfDists[,3] - rfDists[,1], type = "l", ylim = c(-10, 48))
lines(rfDists[,2] - rfDists[,1], type = "l", col = 2)
lines(rfDists[,4] - rfDists[,1], type = "l", col = "blue")
abline(h = 0, lty = 2)
legend(legend = c("neighbor-joining", "upgma", "parsimony"), fill = c(1,2,"blue"), x = "topright"); title("difference in RF-dists from MCC Tree RF-dist")
