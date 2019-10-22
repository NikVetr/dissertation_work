library(abind)
library(parallel)
library(foreach)
library(doParallel)
library(ape)
library(coda)
library(rwty)
library(phangorn)
library(lattice)
library(latticeExtra)
library(MCMCpack)
library(readtext)
library(StatMatch)
library(adephylo)

mode <- function(x){
  dens <- density(x)
  dens$x[which.max(dens$y)]
}

#patristic distance calculation
traitNumIter <- c(2,4,8,16,32,57,64)
degreeMisspecification <- c("perfect", "interregional", "chimpHuman")
ntraits <- traitNumIter[7]
rateMatrixSpecification <- degreeMisspecification[1]
replicateNum <- 50
fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_", "replicate_", replicateNum)
trees1 <- read.tree(file = paste0("/Users/nikolai/mvBM_sims_all/output/", fileName, "_trees_run_1.trees"))
trees2 <- read.tree(file = paste0("/Users/nikolai/mvBM_sims_all/output/", fileName, "_trees_run_2.trees"))
treesAll <- trees <- c(trees1, trees2)
trees <- trees[seq(1, length(trees), by = 10)] #thin further
trueTree <- read.nexus(file = paste0("/Users/nikolai/mvBM_sims_all/trees/", fileName, "_tree.nex"))

tipNames <- trees[[1]]$tip.label
distTips(trees[[1]], tips = tipNames, method = "patristic", useC = TRUE)
patristicDists <- lapply(1:length(trees), 
                         function(x) distTips(trees[[x]], tips = tipNames, method = "patristic", useC = TRUE))
modalDists <- matrix( nrow = nrow(as.matrix(patristicDists[[1]])), ncol = ncol(as.matrix(patristicDists[[1]])), data = 0)
for(i in 1:(dim(modalDists)[1]-1)){
  print(i)
  for(j in (i+1):dim(modalDists)[2])  
    modalDists[i,j] <- mode(sapply(1:length(patristicDists), function(x) as.matrix(patristicDists[[x]])[i,j]))
}
modalDists <- modalDists + t(modalDists)
rownames(modalDists) <- colnames(modalDists) <- rownames(as.matrix(patristicDists[[1]]))
njTree <- nj(modalDists)
dist.topo(x = trueTree, y = njTree, method = "PH85")


pairNames <- vector()
for(i in 1:length(tipNames)){
  for(j in (i+1):length(tipNames)){
    pairNames <- c(pairNames, (paste(tipNames[i], "-", tipNames[j])))
  }
}

#making trace plots and marginal histograms
png(paste0("pairwisePatristicTraces_HowellsAnalyses_MC2vMC3.png"), width = 10000, height = 10080)
par(mfrow = c(29, 15))
for(i in 1:length(patristicDistsMC2[[1]])){
  index <- i
  distsThruChainMC2 <- sapply(1:length(patristicDistsMC2), function(x) patristicDistsMC2[[x]][index])
  distsThruChain3CMC3 <- sapply(1:length(patristicDists3CMC3), function(x) patristicDists3CMC3[[x]][index])
  distsThruChain5CMC3 <- sapply(1:length(patristicDists5CMC3), function(x) patristicDists5CMC3[[x]][index])
  plot(distsThruChainMC2, x = 1:16000, type = "l", col=rgb(1,0,0,0.3), ylab = "Patristic Distance between Tips", main = pairNames[index])
  lines(distsThruChain3CMC3, x = 1:16000, type = "l", col=rgb(0,1,0,0.3))
  lines(distsThruChain5CMC3, x = 1:16000, type = "l", col=rgb(0,0,1,0.3))
  legend('topleft', c("MC2", "3CMC3", "5CMC3"), col = c(rgb(1,0,0,0.3), rgb(0,1,0,0.3), rgb(0,0,1,0.3)), lwd = 1)
}
dev.off()

png(paste0("pairwisePatristicHistograms_HowellsAnalyses_MC2vMC3.png"), width = 10000, height = 10080)
par(mfrow = c(29, 15))
for(i in 1:length(patristicDistsMC2[[1]])){
  index <- i
  histBreaks <- 100
  distsThruChainMC2 <- sapply(1:length(patristicDistsMC2), function(x) patristicDistsMC2[[x]][index])
  distsThruChain3CMC3 <- sapply(1:length(patristicDists3CMC3), function(x) patristicDists3CMC3[[x]][index])
  distsThruChain5CMC3 <- sapply(1:length(patristicDists5CMC3), function(x) patristicDists5CMC3[[x]][index])
  hist(distsThruChainMC2, col=rgb(1,0,0,0.3), breaks = histBreaks, xlab = "Patristic Distance between Tips", main = pairNames[index])
  hist(distsThruChain3CMC3, breaks = histBreaks, col=rgb(0,1,0,0.3), add = T)
  hist(distsThruChain5CMC3, breaks = histBreaks, col=rgb(0,0,1,0.3), add = T)
  legend('topright', c("MC2", "3CMC3", "5CMC3"), col = c(rgb(1,0,0,0.3), rgb(0,1,0,0.3), rgb(0,0,1,0.3)), lwd = 1)
}
dev.off()


traitNumIter <- c(2,4,8,16,32,57,64)

degreeMisspecification <- c("perfect", "interregional", "chimpHuman")

distMatrixNJDists <- matrix(nrow = 2100, ncol = 6)
distMatrixUPGMADists <- matrix(nrow = 2100, ncol = 6)

rowCounter <- 0
# k iterates over degree of misspecification
for (k in 1:length(degreeMisspecification)) {
  
  for (j in 1:7) { #iterates over trait number
    
    for(i in 1:100){ #iterates over replicates
      rowCounter <- rowCounter + 1
      print(rowCounter)
      
      ntraits <- traitNumIter[j]
      rateMatrixSpecification <- degreeMisspecification[k]
      replicateNum <- i
      
      fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_", "replicate_", replicateNum)
      
      print("reading trees")
      trees1 <- read.tree(file = paste0("/Users/nikolai/mvBM_sims_all/output/", fileName, "_trees_run_1.trees"))
      trees2 <- read.tree(file = paste0("/Users/nikolai/mvBM_sims_all/output/", fileName, "_trees_run_2.trees"))
      trees <- c(trees1, trees2)
      trees <- trees[seq(1, length(trees), by = 10)] #thin further
      
      print("done reading")
      
      #get true tree
      trueTree <- read.nexus(file = paste0("/Users/nikolai/mvBM_sims_all/trees/", fileName, "_tree.nex"))
      
      #get modal distance matrix
      
      
      
      #compute trees
      njTree <- nj(mahDistMat)
      upgmaTree <- upgma(mahDistMat)
      rfDistNJ <- dist.topo(x = trueTree, y = njTree, method = "PH85")
      rfDistUPGMA <- dist.topo(x = trueTree, y = upgmaTree, method = "PH85")
      
      distMatrixNJDists[rowCounter,] <- c(ntraits, "NJ", replicateNum, 1, rfDistNJ, NA)
      distMatrixUPGMADists[rowCounter,] <- c(ntraits, "UPGMA", replicateNum, 1, rfDistNJ, NA)
      # kfDist <- dist.topo(x = trueTree, y = njTree, method = "score")
      
      
    }
  }
}

distMatrixUPGMADists <- as.data.frame(distMatrixUPGMADists)
distMatrixNJDists <- as.data.frame(distMatrixNJDists)
colnames(distMatrixUPGMADists) <- c("ntraits", "rateMatrixSpecification", "replicateNum", "idk", "rfDist", "kfDist")
colnames(distMatrixNJDists) <- c("ntraits", "rateMatrixSpecification", "replicateNum", "idk", "rfDist", "kfDist")













tree <- rtree(n=10)
dists <- distTips(tree, tips = tree$tip.label, method = "patristic", useC = TRUE)
njTree <- nj(dists)
dist.topo(x = tree, y = njTree, method = "PH85")
njDists <- distTips(njTree, tips = njTree$tip.label, method = "patristic", useC = TRUE)
all(njDists == dists)
njTree == tree