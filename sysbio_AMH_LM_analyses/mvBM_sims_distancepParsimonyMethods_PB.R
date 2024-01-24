setwd(dir = "/Volumes/2TB/mvBM_sims_PB")


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
library(mvtnorm)

traitNumIter <- c(2,4,8,16,32,64)

degreeCorrelation <- c("3")
degreeCorrelation <- c("2", "4", "6", "8", "9")

maha <- function(x1, x2, vcvm, squared = FALSE){
  squaredDist <- t(x1-x2) %*% solve(vcvm) %*% (x1-x2)
  if(!squared){
    squaredDist
  } else {
    squaredDist^0.5
  }
} 

mahaMatrix <- function(data, vcvm){
  mat <- matrix (nrow = nrow(data), ncol = nrow(data), 0)
  for(i in 1:nrow(mat)){
    for(j in i:nrow(mat)){
      mat[i,j] <- maha(data[i,], data[j,], vcvm)
    }
  }
  mat <- mat + t(mat)
  rownames(mat) <- colnames(mat) <-rownames(data)
  mat
}

distMatrixNJDists <- matrix(nrow = 3000, ncol = 6)
distMatrixUPGMADists <- matrix(nrow = 3000, ncol = 6)
maximumParsimonyDists <- matrix(nrow = 3000, ncol = 6)

rowCounter <- 0
# k iterates over degree of misspecification
for (k in 1:length(degreeCorrelation)) {
  
  for (j in 1:6) { #iterates over trait number
    
    for(i in 1:100){ #iterates over replicates
      rowCounter <- rowCounter + 1
      print(rowCounter)
      
      ntraits <- traitNumIter[j]
      rateMatrixCorrelation <- degreeCorrelation[k]
      replicateNum <- i
      fileName <- paste0("simulations_", ntraits, "traits_perfectRateMatrix_", rateMatrixCorrelation, "correlationStrength_", "replicate_", replicateNum)
      
      #get true tree
      trueTree <- read.nexus(file = paste0("trees/", fileName, "_tree.nex"))
      
      #get traits
      traits <- strsplit(readtext(file = paste0("data/", fileName, "_traits.nex"))$text, split = "\n")[[1]][8:37]
      traits <- t(sapply(1:30, function (x) strsplit(traits[x], split = " ")[[1]]))
      names <- traits[,1]
      traits <- traits [,-1]
      class(traits) <- "numeric"
      rownames(traits) <- names
      
      #get pooled within-group covariance matrix
      wgCM <- as.matrix(read.table(file = paste0("rates/", fileName, "_rates.tsv")))
      
      #get mahalonibis distance matrix
      mahDistMat <- mahaMatrix(traits, wgCM)

      # #get simulated individuals      
      # numIndiv <- 50
      # numPops <- length(traits[,1])
      # simulIndiv <- array(data = 0, dim = c(numPops,ntraits,numIndiv))
      # for (pops in 1:numPops){
      #   simulIndiv[pops,,] <- rmvnorm(n = numIndiv, mean = traits[pops,], sigma = wgCM)
      #   
      # }
      # 
      # #get discrete traits via divergence coding
      # divergenceCodedTraits <- matrix(data = 0, nrow = numPops, ncol = ntraits)
      # 
      # for(trait in 1:ntraits){
      #   divergenceMatrix <- matrix(data = 0, nrow = numPops, ncol = numPops)
      #     for(bar in 1:(numPops-1)){
      #     for(foo in (bar+1):numPops){
      #       significance <- t.test(simulIndiv[foo,trait,], simulIndiv[bar,trait,])$p.value < 0.05
      #       if(significance){
      #         if(mean(simulIndiv[bar,1,]) > mean(simulIndiv[foo,1,])){
      #           divergenceMatrix[foo, bar] <- 1
      #         } else {
      #           divergenceMatrix[foo, bar] <- -1
      #         }
      #       } else {
      #         divergenceMatrix[foo, bar] <- 0
      #       }
      #     }
      #   }
      #   divergenceMatrix <- divergenceMatrix - t(divergenceMatrix)
      #   divergenceCodedTraits[,trait] <- sapply(1:numPops, function(x) sum(divergenceMatrix[,x]))
      # }
      # rownames(divergenceCodedTraits) <- rownames(traits)
      # 
      # divergenceCodedTraitsPhyDat <- phyDat(divergenceCodedTraits[1:5,], type = "USER", levels = (-29:29))
      # pratchet(data = divergenceCodedTraitsPhyDat)
      # parsimony(divergenceCodedTraitsPhyDat, tree = rtree(n = length(rownames(traits)), tip.label = rownames(traits)))
      
      #get discrete traits via gap coding
      # trait <- simulIndiv[,1,]
      # rownames(trait) <- rownames(traits)
      # traitMinMax <- cbind(sapply(1:numPops, function(x) min(trait[x,])),sapply(1:numPops, function(x) max(trait[x,])))
      # rownames(traitMinMax) <- rownames(traits)
      # groups <- rep(0, numPops)
      #nevermind, gap coding doesn't work, within-group differences too much greater than between group differences
      
      #compute trees
      njTree <- nj(mahDistMat)
      upgmaTree <- upgma(mahDistMat)
      rfDistNJ <- dist.topo(x = trueTree, y = njTree, method = "PH85")
      rfDistUPGMA <- dist.topo(x = trueTree, y = upgmaTree, method = "PH85")
      kfDistNJ <- dist.topo(x = trueTree, y = njTree, method = "score")
      kfDistUPGMA <- dist.topo(x = trueTree, y = upgmaTree, method = "score")
      
      distMatrixNJDists[rowCounter,] <- c(ntraits, rateMatrixCorrelation, "NJ", replicateNum, rfDistNJ, kfDistNJ)
      distMatrixUPGMADists[rowCounter,] <- c(ntraits, rateMatrixCorrelation, "UPGMA", replicateNum, rfDistUPGMA, kfDistUPGMA)
      # kfDist <- dist.topo(x = trueTree, y = njTree, method = "score")

      
    }
  }
}

distMatrixUPGMADists <- as.data.frame(distMatrixUPGMADists)
distMatrixNJDists <- as.data.frame(distMatrixNJDists)
colnames(distMatrixUPGMADists) <- c("ntraits", "rateMatrixCorrelation", "type", "replicateNum", "rfDist", "kfDist")
colnames(distMatrixNJDists) <- c("ntraits", "rateMatrixCorrelation", "type", "replicateNum", "rfDist", "kfDist")

write.table(distMatrixUPGMADists, file = "analyses/distMatrixUPGMADists.txt")
write.table(distMatrixNJDists, file = "analyses/distMatrixNJDists.txt")

read.table(file = "analyses/distMatrixUPGMADists.txt")
read.table(file = "analyses/distMatrixNJDists.txt")

