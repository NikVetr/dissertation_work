setwd(dir = "/Volumes/2TB/mvBM_sims_PB_noiseless")
# setwd(dir = "/Volumes/2TB/mvBM_sims_PB_noisy")
# setwd(dir = "/Volumes/2TB/uvBM_sims_PB_noiseless")
# setwd(dir = "/Volumes/2TB/uvBM_sims_PB_noisy")

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

traitNumIter <- c(2,4,8,16,32,64,128)
noisyness <- strsplit(getwd(), split = "_")[[1]][4]
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

distMatrixNJDists <- matrix(nrow = length(degreeCorrelation)*length(traitNumIter)*100, ncol = 6)
distMatrixUPGMADists <- matrix(nrow = length(degreeCorrelation)*length(traitNumIter)*100, ncol = 6)

rowCounter <- 0
# k iterates over degree of misspecification
for (k in 1:length(degreeCorrelation)) {
  
  for (j in 1:length(traitNumIter)) { #iterates over trait number
    
    for(i in 1:100){ #iterates over replicates
      rowCounter <- rowCounter + 1
      print(rowCounter)
      
      ntraits <- traitNumIter[j]
      rateMatrixCorrelation <- degreeCorrelation[k]
      replicateNum <- i
      fileName <- paste0("simulations_", ntraits, "traits_perfectRateMatrix_", rateMatrixCorrelation, "correlationStrength_", noisyness, "_replicate_", replicateNum)
      
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

      #compute trees
      njTree <- nj(mahDistMat)
      upgmaTree <- upgma(mahDistMat)
      rfDistNJ <- dist.topo(x = trueTree, y = njTree, method = "PH85")
      rfDistUPGMA <- dist.topo(x = trueTree, y = upgmaTree, method = "PH85")
      kfDistNJ <- dist.topo(x = trueTree, y = njTree, method = "score")
      kfDistUPGMA <- dist.topo(x = trueTree, y = upgmaTree, method = "score")
      
      distMatrixNJDists[rowCounter,] <- c(ntraits, rateMatrixCorrelation, "NJ", replicateNum, rfDistNJ, kfDistNJ)
      distMatrixUPGMADists[rowCounter,] <- c(ntraits, rateMatrixCorrelation, "UPGMA", replicateNum, rfDistUPGMA, kfDistUPGMA)
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



