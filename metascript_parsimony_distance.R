setwd(dir = "/Volumes/1TB/Harvati/")
replicateNum <- 
ntraits <- 

alreadySimulated <- T
meanValueCodingDo <- T
numStartingTrees <- 50 #number of independent random starting trees
unchangedRun <- 500 #how long to keep the ratchet running without any improvement
traceOutput <- 0 #should the parsimony ratchet output stuff to screen?
l <- 1

library(abind)
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

convNum2Str <- function(nums, key){
  sapply(1:length(nums), function(x) key[nums[x]])
}

wagnerCostMatrix <- function(dim){
  mat <- matrix(0, nrow = dim, ncol = dim)
  for(i in 1:(dim-1)){
    for(j in (i+1):dim){
      mat[i,j] <- abs(i-j)
    }
  }
  mat + t(mat)
}

fileName <- paste0("simulations_", ntraits, "traits_", "replicate_", replicateNum)

#get true tree
trueTree <- read.nexus(file = paste0("trees/", fileName, "_tree.nex"))
tips <- trueTree$tip.label
numPops <- length(tips)

for(type in 1:2){ #type == 1 : raw, type == 2 : PCs
  if(!alreadySimulated){
    #get traits
    traits <- strsplit(readLines(paste0("data/", fileName, "_traits.nex")), split = "\t")
    traits <- traits[8:(length(traits)-2)]
    rwn <- sapply(1:length(traits), function(x) strsplit(traits[[x]], split = " ")[[1]][1])
    traits <- t(as.matrix(sapply(1:length(traits), function(x) as.numeric(strsplit(traits[[x]], split = " ")[[1]][-1]))))
    rownames(traits) <- rwn
    
    #get pooled within-group covariance matrix
    wgCM <- cov <- as.matrix(read.table(paste0("rates/", fileName, "_rates.tsv")))
    
    #apply PC transform if appropriate
    if(type == 2){
      V <- eigen(wgCM)$vectors
      L <- eigen(wgCM)$values
      VarExpl <- round(L / sum(L), 4) * 100
      PC_Scores <- t(V) %*% t(traits) 
      traits <- t(PC_Scores)
    }

    #get simulated individuals      
    numIndiv <- 41 #average sample size in dataset -- 1097 / 27 -- but gives all 0s
    numPops <- length(traits[,1])
    simulIndiv <- array(data = 0, dim = c(numPops,ntraits,numIndiv))
    for (pops in 1:numPops){
      simulIndiv[pops,,] <- t(rmvnorm(n = numIndiv, mean = traits[pops,], sigma = wgCM))
    }
    
    #get discrete traits via divergence coding
    divergenceCodedTraits <- matrix(data = 0, nrow = numPops, ncol = ntraits)
    rownames(divergenceCodedTraits) <- rownames(traits)
    pval <- 0.05
    
    for(trait in 1:ntraits){
      divergenceMatrix <- matrix(data = 0, nrow = numPops, ncol = numPops)
      for(poprow in 2:numPops){
        for(popcol in 1:(poprow - 1)){
          poprowindiv <- simulIndiv[poprow,trait,]
          popcolindiv <- simulIndiv[popcol,trait,]
          significance <- t.test(poprowindiv, popcolindiv)$p.value < pval
          if(significance){
            if(mean(popcolindiv) > mean(poprowindiv)){
              divergenceMatrix[poprow, popcol] <- 1
            } else {
              divergenceMatrix[poprow, popcol] <- -1
            }
          } else {
            divergenceMatrix[poprow, popcol] <- 0
          }
        }
      }
      divergenceMatrix <- divergenceMatrix - t(divergenceMatrix)
      divergenceCodedTraits[,trait] <- sapply(1:numPops, function(x) sum(divergenceMatrix[,x]))
    }
    
    for(foo in 1:length(divergenceCodedTraits[1,])){
      divergenceCodedTraits[,foo] <- divergenceCodedTraits[,foo] + abs(min(divergenceCodedTraits[,foo])) + 1
    }
  } else {
    traitsPath_divergenceCoding <- paste0("data_discrete/", fileName, ifelse(type == 1, "_raw", "_PCs"), "_DC_traits.tsv")
    divergenceCodedTraits <- read.table(traitsPath_divergenceCoding, header = F, sep = "\t", colClasses = "character")
    rwn <- trimws(divergenceCodedTraits[,1]); divergenceCodedTraits[] <- as.numeric(as.matrix(divergenceCodedTraits)); rownames(divergenceCodedTraits) <- rwn; 
    divergenceCodedTraits <- as.matrix(divergenceCodedTraits[,-1])
    for(foo in 1:length(divergenceCodedTraits[1,])){
      divergenceCodedTraits[,foo] <- divergenceCodedTraits[,foo] + 1
    }
  }
  #conver to phyDat format
  DCT.pD <- phyDat(divergenceCodedTraits, type = "USER", levels = 1:(numPops*2-1))
  
  #use the parsimony ratchet to search for most parsimonious trees a number of times
  mpTrees <- list()
  for(sepStart in 1:numStartingTrees){
    mpTree <- pratchet(DCT.pD, maxit = 1e5, method = "sankoff", k = unchangedRun, start = rtree(n = length(tips), tip.label = sample(tips)), trace = traceOutput, all = T, cost = wagnerCostMatrix(numPops*2-1))
    mpTreeScore <- parsimony(tree = mpTree, data = DCT.pD, method = "sankoff", cost = wagnerCostMatrix(numPops*2-1))
    mpTrees[[sepStart]] <- c(mpTree, mpTreeScore) 
  }
  
  #find unique most parsimonious trees
  mpScores <- sapply(1:length(mpTrees), function(x) mpTrees[[x]][[2]][1])
  mostPars <- which(mpScores == min(mpScores))
  bestTreesRough <- lapply(1:length(mostPars), function(x) mpTrees[[mostPars[x]]][[1]])
  bestTrees <- bestTreesRough[[1]]
  if(length(bestTreesRough) > 1){for(treeGroup in 2:length(bestTreesRough)){bestTrees <- c(bestTrees, bestTreesRough[[treeGroup]])}}
  if(class(bestTrees) != "phylo"){bestTrees <- unique(bestTrees)}
  
  #calculate average RF distance from the true tree
  # if(class(bestTrees) == "multiPhylo"){avgRFdistDC <- mean(sapply(1:length(bestTrees), function(x) dist.topo(unroot(trueTree), bestTrees[[x]])[1]))
  # } else {avgRFdistDC <- dist.topo(unroot(trueTree), bestTrees)[1]}
  
  if(meanValueCodingDo){
    #do mean value coding
    if(!alreadySimulated){
      traitMeans <- sapply(1:length(traits[1,]), function(x) mean(traits[,x]))
      meanBinarizedTraits <- t(sapply(1:length(traits[,1]), function(x) traits[x,] > traitMeans))
      rownames(meanBinarizedTraits) <- rownames(traits)
      meanBinarizedTraits[meanBinarizedTraits] <- "L"
      meanBinarizedTraits[meanBinarizedTraits == "FALSE"] <- "S"
    } else {
      traitsPath_meanValueCoding <- paste0("data_discrete/", fileName, ifelse(type == 1, "_raw", "_PCs"), "_MVC_traits.nex")
      tmpMBT <- strsplit(readLines(traitsPath_meanValueCoding)[8:34], split =  " ")
      meanBinarizedTraits <- t(sapply(1:length(tmpMBT), function(x) tmpMBT[[x]][-1]))
      meanBinarizedTraits[meanBinarizedTraits == "1"] <- "L"; meanBinarizedTraits[meanBinarizedTraits == "0"] <- "S" 
      rownames(meanBinarizedTraits) <- sapply(1:length(tmpMBT), function(x) tmpMBT[[x]][1])
    }
    mBT.pD <- phyDat(meanBinarizedTraits, type = "USER", levels = c("L", "S"))
    mpTrees <- list()
    for(sepStart in 1:numStartingTrees){
      mpTree <- pratchet(mBT.pD, maxit = 1e5, method = "sankoff", k = unchangedRun, start = rtree(n = length(tips), tip.label = sample(tips)), trace = traceOutput, all = T)
      mpTreeScore <- parsimony(tree = mpTree, data = mBT.pD, method = "sankoff")
      mpTrees[[sepStart]] <- c(mpTree, mpTreeScore) 
    }
    mpScores <- sapply(1:length(mpTrees), function(x) mpTrees[[x]][[2]][1])
    mostPars <- which(mpScores == min(mpScores))
    bestTreesRough_MVC <- lapply(1:length(mostPars), function(x) mpTrees[[mostPars[x]]][[1]])
    bestTrees_MVC <- bestTreesRough_MVC[[1]]
    if(length(bestTreesRough_MVC) > 1){for(treeGroup in 2:length(bestTreesRough_MVC)){bestTrees_MVC <- c(bestTrees_MVC, bestTreesRough_MVC[[treeGroup]])}}
    if(class(bestTrees_MVC) != "phylo"){bestTrees_MVC <- unique(bestTrees_MVC)}
    # if(class(bestTrees) == "multiPhylo"){avgRFdistMVC <- mean(sapply(1:length(bestTrees), function(x) dist.topo(unroot(trueTree), bestTrees[[x]])[1]))
    # } else {avgRFdistMVC <- dist.topo(unroot(trueTree), bestTrees)[1]}
  }
  if(type == 1){
    write.tree(bestTrees, file = paste0("output_parsimony/", fileName, "_divergenceCodingParsimony.txt"))
    write.tree(bestTrees_MVC, file = paste0("output_parsimony/", fileName, "_meanValueCodingParsimony.txt"))
  } else if (type == 2){
    write.tree(bestTrees, file = paste0("output_parsimony_PCA/", fileName, "_divergenceCodingParsimony.txt"))
    write.tree(bestTrees_MVC, file = paste0("output_parsimony_PCA/", fileName, "_meanValueCodingParsimony.txt"))
  }
}

##############################
### compute distance trees ###
##############################

#get rate matrix
wgCM <- cov <- as.matrix(read.table(paste0("rates/", fileName, "_rates.tsv")))

#get traits
traits <- strsplit(readLines(paste0("data/", fileName, "_traits.nex")), split = "\t")
traits <- traits[8:(length(traits)-2)]
rwn <- sapply(1:length(traits), function(x) strsplit(traits[[x]], split = " ")[[1]][1])
traits <- t(as.matrix(sapply(1:length(traits), function(x) as.numeric(strsplit(traits[[x]], split = " ")[[1]][-1]))))
rownames(traits) <- rwn

mahDistMat <- mahaMatrix(traits, wgCM)
njTree <- nj(mahDistMat)
upgmaTree <- upgma(mahDistMat)
write.tree(njTree, file = paste0("output_distance/", fileName, "_neighborJoining.txt"))
write.tree(upgmaTree, file = paste0("output_distance/", fileName, "_upgma.txt"))

