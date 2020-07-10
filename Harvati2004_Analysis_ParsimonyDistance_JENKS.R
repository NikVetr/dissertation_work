# library(geomorph)
# library(Matrix)
# library(phytools)
library(phangorn)
# library(abind)
# library(ape)
# library(coda)
# library(rwty)
# library(phangorn)
# library(lattice)
# library(latticeExtra)
# library(MCMCpack)
# library(readtext)
# library(StatMatch)
# library(mvtnorm)
# library(classInt)

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

setwd("~")


collapseHomo <- F

fileName <- paste0("harvati_empirical", ifelse(collapseHomo, "_nohomopops", ""))
ntraits <- 46
jenks_exp <- c(2,4,6)
numStartingTrees <- 50 #number of independent random starting trees
unchangedRun <- 500 #how long to keep the ratchet running without any improvement
traceOutput <- 1 #should the parsimony ratchet output stuff to screen?

for(transData in c("raw", "PCs")){
  
  for(jenks_exp_cats in 1:length(jenks_exp)){
    
    print(paste0(transData, " ", jenks_exp[jenks_exp_cats]))
    
    ncats <- as.numeric(strsplit(readLines(paste0("data/", fileName, "_", transData, "_jenks_", jenks_exp[jenks_exp_cats], "_expcats.txt"), warn = F), split = " ")[[1]])
    ncats <- sort(ncats)
    dataFiles <- paste0("data/", fileName, "_", transData, "_jenks_", jenks_exp[jenks_exp_cats], "_expcats_", ncats, "_cats_traits.tsv")
    if(!collapseHomo){
      ds <- lapply(1:length(dataFiles), function(df) trimws(t(sapply(1:27*2 - 1, function(x) (strsplit(readLines(dataFiles[df]), split = "\t"))[[x]]))))
    } else {
      ds <- lapply(1:length(dataFiles), function(df) trimws(t(sapply(1:21*2 - 1, function(x) (strsplit(readLines(dataFiles[df]), split = "\t"))[[x]]))))
    }
    popnames <- ds[[1]][,1]
    ds <- lapply(ds, function(x) x[,-1])
    d <- (do.call(cbind, ds))
    class(d) <- "numeric"
    rownames(d) <- popnames
    
    if(ncol(d) != 46){
      print("error in writing data! doubled the single trait file")
    }
    
    d <- d + 1 #increment so lowest trait is a 1
    
    
    #conver to phyDat format
    d.pD <- phyDat(d, type = "USER", levels = 1:max(d))
  
    #use the parsimony ratchet to search for most parsimonious trees a number of times
    mpTrees <- list()
    for(sepStart in 1:numStartingTrees){
      mpTree <- pratchet(d.pD, maxit = 1e5, method = "sankoff", k = unchangedRun, start = rtree(n = length(popnames), tip.label = sample(popnames)), trace = traceOutput, all = T, cost = wagnerCostMatrix(max(d)))
      mpTreeScore <- parsimony(tree = mpTree, data = d.pD, method = "sankoff", cost = wagnerCostMatrix(max(d)))
      mpTrees[[sepStart]] <- list(mpTree, mpTreeScore) 
    }
    
    #find unique most parsimonious trees
    mpScores <- sapply(1:length(mpTrees), function(x) mpTrees[[x]][[2]][1])
    mostPars <- which(mpScores == min(mpScores))
    bestTreesRough <- lapply(1:length(mostPars), function(x) mpTrees[[mostPars[x]]][[1]])
    bestTrees <- bestTreesRough[[1]]
    if(length(bestTreesRough) > 1){for(treeGroup in 2:length(bestTreesRough)){bestTrees <- c(bestTrees, bestTreesRough[[treeGroup]])}}
    if(class(bestTrees) != "phylo"){bestTrees <- unique(bestTrees)}
    
    write.tree(bestTrees, paste0("output/maximumParsimonyTree_", fileName, "_", transData, "_jenks_", jenks_exp[jenks_exp_cats], "_expcats", ifelse(collapseHomo, "_nohomopops", ""), ".txt"))
    
  }
}

#write the distance matrices to file too
collapseHomo <- F
fileName <- paste0("harvati_empirical", ifelse(collapseHomo, "_nohomopops", ""))
ntraits <- 46
jenks_exp <- c(2,4,6)
eucl <- F
univSlide <- T
for(transData in c("raw", "PCs")){
  
  for(jenks_exp_cats in 1:length(jenks_exp)){
    
    print(paste0(transData, " ", jenks_exp[jenks_exp_cats]))
    
    ncats <- as.numeric(strsplit(readLines(paste0("data/", fileName, "_", transData, "_jenks_", jenks_exp[jenks_exp_cats], "_expcats.txt"), warn = F), split = " ")[[1]])
    ncats <- sort(ncats)
    dataFiles <- paste0("data/", fileName, "_", transData, "_jenks_", jenks_exp[jenks_exp_cats], "_expcats_", ncats, "_cats_traits.tsv")
    if(!collapseHomo){
      ds <- lapply(1:length(dataFiles), function(df) trimws(t(sapply(1:27*2 - 1, function(x) (strsplit(readLines(dataFiles[df]), split = "\t"))[[x]]))))
    } else {
      ds <- lapply(1:length(dataFiles), function(df) trimws(t(sapply(1:21*2 - 1, function(x) (strsplit(readLines(dataFiles[df]), split = "\t"))[[x]]))))
    }
    popnames <- ds[[1]][,1]
    ds <- lapply(ds, function(x) x[,-1])
    d <- (do.call(cbind, ds))
    class(d) <- "numeric"
    rownames(d) <- popnames
    
    if(ncol(d) != 46){
      print("error in writing data! doubled the single trait file")
    }
    
    d <- d + 1 #increment so lowest trait is a 1
    
    # d.pD <- phyDat(d, type = "USER", levels = 1:max(d))
    if(eucl){distmat <- stats::dist(d)}
    if(univSlide){
      distmat_univ <- lapply(1:length(d[1,]), function(tr) stats::dist(d[,tr]))
      distmat <- distmat_univ[[1]]
      for(i in 2:length(distmat_univ)){
        distmat <- distmat + distmat_univ[[i]]
      }
    }
    
    writeDist(distmat, file = paste0("output/jenksDistances_", fileName, "_", transData, "_jenks_", jenks_exp[jenks_exp_cats], "_expcats", ifelse(collapseHomo, "_nohomopops", ""), ".dist"))
  }
}
