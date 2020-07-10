setwd(dir = "/Volumes/1TB/Harvati/")
replicateNum <- 
ntraits <- 
expCat <- 
transData <- 

numStartingTrees <- 50 #number of independent random starting trees
unchangedRun <- 500 #how long to keep the ratchet running without any improvement
traceOutput <- 0 #should the parsimony ratchet output stuff to screen?

library(phangorn)

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

ncats <- as.numeric(strsplit(readLines(paste0("data_discrete/", fileName, "_", transData, "_jenks_", expCat, "_expcats.txt"), warn = F), split = " ")[[1]])
ncats <- sort(ncats)

dataFiles <- paste0("data_discrete/", fileName, "_", transData, "_jenks_", expCat, "_expcats_", ncats, "_cats_traits.tsv")

ds <- lapply(1:length(dataFiles), function(df) trimws(t(sapply(1:27*2 - 1, function(x) (strsplit(readLines(dataFiles[df]), split = "\t"))[[x]]))))

popnames <- ds[[1]][,1]
ds <- lapply(ds, function(x) x[,-1])
d <- (do.call(cbind, ds))
class(d) <- "numeric"
rownames(d) <- popnames

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

write.tree(bestTrees, paste0("output_discrete_jenks/maximumParsimonyTree_", fileName, "_", transData, "_jenks_", expCat, "_expcats.txt"))

q()