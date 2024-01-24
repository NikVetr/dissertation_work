setwd(dir = "/Volumes/2TB/mvBM_sims_Howells")
# setwd(dir = "/Volumes/2TB/uvBM_sims_Howells")

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

traitNumIter <- c(2,4,8,16,32,57)
degreeMisspecification <- c("perfect", "interregional", "chimpHuman")
dir.create("analyses/MCC")


cl <- makeCluster(8, outfile="")
registerDoParallel(cl)
getDoParWorkers()

foreach(m=1:100, .packages = c("ape", "coda", "rwty", "phangorn", "lattice", "latticeExtra", "MCMCpack", "readtext", "StatMatch", "mvtnorm")) %dopar% {
  
  mccDists <- matrix(nrow = length(degreeMisspecification)*length(traitNumIter), ncol = 6)
  rowCounter <- 0
  
  for(i in m:m){ #iterates over replicates
    
    # k iterates over degree of misspecification
    for (k in 1:length(degreeMisspecification)) {
      
      for (j in 1:length(traitNumIter)) { #iterates over trait number
        
        rowCounter <- rowCounter + 1
        ntraits <- traitNumIter[j]
        rateMatrixMisspecification <- degreeMisspecification[k]
        replicateNum <- i
        fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixMisspecification, "RateMatrix_replicate_", replicateNum)
        
        #get true tree
        trueTree <- read.nexus(file = paste0("trees/", fileName, "_tree.nex"))
        
        #get mcmc output
        trees1 <- read.tree(file = paste0("output/", fileName, "_trees_run_1_thinned.trees"))
        trees2 <- read.tree(file = paste0("output/", fileName, "_trees_run_2_thinned.trees"))
        trees <- c(trees1, trees2)
        
        #compute summary trees
        mccTree <- maxCladeCred(trees)
        rfDistMCC <- dist.topo(x = trueTree, y = mccTree, method = "PH85")
        kfDistMCC <- dist.topo(x = trueTree, y = mccTree, method = "score")
        
        mccDists[rowCounter,] <- c(ntraits, rateMatrixMisspecification, "MCC", replicateNum, rfDistMCC, kfDistMCC)
        cat(paste0("iter",m,", ",ntraits," traits, ", rateMatrixMisspecification, " corr , ", rfDistMCC, " RF, ", kfDistMCC, " KF\n"))
      }
    }
  }
  mccDists <- as.data.frame(mccDists)
  colnames(mccDists) <- c("ntraits", "rateMatrixMisspecification", "type", "replicateNum", "rfDist", "kfDist")
  write.table(mccDists, file = paste0("analyses/MCC/", i,".txt"))
}

mccDists <- as.data.frame(matrix(nrow = length(degreeMisspecification)*length(traitNumIter)*100, ncol = 6))
colnames(mccDists) <- c("ntraits", "rateMatrixMisspecification", "type", "replicateNum", "rfDist", "kfDist")
for(i in 1:100){
  print(i)
  # print(c(((i-1)*length(degreeMisspecification)*length(traitNumIter)+1),((i*length(degreeMisspecification)*length(traitNumIter)))))
  mccDists[((i-1)*length(degreeMisspecification)*length(traitNumIter)+1):((i*length(degreeMisspecification)*length(traitNumIter))),1:ncol(mccDists)] <- read.table(file = paste0("analyses/MCC/",i,".txt"))
}
write.table(mccDists, file = paste0("analyses/mccDists.txt"))
read.table(file = "analyses/mccDists.txt")[read.table(file = "analyses/mccDists.txt")[,1] == 57,]
read.table(file = "analyses/mccDists.txt")


