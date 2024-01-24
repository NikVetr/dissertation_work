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
dir.create("analyses/MCC")


cl <- makeCluster(8, outfile="")
registerDoParallel(cl)
getDoParWorkers()

foreach(m=1:100, .packages = c("ape", "coda", "rwty", "phangorn", "lattice", "latticeExtra", "MCMCpack", "readtext", "StatMatch", "mvtnorm")) %dopar% {

  mccDists <- matrix(nrow = length(degreeCorrelation)*length(traitNumIter), ncol = 6)
  rowCounter <- 0
  
  for(i in m:m){ #iterates over replicates
  
  # k iterates over degree of misspecification
  for (k in 1:length(degreeCorrelation)) {
    
    for (j in 1:length(traitNumIter)) { #iterates over trait number
      
        rowCounter <- rowCounter + 1
        ntraits <- traitNumIter[j]
        rateMatrixCorrelation <- degreeCorrelation[k]
        replicateNum <- i
        fileName <- paste0("simulations_", ntraits, "traits_perfectRateMatrix_", rateMatrixCorrelation, "correlationStrength_", noisyness, "_replicate_", replicateNum)
        
        #get true tree
        trueTree <- read.nexus(file = paste0("trees/", fileName, "_tree.nex"))
        
        #get mcmc output
        trees1 <- read.tree(file = paste0("output/", fileName, "_trees_run_1.trees"))
        trees2 <- read.tree(file = paste0("output/", fileName, "_trees_run_2.trees"))
        trees <- c(trees1, trees2)
        
        #compute summary trees
        mccTree <- maxCladeCred(trees)
        rfDistMCC <- dist.topo(x = trueTree, y = mccTree, method = "PH85")
        kfDistMCC <- dist.topo(x = trueTree, y = mccTree, method = "score")
  
        mccDists[rowCounter,] <- c(ntraits, rateMatrixCorrelation, "MCC", replicateNum, rfDistMCC, kfDistMCC)
        cat(paste0("iter",m,", ",ntraits," traits, ", rateMatrixCorrelation, " corr , ", rfDistMCC, " RF, ", kfDistMCC, " KF\n"))
      }
    }
  }
  mccDists <- as.data.frame(mccDists)
  colnames(mccDists) <- c("ntraits", "rateMatrixCorrelation", "type", "replicateNum", "rfDist", "kfDist")
  write.table(mccDists, file = paste0("analyses/MCC/", i,".txt"))
}

mccDists <- as.data.frame(matrix(nrow = length(degreeCorrelation)*length(traitNumIter)*100, ncol = 6))
colnames(mccDists) <- c("ntraits", "rateMatrixCorrelation", "type", "replicateNum", "rfDist", "kfDist")
for(i in 1:100){
  print(i)
  # print(c(((i-1)*length(degreeCorrelation)*length(traitNumIter)+1),((i*length(degreeCorrelation)*length(traitNumIter)))))
  mccDists[((i-1)*length(degreeCorrelation)*length(traitNumIter)+1):((i*length(degreeCorrelation)*length(traitNumIter))),1:ncol(mccDists)] <- read.table(file = paste0("analyses/MCC/",i,".txt"))
}
write.table(mccDists, file = paste0("analyses/mccDists.txt"))
read.table(file = "analyses/mccDists.txt")



