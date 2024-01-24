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

windowCodes <- seq(from = 1, to = 100, by = 1)

cl <- makeCluster(8, outfile="")
registerDoParallel(cl)
getDoParWorkers()

foreach(m=1:100, .packages = c("ape", "coda", "rwty", "phangorn", "lattice", "latticeExtra")) %dopar% {
  
  PairwiseDists <- array(dim = c(21, 5, 1000))
  
  traitNumIter <- c(2,4,8,16,32,64)
  
  degreeCorrelation <- c("2", "4", "3", "6", "8", "9")
  
  
  rowCounter <- 0
  # k iterates over degree of misspecification
  for (k in 1:length(degreeCorrelation)) {
    
    for (j in 1:6) { #iterates over trait number
      
      for(i in m:m){ #iterates over replicates
        
        rowCounter <- rowCounter + 1
        print(paste(round(rowCounter/21, digits = 2),j, i, k, "... ")) #read out percent done
        print(paste0(m, "% done"))
        
        ntraits <- traitNumIter[j]
        rateMatrixCorrelation <- degreeCorrelation[k]
        replicateNum <- i
        
        fileName <- paste0("simulations_", ntraits, "traits_perfectRateMatrix_", rateMatrixCorrelation, "correlationStrength_", "replicate_", replicateNum)
        trueTree <- read.nexus(file = paste0("trees/", fileName, "_tree.nex"))
        #plot(trueTree)
        
        print("reading trees")
        trees1 <- read.tree(file = paste0("output/", fileName, "_trees_run_1.trees"))
        trees2 <- read.tree(file = paste0("output/", fileName, "_trees_run_2.trees"))
        trees <- c(trees1, trees2)
        print("done reading")
        
        samples <- sapply(1:1000, function(x) sample(1:length(trees), size = 2, replace = F))
        RFdists <- sapply(1:1000, function(x) dist.topo(x = trees[[samples[1,x]]], y = trees[[samples[2,x]]], method = "PH85"))
        KFdists <- sapply(1:1000, function(x) dist.topo(x = trees[[samples[1,x]]], y = trees[[samples[2,x]]], method = "score"))
        
        PairwiseDists[rowCounter, 1,] <- rep(ntraits, 1000)
        PairwiseDists[rowCounter, 2,] <- rep(rateMatrixCorrelation, 1000)
        PairwiseDists[rowCounter, 3,] <- rep(replicateNum, 1000)
        PairwiseDists[rowCounter, 4,] <- RFdists
        PairwiseDists[rowCounter, 5,] <- KFdists
        
      }
    }
  }
  
  save(x = PairwiseDists, file = paste0("analyses/PairwiseDists", m))
  
}
stopCluster(cl)


#read back in the files
library(abind)
load(paste0("analyses/PairwiseDists", 1))
allDistances <- PairwiseDists
for(i in 2:100){
  print(i)
  load(paste0("analyses/PairwiseDists", i))
  allDistances <- abind(allDistances, PairwiseDists, along = 1)
}
PairwiseDists <- allDistances

# save(x = PairwiseDists, file = paste0("analyses/PairwiseDistsALL"))
load(paste0("analyses/PairwiseDistsALL"))


par(mfrow = c(1,1))
means <- PairwiseDists[,,1]
means[,4] <- sapply(1:2100, function (x) mean(as.numeric(PairwiseDists[x,4,])))
means[,5] <- sapply(1:2100, function (x) mean(as.numeric(PairwiseDists[x,5,])))
means <- as.data.frame(means)
colnames(means) <- c("ntraits", "rateMatrixCorrelation", "replicateNum", "rfDist", "kfDist")

plot(as.numeric(as.character(means$ntraits)), as.numeric(as.character(means$rfDist)), col = means$rateMatrixCorrelation, 
     xlab = "Number of Traits", ylab = "Mean RF-Distance", main = "Mean Pairwise RF-Distance in Joint Posterior")
legend(54, 54,levels(means$rateMatrixCorrelation),col=1:length(means$rateMatrixCorrelation),pch=1)

plot(as.numeric(as.character(means$ntraits)), as.numeric(as.character(means$kfDist)), col = means$rateMatrixCorrelation, 
     xlab = "Number of Traits", ylab = "Mean KF-Distance", main = "Mean Pairwise KF-Distance in Joint Posterior")
legend(54, 3.225,levels(means$rateMatrixCorrelation),col=1:length(means$rateMatrixCorrelation),pch=1)


#summarize mean values
meansSummary <- data.frame(matrix(nrow = 21, ncol = 6))
colnames(meansSummary) <- c("ntraits", "specification", "meanRF", "meanKF", "sdRF", "sdKF")

meansSummary[,1] <- as.numeric(as.character(means[1:21, 1]))
meansSummary[,2] <- means[1:21, 2]
# meansSummary <- rbind(meansSummary, c(0, NA, NA, NA, NA, NA)); meansSummary[22,2] <- "none"

for(i in 1:nrow(meansSummary)){
  meansSummary[i,3] <- mean(as.numeric(as.character(means[means[,1] == meansSummary[i,1] & means[,2] == meansSummary[i,2], 4])))
  meansSummary[i,4] <- mean(as.numeric(as.character(means[means[,1] == meansSummary[i,1] & means[,2] == meansSummary[i,2], 5])))
  meansSummary[i,5] <- sd(as.numeric(as.character(means[means[,1] == meansSummary[i,1] & means[,2] == meansSummary[i,2], 4])))
  meansSummary[i,6] <- sd(as.numeric(as.character(means[means[,1] == meansSummary[i,1] & means[,2] == meansSummary[i,2], 5])))
}
colnames(meansSummary) <- c("ntraits", "specification", "meanRF", "meanKF", "sdRF", "sdKF")

library(ggplot2)
ggplot(meansSummary, aes(x = ntraits, y = meanRF, colour = specification)) +
  geom_point(size = 2) + 
  geom_errorbar(aes(ymin = meanRF - sdRF, ymax = meanRF + sdRF)) + 
  theme_bw()

ggplot(meansSummary, aes(x = ntraits, y = meanKF, colour = specification)) +
  geom_point(size = 2) + 
  geom_errorbar(aes(ymin = meanKF - sdKF, ymax = meanKF + sdKF)) + 
  theme_bw()

