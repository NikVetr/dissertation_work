setwd(dir = "/Users/nikolai/mvBM_sims_all")

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

windowCodes <- seq(from = 1, to = 100, by = 1)

cl <- makeCluster(8, outfile="")
registerDoParallel(cl)
getDoParWorkers()
stopCluster(cl)

foreach(m=1:100, .packages = c("ape", "coda", "rwty", "phangorn", "lattice", "latticeExtra")) %dopar% {

MCMC_Diagnostics_PSRF <- matrix(ncol = 7)
colnames(MCMC_Diagnostics_PSRF) <- c("ntraits", "rateMatrixSpecification", "replicateNum", "TL_Point", "TL_UpperCI", "LL_Point", "LL_UpperCI")
MCMC_Diagnostics_Geweke <- matrix(ncol = 7)
colnames(MCMC_Diagnostics_Geweke) <- c("ntraits", "rateMatrixSpecification", "replicateNum", "gTL1", "gTL2", "gLL1", "gLL2")
MCMC_Diagnostics_ESS <- matrix(ncol = 7)
colnames(MCMC_Diagnostics_ESS) <- c("ntraits", "rateMatrixSpecification", "replicateNum", "essTL1", "essTL2", "essLL1", "essLL2")
# MCMC_Diagnostics_ASDSF <- matrix()
MCMC_Diagnostics_NodalFreqCorr <- matrix(ncol = 5)
colnames(MCMC_Diagnostics_NodalFreqCorr) <- c("ntraits", "rateMatrixSpecification", "replicateNum", "correlations", "numBipartsNotInChain")

Tree_Dist <- matrix(ncol = 5)
colnames(Tree_Dist) <- c("ntraits", "rateMatrixSpecification", "replicateNum", "rfDist", "kfDist")
MJCT_Compatibility <- matrix(ncol = 5)
colnames(MJCT_Compatibility) <- c("ntraits", "rateMatrixSpecification", "replicateNum", "numInconsistencies", "conTreeResolution")

windowCode <- windowCodes[m]

#Also output marginal histograms and trace plots for continuous parameters, split	frequencies	&	presence/absence for trees

traitNumIter <- c(2,4,8,16,32,57,64)

degreeMisspecification <- c("perfect", "interregional", "chimpHuman")


# k iterates over degree of misspecification
for (k in 1:length(degreeMisspecification)) {

for (j in 1:7) {
  
  for(i in windowCodes[m]:(windowCodes[m])){
    
    cat(c(j, i, k, "... \n"))
    
    ntraits <- traitNumIter[j]
    rateMatrixSpecification <- degreeMisspecification[k]
    replicateNum <- i
    
    fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_", "replicate_", replicateNum)
    trueTree <- read.nexus(file = paste0("/Users/nikolai/mvBM_sims_all/trees/", fileName, "_tree.nex"))
    
    # continuous character diagnostics
    log1 <- read.table(file = paste0("/Users/nikolai/mvBM_sims_all/output/", fileName, "_params_run_1.log"), header = T)
    # log1 <- log1[333:1667,] #discard burn-in
    log1 <- log1[seq(1, dim(log1)[1], by = 2),] #thin further
    log1 <- log1[, c("Likelihood", "TL")]
    TL1 <- mcmc(log1$TL, start = 20001, thin = 30)
    LL1 <- mcmc(log1$Likelihood, start = 20001, thin = 30)
    
    log2 <- read.table(file = paste0("/Users/nikolai/mvBM_sims_all/output/", fileName, "_params_run_2.log"), header = T)
    # log2 <- log2[333:1667,] #discard burn-in
    log2 <- log2[seq(1, dim(log2)[1], by = 2),] #thin further
    log2 <- log2[, c("Likelihood", "TL")]
    TL2 <- mcmc(log2$TL, start = 20001, thin = 30)
    LL2 <- mcmc(log2$Likelihood, start = 20001, thin = 30)
    
    TLs <- mcmc.list(TL1, TL2)
    LLs <- mcmc.list(LL1, LL2)
    
    cat("processing continuous data... ")
    
    trueTL <- sum(trueTree$edge.length)
    
    png(paste0("/Users/nikolai/mvBM_sims_all/plots/", fileName, "_TraceMarginals.png"), width = 1920, height = 1080)
    par(mfrow = c(2,2))
    traceplot(TLs, main = "TLs"); abline(h = trueTL, lwd = 2); traceplot(LLs, main = "LLs")
    hist(x = TL1, breaks = 50, col=rgb(1,0,0,0.3), main = "TLs"); hist(x = TL2, breaks = 50, col=rgb(0,0,1,0.3), add = T); abline(v = trueTL, lwd = 2)
    hist(x = LL1, breaks = 50, col=rgb(1,0,0,0.3), main = "LLs"); hist(x = LL2, breaks = 50, col=rgb(0,0,1,0.3), add = T)
    dev.off()
    
    
    TLsPSRF <- gelman.diag(TLs)$psrf
    LLsPSRF <- gelman.diag(LLs)$psrf
    
    MCMC_Diagnostics_PSRF <- rbind(MCMC_Diagnostics_PSRF, c(ntraits, rateMatrixSpecification, replicateNum, TLsPSRF, LLsPSRF))
    
    gTL1 <- geweke.diag(TL1)[1]$z
    gTL2 <- geweke.diag(TL2)[1]$z
    gLL1 <- geweke.diag(LL1)[1]$z
    gLL2 <- geweke.diag(LL2)[1]$z
    
    MCMC_Diagnostics_Geweke <- rbind(MCMC_Diagnostics_Geweke, c(ntraits, rateMatrixSpecification, replicateNum, gTL1, gTL2, gLL1, gLL2))
    
    essTL1 <- effectiveSize(TL1)[[1]]
    essTL2 <- effectiveSize(TL2)[[1]]
    essLL1 <- effectiveSize(LL1)[[1]]
    essLL2 <- effectiveSize(LL2)[[1]]
    
    MCMC_Diagnostics_ESS <- rbind(MCMC_Diagnostics_ESS, c(ntraits, rateMatrixSpecification, replicateNum, essTL1, essTL2, essLL1, essLL2))
    
    cat("reading trees... ")
    
    # tree diagnostics
    trees1 <- read.tree(file = paste0("/Users/nikolai/mvBM_sims_all/output/", fileName, "_trees_run_1.trees"))
    # trees1 <- trees1[333:1667] #discard burn-in
    trees1 <- trees1[seq(1, length(trees1), by = 2)] #thin further
    
    trees2 <- read.tree(file = paste0("/Users/nikolai/mvBM_sims_all/output/", fileName, "_trees_run_2.trees"))
    # trees2 <- trees2[333:1667] #discard burn-in
    trees2 <- trees2[seq(1, length(trees2), by = 2)] #thin further
    
    #     write.tree(trees1, file = paste0("trees1rwty", windowCode))
    #     write.tree(trees2, file = paste0("trees2rwty", windowCode))
    #     write.table(x = rbind(c(NA, NA),log1), file = paste0("log1rwty", windowCode))
    #     write.table(x = rbind(c(NA, NA),log2), file = paste0("log2rwty", windowCode))
    #     
    #     trees1rwty <- load.trees(file = paste0("A:\\mvBM_sims\\trees1rwty", windowCode), type = "newick", gens.per.tree = 40, logfile = paste0("A:\\mvBM_sims\\log1rwty", windowCode))
    #     trees2rwty <- load.trees(file = paste0("A:\\mvBM_sims\\trees2rwty", windowCode), type = "newick", gens.per.tree = 40, logfile = paste0("A:\\mvBM_sims\\log2rwty", windowCode))
    # 
    #     png(paste0("A:\\mvBM_sims\\plots\\", fileName, "_cumulSplitFreqs.png"), width = 1920, height = 1080)
    #     makeplot.splitfreqs.cumulative(list(trees1rwty, trees2rwty), n.clades = 30)
    #     dev.off()
    
    cat("processing trees...\n")
    
    # assessing accurate retrieval of data generating tree    
    trees <- c(trees1, trees2)
    mccTree <- maxCladeCred(trees,rooted=FALSE)
    conTree <- consensus(trees, p=.75)
    trueTree <- read.nexus(file = paste0("/Users/nikolai/mvBM_sims_all/trees/", fileName, "_tree.nex"))
    
    if(conTree$Nnode != 1) {
      numInconsistensies <- sum(is.na(prop.clades(phy = conTree, unroot(trueTree))))} else {
        numInconsistensies <- 0
      }
    
    MJCT_Compatibility <- rbind(MJCT_Compatibility, c(ntraits, rateMatrixSpecification, replicateNum, numInconsistensies, conTree$Nnode))
    
    rfDist <- dist.topo(x = trueTree, y = mccTree, method = "PH85")
    kfDist <- dist.topo(x = trueTree, y = mccTree, method = "score")
    Tree_Dist <- rbind(Tree_Dist, c(ntraits, rateMatrixSpecification, replicateNum, rfDist, kfDist))
    
    cat("calculating prop. clades...\n")
    
    pp1 <- prop.clades(mccTree, trees1, rooted = F) / length(trees1)
    pp2 <- prop.clades(mccTree, trees2, rooted = F) / length(trees2)
    corT1T2 <- cor(pp1, pp2, use = "com")
    numBipartsNotInChain1orChain2 <- sum(is.na(c(pp1, pp2)))
    
    MCMC_Diagnostics_NodalFreqCorr <- rbind(MCMC_Diagnostics_NodalFreqCorr, c(ntraits, rateMatrixSpecification, replicateNum, corT1T2, numBipartsNotInChain1orChain2))
    
    }
  }
  }

  #for single run
  cat("writing data to disk...\n")
  
  write.table(x = MCMC_Diagnostics_PSRF, file = paste0("MCMC_Diagnostics_PSRF", windowCode))
  write.table(x = MCMC_Diagnostics_Geweke, file = paste0("MCMC_Diagnostics_Geweke", windowCode))
  write.table(x = MCMC_Diagnostics_ESS, file = paste0("MCMC_Diagnostics_ESS", windowCode))
  write.table(x = MCMC_Diagnostics_NodalFreqCorr, file = paste0("MCMC_Diagnostics_NodalFreqCorr", windowCode))
  write.table(x = Tree_Dist, file = paste0("Tree_Dist", windowCode))
  write.table(x = MJCT_Compatibility, file = paste0("MJCT_Compatibility", windowCode))
  
}

#read back in the files
MCMC_Diagnostics_PSRF <- read.table(file = "analyses/MCMC_Diagnostics_PSRF1")[-1,]
MCMC_Diagnostics_Geweke <- read.table(file = "analyses/MCMC_Diagnostics_Geweke1")[-1,]
MCMC_Diagnostics_ESS <- read.table(file = "analyses/MCMC_Diagnostics_ESS1")[-1,]
MCMC_Diagnostics_NodalFreqCorr <- read.table(file = "analyses/MCMC_Diagnostics_NodalFreqCorr1")[-1,]
Tree_Dist <- read.table(file = "analyses/Tree_Dist1")[-1,]
MJCT_Compatibility <- read.table(file = "analyses/MJCT_Compatibility1")[-1,]

for (i in 2:100){
  print(i)
  MCMC_Diagnostics_PSRF <- rbind(MCMC_Diagnostics_PSRF, read.table(file = paste0("analyses/MCMC_Diagnostics_PSRF", i))[-1,])
  MCMC_Diagnostics_Geweke <- rbind(MCMC_Diagnostics_Geweke, read.table(file = paste0("analyses/MCMC_Diagnostics_Geweke", i))[-1,])
  MCMC_Diagnostics_ESS <- rbind(MCMC_Diagnostics_ESS, read.table(file = paste0("analyses/MCMC_Diagnostics_ESS", i))[-1,])
  MCMC_Diagnostics_NodalFreqCorr <- rbind(MCMC_Diagnostics_NodalFreqCorr, read.table(file = paste0("analyses/MCMC_Diagnostics_NodalFreqCorr", i))[-1,])
  Tree_Dist <- rbind(Tree_Dist, read.table(file = paste0("analyses/Tree_Dist", i))[-1,])
  MJCT_Compatibility <- rbind(MJCT_Compatibility, read.table(file = paste0("analyses/MJCT_Compatibility", i))[-1,])
}

###################################################################
#### OK #### NOW #### THAT #### WE #### ARE #### ALL #### HERE ####
###################################################################


panel.hanoi <- function(x, y, horizontal, breaks="Sturges", ...) {  # "Sturges" is hist()'s default
  
  if (horizontal) {
    condvar <- y # conditioning ("independent") variable
    datavar <- x # data ("dependent") variable
  } else {
    condvar <- x
    datavar <- y
  }
  
  conds <- sort(unique(condvar))
  
  # loop through the possible values of the conditioning variable
  for (i in seq_along(conds)) {
    
    h <- hist(datavar[condvar == conds[i]], plot=F, breaks) # use base hist(ogram) function to extract some information
    
    # strip outer counts == 0, and corresponding bins
    brks.cnts <- stripOuterZeros(h$breaks, h$counts)
    brks <- brks.cnts[[1]]
    cnts <- brks.cnts[[2]]
    
    halfrelfs <- (cnts/sum(cnts))/2  # i.e. half of the relative frequency
    center <- i
    
    # All of the variables passed to panel.rec will usually be vectors, and panel.rect will therefore make multiple rectangles.
    if (horizontal) {
      panel.rect(head(brks, -1), center - halfrelfs, tail(brks, -1), center + halfrelfs, ...)
    } else {
      panel.rect(center - halfrelfs, head(brks, -1), center + halfrelfs, tail(brks, -1), ...)
    }
  }
}

# function to strip counts that are all zero on ends of data, along with the corresponding breaks
stripOuterZeros <- function(brks, cnts) { do.call("stripLeftZeros", stripRightZeros(brks, cnts)) }

stripLeftZeros <- function(brks, cnts) {
  if (cnts[1] == 0) {
    stripLeftZeros(brks[-1], cnts[-1])
  } else {
    list(brks, cnts)
  }
}

stripRightZeros <- function(brks, cnts) {
  len <- length(cnts)
  if (cnts[len] ==0) {
    stripRightZeros(brks[-(len+1)], cnts[-len])
  } else {
    list(brks, cnts)
  }
}


Tree_Dist <- as.data.frame(Tree_Dist)
Tree_Dist$ntraits <- as.character(Tree_Dist$ntraits)
Tree_Dist$ntraits <- factor(Tree_Dist$ntraits, levels = c("2", "4", "8", "16", "32", "57", "64"))

bwplot(rfDist ~ ntraits, ylab = "Robinson-Foulds Distance", main = "Comparing inferred MCC tree to Data-Generating Tree", xlab = "Number of Traits Used in Simulation", data=Tree_Dist, pch="|", coef=0, panel=function(...){panel.hanoi(col="orange", ...); panel.bwplot(...)})
bwplot(kfDist ~ ntraits, ylab = "Kuhner-Felsenstein Distance", main = "Comparing inferred MCC tree to Data-Generating Tree", xlab = "Number of Traits Used in Simulation", data=Tree_Dist, pch="|", coef=0, panel=function(...){panel.hanoi(col="orange", ...); panel.bwplot(...)})

MJCT_Compatibility <- as.data.frame(MJCT_Compatibility)
MJCT_Compatibility$ntraits <- as.character(MJCT_Compatibility$ntraits)
MJCT_Compatibility$ntraits <- factor(MJCT_Compatibility$ntraits, levels = c("2", "4", "8", "16", "32", "57", "64"))
MJCT_Compatibility$conTreeResolution <- as.numeric(as.character(MJCT_Compatibility$conTreeResolution))
MJCT_Compatibility$numInconsistencies <- as.numeric(as.character(MJCT_Compatibility$numInconsistencies))
MJCT_Compatibility$numCorrectNodes <- MJCT_Compatibility$conTreeResolution - MJCT_Compatibility$numInconsistencies

key <- list(x = 1.1, y = .92, corner = c(5.5, 0),
            text = list(c("Correctly", "Incorrectly")),
            lines = list(type = c("l", "l"), col = c("green", "red"),
                         pch = 3, lwd = 8, lty = 1))

a <- bwplot(numCorrectNodes ~ ntraits, main = "Consistency with MJ-Consensus Tree -- ALL RATE MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in Strict Majority-Rule Consensus Tree", 
            data=MJCT_Compatibility, ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies ~ ntraits, data=MJCT_Compatibility, ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

a <- bwplot(numCorrectNodes ~ ntraits, main = "Consistency with MJ-Consensus Tree -- PERFECT MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in Strict Majority-Rule Consensus Tree", 
            data=MJCT_Compatibility[MJCT_Compatibility[,2] == "perfect",], ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies ~ ntraits, data=MJCT_Compatibility[MJCT_Compatibility[,2] == "perfect",], ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

a <- bwplot(numCorrectNodes ~ ntraits, main = "Consistency with MJ-Consensus Tree -- INTERREGIONAL MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in Strict Majority-Rule Consensus Tree", 
            data=MJCT_Compatibility[MJCT_Compatibility[,2] == "interregional",], ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies ~ ntraits, data=MJCT_Compatibility[MJCT_Compatibility[,2] == "interregional",], ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

a <- bwplot(numCorrectNodes ~ ntraits, main = "Consistency with MJ-Consensus Tree -- CHIMP-HUMAN MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in Strict Majority-Rule Consensus Tree", 
            data=MJCT_Compatibility[MJCT_Compatibility[,2] == "chimpHuman",], ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies ~ ntraits, data=MJCT_Compatibility[MJCT_Compatibility[,2] == "chimpHuman",], ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

####################################
## ASSESS TOPOLOGICAL CONVERGENCE ##
####################################

windowCodes <- seq(from = 1, to = 100, by = 1)

cl <- makeCluster(8, outfile="")
registerDoParallel(cl)
getDoParWorkers()

foreach(m=1:100, .packages = c("ape", "coda", "rwty", "phangorn", "lattice", "latticeExtra")) %dopar% {
  
  MCMC_Diagnostics_Distances <- array(dim = c(21, 6, 26668))
  
  traitNumIter <- c(2,4,8,16,32,57,64)
  
  degreeMisspecification <- c("perfect", "interregional", "chimpHuman")
  
  
  rowCounter <- 0
  # k iterates over degree of misspecification
  for (k in 1:length(degreeMisspecification)) {
    
    for (j in 1:7) { #iterates over trait number
      
      for(i in m:m){ #iterates over replicates
        
        rowCounter <- rowCounter + 1
        print(paste(round(rowCounter/21, digits = 2),j, i, k, "... ")) #read out percent done
        print(paste0(m, "% done"))
        
        ntraits <- traitNumIter[j]
        rateMatrixSpecification <- degreeMisspecification[k]
        replicateNum <- i
        
        fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_", "replicate_", replicateNum)
        trueTree <- read.nexus(file = paste0("/Users/nikolai/mvBM_sims_all/trees/", fileName, "_tree.nex"))
        #plot(trueTree)
        
        print("reading trees")
        trees1 <- read.tree(file = paste0("/Users/nikolai/mvBM_sims_all/output/", fileName, "_trees_run_1.trees"))
        trees2 <- read.tree(file = paste0("/Users/nikolai/mvBM_sims_all/output/", fileName, "_trees_run_2.trees"))
        print("done reading")
        
        RFdists1 <- vector()
        RFdists2 <- vector()
        KFdists1 <- vector()
        KFdists2 <- vector()
        
        for(l in 1:length(trees1)){
          # print(l)
          dist1 <- treedist(trueTree, trees1[[l]])
          RFdists1[l] <- dist1[[1]]
          KFdists1[l] <- dist1[[2]]
          dist2 <- treedist(trueTree, trees2[[l]])
          RFdists2[l] <- dist2[[1]]
          KFdists2[l] <- dist2[[2]]
        }
        
        MCMC_Diagnostics_Distances[rowCounter, 1,] <- rep(ntraits, 26668)
        MCMC_Diagnostics_Distances[rowCounter, 2,] <- rep(rateMatrixSpecification, 26668)
        MCMC_Diagnostics_Distances[rowCounter, 3,] <- rep(replicateNum, 26668)
        MCMC_Diagnostics_Distances[rowCounter, 4,] <- c(rep(1,13334), rep(2,13334))
        MCMC_Diagnostics_Distances[rowCounter, 5,] <- c(RFdists1, RFdists2)
        MCMC_Diagnostics_Distances[rowCounter, 6,] <- c(KFdists1, KFdists2)
        
        RFdists1 <- mcmc(RFdists1, start = 1)
        RFdists2 <- mcmc(RFdists2, start = 1)
        KFdists1 <- mcmc(KFdists1, start = 1)
        KFdists2 <- mcmc(KFdists2, start = 1)
  
        RFs <- mcmc.list(RFdists1, RFdists2)
        KFs <- mcmc.list(KFdists1, KFdists2)
        
        png(paste0("/Users/nikolai/mvBM_sims_all/plots/treeDist/", fileName, "_TreeDistTraceMarginals.png"), width = 1920, height = 1080)
        par(mfrow = c(2,2))
        traceplot(RFs, main = "RFs"); traceplot(KFs, main = "KFs")
        hist(x = RFdists1, breaks = 50, col=rgb(1,0,0,0.3), main = "RFs from True Tree"); hist(x = RFdists2, breaks = 50, col=rgb(0,0,1,0.3), add = T)
        hist(x = KFdists1, breaks = 50, col=rgb(1,0,0,0.3), main = "KFs from True Tree"); hist(x = KFdists2, breaks = 50, col=rgb(0,0,1,0.3), add = T)
        dev.off()
        
      }
    }
  }
  
  save(x = MCMC_Diagnostics_Distances, file = paste0("analyses/MCMC_Diagnostics_Distances", m))
  
}
stopCluster()



#DOING THE ANALYTIC PRIOR FOR ALL THE REPLICATES
windowCodes <- seq(from = 1, to = 100, by = 1)

cl <- makeCluster(8, outfile="")
registerDoParallel(cl)
getDoParWorkers()

foreach(m=1:100, .packages = c("ape", "coda", "rwty", "phangorn", "lattice", "latticeExtra", "MCMCpack")) %dopar% {
  
  gc()
  
  MCMC_Diagnostics_Distances_prior <- array(dim = c(1, 6, 26668))
  
  traitNumIter <- c(2,4,8,16,32,57,64)
  
  degreeMisspecification <- "perfect"
  
  
  rowCounter <- 0
  # k iterates over degree of misspecification
  for (k in 1:length(degreeMisspecification)) {
    
    for (j in 1:1) { #iterates over trait number
      
      for(i in m:m){ #iterates over replicates
        
        rowCounter <- rowCounter + 1
        print(paste(round(rowCounter/7, digits = 2),j, i, k, "... ")) #read out percent done
        print(paste0(m, "% done"))
        
        ntraits <- traitNumIter[j]
        rateMatrixSpecification <- degreeMisspecification[k]
        replicateNum <- i
        
        fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_", "replicate_", replicateNum)
        trueTree <- read.nexus(file = paste0("/Users/nikolai/mvBM_sims_all/trees/", fileName, "_tree.nex"))
        #plot(trueTree)
        
        RFdists1 <- vector()
        RFdists2 <- vector()
        KFdists1 <- vector()
        KFdists2 <- vector()
        
        for(l in 1:13334){
          # print(l)
          tree1 <- rtree(n=30)
          tree1$edge.length <- 10^runif(n = 1, min = -3, max = 2) * rdirichlet(n = 1, alpha = rep(1,58))
          tree1$tip.label <- sample(trueTree$tip.label, size = 30, replace = F)
          dist1 <- treedist(trueTree, tree1)
          RFdists1[l] <- dist1[[1]]
          KFdists1[l] <- dist1[[2]]
          tree2 <- rtree(n=30)
          tree2$tip.label <- sample(trueTree$tip.label, size = 30, replace = F)
          tree2$edge.length <- 10^runif(n = 1, min = -3, max = 2) * rdirichlet(n = 1, alpha = rep(1,58))
          dist2 <- treedist(trueTree, tree2)
          RFdists2[l] <- dist2[[1]]
          KFdists2[l] <- dist2[[2]]        }
        
        MCMC_Diagnostics_Distances_prior[rowCounter, 1,] <- rep(ntraits, 26668)
        MCMC_Diagnostics_Distances_prior[rowCounter, 2,] <- rep(rateMatrixSpecification, 26668)
        MCMC_Diagnostics_Distances_prior[rowCounter, 3,] <- rep(replicateNum, 26668)
        MCMC_Diagnostics_Distances_prior[rowCounter, 4,] <- c(rep(1,13334), rep(2,13334))
        MCMC_Diagnostics_Distances_prior[rowCounter, 5,] <- c(RFdists1, RFdists2)
        MCMC_Diagnostics_Distances_prior[rowCounter, 6,] <- c(KFdists1, KFdists2)
        
        
      }
    }
  }
  
  save(x = MCMC_Diagnostics_Distances_prior, file = paste0("analyses/MCMC_Diagnostics_Distances_prior", m))
  
}
stopCluster(cl)
  
# save(x = MCMC_Diagnostics_Distances, file = paste0("analyses/MCMC_Diagnostics_Distances", m))
save(x = MCMC_Diagnostics_Distances_prior, file = "MCMC_Diagnostics_Distances_prior")


stopCluster(cl)

# save(x = MCMC_Diagnostics_Distances, file = "MCMC_Diagnostics_Distances_7TipAnalyses")
load("MCMC_Diagnostics_Distances_7TipAnalyses")
load("MCMC_Diagnostics_Distances_prior")

library(abind)
load(paste0("analyses/MCMC_Diagnostics_Distances", 1))
allDistances <- MCMC_Diagnostics_Distances
for(i in 2:100){
  print(i)
  load(paste0("analyses/MCMC_Diagnostics_Distances", i))
  allDistances <- abind(allDistances, MCMC_Diagnostics_Distances, along = 1)
}
MCMC_Diagnostics_Distances <- allDistances

library(abind)
load(paste0("analyses/MCMC_Diagnostics_Distances_prior", 1))
allDistances <- MCMC_Diagnostics_Distances_prior
for(i in 2:100){
  print(i)
  str(allDistances)
  load(paste0("analyses/MCMC_Diagnostics_Distances_prior", i))
  allDistances <- abind(allDistances, MCMC_Diagnostics_Distances_prior, along = 1)
}
MCMC_Diagnostics_Distances_prior <- allDistances

save(x = MCMC_Diagnostics_Distances_prior, file = "MCMC_Diagnostics_Distances_prior")
load("MCMC_Diagnostics_Distances_7TipAnalyses")



###doing the analytic prior for just 100 of the replicates

MCMC_Diagnostics_Distances_prior <- array(dim = c(100, 6, 26668))

traitNumIter <- c(2,4,8,16,32,57,64)

degreeMisspecification <- "perfect"


rowCounter <- 0
# k iterates over degree of misspecification
for (k in 1:length(degreeMisspecification)) {
  
  for (j in 1:1) { #iterates over trait number
    
    for(i in 1:100){ #iterates over replicates
      
      print(paste0(i, "% done"))
      
      ntraits <- traitNumIter[j]
      rateMatrixSpecification <- degreeMisspecification[k]
      replicateNum <- i
      
      fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_", "replicate_", replicateNum)
      trueTree <- read.nexus(file = paste0("/Users/nikolai/mvBM_sims_all/trees/", fileName, "_tree.nex"))
      #plot(trueTree)
      
      RFdists1 <- vector()
      RFdists2 <- vector()
      KFdists1 <- vector()
      KFdists2 <- vector()
      
      for(l in 1:13334){
        # print(paste0(i, ": ", l/1443))
        tree1 <- rtree(n=30)
        tree1$edge.length <- 10^runif(n = 1, min = -3, max = 2) * rdirichlet(n = 1, alpha = rep(1,58))
        tree1$tip.label <- sample(trueTree$tip.label, size = 30, replace = F)
        dist1 <- treedist(trueTree, tree1)
        RFdists1[l] <- dist1[[1]]
        KFdists1[l] <- dist1[[2]]
        tree2 <- rtree(n=30)
        tree2$tip.label <- sample(trueTree$tip.label, size = 30, replace = F)
        tree2$edge.length <- 10^runif(n = 1, min = -3, max = 2) * rdirichlet(n = 1, alpha = rep(1,58))
        dist2 <- treedist(trueTree, tree2)
        RFdists2[l] <- dist2[[1]]
        KFdists2[l] <- dist2[[2]]        }
      
      MCMC_Diagnostics_Distances_prior[i, 1,] <- rep(ntraits, 26668)
      MCMC_Diagnostics_Distances_prior[i, 2,] <- rep(rateMatrixSpecification, 26668)
      MCMC_Diagnostics_Distances_prior[i, 3,] <- rep(replicateNum, 26668)
      MCMC_Diagnostics_Distances_prior[i, 4,] <- c(rep(1,13334), rep(2,13334))
      MCMC_Diagnostics_Distances_prior[i, 5,] <- c(RFdists1, RFdists2)
      MCMC_Diagnostics_Distances_prior[i, 6,] <- c(KFdists1, KFdists2)
      
      
    }
  }
}

# save(x = MCMC_Diagnostics_Distances_prior, file = "analyses/MCMC_Diagnostics_Distances_prior")
# save(x = MCMC_Diagnostics_Distances, file = "analyses/MCMC_Diagnostics_Distances_All_inclPrior")




load("MCMC_Diagnostics_Distances_All")
load("analyses/MCMC_Diagnostics_Distances_prior")
for(i in 1:100){ print(i); MCMC_Diagnostics_Distances_prior[i,1,] <- 0}
for(i in 1:100){ print(i); MCMC_Diagnostics_Distances_prior[i,2,] <- "none"}
MCMC_Diagnostics_Distances <- abind(MCMC_Diagnostics_Distances, MCMC_Diagnostics_Distances_prior, along = 1)























load("analyses/MCMC_Diagnostics_Distances_All_inclPrior")
MCMC_Diagnostics_Distances <- MCMC_Diagnostics_Distances[,,seq(0, length(MCMC_Diagnostics_Distances[1,1,]), 5)]
# save(x = MCMC_Diagnostics_Distances, file = "analyses/MCMC_Diagnostics_Distances_All_inclPrior_thinned")
load("analyses/MCMC_Diagnostics_Distances_All_inclPrior_thinned")


par(mfrow = c(1,1))
means <- MCMC_Diagnostics_Distances[,,1]
means[,5] <- sapply(1:2200, function (x) mean(as.numeric(MCMC_Diagnostics_Distances[x,5,])))
means[,6] <- sapply(1:2200, function (x) mean(as.numeric(MCMC_Diagnostics_Distances[x,6,])))
means <- as.data.frame(means)
colnames(means) <- c("ntraits", "rateMatrixSpecification", "replicateNum", "idk", "rfDist", "kfDist")

means <- rbind(means, distMatrixNJDists)
means <- rbind(means, distMatrixUPGMADists)

plot(as.numeric(as.character(means$ntraits)), as.numeric(as.character(means$rfDist)), col = means$rateMatrixSpecification, 
     xlab = "Number of Traits", ylab = "Mean RF-Distance", main = "Mean RF-Dist from True, Data-Generating Tree")
legend(54, 55,levels(means$rateMatrixSpecification),col=1:length(means$rateMatrixSpecification),pch=1)

plot(as.numeric(as.character(means$ntraits)), as.numeric(as.character(means$kfDist)), col = means$rateMatrixSpecification, 
     xlab = "Number of Traits", ylab = "Mean KF-Distance", main = "Mean KF-Dist from True, Data-Generating Tree")
legend(54, 3.1,levels(means$rateMatrixSpecification),col=1:length(means$rateMatrixSpecification),pch=1)


#summarize mean values
meansSummary <- data.frame(matrix(nrow = 21, ncol = 6))
meansSummary[,1] <- as.numeric(as.character(means[1:21, 1]))
meansSummary[,2] <- means[1:21, 2]
meansSummary <- rbind(meansSummary, c(0, NA, NA, NA, NA, NA)); meansSummary[22,2] <- "none"
distanceMeansSummary <- data.frame(matrix(nrow = 14, ncol = 6))
distanceMeansSummary[,1] <- as.numeric(as.character(means[1:14, 1]))
distanceMeansSummary[,2] <- c(rep("NJ", 7), rep("UPGMA", 7))
meansSummary <- rbind(meansSummary, distanceMeansSummary)

for(i in 1:nrow(meansSummary)){
  meansSummary[i,3] <- mean(as.numeric(as.character(means[means[,1] == meansSummary[i,1] & means[,2] == meansSummary[i,2], 5])))
  meansSummary[i,4] <- mean(as.numeric(as.character(means[means[,1] == meansSummary[i,1] & means[,2] == meansSummary[i,2], 6])))
  meansSummary[i,5] <- sd(as.numeric(as.character(means[means[,1] == meansSummary[i,1] & means[,2] == meansSummary[i,2], 5])))
  meansSummary[i,6] <- sd(as.numeric(as.character(means[means[,1] == meansSummary[i,1] & means[,2] == meansSummary[i,2], 6])))
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

1+1




sds <- MCMC_Diagnostics_Distances[,,1]
sds[,5] <- sapply(1:2200, function (x) sd(as.numeric(MCMC_Diagnostics_Distances[x,5,])))
sds[,6] <- sapply(1:2200, function (x) sd(as.numeric(MCMC_Diagnostics_Distances[x,6,])))
sds <- as.data.frame(sds)
colnames(sds) <- c("ntraits", "rateMatrixSpecification", "replicateNum", "idk", "rfDist", "kfDist")
plot(as.numeric(as.character(sds$ntraits)), as.numeric(as.character(sds$rfDist)), col = sds$rateMatrixSpecification, 
     xlab = "Number of Traits", ylab = "Standard Deviation of RF-Distances", main = "Standard Deviation of RF-Dists from True, Data-Generating Tree")
legend(54, 2.1,levels(sds$rateMatrixSpecification),col=1:length(sds$rateMatrixSpecification),pch=1)

plot(as.numeric(as.character(sds$ntraits)), as.numeric(as.character(sds$kfDist)), col = sds$rateMatrixSpecification, 
     xlab = "Number of Traits", ylab = "Standard Deviation of KF-Distances", main = "Standard Deviation of KF-Dists from True, Data-Generating Tree")
legend(54, 3.2,levels(sds$rateMatrixSpecification),col=1:length(sds$rateMatrixSpecification),pch=1)

#summarize sd values
sdsSummary <- data.frame(matrix(nrow = 21, ncol = 6))
sdsSummary[,1] <- as.numeric(as.character(sds[1:21, 1]))
sdsSummary[,2] <- sds[1:21, 2]
sdsSummary <- rbind(sdsSummary, c(0, NA, NA, NA, NA, NA)); sdsSummary[22,2] <- "none"

for(i in 1:nrow(sdsSummary)){
  sdsSummary[i,3] <- mean(as.numeric(as.character(sds[sds[,1] == sdsSummary[i,1] & sds[,2] == sdsSummary[i,2], 5])))
  sdsSummary[i,4] <- mean(as.numeric(as.character(sds[sds[,1] == sdsSummary[i,1] & sds[,2] == sdsSummary[i,2], 6])))
  sdsSummary[i,5] <- sd(as.numeric(as.character(sds[sds[,1] == sdsSummary[i,1] & sds[,2] == sdsSummary[i,2], 5])))
  sdsSummary[i,6] <- sd(as.numeric(as.character(sds[sds[,1] == sdsSummary[i,1] & sds[,2] == sdsSummary[i,2], 6])))
}
colnames(sdsSummary) <- c("ntraits", "specification", "sdRF", "sdKF", "sdRF", "sdKF")

library(ggplot2)
ggplot(sdsSummary, aes(x = ntraits, y = sdRF, colour = specification)) +
  geom_point(size = 2) + 
  geom_errorbar(aes(ymin = sdRF - sdRF, ymax = sdRF + sdRF)) + 
  theme_bw()

ggplot(sdsSummary, aes(x = ntraits, y = sdKF, colour = specification)) +
  geom_point(size = 2) + 
  geom_errorbar(aes(ymin = sdKF - sdKF, ymax = sdKF + sdKF)) + 
  theme_bw()

1+1








colnames(MCMC_Diagnostics_Distances) <- c("ntraits", "rateMatrixSpecification", "replicateNum", "runNum", "rfDist", "kfDist")
KDEsmoothingBandwidth <- 1
plot(density(as.numeric(MCMC_Diagnostics_Distances[1,"rfDist", 1:1335]), bw = KDEsmoothingBandwidth), xlim = c(0,15), col = rgb(1,0,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[1,"rfDist", 2251:26668]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))

# RF Distances


###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################

# KF Distances
#rowstarts add to 90, need to add to 1800
# perfect rate matrix specification
par(mfrow = c(3,1))
KDEsmoothingBandwidth <- .01
rowStart <- 1
plot(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"kfDist", 1:26668]), bw = KDEsmoothingBandwidth), 
     xlim = c(0,6), ylim = c(0,20), col = rgb(1,0,0), main = "Perfectly Specified Rate Matrix Condition -- KF Distance KDEs")
legend('topright', title = "number of traits", c("2", "4", "8", "16", "32", "64"), col = c(rgb(1,0,0), rgb(1,.5,0), rgb(1,.75,.25), rgb(0,1,.5), rgb(0,.75,1), rgb(0,0,1)), lwd = 1)
for(i in 1:99){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"kfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"kfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(1,0.5,0))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"kfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(1,0.75,0.25))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"kfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(0,1,0.5))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"kfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(0,0.75,1))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"kfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
}

# interregional rate matrix specification
rowStart <- rowStart + 1
plot(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"kfDist", 1:26668]), bw = KDEsmoothingBandwidth), 
     xlim = c(0,6), ylim = c(0,6), col = rgb(1,0,0), main = "Interregional Rate Matrix Condition -- KF Distance KDEs")
legend('topright', title = "number of traits", c("2", "4", "8", "16", "32", "64"), col = c(rgb(1,0,0), rgb(1,.5,0), rgb(1,.75,.25), rgb(0,1,.5), rgb(0,.75,1), rgb(0,0,1)), lwd = 1)
for(i in 1:99){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"kfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"kfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(1,0.5,0))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"kfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(1,0.75,0.25))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"kfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(0,1,0.5))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"kfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(0,0.75,1))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"kfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
}


# chimp-human rate matrix specification
rowStart <- rowStart + 1
plot(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"kfDist", 1:26668]), bw = KDEsmoothingBandwidth), 
     xlim = c(0,6), ylim = c(0,6), col = rgb(1,0,0), main = "Chimp-Human Rate Matrix Condition -- KF Distance KDEs")
legend('topright', title = "number of traits", c("2", "4", "8", "16", "32", "64"), col = c(rgb(1,0,0), rgb(1,.5,0), rgb(1,.75,.25), rgb(0,1,.5), rgb(0,.75,1), rgb(0,0,1)), lwd = 1)
for(i in 1:99){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"kfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"kfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(1,0.5,0))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"kfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(1,0.75,0.25))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"kfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(0,1,0.5))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"kfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(0,0.75,1))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"kfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
}


#########################
#########################
#########################
#########################
#########################
#########################
#########################

# RF Distances
#rowstarts add to 90, need to add to 1800
# perfect rate matrix specification
par(mfrow = c(3,1))
KDEsmoothingBandwidth <- .01
rowStart <- 1
plot(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"rfDist", 1:26668]), bw = KDEsmoothingBandwidth), 
     xlim = c(0,6), ylim = c(0,30), col = rgb(1,0,0), main = "Perfectly Specified Rate Matrix Condition -- RF Distance KDEs")
legend('topright', title = "number of traits", c("2", "4", "8", "16", "32", "57", "64"), col = c(rgb(1,0,0), rgb(1,.5,0), rgb(1,.75,.25), rgb(0,1,.5), rgb(0,.75,1), rgb(0,0,1)), lwd = 1)
for(i in 1:99){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"rfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"rfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(1,0.5,0))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"rfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(1,0.75,0.25))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"rfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(0,1,0.5))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"rfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(0,0.75,1))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"rfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
}

# interregional rate matrix specification
rowStart <- rowStart + 1
plot(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"rfDist", 1:26668]), bw = KDEsmoothingBandwidth), 
     xlim = c(0,6), ylim = c(0,6), col = rgb(1,0,0), main = "Interregional Rate Matrix Condition -- RF Distance KDEs")
legend('topright', title = "number of traits", c("2", "4", "8", "16", "32", "64"), col = c(rgb(1,0,0), rgb(1,.5,0), rgb(1,.75,.25), rgb(0,1,.5), rgb(0,.75,1), rgb(0,0,1)), lwd = 1)
for(i in 1:99){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"rfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"rfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(1,0.5,0))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"rfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(1,0.75,0.25))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"rfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(0,1,0.5))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"rfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(0,0.75,1))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"rfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
}


# chimp-human rate matrix specification
rowStart <- rowStart + 1
plot(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"rfDist", 1:26668]), bw = KDEsmoothingBandwidth), 
     xlim = c(0,6), ylim = c(0,6), col = rgb(1,0,0), main = "Chimp-Human Rate Matrix Condition -- RF Distance KDEs")
legend('topright', title = "number of traits", c("2", "4", "8", "16", "32", "64"), col = c(rgb(1,0,0), rgb(1,.5,0), rgb(1,.75,.25), rgb(0,1,.5), rgb(0,.75,1), rgb(0,0,1)), lwd = 1)
for(i in 1:99){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"rfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"rfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(1,0.5,0))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"rfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(1,0.75,0.25))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"rfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(0,1,0.5))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"rfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(0,0.75,1))
}

for(i in 1:100){
  rowStart <- rowStart+1
  lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart,"rfDist", 1:26668]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
}
