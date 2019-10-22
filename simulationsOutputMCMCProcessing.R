setwd(dir = "A:\\mvBM_sims")

library(ape)
library(coda)
library(rwty)
library(phangorn)
library(lattice)
library(latticeExtra)

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

windowCode <- "_1_160"

#Also output marginal histograms and trace plots for continuous parameters, split	frequencies	&	presence/absence for trees

traitNumIter <- c(2,4,8,16,32,64)

degreeMisspecification <- c("perfect", "interregional", "chimpHuman")


# k iterates over degree of misspecification
for (k in 1:length(degreeMisspecification)) {

for (j in 1:6) {
  
  for(i in 1:5){
    
    cat(c(j, i, k, "... "))
    
    ntraits <- traitNumIter[j]
    rateMatrixSpecification <- degreeMisspecification[k]
    replicateNum <- i
    
    fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_", "replicate_", replicateNum)
    trueTree <- read.nexus(file = paste0("A:\\mvBM_sims\\trees\\", fileName, "_tree.nex"))
    
    # continuous character diagnostics
    log1 <- read.table(file = paste0("A:\\mvBM_sims\\output\\", fileName, "_params_run_1.log"), header = T)
    log1 <- log1[251:2500,] #discard burn-in
    # log1 <- log1[seq(0, 4000, by = 2),] #thin further
    log1 <- log1[, c("Likelihood", "TL")]
    TL1 <- mcmc(log1$TL, start = 1001, thin = 20)
    LL1 <- mcmc(log1$Likelihood, start = 1001, thin = 20)
    
    log2 <- read.table(file = paste0("A:\\mvBM_sims\\output\\", fileName, "_params_run_2.log"), header = T)
    log2 <- log2[251:2500,] #discard burn-in
    # log2 <- log2[seq(0, 4000, by = 2),] #thin further
    log2 <- log2[, c("Likelihood", "TL")]
    TL2 <- mcmc(log2$TL, start = 1001, thin = 20)
    LL2 <- mcmc(log2$Likelihood, start = 1001, thin = 20)
    
    TLs <- mcmc.list(TL1, TL2)
    LLs <- mcmc.list(LL1, LL2)
    
    cat("processing continuous data... ")
    
    trueTL <- sum(trueTree$edge.length)
    
    png(paste0("A:\\mvBM_sims\\plots\\", fileName, "_TraceMarginals.png"), width = 1920, height = 1080)
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
    trees1 <- read.tree(file = paste0("A:\\mvBM_sims\\output\\", fileName, "_trees_run_1.trees"))
    trees1 <- trees1[251:2500] #discard burn-in
    # trees1 <- trees1[seq(0, 4000, by = 2)] #thin further
    
    trees2 <- read.tree(file = paste0("A:\\mvBM_sims\\output\\", fileName, "_trees_run_2.trees"))
    trees2 <- trees2[251:2500] #discard burn-in
    # trees2 <- trees2[seq(0, 4000, by = 2)] #thin further
    
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
    conTree <- consensus(trees, p=.5)
    trueTree <- read.nexus(file = paste0("A:\\mvBM_sims\\trees\\", fileName, "_tree.nex"))
    
    if(conTree$Nnode != 1) {
      numInconsistensies <- sum(is.na(prop.clades(phy = conTree, trueTree)))} else {
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

windowCode <- "_1_180"
write.table(x = MCMC_Diagnostics_PSRF[-1,], file = paste0("MCMC_Diagnostics_PSRF", windowCode))
write.table(x = MCMC_Diagnostics_Geweke[-1,], file = paste0("MCMC_Diagnostics_Geweke", windowCode))
write.table(x = MCMC_Diagnostics_ESS[-1,], file = paste0("MCMC_Diagnostics_ESS", windowCode))
write.table(x = MCMC_Diagnostics_NodalFreqCorr[-1,], file = paste0("MCMC_Diagnostics_NodalFreqCorr", windowCode))
write.table(x = Tree_Dist[-1,], file = paste0("Tree_Dist", windowCode))
write.table(x = MJCT_Compatibility[-1,], file = paste0("MJCT_Compatibility", windowCode))

MCMC_Diagnostics_PSRF_1_180 <- read.table(file = "MCMC_Diagnostics_PSRF_1_180")
MCMC_Diagnostics_Geweke_1_180 <- read.table(file = "MCMC_Diagnostics_Geweke_1_180")
MCMC_Diagnostics_ESS_1_180 <- read.table(file = "MCMC_Diagnostics_ESS_1_180")
MCMC_Diagnostics_NodalFreqCorr_1_180 <- read.table(file = "MCMC_Diagnostics_NodalFreqCorr_1_180")
Tree_Dist_1_180 <- read.table(file = "Tree_Dist_1_180")
MJCT_Compatibility_1_180 <- read.table(file = "MJCT_Compatibility_1_180")

MCMC_Diagnostics_PSRF_1_80100 <- read.table(file = "MCMC_Diagnostics_PSRF_1_80100")
MCMC_Diagnostics_Geweke_1_80100 <- read.table(file = "MCMC_Diagnostics_Geweke_1_80100")
MCMC_Diagnostics_ESS_1_80100 <- read.table(file = "MCMC_Diagnostics_ESS_1_80100")
MCMC_Diagnostics_NodalFreqCorr_1_80100 <- read.table(file = "MCMC_Diagnostics_NodalFreqCorr_1_80100")
Tree_Dist_1_80100 <- read.table(file = "Tree_Dist_1_80100")
MJCT_Compatibility_1_80100 <- read.table(file = "MJCT_Compatibility_1_80100")

MCMC_Diagnostics_PSRF_23 <- read.table(file = "MCMC_Diagnostics_PSRF_23")
MCMC_Diagnostics_Geweke_23 <- read.table(file = "MCMC_Diagnostics_Geweke_23")
MCMC_Diagnostics_ESS_23 <- read.table(file = "MCMC_Diagnostics_ESS_23")
MCMC_Diagnostics_NodalFreqCorr_23 <- read.table(file = "MCMC_Diagnostics_NodalFreqCorr_23")
Tree_Dist_23 <- read.table(file = "Tree_Dist_23")
MJCT_Compatibility_23 <- read.table(file = "MJCT_Compatibility_23")

MCMC_Diagnostics_PSRF_46 <- read.table(file = "MCMC_Diagnostics_PSRF_46")
MCMC_Diagnostics_Geweke_46 <- read.table(file = "MCMC_Diagnostics_Geweke_46")
MCMC_Diagnostics_ESS_46 <- read.table(file = "MCMC_Diagnostics_ESS_46")
MCMC_Diagnostics_NodalFreqCorr_46 <- read.table(file = "MCMC_Diagnostics_NodalFreqCorr_46")
Tree_Dist_46 <- read.table(file = "Tree_Dist_46")
MJCT_Compatibility_46 <- read.table(file = "MJCT_Compatibility_46")

MCMC_Diagnostics_PSRF <- rbind(MCMC_Diagnostics_PSRF_1_180, MCMC_Diagnostics_PSRF_1_80100, MCMC_Diagnostics_PSRF_23, MCMC_Diagnostics_PSRF_46)
MCMC_Diagnostics_Geweke <- rbind(MCMC_Diagnostics_Geweke_1_180, MCMC_Diagnostics_Geweke_1_80100, MCMC_Diagnostics_Geweke_23, MCMC_Diagnostics_Geweke_46)
MCMC_Diagnostics_ESS <- rbind(MCMC_Diagnostics_ESS_1_180, MCMC_Diagnostics_ESS_1_80100, MCMC_Diagnostics_ESS_23, MCMC_Diagnostics_ESS_46)
MCMC_Diagnostics_NodalFreqCorr <- rbind(MCMC_Diagnostics_NodalFreqCorr_1_180, MCMC_Diagnostics_NodalFreqCorr_1_80100, MCMC_Diagnostics_NodalFreqCorr_23, MCMC_Diagnostics_NodalFreqCorr_46)
Tree_Dist <- rbind(Tree_Dist_1_180, Tree_Dist_1_80100, Tree_Dist_23, Tree_Dist_46)
MJCT_Compatibility <- rbind(MJCT_Compatibility_1_180, MJCT_Compatibility_1_80100, MJCT_Compatibility_23, MJCT_Compatibility_46)

###################################################################
#### OK #### NOW #### THAT #### WE #### ARE #### ALL #### HERE ####
###################################################################

MCMC_Diagnostics_PSRF <- MCMC_Diagnostics_PSRF[-1,]
MCMC_Diagnostics_Geweke <- MCMC_Diagnostics_Geweke[-1,] 
MCMC_Diagnostics_ESS <- MCMC_Diagnostics_ESS[-1,] 
MCMC_Diagnostics_NodalFreqCorr <- MCMC_Diagnostics_NodalFreqCorr[-1,]
Tree_Dist <- Tree_Dist[-1,] 
MJCT_Compatibility <- MJCT_Compatibility[-1,]
  
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
Tree_Dist$ntraits <- factor(Tree_Dist$ntraits, levels = c("2", "4", "8", "16", "32", "64"))

bwplot(rfDist ~ ntraits, ylab = "Robinson-Foulds Distance", main = "Comparing inferred MCC tree to Data-Generating Tree", xlab = "Number of Traits Used in Simulation", data=Tree_Dist, pch="|", coef=0, panel=function(...){panel.hanoi(col="orange", ...); panel.bwplot(...)})
bwplot(kfDist ~ ntraits, ylab = "Kuhner-Felsenstein Distance", main = "Comparing inferred MCC tree to Data-Generating Tree", xlab = "Number of Traits Used in Simulation", data=Tree_Dist, pch="|", coef=0, panel=function(...){panel.hanoi(col="orange", ...); panel.bwplot(...)})

MJCT_Compatibility <- as.data.frame(MJCT_Compatibility)
MJCT_Compatibility$ntraits <- as.character(MJCT_Compatibility$ntraits)
MJCT_Compatibility$ntraits <- factor(MJCT_Compatibility$ntraits, levels = c("2", "4", "8", "16", "32", "57"))
MJCT_Compatibility$conTreeResolution <- as.numeric(as.character(MJCT_Compatibility$conTreeResolution))
MJCT_Compatibility$numInconsistencies <- as.numeric(as.character(MJCT_Compatibility$numInconsistencies))
MJCT_Compatibility$numCorrectNodes <- MJCT_Compatibility$conTreeResolution - MJCT_Compatibility$numInconsistencies

key <- list(x = .97, y = .92, corner = c(5.5, 0),
            text = list(c("Correctly", "Incorrectly")),
            lines = list(type = c("l", "l"), col = c("green", "red"),
                         pch = 3, lwd = 8, lty = 1))

a <- bwplot(numCorrectNodes ~ ntraits, main = "Consistency with MJ-Consensus Tree", xlab = "Number of Traits Used in Simulation", key = key, ylab = "Number of Nodes Identified in Strict Majority-Rule Consensus Tree", data=MJCT_Compatibility, ylim = c(0,22), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies ~ ntraits, data=MJCT_Compatibility, ylim = c(0,22), pch="|", coef=0, panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

####################################
## ASSESS TOPOLOGICAL CONVERGENCE ##
####################################

MCMC_Diagnostics_Distances <- array(dim = c(90, 6, 4500))

traitNumIter <- c(2,4,8,16,32,64)

degreeMisspecification <- c("perfect", "interregional", "chimpHuman")

rowCounter <- 0
# k iterates over degree of misspecification
for (k in 1:length(degreeMisspecification)) {
  
  for (j in 1:6) {
    
    for(i in 1:5){
      
      rowCounter <- rowCounter + 1
      print(paste(rowCounter/.9,j, i, k, "... "))
      
      ntraits <- traitNumIter[j]
      rateMatrixSpecification <- degreeMisspecification[k]
      replicateNum <- i
      
      fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_", "replicate_", replicateNum)
      trueTree <- read.nexus(file = paste0("A:\\mvBM_sims\\trees\\", fileName, "_tree.nex"))
      #plot(trueTree)
      
      trees1 <- read.tree(file = paste0("A:\\mvBM_sims\\output\\", fileName, "_trees_run_1.trees"))
      trees1 <- trees1[251:2500] #discard burn-in
      trees2 <- read.tree(file = paste0("A:\\mvBM_sims\\output\\", fileName, "_trees_run_2.trees"))
      trees2 <- trees2[251:2500] #discard burn-in
      
      RFdists1 <- vector()
      RFdists2 <- vector()
      KFdists1 <- vector()
      KFdists2 <- vector()
      
      for(l in 1:length(trees1)){
        dist1 <- treedist(trueTree, trees1[[l]])
        RFdists1[l] <- dist1[[1]]
        KFdists1[l] <- dist1[[2]]
        dist2 <- treedist(trueTree, trees2[[l]])
        RFdists2[l] <- dist2[[1]]
        KFdists2[l] <- dist2[[2]]
      }
      
      MCMC_Diagnostics_Distances[rowCounter, 1,] <- rep(ntraits, 4500)
      MCMC_Diagnostics_Distances[rowCounter, 2,] <- rep(rateMatrixSpecification, 4500)
      MCMC_Diagnostics_Distances[rowCounter, 3,] <- rep(replicateNum, 4500)
      MCMC_Diagnostics_Distances[rowCounter, 4,] <- c(rep(1,2250), rep(2,2250))
      MCMC_Diagnostics_Distances[rowCounter, 5,] <- c(RFdists1, RFdists2)
      MCMC_Diagnostics_Distances[rowCounter, 6,] <- c(KFdists1, KFdists2)
      
      RFdists1 <- mcmc(RFdists1, start = 1)
      RFdists2 <- mcmc(RFdists2, start = 1)
      KFdists1 <- mcmc(KFdists1, start = 1)
      KFdists2 <- mcmc(KFdists2, start = 1)

      RFs <- mcmc.list(RFdists1, RFdists2)
      KFs <- mcmc.list(KFdists1, KFdists2)
      
      png(paste0("A:\\mvBM_sims\\plots\\TreeDist\\", fileName, "_TreeDistTraceMarginals.png"), width = 1920, height = 1080)
      par(mfrow = c(2,2))
      traceplot(RFs, main = "RFs"); traceplot(KFs, main = "KFs")
      hist(x = RFdists1, breaks = 50, col=rgb(1,0,0,0.3), main = "RFs from True Tree"); hist(x = RFdists2, breaks = 50, col=rgb(0,0,1,0.3), add = T)
      hist(x = KFdists1, breaks = 50, col=rgb(1,0,0,0.3), main = "KFs from True Tree"); hist(x = KFdists2, breaks = 50, col=rgb(0,0,1,0.3), add = T)
      dev.off()
      
    }
  }
}

par(mfrow = c(1,1))
means <- MCMC_Diagnostics_Distances[,,1]
means[,5] <- sapply(1:90, function (x) mean(as.numeric(MCMC_Diagnostics_Distances[x,5,])))
means[,6] <- sapply(1:90, function (x) mean(as.numeric(MCMC_Diagnostics_Distances[x,6,])))
means <- as.data.frame(means)
plot(as.numeric(as.character(means$ntraits)), as.numeric(as.character(means$rfDist)), col = means$rateMatrixSpecification, 
     xlab = "Number of Traits", ylab = "Mean RF-Distance", main = "Mean RF-Dist from True, Data-Generating Tree")
legend(50, 53,unique(means$rateMatrixSpecification),col=1:length(means$rateMatrixSpecification),pch=1)

plot(as.numeric(as.character(means$ntraits)), as.numeric(as.character(means$kfDist)), col = means$rateMatrixSpecification, 
     xlab = "Number of Traits", ylab = "Mean KF-Distance", main = "Mean KF-Dist from True, Data-Generating Tree")
legend(50, 3,unique(means$rateMatrixSpecification),col=1:length(means$rateMatrixSpecification),pch=1)


sds <- MCMC_Diagnostics_Distances[,,1]
sds[,5] <- sapply(1:90, function (x) sd(as.numeric(MCMC_Diagnostics_Distances[x,5,])))
sds[,6] <- sapply(1:90, function (x) sd(as.numeric(MCMC_Diagnostics_Distances[x,6,])))
sds <- as.data.frame(sds)
plot(as.numeric(as.character(sds$ntraits)), as.numeric(as.character(sds$rfDist)), col = sds$rateMatrixSpecification, 
     xlab = "Number of Traits", ylab = "Standard Deviation of RF-Distances", main = "Standard Deviation of RF-Dists from True, Data-Generating Tree")
legend(50, 53,unique(sds$rateMatrixSpecification),col=1:length(sds$rateMatrixSpecification),pch=1)

plot(as.numeric(as.character(sds$ntraits)), as.numeric(as.character(sds$kfDist)), col = sds$rateMatrixSpecification, 
     xlab = "Number of Traits", ylab = "Standard Deviation of KF-Distances", main = "Standard Deviation of KF-Dists from True, Data-Generating Tree")
legend(50, 3,unique(sds$rateMatrixSpecification),col=1:length(sds$rateMatrixSpecification),pch=1)



colnames(MCMC_Diagnostics_Distances) <- c("ntraits", "rateMatrixSpecification", "replicateNum", "runNum", "rfDist", "kfDist")

plot(density(as.numeric(MCMC_Diagnostics_Distances[1,"rfDist", 1:2250]), bw = KDEsmoothingBandwidth), xlim = c(15,60), col = rgb(1,0,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[1,"rfDist", 2251:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))

# RF Distances

# perfect rate matrix specification
par(mfrow = c(3,1))
KDEsmoothingBandwidth <- 1
rowStart <- 0
plot(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), 
     xlim = c(15,60), col = rgb(1,0,0), main = "Perfectly Specified Rate Matrix Condition -- RF Distance KDEs")
legend('topright', title = "number of traits", c("2", "4", "8", "16", "32", "64"), col = c(rgb(1,0,0), rgb(1,.5,0), rgb(1,.75,.25), rgb(0,1,.5), rgb(0,.75,1), rgb(0,0,1)), lwd = 1)
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))

rowStart <- rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))

rowStart <- rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))

rowStart <- rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))

rowStart <-  rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))

rowStart <- rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))

# interregional rate matrix specification
rowStart <- rowStart + 5
plot(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), 
     xlim = c(15,60), col = rgb(1,0,0), main = "Interregional Rate Matrix Condition -- RF Distance KDEs")
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))

rowStart <- rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))

rowStart <- rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))

rowStart <- rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))

rowStart <-  rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))

rowStart <- rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))


# Chimp-Human rate matrix specification
rowStart <- rowStart + 5
plot(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), 
     xlim = c(15,60), col = rgb(1,0,0), main = "Chimp-Human Rate Matrix Condition -- RF Distance KDEs")
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))

rowStart <- rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))

rowStart <- rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))

rowStart <- rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))

rowStart <-  rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))

rowStart <- rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"rfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))


###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################

# KF Distances

# perfect rate matrix specification
par(mfrow = c(3,1))
KDEsmoothingBandwidth <- .1
rowStart <- 0
plot(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), 
     xlim = c(0,4), ylim = c(0,4), col = rgb(1,0,0), main = "Perfectly Specified Rate Matrix Condition -- KF Distance KDEs")
legend('topright', title = "number of traits", c("2", "4", "8", "16", "32", "64"), col = c(rgb(1,0,0), rgb(1,.5,0), rgb(1,.75,.25), rgb(0,1,.5), rgb(0,.75,1), rgb(0,0,1)), lwd = 1)
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))

rowStart <- rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))

rowStart <- rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))

rowStart <- rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))

rowStart <-  rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))

rowStart <- rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))

# interregional rate matrix specification
rowStart <- rowStart + 5
plot(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), 
     xlim = c(0,4), ylim = c(0,4), col = rgb(1,0,0), main = "Interregional Rate Matrix Condition -- KF Distance KDEs")
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))

rowStart <- rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))

rowStart <- rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))

rowStart <- rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))

rowStart <-  rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))

rowStart <- rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))


# Chimp-Human rate matrix specification
rowStart <- rowStart + 5
plot(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), 
     xlim = c(0,4), ylim = c(0,4), col = rgb(1,0,0), main = "Chimp-Human Rate Matrix Condition -- KF Distance KDEs")
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,0,0))

rowStart <- rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.5,0))

rowStart <- rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(1,.75,.25))

rowStart <- rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,1,.5))

rowStart <-  rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,.75,1))

rowStart <- rowStart + 5
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 1,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 2,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 3,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 4,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))
lines(density(as.numeric(MCMC_Diagnostics_Distances[rowStart + 5,"kfDist", 1:4500]), bw = KDEsmoothingBandwidth), col = rgb(0,0,1))

