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

badTree <- function(tree) { any(((tree$edge[,1] == 0) + (tree$edge[,2] == 0)) == 2) }

foreach(m=1:100, .packages = c("ape", "coda", "rwty", "phangorn", "lattice", "latticeExtra")) %dopar% {
  
  CT_Compatibility <- matrix(ncol = 15)
  colnames(CT_Compatibility) <- c("ntraits", "rateMatrixSpecification", "replicateNum"
                                    , "numInconsistencies4", "conTreeResolution4", "numInconsistencies5", "conTreeResolution5"
                                    , "numInconsistencies6", "conTreeResolution6", "numInconsistencies7", "conTreeResolution7"
                                    , "numInconsistencies8", "conTreeResolution8", "numInconsistencies9", "conTreeResolution9")
  
  CT_Dists <- matrix(ncol = 9)
  colnames(CT_Dists) <- c("ntraits", "rateMatrixSpecification", "replicateNum", "rfDist4", "rfDist5", "rfDist6", "rfDist7", "rfDist8", "rfDist9")
  
  
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
       
         cat("reading trees... ")
        
        # tree diagnostics
        trees1 <- read.tree(file = paste0("/Users/nikolai/mvBM_sims_all/output/", fileName, "_trees_run_1.trees"))
        trees1 <- trees1[seq(1, length(trees1), by = 5)] #thin further
        write.tree(phy = trees1, file = paste0("/Users/nikolai/mvBM_sims_all/output/", fileName, "_trees_run_1_thinned.trees"))
        
        trees2 <- read.tree(file = paste0("/Users/nikolai/mvBM_sims_all/output/", fileName, "_trees_run_2.trees"))
        trees2 <- trees2[seq(1, length(trees2), by = 5)] #thin further
        write.tree(phy = trees2, file = paste0("/Users/nikolai/mvBM_sims_all/output/", fileName, "_trees_run_2_thinned.trees"))
        
        cat("processing trees...\n")
        
        # assessing accurate retrieval of data generating tree    
        trees <- c(trees1, trees2)
        
        conTreeProbs <- seq(.4, .9, by = .1)
        conTree <- lapply(1:length(conTreeProbs), function(x) consensus(trees, p=conTreeProbs[x]))
        
        numInconsistensies <- rep(0, length(conTreeProbs))
        for(ctp in 1:length(conTreeProbs)){
          if(conTree[[ctp]]$Nnode != 1) {
            if(!badTree(conTree[[ctp]])){
             numInconsistensies[ctp] <- sum(is.na(prop.clades(phy = conTree[[ctp]], unroot(trueTree))))
              } else {numInconsistensies[ctp] <- NA}
            }
            else {
              numInconsistensies[ctp] <- 0
            }
        }
        
        CT_Compatibility <- rbind(CT_Compatibility, c(ntraits, rateMatrixSpecification, replicateNum, 
                                                          numInconsistensies[1], conTree[[1]]$Nnode,
                                                          numInconsistensies[2], conTree[[2]]$Nnode,
                                                          numInconsistensies[3], conTree[[3]]$Nnode,
                                                          numInconsistensies[4], conTree[[4]]$Nnode,
                                                          numInconsistensies[5], conTree[[5]]$Nnode,
                                                          numInconsistensies[6], conTree[[6]]$Nnode))
        
        rfDist <- rep(0, length(conTreeProbs))
        for(ctp in 1:length(conTreeProbs)){
          if(!badTree(conTree[[ctp]])){
            rfDist[ctp] <- RF.dist(tree1 = trueTree, tree2 = conTree[[ctp]])
          } else {
            rfDist[ctp] <- NA
          }
        }
       
        CT_Dists <- rbind(CT_Dists, c(ntraits, rateMatrixSpecification, replicateNum, rfDist[1], rfDist[2], rfDist[3], rfDist[4], rfDist[5], rfDist[6]))
        
        cat("calculating prop. clades...\n")
        
      }
    }
  }
  
  #for single run
  cat("writing data to disk...\n")
  
  write.table(x = CT_Dists, file = paste0("analyses/CT_Dists", windowCode))
  write.table(x = CT_Compatibility, file = paste0("analyses/CT_Compatibility", windowCode))
  
}

#read back in the files
CT_Dists <- read.table(file = "analyses/CT_Dists1")[-1,]
CT_Compatibility <- read.table(file = "analyses/CT_Compatibility1")[-1,]

for (i in 2:100){
  print(i)
  CT_Dists <- rbind(CT_Dists, read.table(file = paste0("analyses/CT_Dists", i))[-1,])
  CT_Compatibility <- rbind(CT_Compatibility, read.table(file = paste0("analyses/CT_Compatibility", i))[-1,])
}


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


CT_Dist <- as.data.frame(CT_Dists)
CT_Dists$ntraits <- as.character(CT_Dists$ntraits)
CT_Dists$ntraits <- factor(CT_Dists$ntraits, levels = c("2", "4", "8", "16", "32", "57", "64"))

bwplot(rfDist4 ~ ntraits, ylim = c(0,54), ylab = "Robinson-Foulds Distance", main = "Comparing Inferred Consensus Tree to Data-Generating Tree (p = 0.4)", xlab = "Number of Traits Used in Simulation", data=CT_Dists, pch="|", coef=0, panel=function(...){panel.hanoi(col="orange", ...); panel.bwplot(...)})
bwplot(rfDist5 ~ ntraits, ylim = c(0,54), ylab = "Robinson-Foulds Distance", main = "Comparing Inferred Consensus Tree to Data-Generating Tree (p = 0.5)", xlab = "Number of Traits Used in Simulation", data=CT_Dists, pch="|", coef=0, panel=function(...){panel.hanoi(col="orange", ...); panel.bwplot(...)})
bwplot(rfDist6 ~ ntraits, ylim = c(0,54), ylab = "Robinson-Foulds Distance", main = "Comparing Inferred Consensus Tree to Data-Generating Tree (p = 0.6)", xlab = "Number of Traits Used in Simulation", data=CT_Dists, pch="|", coef=0, panel=function(...){panel.hanoi(col="orange", ...); panel.bwplot(...)})
bwplot(rfDist7 ~ ntraits, ylim = c(0,54), ylab = "Robinson-Foulds Distance", main = "Comparing Inferred Consensus Tree to Data-Generating Tree (p = 0.7)", xlab = "Number of Traits Used in Simulation", data=CT_Dists, pch="|", coef=0, panel=function(...){panel.hanoi(col="orange", ...); panel.bwplot(...)})
bwplot(rfDist8 ~ ntraits, ylim = c(0,54), ylab = "Robinson-Foulds Distance", main = "Comparing Inferred Consensus Tree to Data-Generating Tree (p = 0.8)", xlab = "Number of Traits Used in Simulation", data=CT_Dists, pch="|", coef=0, panel=function(...){panel.hanoi(col="orange", ...); panel.bwplot(...)})
bwplot(rfDist9 ~ ntraits, ylim = c(0,54), ylab = "Robinson-Foulds Distance", main = "Comparing Inferred Consensus Tree to Data-Generating Tree (p = 0.9)", xlab = "Number of Traits Used in Simulation", data=CT_Dists, pch="|", coef=0, panel=function(...){panel.hanoi(col="orange", ...); panel.bwplot(...)})

CT_Compatibility <- as.data.frame(CT_Compatibility)
CT_Compatibility$ntraits <- as.character(CT_Compatibility$ntraits)
CT_Compatibility$ntraits <- factor(CT_Compatibility$ntraits, levels = c("2", "4", "8", "16", "32", "57", "64"))
CT_Compatibility$conTreeResolution4 <- as.numeric(as.character(CT_Compatibility$conTreeResolution4))
CT_Compatibility$conTreeResolution5 <- as.numeric(as.character(CT_Compatibility$conTreeResolution5))
CT_Compatibility$conTreeResolution6 <- as.numeric(as.character(CT_Compatibility$conTreeResolution6))
CT_Compatibility$conTreeResolution7 <- as.numeric(as.character(CT_Compatibility$conTreeResolution7))
CT_Compatibility$conTreeResolution8 <- as.numeric(as.character(CT_Compatibility$conTreeResolution8))
CT_Compatibility$conTreeResolution9 <- as.numeric(as.character(CT_Compatibility$conTreeResolution9))
CT_Compatibility$numInconsistencies4 <- as.numeric(as.character(CT_Compatibility$numInconsistencies4))
CT_Compatibility$numInconsistencies5 <- as.numeric(as.character(CT_Compatibility$numInconsistencies5))
CT_Compatibility$numInconsistencies6 <- as.numeric(as.character(CT_Compatibility$numInconsistencies6))
CT_Compatibility$numInconsistencies7 <- as.numeric(as.character(CT_Compatibility$numInconsistencies7))
CT_Compatibility$numInconsistencies8 <- as.numeric(as.character(CT_Compatibility$numInconsistencies8))
CT_Compatibility$numInconsistencies9 <- as.numeric(as.character(CT_Compatibility$numInconsistencies9))
CT_Compatibility$numCorrectNodes4 <- CT_Compatibility$conTreeResolution4 - CT_Compatibility$numInconsistencies4
CT_Compatibility$numCorrectNodes5 <- CT_Compatibility$conTreeResolution5 - CT_Compatibility$numInconsistencies5
CT_Compatibility$numCorrectNodes6 <- CT_Compatibility$conTreeResolution6 - CT_Compatibility$numInconsistencies6
CT_Compatibility$numCorrectNodes7 <- CT_Compatibility$conTreeResolution7 - CT_Compatibility$numInconsistencies7
CT_Compatibility$numCorrectNodes8 <- CT_Compatibility$conTreeResolution8 - CT_Compatibility$numInconsistencies8
CT_Compatibility$numCorrectNodes9 <- CT_Compatibility$conTreeResolution9 - CT_Compatibility$numInconsistencies9

key <- list(x = 1.1, y = .92, corner = c(5.5, 0),
            text = list(c("Correctly", "Incorrectly")),
            lines = list(type = c("l", "l"), col = c("green", "red"),
                         pch = 3, lwd = 8, lty = 1))

#.4
a <- bwplot(numCorrectNodes4 ~ ntraits, main = "Consistency of 0.4 Consensus Tree with True, Data-Generating Tree\nALL RATE MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in 0.4 Majority-Rule Consensus Tree", 
            data=CT_Compatibility, ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies4 ~ ntraits, data=CT_Compatibility, ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

a <- bwplot(numCorrectNodes4 ~ ntraits, main = "Consistency of 0.4 Consensus Tree with True, Data-Generating Tree\nPERFECT MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in 0.4 Majority-Rule Consensus Tree", 
            data=CT_Compatibility[CT_Compatibility[,2] == "perfect",], ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies4 ~ ntraits, data=CT_Compatibility[CT_Compatibility[,2] == "perfect",], ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

a <- bwplot(numCorrectNodes4 ~ ntraits, main = "Consistency of 0.4 Consensus Tree with True, Data-Generating Tree\nINTERREGIONAL MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in 0.4 Majority-Rule Consensus Tree", 
            data=CT_Compatibility[CT_Compatibility[,2] == "interregional",], ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies4 ~ ntraits, data=CT_Compatibility[CT_Compatibility[,2] == "interregional",], ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

a <- bwplot(numCorrectNodes4 ~ ntraits, main = "Consistency of 0.4 Consensus Tree with True, Data-Generating Tree\nCHIMP-HUMAN MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in 0.4 Majority-Rule Consensus Tree", 
            data=CT_Compatibility[CT_Compatibility[,2] == "chimpHuman",], ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies4 ~ ntraits, data=CT_Compatibility[CT_Compatibility[,2] == "chimpHuman",], ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

#.5
a <- bwplot(numCorrectNodes5 ~ ntraits, main = "Consistency of 0.5 Consensus Tree with True, Data-Generating Tree\nALL RATE MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in 0.5 Majority-Rule Consensus Tree", 
            data=CT_Compatibility, ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies5 ~ ntraits, data=CT_Compatibility, ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

a <- bwplot(numCorrectNodes5 ~ ntraits, main = "Consistency of 0.5 Consensus Tree with True, Data-Generating Tree\nPERFECT MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in 0.5 Majority-Rule Consensus Tree", 
            data=CT_Compatibility[CT_Compatibility[,2] == "perfect",], ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies5 ~ ntraits, data=CT_Compatibility[CT_Compatibility[,2] == "perfect",], ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

a <- bwplot(numCorrectNodes5 ~ ntraits, main = "Consistency of 0.5 Consensus Tree with True, Data-Generating Tree\nINTERREGIONAL MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in 0.5 Majority-Rule Consensus Tree", 
            data=CT_Compatibility[CT_Compatibility[,2] == "interregional",], ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies5 ~ ntraits, data=CT_Compatibility[CT_Compatibility[,2] == "interregional",], ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

a <- bwplot(numCorrectNodes5 ~ ntraits, main = "Consistency of 0.5 Consensus Tree with True, Data-Generating Tree\nCHIMP-HUMAN MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in 0.5 Majority-Rule Consensus Tree", 
            data=CT_Compatibility[CT_Compatibility[,2] == "chimpHuman",], ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies5 ~ ntraits, data=CT_Compatibility[CT_Compatibility[,2] == "chimpHuman",], ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

#.6
a <- bwplot(numCorrectNodes6 ~ ntraits, main = "Consistency of 0.6 Consensus Tree with True, Data-Generating Tree\nALL RATE MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in 0.6 Majority-Rule Consensus Tree", 
            data=CT_Compatibility, ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies6 ~ ntraits, data=CT_Compatibility, ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

a <- bwplot(numCorrectNodes6 ~ ntraits, main = "Consistency of 0.6 Consensus Tree with True, Data-Generating Tree\nPERFECT MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in 0.6 Majority-Rule Consensus Tree", 
            data=CT_Compatibility[CT_Compatibility[,2] == "perfect",], ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies6 ~ ntraits, data=CT_Compatibility[CT_Compatibility[,2] == "perfect",], ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

a <- bwplot(numCorrectNodes6 ~ ntraits, main = "Consistency of 0.6 Consensus Tree with True, Data-Generating Tree\nINTERREGIONAL MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in 0.6 Majority-Rule Consensus Tree", 
            data=CT_Compatibility[CT_Compatibility[,2] == "interregional",], ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies6 ~ ntraits, data=CT_Compatibility[CT_Compatibility[,2] == "interregional",], ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

a <- bwplot(numCorrectNodes6 ~ ntraits, main = "Consistency of 0.6 Consensus Tree with True, Data-Generating Tree\nCHIMP-HUMAN MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in 0.6 Majority-Rule Consensus Tree", 
            data=CT_Compatibility[CT_Compatibility[,2] == "chimpHuman",], ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies6 ~ ntraits, data=CT_Compatibility[CT_Compatibility[,2] == "chimpHuman",], ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

#.7
a <- bwplot(numCorrectNodes7 ~ ntraits, main = "Consistency of 0.7 Consensus Tree with True, Data-Generating Tree\nALL RATE MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in 0.7 Majority-Rule Consensus Tree", 
            data=CT_Compatibility, ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies7 ~ ntraits, data=CT_Compatibility, ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

a <- bwplot(numCorrectNodes7 ~ ntraits, main = "Consistency of 0.7 Consensus Tree with True, Data-Generating Tree\nPERFECT MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in 0.7 Majority-Rule Consensus Tree", 
            data=CT_Compatibility[CT_Compatibility[,2] == "perfect",], ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies7 ~ ntraits, data=CT_Compatibility[CT_Compatibility[,2] == "perfect",], ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

a <- bwplot(numCorrectNodes7 ~ ntraits, main = "Consistency of 0.7 Consensus Tree with True, Data-Generating Tree\nINTERREGIONAL MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in 0.7 Majority-Rule Consensus Tree", 
            data=CT_Compatibility[CT_Compatibility[,2] == "interregional",], ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies7 ~ ntraits, data=CT_Compatibility[CT_Compatibility[,2] == "interregional",], ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

a <- bwplot(numCorrectNodes7 ~ ntraits, main = "Consistency of 0.7 Consensus Tree with True, Data-Generating Tree\nCHIMP-HUMAN MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in 0.7 Majority-Rule Consensus Tree", 
            data=CT_Compatibility[CT_Compatibility[,2] == "chimpHuman",], ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies7 ~ ntraits, data=CT_Compatibility[CT_Compatibility[,2] == "chimpHuman",], ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

#.8
a <- bwplot(numCorrectNodes8 ~ ntraits, main = "Consistency of 0.8 Consensus Tree with True, Data-Generating Tree\nALL RATE MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in 0.8 Majority-Rule Consensus Tree", 
            data=CT_Compatibility, ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies8 ~ ntraits, data=CT_Compatibility, ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

a <- bwplot(numCorrectNodes8 ~ ntraits, main = "Consistency of 0.8 Consensus Tree with True, Data-Generating Tree\nPERFECT MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in 0.8 Majority-Rule Consensus Tree", 
            data=CT_Compatibility[CT_Compatibility[,2] == "perfect",], ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies8 ~ ntraits, data=CT_Compatibility[CT_Compatibility[,2] == "perfect",], ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

a <- bwplot(numCorrectNodes8 ~ ntraits, main = "Consistency of 0.8 Consensus Tree with True, Data-Generating Tree\nINTERREGIONAL MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in 0.8 Majority-Rule Consensus Tree", 
            data=CT_Compatibility[CT_Compatibility[,2] == "interregional",], ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies8 ~ ntraits, data=CT_Compatibility[CT_Compatibility[,2] == "interregional",], ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

a <- bwplot(numCorrectNodes8 ~ ntraits, main = "Consistency of 0.8 Consensus Tree with True, Data-Generating Tree\nCHIMP-HUMAN MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in 0.8 Majority-Rule Consensus Tree", 
            data=CT_Compatibility[CT_Compatibility[,2] == "chimpHuman",], ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies8 ~ ntraits, data=CT_Compatibility[CT_Compatibility[,2] == "chimpHuman",], ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

#.9
a <- bwplot(numCorrectNodes9 ~ ntraits, main = "Consistency of 0.9 Consensus Tree with True, Data-Generating Tree\nALL RATE MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in 0.9 Majority-Rule Consensus Tree", 
            data=CT_Compatibility, ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies9 ~ ntraits, data=CT_Compatibility, ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

a <- bwplot(numCorrectNodes9 ~ ntraits, main = "Consistency of 0.9 Consensus Tree with True, Data-Generating Tree\nPERFECT MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in 0.9 Majority-Rule Consensus Tree", 
            data=CT_Compatibility[CT_Compatibility[,2] == "perfect",], ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies9 ~ ntraits, data=CT_Compatibility[CT_Compatibility[,2] == "perfect",], ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

a <- bwplot(numCorrectNodes9 ~ ntraits, main = "Consistency of 0.9 Consensus Tree with True, Data-Generating Tree\nINTERREGIONAL MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in 0.9 Majority-Rule Consensus Tree", 
            data=CT_Compatibility[CT_Compatibility[,2] == "interregional",], ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies9 ~ ntraits, data=CT_Compatibility[CT_Compatibility[,2] == "interregional",], ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

a <- bwplot(numCorrectNodes9 ~ ntraits, main = "Consistency of 0.9 Consensus Tree with True, Data-Generating Tree\nCHIMP-HUMAN MATRICES", xlab = "Number of Traits Used in Simulation", 
            key = key, ylab = "Number of Nodes Identified in 0.9 Majority-Rule Consensus Tree", 
            data=CT_Compatibility[CT_Compatibility[,2] == "chimpHuman",], ylim = c(0,25), pch="|", coef=0, panel=function(...){panel.hanoi(col="green", alpha = .7, ...); panel.bwplot(...)})
b <- bwplot(numInconsistencies9 ~ ntraits, data=CT_Compatibility[CT_Compatibility[,2] == "chimpHuman",], ylim = c(0,30), pch="|", coef=0, 
            panel=function(...){panel.hanoi(col="red", alpha = .7, ...)})
a + as.layer(b)

