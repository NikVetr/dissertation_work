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


# let's try the parsimony ratchet in phanghorn
setwd(dir = "/Volumes/2TB/mvBM_sims_PB_noiseless")
# setwd(dir = "/Volumes/2TB/mvBM_sims_PB_noisy")
# setwd(dir = "/Volumes/2TB/uvBM_sims_PB_noiseless")
# setwd(dir = "/Volumes/2TB/uvBM_sims_PB_noisy")
# dir.create("parsRatchetRF")


cl <- makeCluster(2, outfile="")
registerDoParallel(cl)
getDoParWorkers()
# stopCluster(cl)

foreach(m=1:100, .packages = c("ape", "coda", "rwty", "phangorn", "lattice", "latticeExtra", "MCMCpack", "readtext", "StatMatch", "mvtnorm")) %dopar% {
  
  traitNumIter <- c(2,4,8,16,32,64,128)
  degreeMisspecification <- c("perfect", "interregional", "chimpHuman")
  degreeCorrelation <- c("2", "3", "4", "6", "8", "9")
  noisyness <- "noiseless"
  # noisyness <- "noisy"
  meanValueCodingDo <- F; avgRFdistMVC <- NA
  numStartingTrees <- 50 #number of independent random starting trees
  unchangedRun <- 400 #how long to keep the ratchet running without any improvement
  traceOutput <- 0 #should the parsimony ratchet output stuff to screen?
  counter <- 0
  l <- 1
  
  mpRFdists <- matrix(ncol = 5)[-1,]
  colnames(mpRFdists) <- c("ntraits", "rateMatrixCorrelation", "replicateNum", "meanRF.divergenceCoding", "meanRF.meanValueCoding")
  
  # for (k in 1:1) {
  # for (k in 1:length(degreeMisspecification)) {
  for (k in 1:length(degreeCorrelation)) {
    
    for (j in 1:length(traitNumIter)) { #iterates over trait number
      
      for(i in m:m){ #iterates over replicates
        
        print(c(k,j,i)); counter <- counter + 1
        
        ntraits <- traitNumIter[j]
        # rateMatrixSpecification <- degreeMisspecification[k]
        rateMatrixCorrelation <- degreeCorrelation[k]
        replicateNum <- i
        # fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_", "replicate_", replicateNum)
        fileName <- paste0("simulations_", ntraits, "traits_perfectRateMatrix_", rateMatrixCorrelation, "correlationStrength_", noisyness, "_replicate_", replicateNum)
        
        #get true tree
        trueTree <- read.nexus(file = paste0("trees/", fileName, "_tree.nex"))
        tips <- trueTree$tip.label
        
        #get traits
        traits <- strsplit(readtext(file = paste0("data/", fileName, "_traits.nex"))$text, split = "\n")[[1]][8:37]
        traits <- t(sapply(1:30, function (x) strsplit(traits[x], split = " ")[[1]]))
        names <- traits[,1]
        traits <- traits [,-1]
        class(traits) <- "numeric"
        rownames(traits) <- names
        
        #get within group covariance matrix
        wgCM <- as.matrix(read.table(file = paste0("rates/", fileName, "_rates.tsv")))
        
        #do divergence coding
        #get simulated individuals
        numIndiv <- 50
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
        if(class(bestTrees) == "multiPhylo"){avgRFdistDC <- mean(sapply(1:length(bestTrees), function(x) dist.topo(unroot(trueTree), bestTrees[[x]])[1]))
        } else {avgRFdistDC <- dist.topo(unroot(trueTree), bestTrees)[1]}
        
        if(meanValueCodingDo){
          #do mean value coding
          traitMeans <- sapply(1:length(traits[1,]), function(x) mean(traits[,x]))
          meanBinarizedTraits <- t(sapply(1:length(traits[,1]), function(x) traits[x,] > traitMeans))
          rownames(meanBinarizedTraits) <- rownames(traits)
          meanBinarizedTraits[meanBinarizedTraits] <- "L"
          meanBinarizedTraits[meanBinarizedTraits == "FALSE"] <- "S"
          mBT.pD <- phyDat(meanBinarizedTraits, type = "USER", levels = c("L", "S"))
          mpTrees <- list()
          for(sepStart in 1:numStartingTrees){
            mpTree <- pratchet(mBT.pD, maxit = 1e5, method = "sankoff", k = unchangedRun, start = rtree(n = length(tips), tip.label = sample(tips)), trace = traceOutput, all = T)
            mpTreeScore <- parsimony(tree = mpTree, data = mBT.pD, method = "sankoff")
            mpTrees[[sepStart]] <- c(mpTree, mpTreeScore) 
          }
          mpScores <- sapply(1:length(mpTrees), function(x) mpTrees[[x]][[2]][1])
          mostPars <- which(mpScores == min(mpScores))
          bestTreesRough <- lapply(1:length(mostPars), function(x) mpTrees[[mostPars[x]]][[1]])
          bestTrees <- bestTreesRough[[1]]
          if(length(bestTreesRough) > 1){for(treeGroup in 2:length(bestTreesRough)){bestTrees <- c(bestTrees, bestTreesRough[[treeGroup]])}}
          if(class(bestTrees) != "phylo"){bestTrees <- unique(bestTrees)}
          if(class(bestTrees) == "multiPhylo"){avgRFdistMVC <- mean(sapply(1:length(bestTrees), function(x) dist.topo(unroot(trueTree), bestTrees[[x]])[1]))
          } else {avgRFdistMVC <- dist.topo(unroot(trueTree), bestTrees)[1]}
        }
        #write to table
        mpRFdists <- rbind(mpRFdists, c(ntraits, rateMatrixCorrelation, replicateNum, avgRFdistDC, avgRFdistMVC))
        print(paste0(ntraits, "traits, ", rateMatrixCorrelation, " corr, ", avgRFdistDC, " divergence RF, ", avgRFdistMVC, " mean value RF"))
        
      }
    }
  }
  
  #for single run
  cat("writing data to disk...\n")
  
  write.table(x = mpRFdists, file = paste0("parsRatchetRF/mpRFdists", m, ".txt"))
  
}

setwd(dir = "/Volumes/2TB/mvBM_sims_PB/parsimony")
mpRFdists <- matrix(ncol = 5)[-1,]
colnames(mpRFdists) <- c("ntraits", "rateMatrixCorrelation", "replicateNum", "meanRF.divergenceCoding", "meanRF.meanValueCoding")
for(m in 1:100){
  mpRFdists <- rbind(mpRFdists, read.table(file = paste0("parsRatchetRF/mpRFdists", m, ".txt")))
}



#let's do it in PAUP* too

# setwd(dir = "/Volumes/2TB/mvBM_sims_PB")
# dir.create("parsimony")
setwd(dir = "/Volumes/2TB/mvBM_sims_PB/parsimony")
# dir.create("traits")
# dir.create("scripts")
# dir.create("MPtrees")

traitNumIter <- c(2,4,8,16,32,64)

degreeMisspecification <- c("perfect", "interregional", "chimpHuman")
degreeCorrelation <- c("2", "3", "4", "6", "8", "9")


i <- j <- k <- 2

# k iterates over degree of misspecification
counter <- 0
# for (k in 1:length(degreeMisspecification)) {
# for (k in 1:length(degreeCorrelation)) {
l <- 1 #this represents which replicate trait sim we are on
for (k in 1:1) {
    
  for (j in 1:length(traitNumIter)) { #iterates over trait number
    
    for(i in 1:2){ #iterates over replicates

      print(c(k,j,i))
      
      ntraits <- traitNumIter[j]
      # rateMatrixSpecification <- degreeMisspecification[k]
      rateMatrixCorrelation <- degreeCorrelation[k]
      replicateNum <- i
      # fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_", "replicate_", replicateNum)
      fileName <- paste0("simulations_", ntraits, "traits_perfectRateMatrix_", rateMatrixCorrelation, "correlationStrength_", "replicate_", replicateNum)
      
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
      
      redo <- TRUE
      stuck <- 0
      while(redo){
      
        if(stuck > 0){
          print("stuck ")
        }
        
        #get simulated individuals      
        numIndiv <- 30
        numPops <- length(traits[,1])
        simulIndiv <- array(data = 0, dim = c(numPops,ntraits,numIndiv))
        for (pops in 1:numPops){
          simulIndiv[pops,,] <- rmvnorm(n = numIndiv, mean = traits[pops,], sigma = wgCM)
          
        }
  
        #get discrete traits via divergence coding
        divergenceCodedTraits <- matrix(data = 0, nrow = numPops, ncol = ntraits)
        pval <- 0.01 
        
        for(trait in 1:ntraits){
          divergenceMatrix <- matrix(data = 0, nrow = numPops, ncol = numPops)
          for(bar in 1:(numPops-1)){
            for(foo in (bar+1):numPops){
              significance <- t.test(simulIndiv[foo,trait,], simulIndiv[bar,trait,])$p.value < pval
              if(significance){
                if(mean(simulIndiv[bar,1,]) > mean(simulIndiv[foo,1,])){
                  divergenceMatrix[foo, bar] <- 1
                } else {
                  divergenceMatrix[foo, bar] <- -1
                }
              } else {
                divergenceMatrix[foo, bar] <- 0
              }
            }
          }
          divergenceMatrix <- divergenceMatrix - t(divergenceMatrix)
          divergenceCodedTraits[,trait] <- sapply(1:numPops, function(x) sum(divergenceMatrix[,x]))
        }
        rownames(divergenceCodedTraits) <- rownames(traits)
        
        diffs <- sapply(1:length(divergenceCodedTraits[1,]), function (x) max(divergenceCodedTraits[,x] - min(divergenceCodedTraits[,x])))
        print(max(diffs))
        tooMany <- any(diffs > 54)
        
        if(!tooMany){
          # print("ok")
          redo <- FALSE
        } else {
          stuck <- stuck + 1
        }
        
      }
      
      codes <- "1234567890qwertyuiopasdfghjklzxcvbnm!@#$%^&*_+-=':<>/|"
      codes <- strsplit(x = codes, split = "")[[1]]
      
      for(foo in 1:length(divergenceCodedTraits[1,])){
        divergenceCodedTraits[,foo] <- divergenceCodedTraits[,foo] + abs(min(divergenceCodedTraits[,foo])) + 1
      }
      
      # print(max(divergenceCodedTraits))
    
      for(bar in 1:length(divergenceCodedTraits[1,])){
        divergenceCodedTraits[,bar] <- convNum2Str(nums = as.numeric(divergenceCodedTraits[,bar]), key = codes)
        }
      
      traitsPath <- paste0("traits/", fileName, "_traits_", i, ".nex")
      
      #write traits
      sink(traitsPath, append=F)
      cat("#NEXUS\n\nbegin data;\n\tdimensions ntax =", numPops, "nchar =", ntraits, ";\n\tFormat symbols=\"1234567890qwertyuiopasdfghjklzxcvbnm!@#$%^&*_+-=':<>/|\"; \n\tmatrix\n\n")
      for(boo in 1:length(divergenceCodedTraits[,1])) {cat(rownames(divergenceCodedTraits)[boo], paste0(divergenceCodedTraits[boo,], collapse = ""), "\n")}
      cat(";\nEnd;")
      sink()
      
      scriptPath <- paste0("scripts/", fileName, "_script_sim_", l, ".txt")
      treesPathTBR <- paste0("MPtrees/", fileName, "_trees_tbr_sim_", l, ".tre")
      treesPathNNI <- paste0("MPtrees/", fileName, "_trees_nni_sim_", l, ".tre")
      treesPathSPR <- paste0("MPtrees/", fileName, "_trees_spr_sim_", l, ".tre")
      
      #write scripts
      sink(scriptPath, append = F)
      # cat("paup\n") 
      cat("execute", paste0(getwd(), "/", traitsPath), ";\n")
      # cat("execute", traitsPath, "\n")
      cat("typeset *myTypesetName = ord:", paste(1:ntraits, collapse = " "), ";\n")
      cat("set maxtrees = 5000 increase=auto;\n")
      cat("set criterion=parsimony;\nhsearch swap=tbr;\n")
      cat("savetrees file=", paste0(getwd(), "/", treesPathTBR), "brlens;\n")
      cat("set criterion=parsimony;\nhsearch swap=nni;\n")
      cat("savetrees file=", paste0(getwd(), "/", treesPathNNI), "brlens;\n")
      cat("set criterion=parsimony;\nhsearch swap=spr;\n")
      cat("savetrees file=", paste0(getwd(), "/", treesPathSPR), "brlens;\n")
      cat("q;")
      sink()
      
      # add script path to list
      if(counter ==0){
        sink("job.txt", append = F)
        cat(paste0(scriptPath, "\n"))
        sink()
        } else {
        sink("job.txt", append = T)
        cat(paste0(scriptPath, "\n"))
        sink()
        }
      counter <- counter + 1
      
      
    }
  }
}







MPtreesTBR <- read.nexus("C:/Users/Nikolai/parsimony/MPtrees/simulations_2traits_perfectRateMatrix_replicate_2_trees_tbr_rep_2.tre")
MPtreesNNI <- read.nexus("C:/Users/Nikolai/parsimony/MPtrees/simulations_2traits_perfectRateMatrix_replicate_2_trees_nni_rep_2.tre")
sum(MPtreesNNI[[1]]$edge.length)
sum(MPtreesTBR[[1]]$edge.length)
