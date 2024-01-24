setwd("/Users/nikolai/mvBM_sims_PB_extra")

# CREATING REVBAYES REPLICATE SCRIPTS

library(ape)
library(coda)
library(rwty)
library(phangorn)
library(lattice)
library(latticeExtra)
library(MCMCpack)
library(mvMORPH)
library(phytools)
library(TESS)

##########################################
############ USEFUL FUNCTIONS ############
##########################################

#write means to file function
meansNexus <- function (characterMatrix, outputFilePath, traitsOnColumns = T){
  if (traitsOnColumns) {ntax = dim(characterMatrix)[1]; nchar = dim(characterMatrix)[2]}
  else {ntax = dim(characterMatrix)[2]; nchar = dim(characterMatrix)[1]}
  sink(outputFilePath, append=F)
  cat("#NEXUS\n\nBegin data;\n\tDimensions ntax=", ntax, "nchar=", nchar, ";\n\tFormat datatype=Continuous missing=?;\n\tMatrix\n\n")
  for(i in 1:length(characterMatrix[,1])) {cat(rownames(characterMatrix)[i], characterMatrix[i,], "\n")}
  cat(";\nEnd;")
  sink()
}

ratesPhy <- function (rateMatrix, outputFilePath) {
  if(!isSymmetric(rateMatrix)) {stop("Error: This matrix is not symmetric.")}
  sink(outputFilePath, append=F)
  cat(dim(rateMatrix)[1],"\n")
  for(i in 1:length(rateMatrix[,1])) {
    cat(rownames(rateMatrix)[i], paste0(as.vector(rateMatrix[i, (1:(length(rateMatrix[i,])-1))])))
    cat(" ")
    cat(rateMatrix[i, (length(rateMatrix[i,]))], "\n", sep = "")
  }
  sink()
}

ratesTSV <- function (rateMatrix, outputFilePath) {
  if(!isSymmetric(rateMatrix)) {stop("Error: This matrix is not symmetric.")}
  sink(outputFilePath, append=F)
  for(i in 1:length(rateMatrix[,1])) {
    cat(paste0(as.vector(rateMatrix[i, (1:(length(rateMatrix[i,])-1))])), sep = "\t")
    cat("\t")
    cat(rateMatrix[i, (length(rateMatrix[i,]))], "\n\n", sep = "")
  }
  sink()
}

###deletes all files in target directories###
# unlink("trees/*")
# unlink("data/*")
# unlink("scripts/*")
# unlink("rates/*")




## preparing the files for noisy analysis ##

traitNumIter <- c(128)

degreeMisspecification <- "perfect"

degreeCorrelation <- c(0.2, 0.4, 0.6, 0.8, 0.99)

nrun <- 2

counter <- 0
# k iterates over degree of misspecification
for (k in 1:length(degreeCorrelation)) {
  
  # i iterates over replicates
  for(i in 1:100){
    
    # j iterates over trait numbers
    for (j in 1:length(traitNumIter)) {
      
      
      print(c(k,i,j)); counter <- counter + 1
      
      ntraits <- traitNumIter[j]
      
      cors <- matrix(nrow = ntraits, ncol = ntraits, data = degreeCorrelation[k])
      diag(cors) <- 1
      cov <- cors
      
      rateMatrixSpecification <- "perfect"
      
      replicateNum <- i
      
      fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_", tail(strsplit(as.character(degreeCorrelation[k]), split = "")[[1]], n = 1), "correlationStrength_", "noisy_", "replicate_", replicateNum)
      
      ratesPathTSV <- paste0("rates/", fileName, "_rates.tsv")
      ratesTSV(round(cov, 4), ratesPathTSV)
      
      # simulate the traits and write to file
      traitsPath <- paste0("data/", fileName, "_traits.nex")
      
      #use a pbtree
      tree <- tess.sim.taxa(n = 1, nTaxa = 30, lambda = 1, mu = 0, samplingProbability = 1, max = 1e6)[[1]]
      tree$tip.label <- sample(x = tree$tip.label, size = length(tree$tip.label), replace = F)
      # plot(tree)
      
      traits <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = dim(cov)[1], sigma=round(cov, 4), mu= rep(0, times = dim(cov)[1])))
      #wiggle those wiggley traits
      traits <- traits + rnorm(n = ntraits*length(tree$tip.label), sd = 0.5^0.5, mean = 0)
      
      meansNexus(traits, traitsPath)
      treePath <- paste0("trees/", fileName, "_tree.nex")
      write.nexus(tree, file = treePath, translate = T)
      
      # write the master scripts
      for (l in 1:nrun) {
        scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
        sink(scriptPath, append = F)
        
        cat("## Analysis of Simulated Data:", ntraits, "traits,", rateMatrixSpecification, "rate matrix specification, replication number", replicateNum, "##\n\n")
        cat("replicateNum <-", replicateNum, "\n")
        cat("runNum <-", l, "\n")
        cat(paste0("rateMatrixSpecification <- \"", rateMatrixSpecification, "\"\n\n"))
        cat("mi = 0\nseed(runNum)\n\n")
        cat("corrStr =", tail(strsplit(as.character(degreeCorrelation[k]), split = "")[[1]], n = 1), "\n\n")
        
        cat("# read in the data\n")
        cat(paste0("means <- readContinuousCharacterData(\"data/", fileName, "_traits.nex\")\n\n"))
        
        cat("# read in the rate matrix\n")
        cat(paste0("rateMatrix = readMatrix(\"rates/", fileName, "_rates.tsv\")\n\n"))
        
        cat("# get some useful variables from the data\nnumTips = means.ntaxa()\nnumBranches = 2 * numTips - 2\nntraits <- means.nchar()\ntaxa = means.taxa()\nnumNodes = numTips * 2 - 1\n\n")
        
        cat("# specify tree topology with BLs prior and moves\n")
        cat(paste0("speciation <- 1\n"))
        cat(paste0("extinction <- 0\n"))
        cat(paste0("sampFraction <- 1\n"))
        cat(paste0("rootAge <- ", max(nodeHeights(tree = tree)),"\n"))
        cat("source(\"simulationsTreeTopology.Rev\")\n\n")
        
        cat("# put it all together\nmvBNdata ~ dnPhyloMultivariateBrownianREML(phylogeny, rep(1.0, numBranches), rateMatrix, nSites=ntraits, method = \"transform\")\nmvBNdata.clamp(means)\n\n")
        
        cat("# specify monitors\nsource(\"simulationsMonitors.Rev\")\n\n")
        
        cat("# define model\nmymodel = model(phylogeny)\n\n")
        
        cat("# run MCMC\nsource(\"simulationsMCMC.Rev\")\n\n")
        
        cat("# analyze mcmc output\nsource(\"simulationsAnalysis.Rev\")\n\n")
        
        cat("q()")
        
        sink()
        
        # add master script name to list
        
        jobName <- "job_pb_noisy_128.txt"
        if(counter == 1){
          scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
          sink(jobName, append = F)
          cat(scriptPath, "\n", sep = "")
          sink()
        } else {
          scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
          sink(jobName, append = T)
          cat(scriptPath, "\n", sep = "")
          sink()
        }
      }
      
    }
    
  }
  
}


## preparing the files for noiseless analysis under mvBM ##

traitNumIter <- c(2, 4, 8, 16, 32, 64, 128)

degreeMisspecification <- "perfect"

degreeCorrelation <- c(0.2, 0.4, 0.6, 0.8, 0.99)

nrun <- 2

counter <- 0
# k iterates over degree of misspecification
for (k in 1:length(degreeCorrelation)) {
  
  # i iterates over replicates
  for(i in 1:100){
    
    # j iterates over trait numbers
    for (j in 1:length(traitNumIter)) {
      
      
      print(c(k,i,j)); counter <- counter + 1
      
      ntraits <- traitNumIter[j]
      
      cors <- matrix(nrow = ntraits, ncol = ntraits, data = degreeCorrelation[k])
      diag(cors) <- 1
      cov <- cors
      
      rateMatrixSpecification <- "perfect"
      
      replicateNum <- i
      
      fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_", tail(strsplit(as.character(degreeCorrelation[k]), split = "")[[1]], n = 1), "correlationStrength_", "noiseless_", "replicate_", replicateNum)
      
      ratesPathTSV <- paste0("rates/", fileName, "_rates.tsv")
      ratesTSV(round(cov, 4), ratesPathTSV)
      
      # simulate the traits and write to file
      traitsPath <- paste0("data/", fileName, "_traits.nex")
      
      #use a pbtree
      tree <- tess.sim.taxa(n = 1, nTaxa = 30, lambda = 1, mu = 0, samplingProbability = 1, max = 1e6)[[1]]
      tree$tip.label <- sample(x = tree$tip.label, size = length(tree$tip.label), replace = F)
      # plot(tree)
      
      traits <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = dim(cov)[1], sigma=round(cov, 4), mu= rep(0, times = dim(cov)[1])))
      
      meansNexus(traits, traitsPath)
      treePath <- paste0("trees/", fileName, "_tree.nex")
      write.nexus(tree, file = treePath, translate = T)
      
      # write the master scripts
      for (l in 1:nrun) {
        scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
        sink(scriptPath, append = F)
        
        cat("## Analysis of Simulated Data:", ntraits, "traits,", rateMatrixSpecification, "rate matrix specification, replication number", replicateNum, "##\n\n")
        cat("replicateNum <-", replicateNum, "\n")
        cat("runNum <-", l, "\n")
        cat(paste0("rateMatrixSpecification <- \"", rateMatrixSpecification, "\"\n\n"))
        cat("mi = 0\nseed(runNum)\n\n")
        cat("corrStr =", tail(strsplit(as.character(degreeCorrelation[k]), split = "")[[1]], n = 1), "\n\n")
        
        cat("# read in the data\n")
        cat(paste0("means <- readContinuousCharacterData(\"data/", fileName, "_traits.nex\")\n\n"))
        
        cat("# read in the rate matrix\n")
        cat(paste0("rateMatrix = readMatrix(\"rates/", fileName, "_rates.tsv\")\n\n"))
        
        cat("# get some useful variables from the data\nnumTips = means.ntaxa()\nnumBranches = 2 * numTips - 2\nntraits <- means.nchar()\ntaxa = means.taxa()\nnumNodes = numTips * 2 - 1\n\n")
        
        cat("# specify tree topology with BLs prior and moves\n")
        cat(paste0("speciation <- 1\n"))
        cat(paste0("extinction <- 0\n"))
        cat(paste0("sampFraction <- 1\n"))
        cat(paste0("rootAge <- ", max(nodeHeights(tree = tree)),"\n"))
        cat("source(\"simulationsTreeTopology.Rev\")\n\n")
        
        cat("# put it all together\nmvBNdata ~ dnPhyloMultivariateBrownianREML(phylogeny, rep(1.0, numBranches), rateMatrix, nSites=ntraits, method = \"transform\")\nmvBNdata.clamp(means)\n\n")
        
        cat("# specify monitors\nsource(\"simulationsMonitors.Rev\")\n\n")
        
        cat("# define model\nmymodel = model(phylogeny)\n\n")
        
        cat("# run MCMC\nsource(\"simulationsMCMC.Rev\")\n\n")
        
        cat("# analyze mcmc output\nsource(\"simulationsAnalysis.Rev\")\n\n")
        
        cat("q()")
        
        sink()
        
        # add master script name to list
        
        jobName <- "job_pb_noiseless_all.txt"
        if(counter == 1){
          scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
          sink(jobName, append = F)
          cat(scriptPath, "\n", sep = "")
          sink()
        } else {
          scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
          sink(jobName, append = T)
          cat(scriptPath, "\n", sep = "")
          sink()
        }
      }
      
    }
    
  }
  
}




## preparing the files for noiseless analysis under uvBM misspecification ##

traitNumIter <- c(128)

degreeMisspecification <- "perfect"

degreeCorrelation <- c(0.2, 0.4, 0.6, 0.8, 0.99)

nrun <- 2

counter <- 0

# k iterates over degree of misspecification
for (k in 1:length(degreeCorrelation)) {
  
  # i iterates over replicates
  for(i in 1:100){
    
    # j iterates over trait numbers
    for (j in 1:length(traitNumIter)) {
      
      
      print(c(k,i,j)); counter <- counter + 1
      
      ntraits <- traitNumIter[j]
      
      cors <- matrix(nrow = ntraits, ncol = ntraits, data = degreeCorrelation[k])
      diag(cors) <- 1
      cov <- cors
      
      rateMatrixSpecification <- "perfect"
      
      replicateNum <- i
      
      fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_", tail(strsplit(as.character(degreeCorrelation[k]), split = "")[[1]], n = 1), "correlationStrength_", "noiseless_uvBM_", "replicate_", replicateNum)
      
      ratesPathTSV <- paste0("rates/", fileName, "_rates.tsv")
      ratesTSV(round(cov, 4), ratesPathTSV)
      
      # simulate the traits and write to file
      traitsPath <- paste0("data/", fileName, "_traits.nex")
      
      #use a pbtree
      tree <- tess.sim.taxa(n = 1, nTaxa = 30, lambda = 1, mu = 0, samplingProbability = 1, max = 1e6)[[1]]
      tree$tip.label <- sample(x = tree$tip.label, size = length(tree$tip.label), replace = F)
      # plot(tree)
      
      traits <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = dim(cov)[1], sigma=round(cov, 4), mu= rep(0, times = dim(cov)[1])))
      
      meansNexus(traits, traitsPath)
      treePath <- paste0("trees/", fileName, "_tree.nex")
      write.nexus(tree, file = treePath, translate = T)
      
      # write the master scripts
      for (l in 1:nrun) {
        scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
        sink(scriptPath, append = F)
        
        cat("## Analysis of Simulated Data:", ntraits, "traits,", rateMatrixSpecification, "rate matrix specification, replication number", replicateNum, "##\n\n")
        cat("replicateNum <-", replicateNum, "\n")
        cat("runNum <-", l, "\n")
        cat(paste0("rateMatrixSpecification <- \"", rateMatrixSpecification, "\"\n\n"))
        cat("mi = 0\nseed(runNum)\n\n")
        cat("corrStr =", tail(strsplit(as.character(degreeCorrelation[k]), split = "")[[1]], n = 1), "\n\n")
        
        cat("# read in the data\n")
        cat(paste0("means <- readContinuousCharacterData(\"data/", fileName, "_traits.nex\")\n\n"))
        
        cat("# read in the rate matrix\n")
        cat(paste0("rate = 1.0\n\n"))
        
        cat("# get some useful variables from the data\nnumTips = means.ntaxa()\nnumBranches = 2 * numTips - 2\nntraits <- means.nchar()\ntaxa = means.taxa()\nnumNodes = numTips * 2 - 1\n\n")
        
        cat("# specify tree topology with BLs prior and moves\n")
        cat(paste0("speciation <- 1\n"))
        cat(paste0("extinction <- 0\n"))
        cat(paste0("sampFraction <- 1\n"))
        cat(paste0("rootAge <- ", max(nodeHeights(tree = tree)),"\n"))
        cat("source(\"simulationsTreeTopology.Rev\")\n\n")
        
        cat("# put it all together\nmvBNdata ~ dnPhyloBrownianREML(phylogeny, rep(1.0, numBranches), rep(rate, ntraits), nSites=ntraits)\nmvBNdata.clamp(means)\n\n")
        
        cat("# specify monitors\nsource(\"simulationsMonitors.Rev\")\n\n")
        
        cat("# define model\nmymodel = model(phylogeny)\n\n")
        
        cat("# run MCMC\nsource(\"simulationsMCMC.Rev\")\n\n")
        
        cat("# analyze mcmc output\nsource(\"simulationsAnalysis.Rev\")\n\n")
        
        cat("q()")
        
        sink()
        
        # add master script name to list
        
        jobName <- "job_uv_noiseless_128.txt"
        if(counter == 1){
          scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
          sink(jobName, append = F)
          cat(scriptPath, "\n", sep = "")
          sink()
        } else {
          scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
          sink(jobName, append = T)
          cat(scriptPath, "\n", sep = "")
          sink()
        }
      }
      
    }
    
  }
  
}




## preparing the files for noisy analysis under uvBM misspecification ##

traitNumIter <- c(2, 4, 8, 16, 32, 64, 128)

degreeMisspecification <- "perfect"

degreeCorrelation <- c(0.2, 0.4, 0.6, 0.8, 0.99)

nrun <- 2

counter <- 0

# k iterates over degree of misspecification
for (k in 1:length(degreeCorrelation)) {
  
  # i iterates over replicates
  for(i in 1:100){
    
    # j iterates over trait numbers
    for (j in 1:length(traitNumIter)) {
      
      
      print(c(k,i,j)); counter <- counter + 1
      
      ntraits <- traitNumIter[j]
      
      cors <- matrix(nrow = ntraits, ncol = ntraits, data = degreeCorrelation[k])
      diag(cors) <- 1
      cov <- cors
      
      rateMatrixSpecification <- "perfect"
      
      replicateNum <- i
      
      fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_", tail(strsplit(as.character(degreeCorrelation[k]), split = "")[[1]], n = 1), "correlationStrength_", "noisy_uvBM_", "replicate_", replicateNum)
      
      ratesPathTSV <- paste0("rates/", fileName, "_rates.tsv")
      ratesTSV(round(cov, 4), ratesPathTSV)
      
      # simulate the traits and write to file
      traitsPath <- paste0("data/", fileName, "_traits.nex")
      
      #use a pbtree
      tree <- tess.sim.taxa(n = 1, nTaxa = 30, lambda = 1, mu = 0, samplingProbability = 1, max = 1e6)[[1]]
      tree$tip.label <- sample(x = tree$tip.label, size = length(tree$tip.label), replace = F)
      # plot(tree)
      
      traits <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = dim(cov)[1], sigma=round(cov, 4), mu= rep(0, times = dim(cov)[1])))
      #wiggle those wiggley traits
      traits <- traits + rnorm(n = ntraits*length(tree$tip.label), sd = 0.5^0.5, mean = 0)
      
      meansNexus(traits, traitsPath)
      treePath <- paste0("trees/", fileName, "_tree.nex")
      write.nexus(tree, file = treePath, translate = T)
      
      # write the master scripts
      for (l in 1:nrun) {
        scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
        sink(scriptPath, append = F)
        
        cat("## Analysis of Simulated Data:", ntraits, "traits,", rateMatrixSpecification, "rate matrix specification, replication number", replicateNum, "##\n\n")
        cat("replicateNum <-", replicateNum, "\n")
        cat("runNum <-", l, "\n")
        cat(paste0("rateMatrixSpecification <- \"", rateMatrixSpecification, "\"\n\n"))
        cat("mi = 0\nseed(runNum)\n\n")
        cat("corrStr =", tail(strsplit(as.character(degreeCorrelation[k]), split = "")[[1]], n = 1), "\n\n")
        
        cat("# read in the data\n")
        cat(paste0("means <- readContinuousCharacterData(\"data/", fileName, "_traits.nex\")\n\n"))
        
        cat("# read in the rate matrix\n")
        cat(paste0("rate = 1.0\n\n"))
        
        cat("# get some useful variables from the data\nnumTips = means.ntaxa()\nnumBranches = 2 * numTips - 2\nntraits <- means.nchar()\ntaxa = means.taxa()\nnumNodes = numTips * 2 - 1\n\n")
        
        cat("# specify tree topology with BLs prior and moves\n")
        cat(paste0("speciation <- 1\n"))
        cat(paste0("extinction <- 0\n"))
        cat(paste0("sampFraction <- 1\n"))
        cat(paste0("rootAge <- ", max(nodeHeights(tree = tree)),"\n"))
        cat("source(\"simulationsTreeTopology.Rev\")\n\n")
        
        cat("# put it all together\nmvBNdata ~ dnPhyloBrownianREML(phylogeny, rep(1.0, numBranches), rep(rate, ntraits), nSites=ntraits)\nmvBNdata.clamp(means)\n\n")
        
        cat("# specify monitors\nsource(\"simulationsMonitors.Rev\")\n\n")
        
        cat("# define model\nmymodel = model(phylogeny)\n\n")
        
        cat("# run MCMC\nsource(\"simulationsMCMC.Rev\")\n\n")
        
        cat("# analyze mcmc output\nsource(\"simulationsAnalysis.Rev\")\n\n")
        
        cat("q()")
        
        sink()
        
        # add master script name to list
        
        jobName <- "job_uv_noisy_all.txt"
        if(counter == 1){
          scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
          sink(jobName, append = F)
          cat(scriptPath, "\n", sep = "")
          sink()
        } else {
          scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
          sink(jobName, append = T)
          cat(scriptPath, "\n", sep = "")
          sink()
        }
      }
      
    }
    
  }
  
}



