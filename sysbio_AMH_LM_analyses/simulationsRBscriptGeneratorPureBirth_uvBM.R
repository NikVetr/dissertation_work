# dir.create("uvBM_sims_PB")
setwd("/Users/nikolai/uvBM_sims_PB")
setwd("/Volumes/2TB/500/uvBM_sims_PB_noiseless_500/")

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


## preparing the files for noisy analysis ##
noisyness <- "noisy"
setwd("/Volumes/2TB/uvBM_sims_PB_noisy_8tips/")
setwd("/Volumes/2TB/500/uvBM_sims_PB_noisy_500/")


## preparing the files for noiseless analysis ##
noisyness <- "noiseless"
setwd("/Volumes/2TB/uvBM_sims_PB_noiseless_8tips/")
setwd("/Volumes/2TB/500/uvBM_sims_PB_noiseless_500/")

#remove old files if necessary
file.remove(paste0("trees/", list.files("trees/")))
dir.create("trees")
file.remove(paste0("data/", list.files("data/")))
dir.create("data")
file.remove(paste0("rates/", list.files("rates/")))
dir.create("rates")
file.remove(paste0("scripts/", list.files("scripts/")))
dir.create("scripts")

traitNumIter <- c(2, 4, 8, 16, 32, 64, 128)

degreeMisspecification <- "perfect"

degreeCorrelation <- c(0.2, 0.4, 0.6, 0.8, 0.99)

ntips <- 30

nrun <- 2

counter <- 0
# k iterates over degree of misspecification
for (k in 1:length(degreeCorrelation)) {
  
  # i iterates over replicates
  for(i in 1:500){
    
    # j iterates over trait numbers
    for (j in 1:length(traitNumIter)) {
      
      
      print(c(k,i,j)); counter <- counter + 1
      
      ntraits <- traitNumIter[j]
      
      cors <- matrix(nrow = ntraits, ncol = ntraits, data = degreeCorrelation[k])
      diag(cors) <- 1
      cov <- cors
      
      rateMatrixSpecification <- "perfect"
      
      replicateNum <- i
      
      fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_", tail(strsplit(as.character(degreeCorrelation[k]), split = "")[[1]], n = 1), "correlationStrength_", noisyness, "_replicate_", replicateNum)
      
      ratesPathTSV <- paste0("rates/", fileName, "_rates.tsv")
      ratesTSV(round(cov, 4), ratesPathTSV)
      
      # simulate the traits and write to file
      traitsPath <- paste0("data/", fileName, "_traits.nex")
      
      #use a pbtree
      tree <- tess.sim.taxa(n = 1, nTaxa = ntips, lambda = 1, mu = 0, samplingProbability = 1, max = 1e6)[[1]]
      tree$tip.label <- sample(x = tree$tip.label, size = length(tree$tip.label), replace = F)
      # plot(tree)
      
      traits <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = dim(cov)[1], sigma=round(cov, 4), mu= rep(0, times = dim(cov)[1])))
      
      #wiggle those wiggley traits
      if(noisyness == "noisy"){
        traits <- traits + rnorm(n = ntraits*length(tree$tip.label), sd = 0.5^0.5, mean = 0)
      }
      
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
        cat("corrStr =", tail(strsplit(as.character(degreeCorrelation[k]), split = "")[[1]], n = 1), "\n")
        cat(paste0("noisyness = \"", noisyness, "\"\n\n"))
        
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
        
        jobName <- "job_uv.txt"
        if(counter == 1 & l == 1){
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



