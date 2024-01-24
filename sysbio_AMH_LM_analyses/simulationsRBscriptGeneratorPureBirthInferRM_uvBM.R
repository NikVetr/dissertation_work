setwd("/Volumes/2TB/uvBM_sims_PB_InferRM_8tips")

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
library(rethinking)

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

#remove old files if necessary
file.remove(paste0("trees/", list.files("trees/")))
dir.create("trees")
file.remove(paste0("data/", list.files("data/")))
dir.create("data")
file.remove(paste0("rates/", list.files("rates/")))
dir.create("rates")
file.remove(paste0("scripts/", list.files("scripts/")))
dir.create("scripts")


## preparing the files for analysis ##

traitNumIter <- c(2, 4, 8, 16, 32, 64, 128)
traitNumIter <- c(2, 4, 8, 16, 32)

degreeMisspecification <- "inferred"

etas <- c(0.001, 0.01, 0.1, 1, 10)
etas <- c(1)

nrun <- 2
nreps <- 100

ntips <- 8

counter <- 0
# k iterates over degree of misspecification
for (k in 1:length(etas)) {
  
  # i iterates over replicates
  for(i in 1:nreps){
    
    # j iterates over trait numbers
    for (j in 1:length(traitNumIter)) {
      
      
      print(c(k,i,j)); counter <- counter + 1
      
      ntraits <- traitNumIter[j]
      
      avgRate <- runif(1,1,100)
      totalRate <- avgRate*ntraits
      vars <- as.vector(rdirichlet(1, rep(1, ntraits))) * totalRate
      eta <- etas[k]
      cors <- rlkjcorr(1, ntraits, eta = eta)
      # cors <- rcorvine(ntraits, eta = .001)
      # hist(abs(cors[upper.tri(cors)]))
      cov <- diag(sqrt(vars)) %*% cors %*% diag(sqrt(vars))
      
      rateMatrixSpecification <- "inferred"
      
      replicateNum <- i
      
      fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_eta1correlationStrength_", "replicate_", replicateNum)
      
      ratesPathTSV <- paste0("rates/", fileName, "_rates.tsv")
      ratesTSV(round(cov, 6), ratesPathTSV)
      
      # simulate the traits and write to file
      traitsPath <- paste0("data/", fileName, "_traits.nex")
      
      #use a pbtree
      tree <- tess.sim.taxa(n = 1, nTaxa = ntips, lambda = 1, mu = 0, samplingProbability = 1, max = 1e6)[[1]]
      tree$tip.label <- sample(x = tree$tip.label, size = length(tree$tip.label), replace = F)
      # plot(tree)
      
      traits <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = dim(cov)[1], sigma=round(cov, 6), mu= rep(0, times = dim(cov)[1])))
      
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
        cat("eta =", eta, "\n\n")
        
        cat("# read in the data\n")
        cat(paste0("means <- readContinuousCharacterData(\"data/", fileName, "_traits.nex\")\n\n"))
        
        cat("# get some useful variables from the data\nnumTips = means.ntaxa()\nnumBranches = 2 * numTips - 2\nntraits <- means.nchar()\ntaxa = means.taxa()\nnumNodes = numTips * 2 - 1\n\n")
        
        cat("# construct the rates prior\n")
        cat("proportional_var ~ dnDirichlet(rep(1, ntraits))\naverage_rate ~ dnUniform(1,100)\nrelative_var := proportional_var * ntraits * average_rate\nrelative_sd  := relative_var ^ 0.5\n")

        cat("moves[++mi] = mvBetaSimplex(proportional_var, tune = true, weight=30.0)\n")
        cat("moves[++mi] = mvScale(average_rate, tune = true, weight = 5)\nmoves[++mi] = mvSlide(average_rate, tune = true, weight = 5)\n\n")
        
        cat("# specify tree topology with BLs prior and moves\n")
        cat(paste0("speciation <- 1\n"))
        cat(paste0("extinction <- 0\n"))
        cat(paste0("sampFraction <- 1\n"))
        cat(paste0("rootAge <- ", max(nodeHeights(tree = tree)),"\n"))
        cat("source(\"simulationsTreeTopology.Rev\")\n\n")
        
        cat("# put it all together\nmvBNdata ~ dnPhyloBrownianREML(phylogeny, rep(1.0, numBranches), relative_var, nSites=ntraits)\nmvBNdata.clamp(means)\n\n")
        
        cat("# specify monitors\nsource(\"simulationsMonitors.Rev\")\n\n")
        
        cat("# define model\nmymodel = model(phylogeny)\n\n")
        
        cat("# run MCMC\nsource(\"simulationsMCMC.Rev\")\n\n")
        
        cat("# analyze mcmc output\nsource(\"simulationsAnalysis.Rev\")\n\n")
        
        cat("q()")
        
        sink()
        
        # add master script name to list
        
        jobName <- "job_pb.txt"
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

#side script to make seperate job files
njobs <- 6
newCounter <- 0
scriptsPerJob <- nreps * length(traitNumIter) * length(etas) * nrun / njobs
for (k in 1:length(etas)) {
  # i iterates over replicates
  for(i in 1:nreps){
    # j iterates over trait numbers
    for (j in 1:length(traitNumIter)) {
      
      ntraits <- traitNumIter[j]
      rateMatrixSpecification <- "inferred"
      replicateNum <- i
      fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_eta1correlationStrength_", "replicate_", replicateNum)
      
      # write the master scripts
      for (l in 1:nrun) {
        newCounter <- newCounter + 1
        if(newCounter == 0){
          jobNames <- paste0("job_pb", 1:njobs, ".txt")
          for(dels in 1:length(jobNames)){
            if (file.exists(jobNames[dels])) file.remove(jobNames[dels])
          }
        }
        if(newCounter != (counter*2)){jobNum <- newCounter %/% scriptsPerJob + 1}else{jobNum <- njobs}
        # print(paste(jobNum, fileName))
        jobName <- paste0("job_pb", jobNum, ".txt")
        scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
        sink(jobName, append = T)
        cat(scriptPath, "\n", sep = "")
        sink()
      }
      
    }
    
  }
  
}
