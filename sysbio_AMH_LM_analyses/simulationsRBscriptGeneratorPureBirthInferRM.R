setwd("/Volumes/2TB/mvBM_sims_PB_inferRM_8tips/")

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
        
        cat("# construct the rate matrix prior\n")
        cat("proportional_var ~ dnDirichlet(rep(1, ntraits))\naverage_rate ~ dnUniform(1,100)\nrelative_var := proportional_var * ntraits * average_rate\nrelative_sd  := relative_var ^ 0.5\n")
        cat(paste0("partial_correlation_matrix ~ dnLKJPartial(eta, ntraits)\n"))
        cat("correlation_matrix := fnPartialToCorr(partial_correlation_matrix)\nrateMatrix := fnDecompVarCovar(relative_sd, correlation_matrix)\ncorrelation := fnUpperTriangle(correlation_matrix)\n\n")
        
        cat("# specify moves on the rate matrix\n")
        cat("moves[++mi] = mvCorrelationMatrixRandomWalk(partial_correlation_matrix, sigma=0.1, weight=150.0, tune = true, tuneTarget=0.234)\n")
        cat("moves[++mi] = mvBetaSimplex(proportional_var, tune = true, weight=30.0)\n")
        cat("for(i in 1:(ntraits - 1)) {\n\tfor(j in (i+1):ntraits) {\n\t\tmoves[++mi] = mvCorrelationMatrixSpecificElementBeta(partial_correlation_matrix, i, j, alpha=10.0, weight=0.25, tune = true, tuneTarget=0.234)\n\t}\n}\n")
        cat("moves[++mi] = mvScale(average_rate, tune = true, weight = 5)\nmoves[++mi] = mvSlide(average_rate, tune = true, weight = 5)\n\n")
        
        
        cat("# specify tree topology with BLs prior and moves\n")
        cat(paste0("speciation <- 1\n"))
        cat(paste0("extinction <- 0\n"))
        cat(paste0("sampFraction <- 1\n"))
        cat(paste0("rootAge <- ", max(nodeHeights(tree = tree)),"\n"))
        cat("source(\"simulationsTreeTopology.Rev\")\n\n")
        
        cat("# put it all together\nmvBNdata ~ dnPhyloMultivariateBrownianREML(phylogeny, rep(1.0, numBranches), rateMatrix, nSites=ntraits)\nmvBNdata.clamp(means)\n\n")
        
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
    print(i)
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
        print(paste(jobNum, fileName))
        jobName <- paste0("job_pb", jobNum, ".txt")
        scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
        sink(jobName, append = T)
        cat(scriptPath, "\n", sep = "")
        sink()
      }
      
    }
    
  }
  
}

################################################
##
##
##
##
##
##
#modifying the script to do modules of 32 traits
## preparing the files for analysis ##
##
##
##
##
##
##
################################################


traitNumIter <- c(64, 128)

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
      
      if(ntraits > 32){
        cov <- list()
        for(mods in 1:(ntraits/32)){
          avgRate <- runif(1,1,100)
          totalRate <- avgRate*32
          vars <- as.vector(rdirichlet(1, rep(1, 32))) * totalRate
          eta <- etas[k]
          cors <- rlkjcorr(1, 32, eta = eta)
          # cors <- rcorvine(ntraits, eta = .001)

          cov[[mods]] <- diag(sqrt(vars)) %*% cors %*% diag(sqrt(vars))
        }
      }
      
      ntraits <- traitNumIter[j]
      
      rateMatrixSpecification <- "inferred"
      
      replicateNum <- i
      
      fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_eta1correlationStrength_", "replicate_", replicateNum)
      
      if(ntraits > 32){
        for(mods in 1:(ntraits/32)){
        ratesPathTSV <- paste0("rates/", fileName, "_mod", mods, "_rates.tsv")
        ratesTSV(round(cov[[mods]], 6), ratesPathTSV)
        }
      }

      
      #use a pbtree
      tree <- tess.sim.taxa(n = 1, nTaxa = ntips, lambda = 1, mu = 0, samplingProbability = 1, max = 1e6)[[1]]
      tree$tip.label <- sample(x = tree$tip.label, size = length(tree$tip.label), replace = F)
      # plot(tree)
      
      # simulate the traits and write to file
      if(ntraits > 32){
        traits <- list()
        for(mods in 1:(ntraits/32)){
        traits[[mods]] <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = dim(cov[[mods]])[1], sigma=round(cov[[mods]], 6), mu= rep(0, times = dim(cov[[mods]])[1])))
        traitsPath <- paste0("data/", fileName, "_mod", mods, "_traits.nex")
        meansNexus(traits[[mods]], traitsPath)
        }
      }
      
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
        cat("ntraitsTotal =", ntraits, "\n\n")
        
        
        cat("# read in the data\n")
        if(ntraits > 32){
          for(mods in 1:(ntraits/32)){
            cat(paste0("means", mods, " <- readContinuousCharacterData(\"data/", fileName, "_mod", mods, "_traits.nex\")\n"))
          }
        }
        
        cat(paste0("\n# get some useful variables from the data\n"))
        cat(paste0("numTips = means1.ntaxa()\n"))
        cat(paste0("numBranches = 2 * numTips - 2\n"))
        cat(paste0("ntraits <- means1.nchar()\n"))
        cat(paste0("taxa = means1.taxa()\n"))
        cat(paste0("numNodes = numTips * 2 - 1\n\n"))
        
        if(ntraits > 32){
          for(mods in 1:(ntraits/32)){
          cat("# construct the rate matrix prior", mods, "\n")
          cat(paste0("proportional_var", mods, " ~ dnDirichlet(rep(1, ntraits))\n"))
          cat(paste0("average_rate", mods, " ~ dnUniform(1,100)\n"))
          cat(paste0("relative_var", mods, " := proportional_var", mods, " * ntraits * average_rate", mods, "\n"))
          cat(paste0("relative_sd", mods, "  := relative_var", mods, " ^ 0.5\n"))
          cat(paste0("partial_correlation_matrix", mods, " ~ dnLKJPartial(eta, ntraits)\n"))
          cat(paste0("correlation_matrix", mods, " := fnPartialToCorr(partial_correlation_matrix", mods, ")\n"))
          cat(paste0("rateMatrix", mods, " := fnDecompVarCovar(relative_sd", mods, ", correlation_matrix", mods, ")\n"))
          cat(paste0("correlation", mods, " := fnUpperTriangle(correlation_matrix", mods, ")\n\n"))
          
          cat(paste0("# specify moves on the rate matrix\n"))
          cat(paste0("moves[++mi] = mvCorrelationMatrixRandomWalk(partial_correlation_matrix", mods, ", sigma=0.1, weight=150.0, tune = true, tuneTarget=0.234)\n"))
          cat(paste0("moves[++mi] = mvBetaSimplex(proportional_var", mods, ", tune = true, weight=30.0)\n"))
          cat(paste0("for(i in 1:(ntraits - 1)) {\n\tfor(j in (i+1):ntraits) {\n\t\tmoves[++mi] = mvCorrelationMatrixSpecificElementBeta(partial_correlation_matrix", mods, ", i, j, alpha=10.0, weight=0.25, tune = true, tuneTarget=0.234)\n\t}\n}\n"))
          cat(paste0("moves[++mi] = mvScale(average_rate", mods, ", tune = true, weight = 5)\nmoves[++mi] = mvSlide(average_rate", mods, ", tune = true, weight = 5)\n\n"))
          }
        }
        
        cat("# specify tree topology with BLs prior and moves\n")
        cat(paste0("speciation <- 1\n"))
        cat(paste0("extinction <- 0\n"))
        cat(paste0("sampFraction <- 1\n"))
        cat(paste0("rootAge <- ", max(nodeHeights(tree = tree)),"\n"))
        cat("source(\"simulationsTreeTopology.Rev\")\n\n")
        
        if(ntraits > 32){
          for(mods in 1:(ntraits/32)){
        cat(paste0("# put it all together ", mods, "\nmvBNdata", mods, " ~ dnPhyloMultivariateBrownianREML(phylogeny, rep(1.0, numBranches), rateMatrix", mods, ", nSites=ntraits)\nmvBNdata", mods, ".clamp(means", mods, ")\n\n"))
          }
        }
        
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
traitNumIter <- rev(sort(traitNumIter))
for (j in 1:length(traitNumIter)) {
  # j iterates over trait numbers
  for (k in 1:length(etas)) {
  # k iterates over etas
    for(i in 1:nreps){
      # i iterates over replicates  
      print(i)
      
      ntraits <- traitNumIter[j]
      rateMatrixSpecification <- "inferred"
      replicateNum <- i
      fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_eta1correlationStrength_", "replicate_", replicateNum)
      
      # write the master scripts
      for (l in 1:nrun) {
        if(newCounter == 0){
          jobNames <- paste0("job_pb", 1:njobs, ".txt")
          for(dels in 1:length(jobNames)){
            if (file.exists(jobNames[dels])) file.remove(jobNames[dels])
          }
        }
        newCounter <- newCounter + 1
        if(newCounter != (counter*2)){jobNum <- newCounter %/% scriptsPerJob + 1}else{jobNum <- njobs}
        print(paste(jobNum, fileName))
        jobName <- paste0("job_pb", jobNum, ".txt")
        scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
        sink(jobName, append = T)
        cat(scriptPath, "\n", sep = "")
        sink()
      }
      
    }
    
  }
  
}
