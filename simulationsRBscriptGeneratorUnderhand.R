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
lib

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

## preparing the files ##

traitNumIter <- c(128, 256)
traitNumIter <- c(2, 4, 8, 16, 32, 64)

degreeMisspecification <- "perfect"

degreeCorrelation <- c(0.2, 0.4, 0.6, 0.8, 0.999)

nrun <- 2

# k iterates over degree of misspecification
for (k in 1:length(degreeCorrelation)) {
  
  # i iterates over replicates
  for(i in 16:20){
    
    # j iterates over trait numbers
    for (j in 1:length(traitNumIter)) {
      
      
      print(c(j, k, i))
      
      ntraits <- traitNumIter[j]
      
      cors <- matrix(nrow = ntraits, ncol = ntraits, data = degreeCorrelation[k])
      diag(cors) <- 1
      # sds <- sample(1:20, ntraits, replace = T)
      # cov <- diag(sds) %*% cors %*% diag(sds)
      cov <- cors
      
      rateMatrixSpecification <- "perfect"
      
      replicateNum <- i
      
      fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_", tail(strsplit(as.character(degreeCorrelation[k]), split = "")[[1]], n = 1), "correlationStrength_", "replicate_", replicateNum)
      
      # write the covariance matrices
      # ratesPath <- paste0("A:\\Google Drive\\mvBM_sims_underhand\\rates\\", fileName, "_rates.Phy")
      # ratesPhy(round(cov, 4), ratesPath)
      
      ratesPathTSV <- paste0("A:\\Google Drive\\mvBM_sims_underhand\\rates\\", fileName, "_rates.tsv")
      ratesTSV(round(cov, 4), ratesPathTSV)
      
      # simulate the traits and write to file
      traitsPath <- paste0("A:\\Google Drive\\mvBM_sims_underhand\\data\\", fileName, "_traits.nex")
      
      #use softball tree
      tree <- read.tree(text = "(((((t1:.1, t2:.1):.1,(t3:.1,t4:.1):.1):.1,((t5:.1, t6:.1):.1,(t7:.1,t8:.1):.1):.1):.1,(((t9:.1, t10:.1):.1,(t11:.1,t12:.1):.1):.1,((t13:.1, t14:.1):.1,(t15:.1,t16:.1):.1):.1):.1):.1,((((t17:.1, t18:.1):.1,(t19:.1,t20:.1):.1):.1,((t21:.1, t22:.1):.1,(t23:.1,t24:.1):.1):.1):.1,(((t25:.1, t26:.1):.1,(t27:.1,t28:.1):.1):.1,(t29:.1, t30:.1):.2):.1):.1):.1;")
      tree$tip.label <- sample(x = tree$tip.label, size = length(tree$tip.label), replace = F)
      # plot(tree)
      
      traits <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = dim(cov)[1], sigma=round(cov, 4), mu= rep(0, times = dim(cov)[1])))
      meansNexus(traits, traitsPath)
      treePath <- paste0("A:\\Google Drive\\mvBM_sims_underhand\\trees\\", fileName, "_tree.nex")
      write.nexus(tree, file = treePath, translate = T)
      # treePathF <- paste0("A:\\Google Drive\\mvBM_sims_underhand\\trees\\translateEqualsFalse_", fileName, "_tree.nex")
      # write.nexus(tree, file = treePathF, translate = F)
      
      # write the master scripts
      for (l in 1:nrun) {
        scriptPath <- paste0("A:\\Google Drive\\mvBM_sims_underhand\\scripts\\", fileName, "_script_run_", l, ".Rev")
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
        
        cat("# get some useful variables from the data\nnumTips = means.ntaxa()\nnumBranches = 2 * numTips - 3\nntraits <- means.nchar()\ntaxa = means.taxa()\nnumNodes = numTips * 2 - 1\n\n")
        
        cat("# specify tree topology prior and moves\nsource(\"simulationsTreeTopology.Rev\")\n\n")
        
        cat("# specify branch length priors and moves\nsource(\"simulationsBranchLengths.Rev\")\n\n")
        
        cat("# assemble the tree\nphylogeny := treeAssembly(topology, br_lens)\n\n")
        
        cat("# put it all together\nmvBNdata ~ dnPhyloMultivariateBrownianREML(phylogeny, rep(1.0, numBranches), rateMatrix, nSites=ntraits)\nmvBNdata.clamp(means)\n\n")
        
        cat("# specify monitors\nsource(\"simulationsMonitors.Rev\")\n\n")
        
        cat("# define model\nmymodel = model(phylogeny)\n\n")
        
        cat("# run MCMC\nsource(\"simulationsMCMC.Rev\")\n\n")
        
        cat("# analyze mcmc output\nsource(\"simulationsAnalysis.Rev\")\n\n")
        
        cat("q()")
        
        sink()
        
        # add master script name to list
        
        scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
        sink("A:\\Google Drive\\mvBM_sims_underhand\\job_highCorr_extra.txt", append = T)
        cat(scriptPath, "\n")
        sink()
        
      }
      
    }
    
  }
  
}


## preparing the files for noisy analysis ##

traitNumIter <- c(128, 256)
traitNumIter <- c(2, 4, 8, 16, 32, 64)

degreeMisspecification <- "perfect"

degreeCorrelation <- c(0.2, 0.4, 0.6, 0.8, 0.999)

nrun <- 2

# k iterates over degree of misspecification
for (k in 1:length(degreeCorrelation)) {
  
  # i iterates over replicates
  for(i in 211:220){
    
    # j iterates over trait numbers
    for (j in 1:length(traitNumIter)) {
      
      
      print(c(j, k, i))
      
      ntraits <- traitNumIter[j]
      
      cors <- matrix(nrow = ntraits, ncol = ntraits, data = degreeCorrelation[k])
      diag(cors) <- 1
      # sds <- sample(1:20, ntraits, replace = T)
      # cov <- diag(sds) %*% cors %*% diag(sds)
      cov <- cors
      
      rateMatrixSpecification <- "perfect"
      
      replicateNum <- i
      
      fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_", tail(strsplit(as.character(degreeCorrelation[k]), split = "")[[1]], n = 1), "correlationStrength_", "replicate_", replicateNum)
      
      # write the covariance matrices
      # ratesPath <- paste0("A:\\Google Drive\\mvBM_sims_underhand\\rates\\", fileName, "_rates.Phy")
      # ratesPhy(round(cov, 4), ratesPath)
      
      ratesPathTSV <- paste0("A:\\Google Drive\\mvBM_sims_underhand\\rates2\\", fileName, "_rates.tsv")
      ratesTSV(round(cov, 4), ratesPathTSV)
      
      # simulate the traits and write to file
      traitsPath <- paste0("A:\\Google Drive\\mvBM_sims_underhand\\data2\\", fileName, "_traits.nex")
      
      #use softball tree
      tree <- read.tree(text = "(((((t1:.1, t2:.1):.1,(t3:.1,t4:.1):.1):.1,((t5:.1, t6:.1):.1,(t7:.1,t8:.1):.1):.1):.1,(((t9:.1, t10:.1):.1,(t11:.1,t12:.1):.1):.1,((t13:.1, t14:.1):.1,(t15:.1,t16:.1):.1):.1):.1):.1,((((t17:.1, t18:.1):.1,(t19:.1,t20:.1):.1):.1,((t21:.1, t22:.1):.1,(t23:.1,t24:.1):.1):.1):.1,(((t25:.1, t26:.1):.1,(t27:.1,t28:.1):.1):.1,(t29:.1, t30:.1):.2):.1):.1):.1;")
      tree$tip.label <- sample(x = tree$tip.label, size = length(tree$tip.label), replace = F)
      # plot(tree)
      
      traits <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = dim(cov)[1], sigma=round(cov, 4), mu= rep(0, times = dim(cov)[1])))
      traits <- traits + rnorm(n = ntraits*length(tree$tip.label), sd = 0.02, mean = 0)
      
      meansNexus(traits, traitsPath)
      treePath <- paste0("A:\\Google Drive\\mvBM_sims_underhand\\trees2\\", fileName, "_tree.nex")
      write.nexus(tree, file = treePath, translate = T)
      # treePathF <- paste0("A:\\Google Drive\\mvBM_sims_underhand\\trees\\translateEqualsFalse_", fileName, "_tree.nex")
      # write.nexus(tree, file = treePathF, translate = F)
      
      # write the master scripts
      for (l in 1:nrun) {
        scriptPath <- paste0("A:\\Google Drive\\mvBM_sims_underhand\\scripts2\\", fileName, "_script_run_", l, ".Rev")
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
        
        cat("# get some useful variables from the data\nnumTips = means.ntaxa()\nnumBranches = 2 * numTips - 3\nntraits <- means.nchar()\ntaxa = means.taxa()\nnumNodes = numTips * 2 - 1\n\n")
        
        cat("# specify tree topology prior and moves\nsource(\"simulationsTreeTopology.Rev\")\n\n")
        
        cat("# specify branch length priors and moves\nsource(\"simulationsBranchLengths.Rev\")\n\n")
        
        cat("# assemble the tree\nphylogeny := treeAssembly(topology, br_lens)\n\n")
        
        cat("# put it all together\nmvBNdata ~ dnPhyloMultivariateBrownianREML(phylogeny, rep(1.0, numBranches), rateMatrix, nSites=ntraits)\nmvBNdata.clamp(means)\n\n")
        
        cat("# specify monitors\nsource(\"simulationsMonitors.Rev\")\n\n")
        
        cat("# define model\nmymodel = model(phylogeny)\n\n")
        
        cat("# run MCMC\nsource(\"simulationsMCMC.Rev\")\n\n")
        
        cat("# analyze mcmc output\nsource(\"simulationsAnalysis.Rev\")\n\n")
        
        cat("q()")
        
        sink()
        
        # add master script name to list
        
        scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
        sink("A:\\Google Drive\\mvBM_sims_underhand\\job_noise_bm.txt", append = T)
        cat(scriptPath, "\n")
        sink()
        
      }
      
    }
    
  }
  
}





## preparing the files for perfectly correlated analyses ##

traitNumIter <- c(128, 256)
traitNumIter <- c(2, 4, 8, 16, 32, 64)

degreeMisspecification <- "perfect"

degreeCorrelation <- c(1.0)

nrun <- 2

# k iterates over degree of misspecification
for (k in 1:length(degreeCorrelation)) {
  
  # i iterates over replicates
  for(i in 1:10){
    
    # j iterates over trait numbers
    for (j in 1:length(traitNumIter)) {
      
      
      print(c(j, k, i))
      
      ntraits <- traitNumIter[j]
      
      cors <- matrix(nrow = ntraits, ncol = ntraits, data = degreeCorrelation[k])
      diag(cors) <- 1
      # sds <- sample(1:20, ntraits, replace = T)
      # cov <- diag(sds) %*% cors %*% diag(sds)
      cov <- cors
      
      rateMatrixSpecification <- "perfect"
      
      replicateNum <- i
      
      fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_", tail(strsplit(as.character(degreeCorrelation[k]), split = "")[[1]], n = 1), "correlationStrength_", "replicate_", replicateNum)
      
      # write the covariance matrices
      # ratesPath <- paste0("A:\\Google Drive\\mvBM_sims_underhand\\rates\\", fileName, "_rates.Phy")
      # ratesPhy(round(cov, 4), ratesPath)
      
      ratesPathTSV <- paste0("A:\\Google Drive\\mvBM_sims_underhand\\rates2\\", fileName, "_rates.tsv")
      ratesTSV(round(cov, 4), ratesPathTSV)
      
      # simulate the traits and write to file
      traitsPath <- paste0("A:\\Google Drive\\mvBM_sims_underhand\\data2\\", fileName, "_traits.nex")
      
      #use softball tree
      tree <- read.tree(text = "(((((t1:.1, t2:.1):.1,(t3:.1,t4:.1):.1):.1,((t5:.1, t6:.1):.1,(t7:.1,t8:.1):.1):.1):.1,(((t9:.1, t10:.1):.1,(t11:.1,t12:.1):.1):.1,((t13:.1, t14:.1):.1,(t15:.1,t16:.1):.1):.1):.1):.1,((((t17:.1, t18:.1):.1,(t19:.1,t20:.1):.1):.1,((t21:.1, t22:.1):.1,(t23:.1,t24:.1):.1):.1):.1,(((t25:.1, t26:.1):.1,(t27:.1,t28:.1):.1):.1,(t29:.1, t30:.1):.2):.1):.1):.1;")
      tree$tip.label <- sample(x = tree$tip.label, size = length(tree$tip.label), replace = F)
      # plot(tree)
      
      traits <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = 1, sigma=1, mu= 0))
      traits <- matrix(nrow = nrow(traits), ncol = ntraits, data = rep(traits, ntraits), dimnames = list(rownames(traits)))
      meansNexus(traits, traitsPath)
      treePath <- paste0("A:\\Google Drive\\mvBM_sims_underhand\\trees2\\", fileName, "_tree.nex")
      write.nexus(tree, file = treePath, translate = T)
      # treePathF <- paste0("A:\\Google Drive\\mvBM_sims_underhand\\trees\\translateEqualsFalse_", fileName, "_tree.nex")
      # write.nexus(tree, file = treePathF, translate = F)
      
      # write the master scripts
      for (l in 1:nrun) {
        scriptPath <- paste0("A:\\Google Drive\\mvBM_sims_underhand\\scripts2\\", fileName, "_script_run_", l, ".Rev")
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
        
        cat("# get some useful variables from the data\nnumTips = means.ntaxa()\nnumBranches = 2 * numTips - 3\nntraits <- means.nchar()\ntaxa = means.taxa()\nnumNodes = numTips * 2 - 1\n\n")
        
        cat("# specify tree topology prior and moves\nsource(\"simulationsTreeTopology.Rev\")\n\n")
        
        cat("# specify branch length priors and moves\nsource(\"simulationsBranchLengths.Rev\")\n\n")
        
        cat("# assemble the tree\nphylogeny := treeAssembly(topology, br_lens)\n\n")
        
        cat("# put it all together\nmvBNdata ~ dnPhyloMultivariateBrownianREML(phylogeny, rep(1.0, numBranches), rateMatrix, nSites=ntraits)\nmvBNdata.clamp(means)\n\n")
        
        cat("# specify monitors\nsource(\"simulationsMonitors.Rev\")\n\n")
        
        cat("# define model\nmymodel = model(phylogeny)\n\n")
        
        cat("# run MCMC\nsource(\"simulationsMCMC.Rev\")\n\n")
        
        cat("# analyze mcmc output\nsource(\"simulationsAnalysis.Rev\")\n\n")
        
        cat("q()")
        
        sink()
        
        # add master script name to list
        
        scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
        sink("A:\\Google Drive\\mvBM_sims_underhand\\job_highCorr_extra_extra.txt", append = T)
        cat(scriptPath, "\n")
        sink()
        
      }
      
    }
    
  }
  
}
