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


#################################################
###### Getting Modern Human Var-Cov Matrix ######
#################################################

d.orig <- read.csv("C:\\Users\\Nikolai\\Documents\\data\\Howell.csv")
d <- d.orig
# str(d)
# d[d == 0] <- NA
# d <- d[complete.cases(d),]
men <- d[d$Sex == "M",]
#men <- d
#head(as.matrix(men))
#unique(men$Population)[1]
pops <- sapply(1:length(unique(men$Population)), function (x) as.matrix((t(men[men$Population == unique(men$Population)[x],-(1:4)]))))

#get traits
traits <- list()
for(i in 1:82){ 
  print(i)
  trait <- list()
  trait[[1]] <- pops[[1]][i,]
  for (j in 2:30){
    trait[[j]] <- pops[[j]][i,]
  } 
  traits[[i]] <- trait
}

#make cov matrix
cov <- matrix(nrow = 82, ncol = 82)
for (i in 1:82){
  print(i)
  for(j in 1:82){
    trait1 <- traits[[i]]
    trait2 <- traits[[j]]
    trait1means <- sapply(1:length(trait1), function(x) mean(trait1[[x]]))
    trait2means <- sapply(1:length(trait2), function(x) mean(trait2[[x]]))
    trait1deviations <- sapply(1:length(trait1), function(x) trait1[[x]] - trait1means[x])
    trait2deviations <- sapply(1:length(trait2), function(x) trait2[[x]] - trait2means[x])
    deviationsProducts <- sapply(1:length(trait1), function(x) trait1deviations[[x]] * trait2deviations[[x]])
    sumALL <- 0
    for (k in (1:length(deviationsProducts))){
      sumALL <- sumALL + sum(deviationsProducts[[k]])
    }
    indiv <- sum(sapply(1:length(trait1), function(x) length(trait1[[x]])))
    covariance <- sumALL/(indiv - length(deviationsProducts))
    cov[i,j] <- covariance
  }
}

npop <- length(pops)
#compute matrix of population means
traitsHOWELL <- matrix(nrow=npop, ncol=82)
for(i in 1:npop){
  for(j in 1:82){
    traitsHOWELL[i,j] <- mean(pops[[i]][j,])
  }
}
rownames(traitsHOWELL) <- unique(men$Population)
colnames(traitsHOWELL) <- colnames(men)[-(1:4)]

rownames(cov) <- colnames(cov) <- colnames(men)[-(1:4)]
covHuman <- cov

#taken from roseman and weaver 2004, only BPC -> BPL (due to a typo? it's the basion-prosthion length )
linMeasNames <- c("GOL", "NOL", "BNL", "BBH", "XCB", "XFB", "STB", "ZYB", "AUB", "WCB", "ASB", "BPL", "NPH", "NLH", "OBH", "OBB", "JUB", "NLB", "MAB", "MDH", "MDB", "ZMB", "SSS", "FMB", "NAS", "EKB", "DKB", "NDS", "WNB", "SIS", "IML", "XML", "MLS", "WMH", "FOL", "FRC", "FRS", "FRF", "PAC", "PAS", "PAF", "OCC", "OCS", "OCF", "VRR", "NAR", "SSR", "PRR", "DKR", "ZOR", "FMR", "EKR", "ZMR", "AVR", "DKS", "SOS", "GLS")
linMeasTraitsHuman <- traitsHOWELL[,linMeasNames]
covFULL <- linMeasCovHuman <- covHuman[linMeasNames, linMeasNames]

##########################################
###### Getting Chimp Var-Cov Matrix ######
##########################################


############# DATA PREPARATION #############

d.orig <- read.csv("C:\\Users\\Nikolai\\Documents\\data\\PanSorted.csv")
d <- d.orig
men <- d[d$Sex == "M",]
d <- men
str(d)


head(as.matrix(d))
unique(d$Taxon)[1]
sp <- sapply(1:length(unique(d$Taxon)), function (x) as.matrix((t(d[d$Taxon == unique(d$Taxon)[x],-(1:2)]))))
ntraits <- 27
nspecies <- 4

#get traits
traits <- list()
for(i in 1:ntraits){ 
  print(i)
  trait <- list()
  trait[[1]] <- sp[[1]][i,]
  for (j in 2:nspecies){
    trait[[j]] <- sp[[j]][i,]
  } 
  traits[[i]] <- trait
}

#make cov matrix
cov <- matrix(nrow = ntraits, ncol = ntraits)
for (i in 1:ntraits){
  print(i)
  for(j in 1:ntraits){
    trait1 <- traits[[i]]
    trait2 <- traits[[j]]
    trait1means <- sapply(1:length(trait1), function(x) mean(trait1[[x]]))
    trait2means <- sapply(1:length(trait2), function(x) mean(trait2[[x]]))
    trait1deviations <- sapply(1:length(trait1), function(x) trait1[[x]] - trait1means[x])
    trait2deviations <- sapply(1:length(trait2), function(x) trait2[[x]] - trait2means[x])
    deviationsProducts <- sapply(1:length(trait1), function(x) trait1deviations[[x]] * trait2deviations[[x]])
    sumALL <- 0
    for (k in (1:length(deviationsProducts))){
      sumALL <- sumALL + sum(deviationsProducts[[k]])
    }
    indiv <- sum(sapply(1:length(trait1), function(x) length(trait1[[x]])))
    covariance <- sumALL/(indiv - length(deviationsProducts))
    cov[i,j] <- covariance
  }
}

offDiagCorr <- c(cov2cor(cov))
offDiagCorr <- offDiagCorr[offDiagCorr < 1] 
par(mfrow = c(1,2))
hist(offDiagCorr, main = "")
title("Correlations Between Traits")
hist(diag(cov), main = "")
title("Variances of Traits")

npop <- length(sp)
#compute matrix of Species means
traitsCHIMP <- matrix(nrow=npop, ncol=ntraits)
for(i in 1:npop){
  for(j in 1:ntraits){
    traitsCHIMP[i,j] <- mean(sp[[i]][j,])
  }
}
rownames(traitsCHIMP) <- unique(d$Taxon)
colnames(traitsCHIMP) <- colnames(d)[-(1:2)]

rownames(cov) <- colnames(cov) <- colnames(d)[-(1:2)]
covChimp <- cov


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

extendedOnionExpansion <- function(cov, d, eta){
  sds <- diag(cov)^0.5
  cor <- cov2cor(cov)
  a <- dim(cor)[1]
  beta <- eta + (d-2)/2
  beta <- beta - (1/2)*(a-2)
  for(k in (a+1):d){
    beta <- beta - 1/2
    y = rbeta(n = 1, shape1 = (k-1)/2, shape2 = beta)
    r <- y^0.5
    theta <- rnorm(n = (k-1), mean = 0.5, sd = 1)
    theta <- theta / (sum(theta^2)^0.5)
    w <- r * theta
    L <- t(chol(cor))
    q <- as.vector(L %*% w)
    traitNames <- c(rownames(cor), paste0("trait", k)) 
    cor <- rbind(cor, q)
    cor <- cbind(cor, c(q,1))
    rownames(cor) <- colnames(cor) <- traitNames
  }
  newSDS <- c(sds, sample(x = sds, size = (d-a), replace = T))
  ncov <- diag(newSDS) %*% cor %*% diag(newSDS)
  ncov
}

#degrees of freedom for each imperfection condition 
dfMatrix <- matrix(c(0,0,0,0,0,0,58,67,76,86,98,118,10,16,26,41,62,114), nrow = 6, ncol = 3)

#getting the trees
trees1 <- read.tree(file = "A:\\Google Drive\\outputGD\\howellsLogUnifMC3_3c_Sims1.trees")
trees1 <- trees1[2001:10000] #discard burn-in
# trees1 <- trees1[seq(0, 30000, by = 5)] #thin further

trees2 <- read.tree(file = "A:\\Google Drive\\outputGD\\howellsLogUnifMC3_3c_Sims2.trees")
trees2 <- trees2[2001:10000] #discard burn-in
# trees2 <- trees2[seq(0, 30000, by = 5)] #thin further

trees <- c(trees1, trees2)

## actually preparing the files ##

traitNumIter <- c(2,4,8,16,32,64)

degreeMisspecification <- c("perfect", "interregional", "chimpHuman")

nrun <- 2

# i iterates over replicates
for(i in 51:100){
  
  # k iterates over degree of misspecification
  for (k in 1:length(degreeMisspecification)) {
    
    
    # j iterates over trait numbers
    for (j in 1:length(traitNumIter)) {
      
      
      print(c(j, k, i))
      
      ntraits <- traitNumIter[j]
      if (ntraits < 58){
        traitsToUse <- sample(linMeasNames, size = ntraits, replace = F)
      } else {
        traitsToUse <- linMeasNames
      }
      cov <- covFULL[traitsToUse, traitsToUse]
      rateMatrixSpecification <- degreeMisspecification[k]
      replicateNum <- i
      
      # generate more traits if needed
      if (traitNumIter[j] > 57){
        cov <- extendedOnionExpansion(cov = cov, d = ntraits, eta = 1)
      }
      
      covOrig <- cov
      
      # wiggle traits if needed
      
      if (rateMatrixSpecification != "perfect"){
        dfWish <- dfMatrix[j,k]
        cov <- rwish(dfWish, cov)/dfWish
      }
      
      fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_", "replicate_", replicateNum)
      
      # write the covariance matrices
      ratesPath <- paste0("A:\\Google Drive\\mvBM_sims\\rates\\", fileName, "_rates.Phy")
      ratesPhy(round(covOrig, 4), ratesPath)
      
      ratesPathTSV <- paste0("A:\\Google Drive\\mvBM_sims\\rates\\", fileName, "_rates.tsv")
      ratesTSV(round(covOrig, 4), ratesPathTSV)
      
      # simulate the traits and write to file
      traitsPath <- paste0("A:\\Google Drive\\mvBM_sims\\data\\", fileName, "_traits.nex")
      tree <- trees[[sample(size = 1, x = 1:length(trees))]]
      tree <- root(tree, outgroup = "BUSHMAN", resolve.root = T)
      traits <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = dim(cov)[1], sigma=round(cov, 4), mu= rep(0, times = dim(cov)[1])))
      meansNexus(traits, traitsPath)
      treePath <- paste0("A:\\Google Drive\\mvBM_sims\\trees\\", fileName, "_tree.nex")
      write.nexus(tree, file = treePath, translate = T)
      treePathF <- paste0("A:\\Google Drive\\mvBM_sims\\trees\\translateEqualsFalse_", fileName, "_tree.nex")
      write.nexus(tree, file = treePathF, translate = F)
      
      # write the scripts for the burnin revbayes bit
      for (l in 1:nrun) {
        scriptPath <- paste0("A:\\Google Drive\\mvBM_sims\\scripts\\", fileName, "_burnin_script_run_", l, ".Rev")
        sink(scriptPath, append = F)
        
        cat("## Analysis of Simulated Data:", ntraits, "traits,", rateMatrixSpecification, "rate matrix specification, replication number", replicateNum, "##\n\n")
        cat("replicateNum <-", replicateNum, "\n")
        cat("runNum <-", l, "\n")
        cat(paste0("rateMatrixSpecification <- \"", rateMatrixSpecification, "\"\n\n"))
        cat("mi = 0\nseed(runNum)\n\n")
        
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
        
        # add RB script name to list, this doesn't actually get used anymore
        
        # scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
        # sink("A:\\Google Drive\\mvBM_sims\\jobRB.txt", append = T)
        # cat(scriptPath, "\n")
        # sink()
        
      }
      
      #write the scripts for the R bit, which calls the burnin script, gets tuned parameter values, writes a new RB script, and executes it 
      for (l in 1:nrun) {
        
        scriptPath <- paste0("A:\\Google Drive\\mvBM_sims\\scripts\\", fileName, "_script_run_", l, "_Rscript.R")
        scriptPathBurnin <- paste0("scripts/", fileName, "_burnin_script_run_", l, ".Rev")
        scriptPathMCMCMC <- paste0("scripts/", fileName, "_MCMCMC_script_run_", l, ".Rev")
        
        sink(scriptPath, append = F)
        
        cat("## Analysis of Simulated Data:", ntraits, "traits,", rateMatrixSpecification, "rate matrix specification, replication number", replicateNum, "##\n\n")
        
        cat("screenOutput <- system(command = \"~/repos/revbayes-development/projects/cmake/rb", scriptPathBurnin, "\", intern = T)\n")
        cat("screenOutput <- strsplit(screenOutput, split = \" \")\n")
        cat("screenOutput <- lapply(screenOutput, function (x) {x[!x == \"\"]} )\n\n")
        
        cat("firstWords <- sapply(1:length(screenOutput), function (x) screenOutput[x][[1]][1])\n\n")
        
        cat("slidingTunedParam <- screenOutput[which(firstWords == \"Sliding\")][[1]][which(screenOutput[which(firstWords == \"Sliding\")][[1]] == \"=\") + 1]\n")
        
        cat("betasimplexTunedParam <- screenOutput[which(firstWords == \"BetaSimplex\")][[1]][which(screenOutput[which(firstWords == \"BetaSimplex\")][[1]] == \"=\") + 1]\n")
        
        cat("dirichletsimplexTunedParam <- screenOutput[which(firstWords == \"DirichletSimplex\")][[1]][which(screenOutput[which(firstWords == \"DirichletSimplex\")][[1]] == \"=\") + 1]\n\n")
        
        cat(paste0("rateMatrixSpecification <- \"", rateMatrixSpecification, "\"\n\n"))
        
        cat("fileName <-", paste0("\"", fileName, "\""), "\n\n")
          
          cat(paste0("sink(\"", scriptPathMCMCMC, "\", append = F)\n\n"))
          
          cat("\tcat(\"replicateNum <-\",", replicateNum, ", \"\\n\")\n")
          
          cat("\tcat(\"runNum <-\",", l, ", \"\\n\")\n")
          
          cat("\tcat(paste0(\"rateMatrixSpecification <- \\\"\", rateMatrixSpecification, \"\\\"\\n\\n\"))\n")
          
          cat("\tcat(\"mi = 0\\nseed(runNum)\\n\\n\")\n")
          
          cat("\tcat(\"# read in the data\\n\")\n\n")
          cat("\tcat(paste0(\"means <- readContinuousCharacterData(\\\"data/\", fileName, \"_traits.nex\\\")\\n\\n\"))\n")
          
          cat("\tcat(\"# read in the rate matrix\\n\")\n\n")
          
          cat("\tcat(paste0(\"rateMatrix = readMatrix(\\\"rates/\", fileName, \"_rates.tsv\\\")\\n\\n\"))\n")
          
          cat("\tcat(\"# get some useful variables from the data\\nnumTips = means.ntaxa()\\nnumBranches = 2 * numTips - 3\\nntraits <- means.nchar()\\ntaxa = means.taxa()\\nnumNodes = numTips * 2 - 1\\n\\n\")\n\n")
          
          cat("\tcat(\"# specify tree topology prior and moves\\nsource(\\\"simulationsTreeTopology.Rev\\\")\\n\\n\")\n\n")
          
          cat("\tcat(\"slidingTunedParam <- \", slidingTunedParam, \"\\n\")\n")
          cat("\tcat(\"betasimplexTunedParam <- \", betasimplexTunedParam, \"\\n\")\n")
          cat("\tcat(\"dirichletsimplexTunedParam <- \", dirichletsimplexTunedParam, \"\\n\")\n")
          cat("\tcat(\"# specify branch length priors and moves\\nsource(\\\"simulationsBranchLengthsTuned.Rev\\\")\\n\\n\")\n\n")

          cat("\tcat(\"# assemble the tree\\nphylogeny := treeAssembly(topology, br_lens)\\n\\n\")\n\n")
          
          cat("\tcat(\"# put it all together\\nmvBNdata ~ dnPhyloMultivariateBrownianREML(phylogeny, rep(1.0, numBranches), rateMatrix, nSites=ntraits)\\nmvBNdata.clamp(means)\\n\\n\")\n\n")
          
          cat("\tcat(\"# specify monitors\\nsource(\\\"simulationsMonitors.Rev\\\")\\n\\n\")\n\n")
          
          cat("\tcat(\"# define model\\nmymodel = model(phylogeny)\\n\\n\")\n\n")
          
          cat("\tcat(\"# run MCMC\\nsource(\\\"simulationsMCMCMC.Rev\\\")\\n\\n\")\n\n")
          
          cat("\tcat(\"# analyze mcmc output\\nsource(\\\"simulationsAnalysis.Rev\\\")\\n\\n\")\n\n")
          
          cat("\tcat(\"q()\")\n\n")
          
          cat("sink()\n\n")
          
          cat("system(command = \"~/repos/revbayes-development/projects/cmake/rb", scriptPathMCMCMC, "\", intern = F)\n\n")
          
        cat("q()")
        
        sink()
        # add master script name to list
        
        scriptPath <- paste0("scripts/", fileName, "_script_run_", l, "_RScript.R")
        sink("A:\\Google Drive\\mvBM_sims\\jobR.txt", append = T)
        cat(scriptPath, "\n")
        sink()
        
      }
      
    }
    
  }
  
}
