# CREATING REVBAYES REPLICATE SCRIPTS
setwd("/Volumes/2TB/uvBM_sims_Howells_8tips")
setwd("/Volumes/2TB/500/uvBM_sims_Howells_500/")

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

d.orig <- read.csv("Howell.csv")
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
targetRegionalMBD <- c(0.3016502, 0.5492471, 1.0165738, 1.9375714, 3.8145961, 7.8122877) #(2,4,8,16,32,64)
targetRegionalMBD <- c(0.3016502, 0.5492471, 1.0165738, 1.9375714, 3.8145961, 6.939882, 7.8122877) #incl 57
targetRegionalMBD <- c(0.3016502, 0.5492471, 1.0165738, 1.9375714, 3.8145961, 6.939882, 7.8122877) #(2,4,8,16,32,57)

targetChimpHumanMBD <- c(0.7828651, 1.2172830, 1.9157832, 3.1369000, 5.3419662, 7.6690310, 8.0803367) #(2,4,8,16,32,57,64)
targetChimpHumanMBD <- c(0.7828651, 1.2172830, 1.9157832, 3.1369000, 5.3419662, 7.6690310) #(2,4,8,16,32,57)


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


#########################################################################################
## figure out how many degrees of freedom are needed to draw from Wishart distribution ##
#########################################################################################

traitNumIter <- c(2,4,8,16,32,57)

wishartSamplingIntensity <- 5000

traitNumIter <- c(2,4,8,16,32,57)

dfMatrix <- matrix(nrow = length(traitNumIter), ncol = length(degreeMisspecification))

degreeMisspecification <- c("perfect", "interregional", "chimpHuman")


# k iterates over degree of misspecification
for (k in 1:length(degreeMisspecification)) {
  
  # j iterates over trait numbers
  for (j in 1:length(traitNumIter)) {
    
    print(c(k,j))
    
    
    rateMatrixSpecification <- degreeMisspecification[k]
    ntraits <- traitNumIter[j]
    
    if (ntraits < 58){
      traitsToUse <- sample(linMeasNames, size = ntraits, replace = F)
    } else {
      traitsToUse <- linMeasNames
    }
    
    cov <- covFULL[traitsToUse, traitsToUse]
    
    if (ntraits > 57){
      cov <- extendedOnionExpansion(cov = cov, d = ntraits, eta = 1)
    }
    
    if (rateMatrixSpecification == "interregional"){
      targetMBD <- targetRegionalMBD[j]
    } else if (rateMatrixSpecification == "chimpHuman"){
      targetMBD <- targetChimpHumanMBD[j]
    } 
    
    
    avgDist <- vector()
    dist <- vector()
    dfWish <- ntraits
    for(r in 1:wishartSamplingIntensity){
      sampCov <- rwish(dfWish, cov)/dfWish
      dist[r] <- sum(log(eigen(solve(sampCov) %*% cov)$values)^2)^0.5
    }
    avgDist[1] <- mean(dist)
    
    while(avgDist[length(avgDist)] > targetMBD){
      dfWish <- dfWish + 20
      dist <- vector()  
      for(r in 1:wishartSamplingIntensity){
        sampCov <- rwish(dfWish, cov)/dfWish
        dist[r] <- sum(log(eigen(solve(sampCov) %*% cov)$values)^2)^0.5
      }
      # print(mean(dist))
      # print(sd(dist)/sqrt(length(dist)))
      avgDist <- c(avgDist, mean(dist))
    }
    
    dfWish <- dfWish - 20
    avgDist <- avgDist[-length(avgDist)]
    
    while(avgDist[length(avgDist)] > targetMBD){
      dfWish <- dfWish + 10
      dist <- vector()  
      for(r in 1:wishartSamplingIntensity){
        sampCov <- rwish(dfWish, cov)/dfWish
        dist[r] <- sum(log(eigen(solve(sampCov) %*% cov)$values)^2)^0.5
      }
      #print(mean(dist))
      avgDist <- c(avgDist, mean(dist))
    }
    
    dfWish <- dfWish - 10
    avgDist <- avgDist[-length(avgDist)]
    
    while(avgDist[length(avgDist)] > targetMBD){
      dfWish <- dfWish + 5
      dist <- vector()  
      for(r in 1:wishartSamplingIntensity){
        sampCov <- rwish(dfWish, cov)/dfWish
        dist[r] <- sum(log(eigen(solve(sampCov) %*% cov)$values)^2)^0.5
      }
      #print(mean(dist))
      avgDist <- c(avgDist, mean(dist))
    }
    
    dfWish <- dfWish - 5
    avgDist <- avgDist[-length(avgDist)]
    
    while(avgDist[length(avgDist)] > targetMBD){
      dfWish <- dfWish + 1
      dist <- vector()  
      for(r in 1:wishartSamplingIntensity){
        sampCov <- rwish(dfWish, cov)/dfWish
        dist[r] <- sum(log(eigen(solve(sampCov) %*% cov)$values)^2)^0.5
      }
      #print(mean(dist))
      avgDist <- c(avgDist, mean(dist))
    }
    
    if (which.min(abs(avgDist - targetMBD)) != length(avgDist)) {
      dfWish <- dfWish - 1
    }
    
    dfMatrix[j,k] <- dfWish
    
  }
  
}

##########################################
## figure out df when traitnum includes 57
##########################################

#w/out 57
dfMatrix <- matrix(c(0,0,0,0,0,0,58,67,76,86,98,118,10,16,26,41,62,114), nrow = 7, ncol = 3)

#incl 57
dfMatrix <- matrix(c(0,0,0,0,0,0,0,58,67,76,86,98,113,118,10,16,26,41,62,101,114), nrow = 7, ncol = 3)

#incl 57 but not 64
dfMatrix <- matrix(c(0,0,0,0,0,0,58,67,76,86,98,113,10,16,26,41,62,101), nrow = 6, ncol = 3)


#getting the trees
howellsTrees1 <- read.tree(file = "~/output/howells_nonlog.trees")

howellsTrees2 <- read.tree(file = "~/output/howells_nonlog2.trees")

trees <- c(howellsTrees1, howellsTrees2)


## actually preparing the files ##

traitNumIter <- c(2,4,8,16,32,57)

degreeMisspecification <- c("perfect", "interregional", "chimpHuman")

processMisspecification <- "uvBM"

# ntips <- 8
# setwd("/Volumes/2TB/uvBM_sims_Howells_8tips")
ntips <- 30
setwd("/Volumes/2TB/500/uvBM_sims_Howells_500/")


#remove old files if necessary
file.remove(paste0("trees/", list.files("trees/")))
dir.create("trees")
file.remove(paste0("data/", list.files("data/")))
dir.create("data")
file.remove(paste0("rates/", list.files("rates/")))
dir.create("rates")
file.remove(paste0("scripts/", list.files("scripts/")))
dir.create("scripts")

nrun <- 2
noisyness <- "noiseless"

# i iterates over replicates
for(i in 1:500){
  
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
      if(processMisspecification == "uvBM"){
        covOrigTrue <- covOrig
        covOrig <- diag(diag(covOrig))
      }
      
      ratesPath <- paste0("rates/", fileName, "_rates.Phy")
      ratesPhy(round(covOrig, 4), ratesPath)
      
      ratesPathTSV <- paste0("rates/", fileName, "_rates.tsv")
      ratesTSV(round(covOrig, 4), ratesPathTSV)
      
      # simulate the traits and write to file
      traitsPath <- paste0("data/", fileName, "_traits.nex")
      
      #use howells trees
      tree <- trees[[sample(size = 1, x = 1:length(trees))]]
      if(ntips < 30){
        tipsToKeep <- sample(tree$tip.label, ntips, replace = F)
        tipsToDrop <- setdiff(tree$tip.label, tipsToKeep)
        tree <- drop.tip(phy = tree, tip = tipsToDrop)
        tree <- root(tree, outgroup = tree$tip.label[1], resolve.root = T)
      } else {tree <- root(tree, outgroup = "BUSHMAN", resolve.root = T)}
      
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
        cat(paste0("noisyness = \"", noisyness, "\"\n\n"))
        cat("mi = 0\nseed(runNum)\n\n")
        
        cat("# read in the data\n")
        cat(paste0("means <- readContinuousCharacterData(\"data/", fileName, "_traits.nex\")\n\n"))
        
        cat("# read in the rate matrix\n")
        cat(paste0("rate = readDistanceMatrix(\"rates/", fileName, "_rates.Phy\")\n"))
        cat("rateMatrix = rate.symmetricMatrix()\n\n")
        
        cat("# get some useful variables from the data\nnumTips = means.ntaxa()\nnumBranches = 2 * numTips - 3\nntraits <- means.nchar()\ntaxa = means.taxa()\nnumNodes = numTips * 2 - 1\n\n")
        
        cat("# specify tree topology prior and moves\nsource(\"simulationsTreeTopology.Rev\")\n\n")
        
        cat("# specify branch length priors and moves\nsource(\"simulationsBranchLengths.Rev\")\n\n")
        
        cat("# assemble the tree\nphylogeny := treeAssembly(topology, br_lens)\n\n")
        
        cat("# put it all together\nmvBNdata ~ dnPhyloMultivariateBrownianREML(phylogeny, rep(1.0, numBranches), rateMatrix, nSites=ntraits, method=\"transform\")\nmvBNdata.clamp(means)\n\n")
        
        cat("# specify monitors\nsource(\"simulationsMonitors.Rev\")\n\n")
        
        cat("# define model\nmymodel = model(phylogeny)\n\n")
        
        cat("# run MCMC\nsource(\"simulationsMCMC.Rev\")\n\n")
        
        cat("# analyze mcmc output\nsource(\"simulationsAnalysis.Rev\")\n\n")
        
        cat("q()")
        
        sink()
        
        # add master script name to list
        
        scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
        sink("job.txt", append = T)
        cat(scriptPath, "\n", sep = "")
        sink()
        
      }
      
    }
    
  }
  
}
