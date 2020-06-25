setwd(dir = "/Volumes/1TB/Harvati/")

library(Matrix)
library(matrixcalc)
library(adephylo)
library(phytools)
library(mvMORPH)
library(mvtnorm)
library(mnormt)
library(ggplot2)
library(phangorn)
library(distory)
library(psych)
library(abind)
library(sirt)
library(irtoys)
library(miscTools)
library(bayesm)
library(phyloTop)
library(parallel)
library(doParallel)
library(foreach)
library(rwty)
library(latticeExtra)

convNum2Str <- function(nums, key){
  sapply(1:length(nums), function(x) key[nums[x]])
}

prop.part.df <- function(trees, cutoff = 0.01, bs = T){
  if(class(trees) == "multiPhylo"){
    if(trees[[1]]$Nnode == (length(trees[[1]]$tip.label) - 1)){
      trees <- unroot(trees) #unroot rooted trees
    }
    tipLabs <- trees[[1]]$tip.label
    numTrees <- length(trees)
    if(bs) {
      out <- as.prop.part(bitsplits(trees))
    } else {
      out <- prop.part(trees)
    }
    outList <- as.list.data.frame(out)
    pps <- attributes(outList)$number/numTrees
    props <- data.frame(pps, as.matrix(outList)); colnames(props) <- c("postProbs", "clade")
    props <- props[order(-pps),]
    props <- props[props[,1] > cutoff,]
    rownames(props) <- 1:nrow(props)
    props$cladeNames <- sapply(1:length(props[,1]), function(x) sort(convNum2Str(props$clade[[x]], attributes(outList)$labels)))
    props <- props[,c(1,3)]
    if(!bs) {
      props <- props[-1,]
    }
    allClades <- c(props$cladeNames, lapply(1:length(props$postProbs), function(x) sort(setdiff(tipLabs, props$cladeNames[[x]]))))
    if(any(duplicated(allClades))){
      # print("duplicate clades found")
      dupClades <- allClades[duplicated(allClades)]
      dupCladesComp <- lapply(1:length(dupClades), function(x) sort(setdiff(tipLabs, dupClades[[x]])))
      matchedDupes <- cbind(1:length(dupClades), sapply(1:length(dupClades), function(x) which(sapply(1:length(dupCladesComp), function(y) setequal(dupClades[[x]], dupCladesComp[[y]])))))
      dupes <- matrix(data = 0, nrow = length(matchedDupes[,1])/2, ncol = 2)
      for(i in 1:length(matchedDupes[,1])){
        pair <- sort(matchedDupes[i,])
        if(!any(sapply(1:length(dupes[,1]), function(x) pair == dupes[x,]))){
          dupes[min(which(apply(dupes, 1, sum) == 0)),] <- pair
        }
      }
      for(i in 1:length(dupes[,1])){
        clade1 <- dupClades[dupes[i,1]][[1]]
        clade2 <- dupClades[dupes[i,2]][[1]]
        propsInd1 <- which(sapply(1:length(props[,1]), function(x) setequal(clade1, props$cladeNames[[x]])))
        propsInd2 <- which(sapply(1:length(props[,1]), function(x) setequal(clade2, props$cladeNames[[x]])))
        props[propsInd1,1] <- props[propsInd1,1] + props[propsInd2,1]
        props <- props[-propsInd2,]
      }
      props <- props[order(props[,1], decreasing = T),]
    }
    props
  } else if (class(trees) == "phylo"){
    if(trees$Nnode == (length(trees$tip.label) - 1)){
      trees <- unroot(trees) #unroot rooted trees
    }
    tipLabs <- trees$tip.label
    numTrees <- 1
    if(bs) {
      out <- as.prop.part(bitsplits(trees))
    } else {
      out <- prop.part(trees)
    }
    outList <- as.list.data.frame(out)
    pps <- attributes(outList)$number/numTrees
    props <- data.frame(pps, as.matrix(outList)); colnames(props) <- c("postProbs", "clade")
    props <- props[order(-pps),]
    props <- props[props[,1] > cutoff,]
    rownames(props) <- 1:nrow(props)
    props$cladeNames <- sapply(1:length(props[,1]), function(x) sort(convNum2Str(props$clade[[x]], attributes(outList)$labels)))
    props <- props[,c(1,3)]
    if(!bs) {
      props <- props[-1,]
    }
    props
  }
}

model = "Mk"
windowCodes <- seq(from = 1, to = 100, by = 1)

#change filenames if needed
addNoisyName <- F
if(addNoisyName){
  files <- list.files("output")
  noisyness <- strsplit(getwd(), "_")[[1]][startsWith(strsplit(getwd(), "_")[[1]], "nois")]
  for(i in 1:length(files)){
    if(i%%1000==0){print(i)}
    file.rename(from = paste0("output/", files[i]), to = paste0("output/", gsub("Strength", paste0("Strength_", noisyness), files[i])))
  }
}

cl <- makeCluster(1, outfile="")
registerDoParallel(cl)
getDoParWorkers()

if(!dir.exists(paste0("bipartitionPosteriors_", model))){dir.create(paste0("bipartitionPosteriors_", model))}

foreach(m=1:100, .packages = c("ape", "coda", "rwty", "phangorn", "lattice", "latticeExtra")) %dopar% {
  
  traitNumIter <- (15*(2^(0:4)))*3+1
  traitNumIter <- rev(traitNumIter)
  transforms <- c("raw", "PCs")
  
  for (j in 1:length(traitNumIter)) {
    
    for(i in m:m){
      
      for(transform in transforms){
      # for(i in 1:100){
      
      cat(c(j,i, "... \n"))
      
      ntraits <- traitNumIter[j]
      replicateNum <- i
      
      fileName <- paste0("simulations_", ntraits, "traits_", "replicate_", replicateNum)
      trueTree <- read.nexus(file = paste0("trees/", fileName, "_tree.nex"))
      
      if(model == "Mk"){fileName <- paste0(fileName, "_MVC_", transform)}
      
      
      cat("reading trees... ")
      
      # analysis
      
      #find the name of the output directory
      if(model == "mvBM"){
        outputDirName <- "output/"
      } else if(model == "Mk"){
        outputDirName <- "output_discrete/"
      } else if(model == "orderedCTMC"){
        outputDirName <- "output_discrete_jenks/"
      }
      trees1 <- read.tree(file = paste0(outputDirName, fileName, "_trees_run_1.trees"))
      trees2 <- read.tree(file = paste0(outputDirName, fileName, "_trees_run_2.trees"))
      trees <- c(trees1, trees2)
      
      cat("processing trees...\n")
      biparts <- prop.part.df(trees, 0.0000005)
      biparts$trueState <- 0
      trueClades <- prop.part.df(unroot(trueTree))[,2]
      
      for(clade in 1:length(trueClades)){
        cladeName <- trueClades[[clade]]
        isThere <- sapply(1:length(biparts[,1]), function(x) identical(cladeName, biparts[x,2][[1]]))
        altCladeName <- sort(setdiff(trueTree$tip.label, cladeName))
        orIsThere <- sapply(1:length(biparts[,1]), function(x) identical(altCladeName, biparts[x,2][[1]]))
        if(any(isThere)){
          biparts$trueState[isThere] <- 1
        } else if (any(orIsThere)) {
          biparts$trueState[orIsThere] <- 1
        } else {
          biparts <- rbind(biparts, c(0, paste(cladeName, collapse = " "), 1))
        }
      }
      save(x = biparts, file = paste0(paste0("bipartitionPosteriors_", model), "/", fileName, "_bipartitionProbs"))
    }
  }
}
}

setwd(dir = "/Volumes/2TB/500/full/mvBM_sims_PB_noiseless_500/")
traitNumIter <- c(2,4,8,16,32,64,128)
# traitNumIter <- c(128)
degreeCorrelation <- c("2", "4", "6", "8", "9")
# degreeCorrelation <- c("2")

for (k in 1:length(degreeCorrelation)) { #degree correlation
  
  bipartsReduced<- data.frame()
  
  for (j in 1:length(traitNumIter)) { #num traits
    
    for(i in 1:500){ #replicate
      
      cat(c(k,j,i, "... \n"))
      
      ntraits <- traitNumIter[j]
      rateMatrixCorrelation <- degreeCorrelation[k]
      replicateNum <- i
      
      fileName <- paste0("simulations_", ntraits, "traits_perfectRateMatrix_", rateMatrixCorrelation, "correlationStrength_noiseless_", "replicate_", replicateNum)
      load(file = paste0("bipartitionPosteriors/", fileName, "_bipartitionProbs"))
      # nums[i] <- (mean(biparts[biparts[,1] > 0.9 & biparts[,1] < 1,3]))
      
      #get internality measure
      internality <- sapply(1:length(biparts[,1]), function(x) length(biparts[x,2][[1]]))
      internality[internality > 15] <- 30 - internality[internality > 15]
      biparts$internality <- internality
      
      bipartsReduced <- rbind(bipartsReduced,biparts[biparts[,1] > .005,][-1,])
      
      
    }
  }
  # save(bipartsReduced, file = paste0("bipartitionPosteriors/bipartsReduced_PB_128Traits_noError0.", rateMatrixCorrelation))
  save(bipartsReduced, file = paste0("bipartitionPosteriors/bipartsReduced_PB_allTraits_noError0.", rateMatrixCorrelation))
}

##################################################################################################################

# load(file = paste0("bipartitionPosteriors/bipartsReduced_PB_128Traits_noError0.", correlationStrength))
# load(file = paste0("bipartitionPosteriors/bipartsReduced_PB_64Traits_noError0.", correlationStrength))

#making the manuscript figure

png(filename = "~/mvbm_manuscript/figures/calibrationCurvePBT.png", width = 1200, height = 600)
par(mfrow = c(1,2), mar = c(3,3,3,0.75))
for(i in 1:2){
  if(i == 1){noisyness <- "noiseless"; error <- "noError"}
  if(i == 2){noisyness <- "noisy"; error <- "error"}
  setwd(dir = paste0("/Volumes/2TB/500/full/mvBM_sims_PB_", noisyness, "_500"))
  degreeCorrelation <- c("2", "4", "6", "8", "9")
  correlationStrength <- "2"
  if(!dir.exists("bipartitionPosteriors/exp_obs")){dir.create("bipartitionPosteriors/exp_obs")}
  if(file.exists(paste0("bipartitionPosteriors/exp_obs/", correlationStrength))){
    load(paste0("bipartitionPosteriors/exp_obs/", correlationStrength))} else {
      load(file = paste0("bipartitionPosteriors/bipartsReduced_PB_allTraits_", error, "0.", correlationStrength))
      intervalSize <- 0.1
      breaks <- seq(0,1,intervalSize)
      # for all the rate matrix probs
      expected <- rep(0, length(breaks)-1)
      observed <- rep(0, length(breaks)-1)
      for(i in 1:length(breaks)-1){
        inInterval <- bipartsReduced[bipartsReduced[,1] > breaks[i] & bipartsReduced[,1] < breaks[i+1],]
        expected[i] <- mean(as.numeric(inInterval[,1]))
        observed[i] <- mean(as.numeric(inInterval[,3]))
      }
      expobs <- list(expected, observed)
      save(expobs, file = paste0("bipartitionPosteriors/exp_obs/", correlationStrength))
    }
  expected <- expobs[[1]]; observed <- expobs[[2]]
  plot(expected, observed, ylab = "Observed Frequency", xlab = "Inferred Probability", xlim = c(0,1), ylim = c(0,1),
       # main = paste0("Calibration Curve for Binned (", intervalSize, ") Nodal Posterior Probabilities\nPure-Birth Trees, Meas. Error, Corr = 0.", correlationStrength, ", All Traits Analyses"))
       main = paste0("Calibration Curve for Binned (", intervalSize, ") Nodal Posterior Probabilities\nPure-Birth Trees, ", tools::toTitleCase(noisyness)))
  abline(a = 0, b = 1, col = 2, lwd = 2)
  lines(expected, observed)
  if(noisyness == "noiseless"){
    text(paste0("0.", correlationStrength), x = expected[1] - 0.02, y = expected[1]+0.03, cex = 0.95)
  } else {
    text(paste0("0.", correlationStrength), x = expected[length(expected)]+0.03, y = observed[length(expected)], cex = 0.9)
  }
  for(cs in 2:length(degreeCorrelation)){
    correlationStrength <- degreeCorrelation[cs]
    if(file.exists(paste0("bipartitionPosteriors/exp_obs/", correlationStrength))){
      load(paste0("bipartitionPosteriors/exp_obs/", correlationStrength))} else {
        load(file = paste0("bipartitionPosteriors/bipartsReduced_PB_allTraits_", error, "0.", correlationStrength))
        
        # for all the rate matrix probs
        expected <- rep(0, length(breaks)-1)
        observed <- rep(0, length(breaks)-1)
        for(i in 1:length(breaks)-1){
          inInterval <- bipartsReduced[bipartsReduced[,1] > breaks[i] & bipartsReduced[,1] < breaks[i+1],]
          expected[i] <- mean(as.numeric(inInterval[,1]))
          observed[i] <- mean(as.numeric(inInterval[,3]))
        }
        expobs <- list(expected, observed)
        save(expobs, file = paste0("bipartitionPosteriors/exp_obs/", correlationStrength))
      }
    expected <- expobs[[1]]; observed <- expobs[[2]]
    points(expected, observed)
    lines(expected, observed)
    if(correlationStrength == 9){correlationStrength <- 99}
    if(noisyness == "noiseless"){
      if(correlationStrength == 9){correlationStrength <- 99}
      text(paste0("0.", correlationStrength), x = expected[cs] - 0.02, y = expected[cs]+0.03, cex = 0.95)
    } else {
      text(paste0("0.", correlationStrength), x = expected[length(expected)]+0.035, y = observed[length(expected)], cex = 0.95)
    }
  }
  
  #getting the uvBM analyses in there too
  setwd(paste0("/Volumes/2TB/500/full/uvBM_sims_PB_", noisyness, "_500"))
  if(!dir.exists("bipartitionPosteriors/exp_obs")){dir.create("bipartitionPosteriors/exp_obs")}
  for(cs in 1:length(degreeCorrelation)){
    correlationStrength <- degreeCorrelation[cs]
    if(file.exists(paste0("bipartitionPosteriors/exp_obs/", correlationStrength))){
      load(paste0("bipartitionPosteriors/exp_obs/", correlationStrength))} else {
        load(file = paste0("bipartitionPosteriors/bipartsReduced_PB_allTraits_", error, "0.", correlationStrength))
        expected <- rep(0, length(breaks)-1)
        observed <- rep(0, length(breaks)-1)
        for(i in 1:length(breaks)-1){
          inInterval <- bipartsReduced[bipartsReduced[,1] > breaks[i] & bipartsReduced[,1] < breaks[i+1],]
          expected[i] <- mean(as.numeric(inInterval[,1]))
          observed[i] <- mean(as.numeric(inInterval[,3]))
        }
        expobs <- list(expected, observed)
        save(expobs, file = paste0("bipartitionPosteriors/exp_obs/", correlationStrength))
      }
    expected <- expobs[[1]]; observed <- expobs[[2]]
    lines(expected, observed, lty = 3)
    points(expected, observed, pch = 0)
    if(correlationStrength == 9){correlationStrength <- 99}
    text(paste0("0.", correlationStrength), x = expected[length(expected)]+.035, cex = 0.95,  y = observed[length(expected)], col = "gray50")
  }
  legend(0, 1, c("multivariate inference model","univariate inference model"), cex=1, pch=21:22, lty=c(1,3));
  if(i == 1){text(x = -0.02, y = 1.02, "A)", font = 3)}
  if(i == 2){text(x = -0.02, y = 1.02, "B)", font = 3)}
}
dev.off()
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################

# sliding window style
pointSpacing <- 0.01
intervalSize <- 0.05
breaks <- seq(0+intervalSize/2,1-intervalSize/2,pointSpacing)
expected <- rep(0, length(breaks))
observed <- rep(0, length(breaks))
for(i in 1:length(breaks)){
  inInterval <- bipartsReduced[bipartsReduced[,1] > breaks[i] - intervalSize/2 & bipartsReduced[,1] < breaks[i] + intervalSize/2,]
  expected[i] <- mean(as.numeric(inInterval[,1]))
  observed[i] <- mean(as.numeric(inInterval[,3]))
}
par(mfrow = c(1,1))
plot(expected, observed, ylab = "Observed Frequency", xlab = "Inferred Probability", xlim = c(0,1), ylim = c(0,1), #type = "l",
     main = paste0("Calibration Curve for Binned (", intervalSize, ") Nodal Posterior Probabilities\nPure-Birth Trees, Meas. Error, Corr = 0.", correlationStrength, ", All Traits Analyses"))
abline(a = 0, b = 1, col = 2, lwd = 2)




# looking at node-internality-specific calibration curves

par(mfrow = c(3,5))
intervalSize <- 0.05
breaks <- seq(0,1,intervalSize)
for(j in 2:15){
  expected <- rep(0, length(breaks)-1)
  observed <- rep(0, length(breaks)-1)
  for(i in 1:length(breaks)-1){
    inInterval <- bipartsReduced[bipartsReduced[,1] > breaks[i] & bipartsReduced[,1] < breaks[i+1] & bipartsReduced[,4] == j,]
    expected[i] <- mean(as.numeric(inInterval[,1]))
    observed[i] <- mean(as.numeric(inInterval[,3]))
  }
  plot(expected, observed, ylab = "Observed Frequency", xlab = "Inferred Probability",
       main = paste0("Bin Size: ", intervalSize, ", Perfect RM\nAll Traits, Node Internality: ", j))
  abline(a = 0, b = 1, col = 2, lwd = 2)
}

#test out intuition about averaging
prob1 <- .279
prob2 <- .15
n1 <- 1e6
n2 <- 4e7
(rbinom(n = 1, size = n1, prob = prob1) + rbinom(n = 1, size = n2, prob = prob2)) /(n1+n2)
(prob1*n1+prob2*n2)/(n1+n2)
#intuition success!


#are more internal nodes harder to idenfity?
setwd(dir = "/Users/nikolai/mvBM_sims_underhand")
traitNumIter <- c(2,4,8,16,32,64)
degreeCorrelation <- c("2", "4", "6", "8")

nodalDifficulty <- matrix(nrow = 14, ncol = 100, data = 0)
counter <- 0
for (k in 1:3) {
  
  for (j in 1:6) {
    
    for(i in 1:100){
      counter <- counter + 1
      
      cat(c(k,j,i, "... \n"))
      
      ntraits <- traitNumIter[j]
      rateMatrixCorrelation <- degreeCorrelation[k]
      replicateNum <- i
      
      fileName <- paste0("simulations_", ntraits, "traits_perfectRateMatrix_", rateMatrixCorrelation, "correlationStrength_", "replicate_", replicateNum)
      load(file = paste0("bipartitionPosteriors/", fileName, "_bipartitionProbs"))
      # bipartsAll <- rbind(bipartsAll, biparts)
      bipartsReduced <- biparts[biparts[,1] > .005,]
      internality <- sapply(1:length(bipartsReduced[,1]), function(x) length(bipartsReduced[x,2][[1]]))
      internality[internality > 15] <- 30 - internality[internality > 15]
      bipartsReduced$internality <- internality
      nodalDifficulty[,counter] <- sapply(2:15, function(x) mean(bipartsReduced[bipartsReduced$internality == x & bipartsReduced$trueState == 1, 1]))
    }
  }
}
par(mfrow = c(1,1))
plot(2:15,sapply(1:14, function(x) mean(nodalDifficulty[x,], na.rm = T)), main = "Average Difficulty Inferring Node of Given Internality (64 Traits)",
     xlab = "Node Internality", ylab = "Mean Posterior Probability of Nodes of Given Interality, when they exist")



#size of nodal credible set estimation
credSetTrees <- function(trees, prob = 0.95, returnMAP = FALSE){
  trees <- unroot(trees)
  alreadyCounted <- rep(FALSE, length(trees))
  # if(exists("uniqueTrees")){rm("uniqueTrees")}
  while(any(!alreadyCounted)){
    # print(paste0(sum(alreadyCounted)/length(trees)*100," %"))
    needToSee <- which(!alreadyCounted)
    tree <- trees[[needToSee[1]]]
    compareTo <- trees[needToSee][-1]
    if(length(compareTo) != 0) {
      dists <- RF.dist(tree, compareTo)
      if(any(dists == 0)) {sameTrees <- which(dists == 0)} else {sameTrees <- 0}
      alreadyCounted[needToSee[c(1,sameTrees)]] <- TRUE
    } else {
      sameTrees <- 0; alreadyCounted[needToSee[c(1,sameTrees)]] <- TRUE
    }
    tree$numAppearances <- sum(sameTrees != 0) + 1
    if(exists("uniqueTrees")){uniqueTrees <- c(uniqueTrees, tree)} else {uniqueTrees <- tree}
  }
  uniqueTrees <- uniqueTrees[order(sapply(1:length(uniqueTrees), function(x) uniqueTrees[[x]]$numAppearances), decreasing = T)]
  pps <- sapply(1:length(uniqueTrees), function(x) uniqueTrees[[x]]$numAppearances) / length(trees)
  numTrees <- min(which(cumsum(pps) > prob))
  return(numTrees)
  if(returnMAP){
    return(list(numTrees, length(trees), uniqueTrees[[1]]))
  } else {
    return(c(numTrees, length(trees)))
  }
}

setwd(dir = "/Volumes/2TB/mvBM_sims_PB_noiseless/")
traitNumIter <- c(2,4,8,16,32,64,128)
degreeCorrelation <- c("2", "4", "6", "8", "9")
noisyness <- "noiseless"
credibleSetProb <- 0.95
nodalCredSetProb <- credibleSetProb * 28
counter <- 0
credsetMat <- matrix(nrow = 100*length(traitNumIter)*length(degreeCorrelation), ncol = 5, data = 0)
colnames(credsetMat) <- c("rateMatrixSpec", "ntraits", "replicateNum", "credSetSizeNode", "credSetSizeTree")

for (k in 1:length(degreeCorrelation)) {
  
  for (j in 1:length(traitNumIter)) {
    
    for(i in 1:100){
      counter <- counter + 1
      
      cat(c(k,j,i, "... \n"))
      
      ntraits <- traitNumIter[j]
      rateMatrixCorrelation <- degreeCorrelation[k]
      replicateNum <- i
      
      fileName <- paste0("simulations_", ntraits, "traits_perfectRateMatrix_", rateMatrixCorrelation, "correlationStrength_", noisyness, "_replicate_", replicateNum)
      load(file = paste0("bipartitionPosteriors/", fileName, "_bipartitionProbs"))
      trees1 <- read.tree(file = paste0("output/", fileName, "_trees_run_1.trees"))
      trees2 <- read.tree(file = paste0("output/", fileName, "_trees_run_2.trees"))
      trees <- c(trees1, trees2)
      trees <- trees[round(seq(1, length(trees), length.out = 2000))]
      credSetSizeNode <- min(which(cumsum(as.numeric(biparts[,1])) > nodalCredSetProb))
      credSetSizeTree <- credSetTrees(trees, credibleSetProb)[1]
      credsetMat[counter,] <- c(k,ntraits, replicateNum, credSetSize)
    }
  }
}

degreeMiss <- 1
plot(traitNumIter, sapply(1:length(traitNumIter), function(x) mean(credsetMat[credsetMat[,2] == traitNumIter[x] & credsetMat[,1] == degreeMiss,4])), log = "y",
     xlab = "Number of Traits Used in Analysis", ylab = "Number of Bipartitions in Credible Set", type = "l", col = 1,
     main = paste0(credibleSetProb * 100, "% credible set size\n(on the scale of bipartitions, so really a ", nodalCredSetProb * 100, "% credible set)"))
degreeMiss <- 2
lines(traitNumIter, sapply(1:length(traitNumIter), function(x) mean(credsetMat[credsetMat[,2] == traitNumIter[x] & credsetMat[,1] == degreeMiss,4])),
      xlab = "Number of Traits Used in Analysis", ylab = "Number of Bipartitions in Credible Set", col = 2,
      main = paste0(credibleSetProb * 100, "% credible set size\n(on the scale of bipartitions, so really a ", nodalCredSetProb * 100, "% credible set)"))
degreeMiss <- 3
lines(traitNumIter, sapply(1:length(traitNumIter), function(x) mean(credsetMat[credsetMat[,2] == traitNumIter[x] & credsetMat[,1] == degreeMiss,4])),
      xlab = "Number of Traits Used in Analysis", ylab = "Number of Bipartitions in Credible Set", col = 3,
      main = paste0(credibleSetProb * 100, "% credible set size\n(on the scale of bipartitions, so really a ", nodalCredSetProb * 100, "% credible set)"))
legend(51, 55000,degreeCorrelation,col=1:3,lwd = 2)


## look at p-value/percentile distributions for tree length,
## likelihood, patristic distances, tree balance, and internal and external branch lengths

setwd(dir = "/Users/nikolai/mvBM_sims_underhand")

library(adephylo)
library(readtext)
traitNumIter <- c(2,4,8,16,32,64)
traitNumIter <- rev(traitNumIter)

degreeCorrelation <- c("2", "4", "6", "8")

BMpruneLL <- function(traits, sig, tree){
  LLs <- 0
  for(i in 1:(length(tree$tip.label) - 1)){
    #Find two sister tips
    #preterminal nodes
    ptn <- tree$edge[sapply(1:length(tree$tip.label), function(t) which(tree$edge[,2] == t)),1]
    uniq <- unique(ptn)
    parent <- uniq[which(sapply(1:length(uniq), function (x) sum(uniq[x] == ptn)) > 1)][1]
    sisters <- tree$edge[tree$edge[,1] == parent, 2]
    sisters <- intersect(sisters, 1:(length(tree$tip.label)))
    
    #calculate contrast between two sister nodes
    contrast <- traits[tree$tip.label[sisters[1]],] - traits[tree$tip.label[sisters[2]],]
    
    #calculate BL separating two sister nodes
    BLength <- tree$edge.length[tree$edge[,2] == sisters[1]] + tree$edge.length[tree$edge[,2] == sisters[2]]
    
    #calculate multivariate normal density of contrast along Branch Length
    LLs[i] <- dmvnorm(x = contrast, sigma = BLength*sig, log = T)
    
    if(length(tree$tip.label) > 2){  
      #Fix Tree and Trait Matrix
      tipName = paste(tree$tip.label[sisters], collapse = "+")
      #Fix Trait Matrix
      nodeValues <- (traits[tree$tip.label[sisters[1]],]*tree$edge.length[tree$edge[,2] == sisters[2]] + traits[tree$tip.label[sisters[2]],]*tree$edge.length[tree$edge[,2] == sisters[1]])/BLength
      traits <- rbind(traits, nodeValues)
      rownames(traits)[length(traits[,1])] <- tipName
      traits <- traits[-c(which(rownames(traits) == tree$tip.label[sisters[1]]), which(rownames(traits) == tree$tip.label[sisters[2]])),]
      
      #Fix Tree
      newTip <- list(edge=matrix(c(2,1),1,2),
                     tip.label=tipName,
                     edge.length=(prod(tree$edge.length[tree$edge[,2] == sisters[1]], tree$edge.length[tree$edge[,2] == sisters[2]])/BLength) - tree$edge.length[tree$edge[,2] == sisters[1]],
                     Nnode=1)
      class(newTip)<-"phylo"
      tree <- bind.tree(x = tree, y = newTip, where = sisters[1])
      tree <- drop.tip(tree, sisters[2])
    }
  }
  return(sum(LLs))
}

counter <- 0
percLLs <- matrix(nrow = 2100, ncol = 4, data = 0)
percTLs <- matrix(nrow = 2100, ncol = 4, data = 0)
percColless <- matrix(nrow = 2100, ncol = 4, data = 0)
percPats <- matrix(nrow = 2100, ncol = 3 + (30^2-30)/2, data = 0)
percTermBLs <- matrix(nrow = 2100, ncol = 3 + 30, data = 0)
#
# write.table(x = percLLs, file = "analyses/percLLs.txt")
# write.table(x = percTLs, file = "analyses/percTLs.txt")
# write.table(x = percColless, file = "analyses/percColless.txt")
# write.table(x = percPats, file = "analyses/percPats.txt")
# write.table(x = percTermBLs, file = "analyses/percTermBLs.txt")

for (k in 1:length(degreeCorrelation)) {
  # for (k in 1:1) {
  
  for (j in 1:6) {
    
    for(i in 1:100){
      
      counter <- counter + 1
      
      ntraits <- traitNumIter[j]
      rateMatrixCorrelation <- degreeCorrelation[k]
      replicateNum <- i
      cat(c(ntraits, rateMatrixCorrelation, replicateNum, "... \n"))
      
      percLLs[counter, 1:3] <- percTLs[counter, 1:3] <- percPats[counter, 1:3] <- percTermBLs[counter, 1:3] <- percColless[counter, 1:3] <- c(ntraits, rateMatrixCorrelation, replicateNum)
      
      fileName <- paste0("simulations_", ntraits, "traits_perfectRateMatrix_", rateMatrixCorrelation, "correlationStrength_", "replicate_", replicateNum)
      trueTree <- read.nexus(file = paste0("/Users/nikolai/mvBM_sims_underhand/trees/", fileName, "_tree.nex"))
      
      cat("reading trees... ")
      
      # analysis
      trees1 <- read.tree(file = paste0("/Users/nikolai/mvBM_sims_underhand/output/", fileName, "_trees_run_1.trees"))
      trees2 <- read.tree(file = paste0("/Users/nikolai/mvBM_sims_underhand/output/", fileName, "_trees_run_2.trees"))
      trees <- c(trees1, trees2)
      
      tipNames <- trees[[1]]$tip.label
      distMats <- lapply(1:length(trees), function (x) as.matrix(distTips(trees[[x]]))[sort(tipNames),sort(tipNames)])
      trueDist <- as.matrix(distTips(trueTree))[sort(tipNames),sort(tipNames)]
      
      patDistPercentiles <- trueDist
      patDistPercentiles[1:length(patDistPercentiles)] <- 0
      for(a in 2:length(patDistPercentiles[1,])){
        for(b in 1:(a-1)){
          patDistPercentiles[a,b] <- sum(trueDist[a,b] > sapply(1:length(distMats), function(x) distMats[[x]][a,b]))/length(distMats)
        }
      }
      patDistPercentilesVec <- patDistPercentiles[lower.tri(patDistPercentiles)]
      percPats[counter, 4:length(percPats[1,])] <- patDistPercentilesVec
      
      # plot(patDistPercentiles[patDistPercentiles != 0], trueDist[patDistPercentiles != 0])
      # hist(patDistPercentiles[patDistPercentiles != 0])
      
      #get traits
      traits <- strsplit(readtext(file = paste0("/Users/nikolai/mvBM_sims_underhand/data/", fileName, "_traits.nex"))$text, split = "\n")[[1]][8:37]
      traits <- t(sapply(1:30, function (x) strsplit(traits[x], split = " ")[[1]]))
      names <- traits[,1]
      traits <- traits [,-1]
      class(traits) <- "numeric"
      rownames(traits) <- names
      
      #get rate matrix
      wgCM <- as.matrix(read.table(file = paste0("/Users/nikolai/mvBM_sims_underhand/rates/", fileName, "_rates.tsv")))
      
      #get true tl and ll
      trueLL <- BMpruneLL(tree = trueTree, traits = traits, sig = wgCM)
      trueTL <- sum(trueTree$edge.length)
      
      #get sampled TLs and LLs
      log1 <- read.table(file = paste0("/Users/nikolai/mvBM_sims_underhand/output/", fileName, "_params_run_1.log"), header = T)
      log1 <- log1[seq(1, dim(log1)[1], by = 2),] #thin further
      log1 <- log1[, c("Likelihood", "TL")]
      log2 <- read.table(file = paste0("/Users/nikolai/mvBM_sims_underhand/output/", fileName, "_params_run_2.log"), header = T)
      log2 <- log2[seq(1, dim(log2)[1], by = 2),] #thin further
      log2 <- log2[, c("Likelihood", "TL")]
      log <- rbind(log1, log2)  
      
      percentileLL <- sum(trueLL > log$Likelihood)/length(log$Likelihood)
      percLLs[counter, 4] <- percentileLL
      percentileTL <- sum(trueTL > log$TL)/length(log$TL)
      percTLs[counter, 4] <- percentileTL
      
      #get Colless' Balance/Clade Shape Metric
      trueColless <- colless.phylo(trueTree, normalise = F)
      mcmcColless <- sapply(1:length(trees), function(x) colless.phylo(trees[[x]], normalise = F))
      percentileColless <- sum(trueColless > mcmcColless)/length(mcmcColless)
      percColless[counter, 4] <- percentileColless
      
      
      #get terminal BLs
      trueTerminalBLs <- cbind(trueTree$tip.label, trueTree$edge.length[sapply(1:length(trueTree$tip.label), function(x) which(trueTree$edge[,2] == x))])
      trueTerminalBLs <- trueTerminalBLs[order(trueTerminalBLs[,1]),]
      abcTrueTerminalBLs <- as.numeric(trueTerminalBLs[,2])
      
      abcMCMCTerminalBLs <- matrix(ncol = length(trees[[a]]$tip.label), nrow = length(trees), data = 0)
      for(a in 1:length(trees)){
        mcmcTerminalBLs <- cbind(trees[[a]]$tip.label, trees[[a]]$edge.length[sapply(1:length(trees[[a]]$tip.label), function(x) which(trees[[a]]$edge[,2] == x))])
        mcmcTerminalBLs <- mcmcTerminalBLs[order(mcmcTerminalBLs[,1]),]
        abcMCMCTerminalBLs[a,] <- as.numeric(mcmcTerminalBLs[,2])
      }
      terminalBLPercentiles <- sapply(1:length(abcTrueTerminalBLs), function(x) sum(abcTrueTerminalBLs[x] > abcMCMCTerminalBLs[,x])/length(abcMCMCTerminalBLs[,x]))
      percTermBLs[counter, 4:length(percTermBLs[1,])] <- terminalBLPercentiles
      
    }
  }
}

percLLs <- read.table(file = "analyses/percLLs.txt")
percTLs <- read.table(file = "analyses/percTLs.txt")
percColless <- read.table(file = "analyses/percColless.txt")
percPats <- read.table(file = "analyses/percPats.txt")
percTermBLs <- read.table(file = "analyses/percTermBLs.txt")

#get info on prior distribution
percLLsPrior <- matrix(nrow = 100, ncol = 4, data = 0)
percTLsPrior <- matrix(nrow = 100, ncol = 4, data = 0)
percCollessPrior <- matrix(nrow = 100, ncol = 4, data = 0)
percPatsPrior <- matrix(nrow = 100, ncol = 3 + (30^2-30)/2, data = 0)
percTermBLsPrior <- matrix(nrow = 100, ncol = 3 + 30, data = 0)

counter <- 0
for (k in 1:1) {
  
  for (j in 1:1) {
    
    for(i in 1:100){
      
      counter <- counter + 1
      
      ntraits <- traitNumIter[j]
      rateMatrixCorrelation <- degreeCorrelation[k]
      replicateNum <- i
      cat(c(ntraits, rateMatrixCorrelation, replicateNum, "... \n"))
      
      percLLsPrior[counter, 1:3] <- percTLsPrior[counter, 1:3] <- percPatsPrior[counter, 1:3] <- percTermBLsPrior[counter, 1:3] <- percCollessPrior[counter, 1:3] <- c(ntraits, rateMatrixCorrelation, replicateNum)
      
      fileName <- paste0("simulations_", ntraits, "traits_perfectRateMatrix_", rateMatrixCorrelation, "correlationStrength_", "replicate_", replicateNum)
      trueTree <- read.nexus(file = paste0("/Users/nikolai/mvBM_sims_underhand/trees/", fileName, "_tree.nex"))
      
      # analysis
      trees <- rmtree(N = 2000, n = 30)
      for(l in 1:2000){
        trees[[l]]$edge.length <- 10^runif(n = 1, min = -3, max = 2) * rdirichlet(n = 1, alpha = rep(1,58))
        trees[[l]]$tip.label <- sample(trueTree$tip.label, size = 30, replace = F)
      }
      
      tipNames <- trees[[1]]$tip.label
      distMats <- lapply(1:length(trees), function (x) as.matrix(distTips(trees[[x]]))[sort(tipNames),sort(tipNames)])
      trueDist <- as.matrix(distTips(trueTree))[sort(tipNames),sort(tipNames)]
      
      patDistPercentiles <- trueDist
      patDistPercentiles[1:length(patDistPercentiles)] <- 0
      for(a in 2:length(patDistPercentiles[1,])){
        for(b in 1:(a-1)){
          patDistPercentiles[a,b] <- sum(trueDist[a,b] > sapply(1:length(distMats), function(x) distMats[[x]][a,b]))/length(distMats)
        }
      }
      patDistPercentilesVec <- patDistPercentiles[lower.tri(patDistPercentiles)]
      percPatsPrior[counter, 4:length(percPatsPrior[1,])] <- patDistPercentilesVec
      
      #get traits
      traits <- strsplit(readtext(file = paste0("/Users/nikolai/mvBM_sims_underhand/data/", fileName, "_traits.nex"))$text, split = "\n")[[1]][8:37]
      traits <- t(sapply(1:30, function (x) strsplit(traits[x], split = " ")[[1]]))
      names <- traits[,1]
      traits <- traits [,-1]
      class(traits) <- "numeric"
      rownames(traits) <- names
      
      #get rate matrix
      wgCM <- as.matrix(read.table(file = paste0("/Users/nikolai/mvBM_sims_underhand/rates/", fileName, "_rates.tsv")))
      
      #get true tl and ll
      trueLL <- BMpruneLL(tree = trueTree, traits = traits, sig = wgCM)
      trueTL <- sum(trueTree$edge.length)
      
      #get sampled TLs and LLs
      log <- data.frame(row.names = 1:length(trees))
      
      log$TL <- sapply(1:length(trees), function(x) sum(trees[[x]]$edge.length))
      log$Likelihood <- sapply(1:length(trees[1:100]), function(x) BMpruneLL(tree = trees[[x]], traits = traits, sig = wgCM))
      
      percentileLL <- sum(trueLL > log$Likelihood)/length(log$Likelihood)
      percLLsPrior[counter, 4] <- percentileLL
      percentileTL <- sum(trueTL > log$TL)/length(log$TL)
      percTLsPrior[counter, 4] <- percentileTL
      
      #get Colless' Balance/Clade Shape Metric
      trueColless <- colless.phylo(trueTree, normalise = F)
      mcmcColless <- sapply(1:length(trees), function(x) colless.phylo(trees[[x]], normalise = F))
      percentileColless <- sum(trueColless > mcmcColless)/length(mcmcColless)
      percCollessPrior[counter, 4] <- percentileColless
      
      
      #get terminal BLs
      trueTerminalBLs <- cbind(trueTree$tip.label, trueTree$edge.length[sapply(1:length(trueTree$tip.label), function(x) which(trueTree$edge[,2] == x))])
      trueTerminalBLs <- trueTerminalBLs[order(trueTerminalBLs[,1]),]
      abcTrueTerminalBLs <- as.numeric(trueTerminalBLs[,2])
      
      abcMCMCTerminalBLs <- matrix(ncol = length(trees[[a]]$tip.label), nrow = length(trees), data = 0)
      for(a in 1:length(trees)){
        mcmcTerminalBLs <- cbind(trees[[a]]$tip.label, trees[[a]]$edge.length[sapply(1:length(trees[[a]]$tip.label), function(x) which(trees[[a]]$edge[,2] == x))])
        mcmcTerminalBLs <- mcmcTerminalBLs[order(mcmcTerminalBLs[,1]),]
        abcMCMCTerminalBLs[a,] <- as.numeric(mcmcTerminalBLs[,2])
      }
      terminalBLPercentiles <- sapply(1:length(abcTrueTerminalBLs), function(x) sum(abcTrueTerminalBLs[x] > abcMCMCTerminalBLs[,x])/length(abcMCMCTerminalBLs[,x]))
      percTermBLsPrior[counter, 4:length(percTermBLsPrior[1,])] <- terminalBLPercentiles
      
    }
  }
}

#whoops let's redo the TL
percTLsPrior <- matrix(nrow = 100, ncol = 4, data = 0)
counter <- 0
for (k in 1:1) {
  
  for (j in 1:1) {
    
    for(i in 1:100){
      
      counter <- counter + 1
      
      ntraits <- traitNumIter[j]
      rateMatrixCorrelation <- degreeCorrelation[k]
      replicateNum <- i
      cat(c(ntraits, rateMatrixCorrelation, replicateNum, "... \n"))
      
      percTLsPrior[counter, 1:3] <- c(ntraits, rateMatrixCorrelation, replicateNum)
      
      fileName <- paste0("simulations_", ntraits, "traits_perfectRateMatrix_", rateMatrixCorrelation, "correlationStrength_", "replicate_", replicateNum)
      trueTree <- read.nexus(file = paste0("/Users/nikolai/mvBM_sims_underhand/trees/", fileName, "_tree.nex"))
      
      percentileTL <- sum(trueTL > 10^runif(n = 5334, min = -3, max = 2))/5334
      percTLsPrior[counter, 4] <- percentileTL
    }
  }
}

# making some graphs
percTLsPrior
percCollessPrior
percPatsPrior
percTermBLsPrior

traitNumIter <- c(2,4,8,16,32,64)

png(filename = "plots/PatristicDistances.png", width = 3000, height = 5000, pointsize = 60)
par(mfrow = c(8,3))
plot.new()
hist(as.numeric(percPatsPrior[,4:435]), xlim = c(0,1), breaks = 30, main = "Patristic Distance Distribution\nPrior", freq = F, xlab = "percentile", ylab = "")
plot.new()
for (j in 1:6) {
  
  for (k in 1:length(degreeCorrelation)) {
    
    ntraits <- traitNumIter[j]
    rateMatrixCorrelation <- degreeCorrelation[k]
    hist(as.numeric(percPats[percPats[,1] == ntraits & percPats[,2] == rateMatrixCorrelation, 4:435]),
         xlim = c(0,1), breaks = 100, main = paste0(rateMatrixCorrelation, " RMS, ", ntraits, " traits"), freq = F, xlab = "", ylab = "")
    
  }
  
}
dev.off()


png(filename = "plots/TreeLengths.png", width = 3000, height = 5000, pointsize = 60)
par(mfrow = c(8,3))
plot.new()
hist(as.numeric(percTLsPrior[,4]), xlim = c(0,1), breaks = 5, main = "Tree Length Distribution\nPrior", freq = F, xlab = "percentile", ylab = "")
plot.new()
for (j in 1:6) {
  
  for (k in 1:length(degreeCorrelation)) {
    
    ntraits <- traitNumIter[j]
    rateMatrixCorrelation <- degreeCorrelation[k]
    hist(as.numeric(percTLs[percTLs[,1] == ntraits & percTLs[,2] == rateMatrixCorrelation, 4]),
         xlim = c(0,1), breaks = 5, main = paste0(rateMatrixCorrelation, " RMS, ", ntraits, " traits"), freq = F, xlab = "", ylab = "")
    
  }
  
}
dev.off()

png(filename = "plots/Colless.png", width = 3000, height = 5000, pointsize = 60)
par(mfrow = c(8,3))
plot.new()
hist(as.numeric(percCollessPrior[,4]), xlim = c(0,1), breaks = 10, main = "Colless Distribution\nPrior", freq = F, xlab = "percentile", ylab = "")
plot.new()
for (j in 1:6) {
  
  for (k in 1:length(degreeCorrelation)) {
    
    ntraits <- traitNumIter[j]
    rateMatrixCorrelation <- degreeCorrelation[k]
    hist(as.numeric(percColless[percColless[,1] == ntraits & percColless[,2] == rateMatrixCorrelation, 4]),
         xlim = c(0,1), breaks = 10, main = paste0(rateMatrixCorrelation, " RMS, ", ntraits, " traits"), freq = F, xlab = "", ylab = "")
    
  }
  
}
dev.off()

png(filename = "plots/TerminalBLs.png", width = 3000, height = 5000, pointsize = 60)
par(mfrow = c(8,3))
plot.new()
hist(as.numeric(percTermBLsPrior[,4:33]), xlim = c(0,1), breaks = 30, main = "Terminal BLs Distribution\nPrior", freq = F, xlab = "percentile", ylab = "")
plot.new()
for (j in 1:6) {
  
  for (k in 1:length(degreeCorrelation)) {
    
    ntraits <- traitNumIter[j]
    rateMatrixCorrelation <- degreeCorrelation[k]
    hist(as.numeric(percTermBLs[percTermBLs[,1] == ntraits & percTermBLs[,2] == rateMatrixCorrelation, 4:33]),
         xlim = c(0,1), breaks = 30, main = paste0(rateMatrixCorrelation, " RMS, ", ntraits, " traits"), freq = F, xlab = "", ylab = "")
    
  }
  
}
dev.off()