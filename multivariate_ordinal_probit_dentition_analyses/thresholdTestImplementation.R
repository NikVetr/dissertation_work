#Felsenstein Pruning Algorithm
library(Matrix)
library(matrixcalc)
library(adephylo)
library(phytools)
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
library(mvProbit)
library( miscTools )


# plot(tree); axisPhylo(); tiplabels(); nodelabels(); edgelabels()


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

#Test out Speed in Single Trial
#Simulate Traits
tree1 <- (pbtree(n=30))
tree2 <- (pbtree(n=30))
tree <- consensus(tree1, tree2, p=.6)
tree <- tree1
plot(tree)


d.orig <- read.csv("~/data/Howell.csv")
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
cov <- covFULL <- linMeasCovHuman <- covHuman[linMeasNames, linMeasNames]


tree <- pbtree(n=10)
plot(tree)
#plot(tree); axisPhylo(); tiplabels(); nodelabels(); edgelabels()
nTraits <- 15
sig <- cov2cor(cov)[1:nTraits, 1:nTraits]
meanLiabs <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = nTraits, sigma=sig, mu= rep(0,nTraits)))
populations <- tree$tip.label
numTips <- length(populations)
nObs <- 50

# discrTraits <- rmvnorm(n = nObs, mean = meanLiabs[1,], sigma = sig)
# discrTraits <- discrTraits > 0 
# discrTraits <- matrix(as.numeric(discrTraits), nrow = nrow(discrTraits))
# traitCounts <- apply(X = discrTraits, MARGIN = 2, FUN = sum)
# discrTraitsAll <- list()
# for(i in 1:numTips){
#   discrTraits <- rmvnorm(n = nObs, mean = meanLiabs[i,], sigma = sig)
#   discrTraits <- discrTraits > 0 
#   discrTraits <- matrix(as.numeric(discrTraits), nrow = nrow(discrTraits))
#   discrTraitsAll[[i]] <- discrTraits
# }



###############################

tree <- pbtree(n=5)
plot(tree)
#plot(tree); axisPhylo(); tiplabels(); nodelabels(); edgelabels()
nTraits <- 5
traitsToUse <- sample(linMeasNames, nTraits)
sig <- cov2cor(cov)[traitsToUse, traitsToUse]
meanLiabs <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = nTraits, sigma=sig, mu= rep(0,nTraits)))
populations <- tree$tip.label
numTips <- length(populations)
nObs <- 20

#construct design matrix
xMat <- matrix(data = 0, ncol = numTips, nrow = nObs * numTips)
for (i in 1:numTips){xMat[((nObs*(i-1)+1):(nObs*i)),i] <- rep(1, nObs)}
colnames(xMat) <- paste0("x", 1:numTips)
xMat <- cbind(0, xMat)

# model coefficients
beta <- meanLiabs
beta <- rbind(0, beta)

# covariance matrix of error terms
sigma <- sig

# generate dependent variables
yMatLin <- xMat %*% beta
yMat <- ( yMatLin + rmvnorm( nObs*numTips, sigma = sigma ) ) > 0
colnames( yMat ) <- paste0( "y", 1:nTraits)

# log likelihood values
myData <- as.data.frame( cbind( xMat, yMat ) )
formulaMvProbit <- as.formula(paste0("cbind( ", paste0("y", 1:(nTraits-1), sep = ",", collapse = ""), paste0("y", nTraits), 
") ~ ", paste0("x", 1:(numTips-1), sep = "+", collapse = ""), paste0("x", numTips)))
formulaMvProbit
 
logLikValTrue <- mvProbitLogLik( formula = formulaMvProbit, coef = c( beta ), sigma = sigma, data = myData, algorithm = "GHK" )

#simple met hastings sampler
sigma <- diag(dim(sigma)[1])
beta <- matrix(0, nrow = dim(beta)[1], ncol = dim(beta)[2])
logLikVal <- mvProbitLogLik( formula = formulaMvProbit, coef = c( beta ), sigma = sigma, data = myData, algorithm = "GHK" )
nIter <- 20000
sigmaSamples <- array(dim = c(dim(sigma), nIter))
sigmaSamples[,,1] <- sigma
betaSamples <- array(dim = c(dim(beta), nIter))
betaSamples[,,1] <- beta
for (i in 2:nIter){
  
  
  pSI <- sample(1:dim(sigma)[1], 2, replace = F)
  propSigma <- sigma
  newCorr <- propSigma[pSI[1], pSI[2]] + runif(1, -0.1, 0.1)
  if (newCorr > 1){newCorr <- 2 - newCorr}
  if (newCorr < -1){newCorr <- -2 - newCorr}
  propSigma[pSI[1], pSI[2]] <- propSigma[pSI[2], pSI[1]] <- newCorr
  while(!all(eigen(propSigma)$values >= 0)){
    print("fail")
    propSigma <- sigma
    newCorr <- propSigma[pSI[1], pSI[2]] + runif(1, -0.1, 0.1)
    if (newCorr > 1){newCorr <- 2 - newCorr}
    if (newCorr < -1){newCorr <- -2 - newCorr}
    propSigma[pSI[1], pSI[2]] <- propSigma[pSI[2], pSI[1]] <- newCorr
  }

  
  dispBeta <- rbind(0, matrix(0, nrow = (dim(beta)[1]-1), ncol = dim(beta)[2]))
  dispBeta[sample(2:(numTips+1), 1),] <- rnorm(n = nTraits, mean = 0, sd = .1)
  propBeta <- beta + dispBeta  
  
  if(rbinom(n = 1, size = 1, prob = .5)){
    propBeta <- beta + rbind(0, matrix(rnorm(n = prod(dim(beta)) - dim(beta)[2], mean = 0, sd = .1), nrow = (dim(beta)[1]-1), ncol = dim(beta)[2]))
  }
  
  propLogLikVal <- mvProbitLogLik( formula = formulaMvProbit, coef = c( propBeta ), sigma = propSigma, data = myData, algorithm = "GHK" )
  if(exp(sum(propLogLikVal) - sum(logLikVal)) > runif(1, 0, 1)){
    sigma <- propSigma
    beta <- propBeta
    logLikVal <- propLogLikVal
  }
  
  sigmaSamples[,,i] <- sigma
  betaSamples[,,i] <- beta
  print(paste0("iter ", i, ", loglik: ", sum(logLikVal), ", trueValuesLogLik: ", sum(logLikValTrue)))
}

plot(density(sigmaSamples[2,3,][-(1:1000)][1:5500])); abline(v = sig[2,3], col = 2, lwd = 2)
plot(density(sigmaSamples[1,3,][-(1:1000)][1:5500])); abline(v = sig[1,3], col = 2, lwd = 2)
plot(density(sigmaSamples[1,2,][-(1:1000)][1:5500])); abline(v = sig[1,2], col = 2, lwd = 2)

meanLiabs
plot(density(betaSamples[-1,,][1,1,][-(1:1000)][1:5500])); abline(v = meanLiabs[1,1], col = 2, lwd = 2)
plot(density(betaSamples[-1,,][1,2,][-(1:1000)][1:5500])); abline(v = meanLiabs[1,2], col = 2, lwd = 2)
plot(density(betaSamples[-1,,][1,3,][-(1:1000)][1:5500])); abline(v = meanLiabs[1,3], col = 2, lwd = 2)
plot(density(betaSamples[-1,,][2,1,][-(1:1000)][1:5500])); abline(v = meanLiabs[2,1], col = 2, lwd = 2)
plot(density(betaSamples[-1,,][2,2,][-(1:1000)][1:5500])); abline(v = meanLiabs[2,2], col = 2, lwd = 2)
plot(density(betaSamples[-1,,][2,3,][-(1:1000)][1:5500])); abline(v = meanLiabs[2,3], col = 2, lwd = 2)
plot(density(betaSamples[-1,,][3,1,][-(1:1000)][1:5500])); abline(v = meanLiabs[3,1], col = 2, lwd = 2)
plot(density(betaSamples[-1,,][3,2,][-(1:1000)][1:5500])); abline(v = meanLiabs[3,2], col = 2, lwd = 2)
plot(density(betaSamples[-1,,][3,3,][-(1:1000)][1:5500])); abline(v = meanLiabs[3,3], col = 2, lwd = 2)



## goofing around ##

rownames(traitCounts) <- populations
colnames(traitCounts) <- colnames(sig)
tet2est <- tetrachoric2(dat = discrTraitsAll, method = "Ol")$rho
tetest <- tetrachoric(discrTraitsAll)$rho
tet2est
tetirtoys <- tet(discrTraitsAll)
diag(tetirtoys) <- rep(1, dim(tetirtoys)[1])
sig
tet2est - sig
mean(abs(tet2est - sig))
polychoric2

psEst <- matrix(data = 0.5, nrow = nrow(meanLiabs), ncol = ncol(meanLiabs))
meanLiabsEst <- qnorm(psEst, mean = 0, sd = 1)
rownames(meanLiabsEst) <- populations
table(as.data.frame(discrTraitsAll[,2:3]))
cos(pi / (1 + ((95*138) / (12*55)) ^ 0.5))

#functions to construct a tetrachoric correlation matrix from scratch

#compute correlation coefficients
rTetr <- function(a, b, c, d){
  cos(pi / (1 + ((a*d) / (b*c)) ^ 0.5))
}

tetrCorrMat <- diag(x = 1, nrow = length(discrTraitsAll[1,]))
for (i in 1:length(discrTraitsAll[1,])){
  for(j in (i):length(discrTraitsAll[1,])){
    contTab <- table(as.data.frame(discrTraitsAll[,c(i,j)]))
    #print(contTab)
    tetrCorrMat[i,j] <- rTetr(contTab[1,1], contTab[1,2], contTab[2,1], contTab[2,2])
  }
}

tetrCorrMat <- tetrCorrMat + t(tetrCorrMat) - diag(1, dim(tetrCorrMat)[1])
tetrCorrMat
isSymmetric(tetrCorrMat)
is.positive.definite(tetrCorrMat)
tetrCorrMat <- cov2cor(as.matrix(nearPD(tetrCorrMat)$mat))
tetrCorrMat - sig
tetrCorrMat - tet2est
mean(abs(tetrCorrMat - tet2est))
mean(abs(tetrCorrMat - sig))
mean(abs(tetrCorrMat - tetirtoys))

#let's try to make an MCMC sampler

iterMCMC <- 2e6
numTraits <- length(cov[1,])
sampInt <- 200
up <- 0
down <- 0
same <- 0
MCMCTree <- unroot(rtree(n=numTips)) ##TREE SHOULD BE UNROOTED
MCMCTree$edge.length <- rep(1, (2*numTips-3))
MCMCTree$tip.label <- sample(populations, length(populations), replace = F)
sigmaBM <- sig
MCMCLL <- BMpruneLL(meanLiabsEst, sigmaBM, MCMCTree) + 
  sum(sapply(1:length(MCMCTree$edge.length), function(edge) dexp(x = MCMCTree$edge.length[edge], rate = 2, log = T))) + #exponential priors on BLs
  sum(dbinom(prob = psEst, x = traitCounts, size = nObs, log = T)) 
  # sum(dbeta(x = psEst, shape1 = 2, shape2 = 1, log = T))
  
MCMCTrees <- rmtree(N = iterMCMC/sampInt, n = 3, rooted = F)
testLLs <- vector(length = iterMCMC)
estimatedMeanLiabilities <- array(data = meanLiabsEst, dim = c(dim(meanLiabsEst), 1))

for (i in 1:iterMCMC){
  if(i%%100==0){print(c(i, testLLs[i-1]))}
  prop <- sample(x = 1:30, size = 1)
  edgeBAD <- T
  while(edgeBAD){
    
    if(prop <= 5){testTree <- rNNI(MCMCTree)} else
      if(prop > 5 & prop < 11){testTree <- rSPR(MCMCTree)} else
      if(prop > 10 & prop < 16){testTree <- MCMCTree
        edge <- sample(1:length(testTree$edge.length), size=1)
        testTree$edge.length[edge] <- testTree$edge.length[edge] + rnorm(n = 1, mean = 0, sd = .05)} else
      if(prop==17){testTree <- MCMCTree
        edge <- sample(1:length(testTree$edge.length), size=1)
        testTree$edge.length[edge] <- testTree$edge.length[edge] + rnorm(n = 1, mean = 0, sd = .05)
        edge <- sample(1:length(testTree$edge.length), size=1)
        testTree$edge.length[edge] <- testTree$edge.length[edge] + rnorm(n = 1, mean = 0, sd = .05)} else
      if(prop>17 & prop < 20){testTree <- MCMCTree
        testTree$edge.length <- testTree$edge.length + rnorm(n = 1, mean = 0, sd = .05)} 
      if(prop>19 & prop < 27){testTree <- MCMCTree
        psEstTest <- psEst
        pop <- sample(1:length(psEstTest[,1]), size = 1)
        psEstTest[pop,] <- psEstTest[pop,] + 
          rmvnorm(n = 1, mean = rep(0, length(psEstTest[pop,])), sigma = diag(.002, nrow = length(psEstTest[pop,])))
        psEstTest <- abs(psEstTest)
        tooBig <- which(psEstTest > 1, arr.ind = T)
        psEstTest[tooBig] <- 2 - psEstTest[tooBig]
        meanLiabsEst <- qnorm(psEstTest, mean = 0, sd = 1)
        rownames(meanLiabsEst) <- populations
        }
      if(prop>26 & prop < 31){testTree <- MCMCTree
      psEstTest <- psEst
      psEstTest <- psEstTest + 
        rmvnorm(n = length(psEstTest[,1]), mean = rep(0, length(psEstTest[pop,])), sigma = diag(.002, nrow = length(psEstTest[pop,])))
      psEstTest <- abs(psEstTest)
      tooBig <- which(psEstTest > 1, arr.ind = T)
      psEstTest[tooBig] <- 2 - psEstTest[tooBig]
      meanLiabsEst <- qnorm(psEstTest, mean = 0, sd = 1)
      rownames(meanLiabsEst) <- populations
      } 

      
    if(all(testTree$edge.length > 0)) {edgeBAD <- F}
  }
  testLL <- BMpruneLL(meanLiabsEst, sigmaBM, testTree) + 
    sum(sapply(1:length(testTree$edge.length), function(edge) dexp(x = testTree$edge.length[edge], rate = 2, log = T))) +
    sum(dbinom(prob = psEstTest, x = traitCounts, size = nObs, log = T)) 
  
  testLLs[i] <- testLL
  if(testLL > MCMCLL){
    MCMCLL <- testLL
    MCMCTree <- testTree
    psEst <- psEstTest
    up <- up + 1
  } else if (runif(n = 1, min = 0, max = 1) < exp(testLL - MCMCLL)){
    MCMCLL <- testLL
    MCMCTree <- testTree
    psEst <- psEstTest
    down <- down + 1
    # print("lower")    }
  } else {same <- same + 1}
  if(i %% sampInt == 0){
    MCMCTrees[[i/sampInt]] <- MCMCTree
    estimatedMeanLiabilities <- abind(estimatedMeanLiabilities, meanLiabsEst, along = 3)
  }
  if(i %% 5000 == 0){
    par(mfrow = c(1,2))
    plot(tree); axisPhylo(); title("true, data-generating tree")
    plot(MCMCTree); axisPhylo(); title(paste("Current Tree, iteration:", i))
  }
}
prematurePause <- MCMCTrees[sapply(1:length(MCMCTrees), function(x) MCMCTrees[[x]]$Nnode != 1)]
MCMCTrees_noBurn <- prematurePause[-(1:(length(prematurePause)/5))]
conTree <- consensus(MCMCTrees_noBurn, p=.35)
conTree <- root(conTree, "BUSHMAN", resolve.root = T)
plot(conTree)
title("Inferred Majority-Rule Consensus Tree (no BLs)")
plot(tree)

#mean liability estimate
foo <- matrix(0, nrow = dim(estimatedMeanLiabilities)[1], ncol = dim(estimatedMeanLiabilities)[2]); 
for(i in 1:dim(estimatedMeanLiabilities)[3]){foo <- foo + estimatedMeanLiabilities[,,i]}
meanLiabEstMean <- foo / dim(estimatedMeanLiabilities)[3]
meanLiabs

rangeS <- floor(length(estimatedMeanLiabilities[1,1,])/5)
rangeE <- length(estimatedMeanLiabilities[1,1,])

par(mfrow = c(5,5))
estimatedMeanLiabilities_noBurn <- estimatedMeanLiabilities[,,rangeS:rangeE]
for(i in 1:dim(estimatedMeanLiabilities_noBurn[,,1])[1]){
  for(j in 1:dim(estimatedMeanLiabilities_noBurn[,,1])[2]){
    plot(density(estimatedMeanLiabilities_noBurn[i,j,]), main = ""); abline(v = meanLiabs[i,j])
  }
}

### POSTERIOR PROBABILITIES ###

pp <- prop.clades(conTree, MCMCTrees_noBurn, rooted = F)/length(MCMCTrees_noBurn)

## plot our tree with PPs
plot(conTree)
edgelabels(round(pp, 3), sapply(1:conTree$Nnode + length(conTree$tip.label), 
                                match, conTree$edge[, 2]), bg = "white")
title("Inferred Consensus Tree (no BLs)")
plot(tree); title("true, data-generating tree")

