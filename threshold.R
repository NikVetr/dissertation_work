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
library(rethinking)


############## FUNCTIONS ##############


#here we have a function that accepts the name of a given tip and a tree (in phylo format) in which that tip is found and returns the names of all the branches
#that lead from that tip back to the root
branchesToRoot <- function(tip, tree) {
  bar <- 1 #initialize a counter variable
  branches <- vector() #initialize an empty vector to hold the branches
  branches[bar] <- which(tree$edge[, 2] == tip) #find the row of the branch/topology matrix in which the desired tip terminates
  repeat{ #keep going through the tree until you reach the root, which is a node that only appears in the "originating" column and not the "terminating" column
    bar <- bar + 1 #increase the counter variable
    if(is.na(which(tree$edge[, 2] == tree$edge[branches[bar-1], 1])[1])){
      break #if you try to look at the terminating column and don't find the node that originated the previous branch, you're at the root! Congrats! Exit the repeat() loop
    }
    branches[bar] <- which(tree$edge[, 2] == tree$edge[branches[bar-1],1]) #otherwise keep going through the tree until you hit the root
  }
  return(branches) #and return all the rows numbers of the tip -> root branches (which are also the branch labels)
}


#this function finds which branches connect two tips on a tree (aka the contrast between those branches).
#to do this it runs the branchesToRoot() function on both tips and then it finds the branches they *do not* share going to the root
#since to share a branch they'd need to share its terminating node. It then returns that set of unshared branches, which represents the 
#path from one tip to another
branchConnect <- function(tip1, tip2, tree) {
  a <- branchesToRoot(tip1, tree)
  b <- branchesToRoot(tip2, tree)
  return(c(setdiff(a, b), setdiff(b, a))) #setdiff() gives the elements of one vector that aren't in another vector
}

#can accept both the index and the name of the tip; returns a list containing branches containing that tip and all other tips
branchesToTip <- function(tree, tip) {
  if(is.character(tip)) {return(sapply(tree$tip.label[tree$tip.label != tip], function(x) 
    branchConnect(which(tree$tip.label == tip), which(tree$tip.label == x), tree), simplify = F)
  )}
  if(is.numeric(tip)) {return(sapply(tree$tip.label[-tip], function(x) branchConnect(tip, which(tree$tip.label == x), tree), simplify = F)
  )}
}

#this function returns the shared branch lengths of n-1 tips from some planted tip
plantRootCov <- function (tip, tree) {
  tipShare <- branchesToTip(tree, tip)
  sharedBLs <- matrix(nrow = length(tipShare), ncol = length(tipShare))
  for(i in 1:length(tipShare)){
    for(j in 1:length(tipShare)){
      sharedBLs[i,j] <- (sum(tree$edge.length[intersect(tipShare[[i]],tipShare[[j]])]))
    }
  }
  rownames(sharedBLs) <- colnames(sharedBLs) <- names(tipShare)
  return(sharedBLs)
}

#traits is a matrix of traits, sigmaBM is your mvBM rate matrix, tree is your tree, tip is the tip to plant (which doesn't affect the result)
PlantQuasiLL <- function(traits, sigmaBM, tree, tip=1){
  if(!isTRUE(all.equal(rownames(traits), tree$tip.label))){ #reorder traits if in incorrect order
    traits <- traits[tree$tip.label,]
  }
  traitStack <- c(t(as.matrix(traits[rownames(traits)[-tip],])))
  return(dmvnorm(traitStack, rep(traits[tip,], times=(nrow(traits)-1)), sigma = kronecker(plantRootCov(tip, tree), sigmaBM), log=T))
} #this function returns the quasiLL of all tip traits from some planted tip

CovMatrix <- function(num){ #generate a random covariance matrix with size of num x num
  sigmas <-  (runif(num, min = .5, max = 2))
  Rho <- matrix(runif(num^2, min=-.2, max=.4), nrow=num)
  diag(Rho) <- rep(1,num)
  Rho[lower.tri(Rho)] = t(Rho)[lower.tri(Rho)]
  return((diag(sigmas) %*% Rho %*% diag(sigmas)))
}

CovMatrix2 <- function(num){ #generate a random covariance matrix with size of num x num
  sigmas <-  (runif(num, min = 0, max = 2))
  Rho <- matrix(runif(num^2, min=-.1, max=.1), nrow=num)
  diag(Rho) <- rep(1,num)
  Rho[lower.tri(Rho)] = t(Rho)[lower.tri(Rho)]
  return((diag(sigmas) %*% Rho %*% diag(sigmas)))
}

#simulate data
numTraits <- 5
numTips <- 5
numIndiv <- 20
tree <- rtree(n=numTips)
sig <- as.matrix(nearPD(CovMatrix(numTraits))$mat)
traits <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = numTraits, sigma=sig, mu= rep(0, numTraits)))
data <- lapply(1:length(traits[,1]), function(y)
  sapply(1:length(traits[1,]), function(x) rbinom(prob = pnorm(0, sd = 1, mean = traits[y,x]), n = numIndiv, size = 1))
)


#trait frequencies in each population

#given some frequency of presence/absence, returns a Bayesian estimate of the mean liability from weak beta prior
liabPost <- function(nHave, nTot, shape1 = 1.01, shape2 = 1.01) {
  pseq <- seq(0, 1, by = .00001)
  prior <- dbeta(x = pseq, shape1 = shape1, shape2 = shape2)
  posterior <- dbinom(nHave, nTot, pseq) * prior; posterior <- posterior / sum(posterior)
  freqDistr <- sample(x = pseq, size = 1e6, replace = T, prob = posterior)
  liabs <- qnorm(freqDistr, mean = 0, sd = 1)
  return(c(mean(liabs), sd(liabs)))
}

freqs <- sapply(1:length(data), function (x) colSums(data[[x]]))/numIndiv
freqs[freqs == 0] <- 1/numIndiv
freqs[freqs == 1] <- 1 - 1/numIndiv

infTraits <- -t(qnorm(freqs, mean = 0, sd = 1))
infTraits/traits


#testing it out
#... works ok

iterOptim <- 5000
numTips <- 5
numTraits <- 50
numIndiv <- 50
iterSim <- 250
timesCorrect <- 0
correctTally <- rep(0, iterSim)
RFDists <- vector()
trueTrees <- list()
bestTrees <- list()
for (j in 1:iterSim){
  change <- rep(0,19)
  trueTree <- (rtree(n=numTips))
  
  #   #test long branch attraction
  #   tips <- sample(1:numTips, 2)
  #   branch1 <- which(trueTree$edge[,2] == tips[1])
  #   branch2 <- which(trueTree$edge[,2] == tips[2])
  #   longBranch <- (trueTree$edge.length[branch1] + trueTree$edge.length[branch2]) * 3
  #   trueTree$edge.length[branch1] <- longBranch
  #   trueTree$edge.length[branch2] <- longBranch
  
  bestTree <- unroot(rtree(n=numTips))
  bestTree$edge.length <- bestTree$edge.length * 1 #check for bias in starting tree BLs
  
  #   sigmaBM <- CovMatrix(numTraits)
  #   sigmaBM <- as.matrix(nearPD(sigmaBM)$mat)
  
  sigmaBM <- as.matrix(nearPD(CovMatrix2(numTraits))$mat) 
  
  #while(!is.positive.definite(sigmaBM)){
  #  sigmaBM <- CovMatrix(numTraits)
  #}
  
  traits <- mvSIM(tree = trueTree, nsim = 1, model = "BM1", param = list(ntraits = numTraits, sigma=sigmaBM, mu= rep(0, numTraits)))

  
  data <- lapply(1:length(traits[,1]), function(y)
    sapply(1:length(traits[1,]), function(x) rbinom(prob = pnorm(0, sd = 1, mean = traits[y,x]), n = numIndiv, size = 1))
  )
  
  #trait frequencies in each population
  freqs <- sapply(1:length(data), function (x) colSums(data[[x]]))/numIndiv
  freqs[freqs == 0] <- 1/numIndiv
  freqs[freqs == 1] <- 1 - 1/numIndiv
  
  infTraits <- -t(qnorm(freqs, mean = 0, sd = 1))
  rownames(infTraits) <- rownames(traits)
  traits <- infTraits
  
  #traits <- traits*0
  trueTree <- unroot(trueTree)
  bestLL <- PlantQuasiLL(traits, sigmaBM, bestTree, tip=1)
  trueLL <- PlantQuasiLL(traits, sigmaBM, trueTree, tip=1)
  for (i in 1:iterOptim){
    prop <- sample(x = 1:19, size = 1)
    if(prop <= 6){testTree <- rNNI(bestTree)} else
      if(prop > 6 & prop < 13){testTree <- rSPR(bestTree)} else
        if(prop > 12 & prop < 17){testTree <- bestTree
        edge <- sample(1:length(testTree$edge.length), size=1)
        testTree$edge.length[edge] <- testTree$edge.length[edge]*runif(n = 1, min = .8, max = 1.2)} else
          if(prop==18){testTree <- bestTree
          edge <- sample(1:length(testTree$edge.length), size=1)
          testTree$edge.length[edge] <- testTree$edge.length[edge]*runif(n = 1, min = .8, max = 1.2)
          edge <- sample(1:length(testTree$edge.length), size=1)
          testTree$edge.length[edge] <- testTree$edge.length[edge]*runif(n = 1, min = .8, max = 1.2)} else
            if(prop==19){testTree <- bestTree
            testTree$edge.length <- testTree$edge.length*runif(n = 1, min = .8, max = 1.2)} 
    
    testLL <- PlantQuasiLL(traits, sigmaBM, testTree, tip=1)
    if(testLL > bestLL){
      bestLL <- testLL
      bestTree <- testTree
      change[prop] <- change[prop] + 1
    }
    if(i %% 100 == 0){print(c(i, paste0("trueLL=",trueLL), paste0("bestLL=",bestLL)))}
    if(i %% 1000 == 0){print(paste("Tree Rearrangements = ", sum(change[1:12]), 
                                   "; Branch Length Changes = ", sum(change[13:19])))
      print(paste0("Robinson-Foulds Distance = ", RF.dist(trueTree, bestTree)))}
  }
  if(RF.dist(trueTree, bestTree) == 0) {timesCorrect <- timesCorrect + 1}
  print(c(j, timesCorrect/j))
  correctTally[j] <- timesCorrect/j
  RFDists[j] <- RF.dist(trueTree, bestTree)
  trueTrees[[j]] <- trueTree
  bestTrees[[j]] <- bestTree
}
