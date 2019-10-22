#I forget which packages I used for this so just load these and all the code should run
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
library(tictoc)

### explore on a single branch ###

#sample a covariance matrix
dim <- 10
weGood <- FALSE
while(!weGood){ #silly do-while loop to make sure we get a safe correlation matrix
  cc <- 0
  cor <- matrix(data = cc, nrow = dim, ncol = dim)
  cc <- runif(min = 0.4, max = 0.5, n = dim*(dim - 1)/2) #sample hopefully plausible correlations
  cor[upper.tri(cor)] <- cc
  cor <- cor + t(cor)
  diag(cor) <- 1
  sds <- runif(n = dim, min = 1, max = 10)
  cov <- diag(sds) %*% cor %*% diag(sds)
  
  if(all(eigen(cov)$values > 0)){
    weGood <- TRUE
  }
}

#simulate node 1's traits
start <- rmvnorm(n=1, mean = rep(0, dim), sigma = cov)

#simulate character evolution along branch to node 2
BL <- 10
end <- rmvnorm(n=1, mean = start, sigma = cov*BL)

#calculate likelihood as multivariate normal density
mvdens <- dmvnorm(end, start, sigma = cov*BL, log = T)

#eigentransform the data to PC scores
startEig <- start %*% eigen(cov)$vectors
endEig <- end %*% eigen(cov)$vectors

#calculate the product of univariate normal densities of transformed chars
uvdens <- dmvnorm(endEig, startEig, sigma = diag(eigen(cov)$values)*BL, log = T) #done with mvtnorm density where sigma is a diag matrix
uvdens
vars <- eigen(cov)$values*BL
uvdens <- sum(sapply(1:dim, function(x) dnorm(endEig[x], startEig[x], sd = (vars[x])^0.5, log = T))) #same as above uvdens
mvdens
uvdens
nBL <- 5
vars <- eigen(cov)$values*nBL
sum(sapply(1:dim, function(x) dnorm(endEig[x], startEig[x], sd = (vars[x])^0.5, log = T)))
dmvnorm(end, start, sigma = cov*nBL, log = T)

#sample a new covariance matrix -- does it stil hold across rate matrix variation?
weGood <- FALSE
while(!weGood){ #silly do-while loop to make sure we get a safe correlation matrix
  cc <- 0
  cor <- matrix(data = cc, nrow = dim, ncol = dim)
  cc <- runif(min = 0.4, max = 0.5, n = dim*(dim - 1)/2) #sample hopefully plausible correlations
  cor[upper.tri(cor)] <- cc
  cor <- cor + t(cor)
  diag(cor) <- 1
  sds <- runif(n = dim, min = 1, max = 10)
  cov <- diag(sds) %*% cor %*% diag(sds)
  
  if(all(eigen(cov)$values > 0)){
    weGood <- TRUE
  }
}
#retransform data
startEig <- start %*% eigen(cov)$vectors
endEig <- end %*% eigen(cov)$vectors

#recalculate likelihoods
nBL <- 5
vals <- eigen(cov)$values
vars <- vals*nBL
sum(sapply(1:dim, function(x) dnorm(endEig[x], startEig[x], sd = (vars[x])^0.5, log = T))) #univ likelihood of transf data
dmvnorm(end, start, sigma = cov*nBL, log = T) #multivariate likelihood
#same thing!

#is there a more elegant way to get R to execute a section of code with ctrl-enter than a 1-fold for-loop?
#ctrl-enter will otherwise only do a single statement (or a single line if you have that option toggled)
#and I'm too lazy to always highlight stuff, and need tic() and toc() (or other benchmarking commands of choice) to run asap
numIter <- 1e4
BLs <- runif(n = numIter, min = 10, max = 100)

for(gogogo in 1:1){
  tic()
  for(i in 1:numIter){ #use a for loop cos nested sapply's look silly
    vars <- vals*BLs[i]
    sum(sapply(1:dim, function(x) dnorm(endEig[x], startEig[x], sd = (vars[x])^0.5, log = T))) #univ likelihood of transf data
  }
  toc()
}

for(gogogo in 1:1){
  tic()
  for(i in 1:numIter){
    dmvnorm(end, start, sigma = cov*BLs[i], log = T) #multivariate likelihood
  }
  toc()
}

#explore on a tree
BMpruneLLmvt <- function(traits, sig, tree){ #this function evaluates multivariate normal densities
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

#lets test it out on some data
tree <- pbtree(b=1, n=10)
traits <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = dim, sigma=cov, mu= rep(0,dim)))
BMpruneLLmvt(traits, cov, tree)

#do we get the same result with a diagonal matrix?
vecs <- eigen(cov)$vectors
vals <- eigen(cov)$values
eigTraits <- traits %*% vecs
BMpruneLLmvt(eigTraits, diag(vals), tree) #gives the same result

#lets toss it into the pruning algorithm directly
BMpruneLLuvt <- function(rawTraits, sig, tree){ #this function evaluates univariate normal densities
  vecs <- eigen(cov)$vectors
  vals <- eigen(cov)$values
  traits <- rawTraits %*% vecs
  dim <- dim(sig)[1]
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
    LLs[i] <- sum(sapply(1:dim, function(x) dnorm(contrast[x], 0, sd = (vals[x]*BLength)^0.5, log = T)))
    
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

#these appear to give the same answer?
cov <- matrix(rnorm(100), 10, 10); cov <- cov %*% t(cov); dim <- dim(cov)[1]
tree <- pbtree(b=1, n=10)
traits <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = dim, sigma=cov, mu= rep(0,dim)))
BMpruneLLmvt(traits, cov, tree)
BMpruneLLuvt(traits, cov, tree)
BMpruneLLuvtChol(traits, cov, tree)

#now let's diagonalize via cholesky decomp
BMpruneLLuvtChol <- function(rawTraits, sig, tree){ #this function evaluates univariate normal densities
  U <- chol(sig)
  L <- t(U)
  Li <- solve(L)
  detcov <- det(sig)
  traits <- t(Li %*% t(rawTraits)); rownames(traits) <- rownames(rawTraits)
  dim <- dim(sig)[1]
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
    LLs[i] <- sum(dnorm(contrast / BLength^0.5, mean = 0, sd = 1, log = T)) - log(detcov*(BLength^dim))/2
    
    #if not pretransforming
    # LLs[i] <- sum(dnorm((Li/BLength^0.5) %*% contrast, mean = 0, sd = 1, log = T)) - log(detcov*(BLength^dim))/2
    # LLs[i] <- sum(dnorm(solve(L*BLength^0.5) %*% contrast, mean = 0, sd = 1, log = T)) - log(detcov*(BLength^dim))/2
    # LLs[i] <- sum(dnorm(solve(t(chol(sig*BLength))) %*% contrast, mean = 0, sd = 1, log = T)) - log(detcov*(BLength^dim))/2
    
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


##############################
#benchmark across lotsa trees#
##############################
numTrees <- 10
numTips <- 20
numTraits <- 300

#get a correctly sized covariance matrix
if (numTraits < 100){
weGood <- FALSE
while(!weGood){ #silly do-while loop to make sure we get a safe correlation matrix
  dim <- numTraits
  cc <- 0
  cor <- matrix(data = cc, nrow = dim, ncol = dim)
  cc <- runif(min = 0.4, max = 0.5, n = dim*(dim - 1)/2)
  cor[upper.tri(cor)] <- cc
  cor <- cor + t(cor)
  diag(cor) <- 1
  sds <- runif(n = dim, min = 1, max = 10)
  cov <- diag(sds) %*% cor %*% diag(sds)
  
  if(all(eigen(cov)$values > 0)){
    weGood <- TRUE
  }
}
} else{
  dim <- numTraits
  cc <- 0.3
  cor <- matrix(data = cc, nrow = dim, ncol = dim)
  diag(cor) <- 1
  sds <- runif(n = dim, min = 1, max = 10)
  cov <- diag(sds) %*% cor %*% diag(sds)
  
}

#simulate character data on some tree
tree <- pbtree(b=1, n=numTips)
traits <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = dim, sigma=cov, mu= rep(0,dim)))
BMpruneLLmvt(traits, cov, tree)
BMpruneLLuvt(traits, cov, tree)

#sample a buncha trees for benchmarking
trees <- rmtree(N = numTrees, n = numTips)

#run the mvtnormal benchmark
for(gogogo in 1:1){
  tic()
  likesMVT <- sapply(1:numTrees, function(x) BMpruneLLmvt(traits, cov, trees[[x]]))
  toc()
}

#run the univariate normal benchmark
for(gogogo in 1:1){
  tic()
  likesUVT <- sapply(1:numTrees, function(x) BMpruneLLuvt(traits, cov, trees[[x]]))
  toc()
}

#this fails for a lot of runs, probably b/c of rounding or machine precision or underflow or something idk
likesMVT == likesUVT

#but looking at their difference
mean(likesMVT - likesUVT)
#it is very tiny


