#inference under the phylogenetics multivariate brownian threshold model

#load compute libraries
library(gpuR)
library(Rcpp)
library(tictoc); library(microbenchmark)

#load phylogenetic libraries
library(mvMORPH)
library(TESS)

tree <- tess.sim.taxa(n = 1, nTaxa = 8, lambda = 1, mu = 0, max = 1E3)[[1]]
R <- diag(c(1.5,2.3,1.7)) %*% matrix(c(1,0.5,0.6,0.5,1,0.4,0.6,0.4,1), 3, 3) %*% diag(c(1.5,2.3,1.7))
nt <- 8; R <- cov[1:nt, 1:nt]
R <- diag(nt)
root_state <- rep(0, dim(R)[1])
traits <- mvSIM(tree, model="BM1", param=list(sigma=R, theta=root_state)) 

prunes <- prune(tree)
prunes <- prunePCV(vcv.phylo(tree))
microbenchmark(prune(tree), vcv.phylo(tree), prunePCV(vcv.phylo(tree)))

contrasts <- prunes[[1]] %*% traits[rownames(prunes[[3]])[1:length(tree$tip.label)],]
#BLs_alt <- diag(prunes[[1]] %*% vcv.phylo(tree) %*% t(prunes[[1]]))
phylocv <- vcv.phylo(reroot(tree, node.number = which(tree$tip.label == rownames(traits)[1])))[-1,-1]
dmvnorm(x = c(t(traits[-1,])), mean = rep(traits[1,], length(traits[-1,1])), sigma = kronecker(phylocv,R), log = T)
sum(sapply(1:(length(tree$tip.label)-1), function(x) dmvnorm(contrasts[x,], sigma = R*prunes[[2]][x], log = T)))
if(all(R[upper.tri(R)] == 0))(
  sum(sapply(1:(length(tree$tip.label)-1), function(x) dnorm(contrasts[x,], sd = sqrt(prunes[[2]][x]), log = T)))
)
BMpruneLL(traits = traits, sig = R, tree = tree)
BMpruneLLuvtChol(rawTraits = traits, sig = R, tree = tree)

t_traits <- t(solve(t(chol(R))) %*% t(contrasts))
detR <- det(R)
dim <- dim(R)[1]
CBLs <- prunes[[2]]
sum(sapply(1:(length(tree$tip.label)-1), function(x) sum(dnorm(x = t_traits[x,] / sqrt(CBLs[x]), mean = 0, sd = 1, log = T)) - log(detR*(CBLs[x]^dim))/2))

#quick benchmarking of pruning algorithm
# microbenchmark(
#   dmvnorm(x = c(t(traits[-1,])), mean = rep(traits[1,], length(traits[-1,1])), sigma = kronecker(phylocv,R), log = T),
#   sum(sapply(1:(length(tree$tip.label)-1), function(x) dmvnorm(contrasts[x,], sigma = R*prunes[[2]][x], log = T))),
#   BMpruneLL(traits = traits, sig = R, tree = tree),
#   BMpruneLLuvtChol(rawTraits = traits, sig = R, tree = tree),
#   vcv.phylo(tree),
#   prunes <- prune(tree),
#   contrasts <- prunes[[1]] %*% traits[rownames(prunes[[3]])[1:length(tree$tip.label)],],
#   t_traits <- t(solve(t(chol(R))) %*% t(contrasts)),
#   detR <- det(R),
#   sum(sapply(1:(length(tree$tip.label)-1), function(x) sum(dnorm(x = t_traits[x,] / sqrt(CBLs[x]), mean = 0, sd = 1, log = T)) - log(detR*(CBLs[x]^dim))/2))
# )

#simulate thresholded traits at the population level
n_indiv <- 20
threshold <- rep(0, dim(traits)[2])
traits_indiv <- lapply(1:dim(traits)[1], function(tip) rmvnorm(n = n_indiv, mean = traits[tip,], sigma = R)); names(traits_indiv) <- rownames(traits)
traits_indiv_discr <- lapply(1:dim(traits)[1], function(tip) t(sapply(1:n_indiv, function(indiv) as.numeric(traits_indiv[[tip]][indiv,] > threshold)))); names(traits_indiv_discr) <- names(traits_indiv) 

#identify unique site patterns
traits_indiv_discr_unique <- lapply(1:dim(traits)[1], function(tip) uniqueSP(traits_indiv_discr[[tip]]))

#find probs of these site patterns and multiply by freq
logLLs <- sapply(1:length(traits_indiv_discr_unique), function(tip) sum(log(multi_pmvnorm(mean = traits[tip,], sigma = R, binaries = traits_indiv_discr_unique[[tip]]$SPs)) * traits_indiv_discr_unique[[tip]]$counts))
microbenchmark(sapply(1:length(traits_indiv_discr_unique), 
                      function(tip) sum(log(multi_pmvnorm(mean = traits[tip,], sigma = R, binaries = traits_indiv_discr_unique[[tip]]$SPs)) * traits_indiv_discr_unique[[tip]]$counts))
               , times = 10)

#verify anticipated operation
mean <- c(1,2,3) - 2
cov <- R[1:3,1:3]/20
binaries <- t(as.matrix(c(1,1,1)))

nsamp <- 1E4
Xs <- rmvnorm(nsamp, mean = mean, sigma = cov)
sum(sapply(1:nsamp, function(samp)  all(sapply(1:length(Xs[samp,]), function(tr) ifelse(binaries[1,tr], Xs[samp,tr], -Xs[samp,tr]) > 0)))) / nsamp
multi_pmvnorm(mean = mean, sigma = cov, binaries = binaries)
#OK seems to work


#########################################
### let's try to make an MCMC sampler ###
#########################################

#first let's resimulate data
nTaxa <- 8
tree <- tess.sim.taxa(n = 1, nTaxa = nTaxa, lambda = 1, mu = 0, max = 1E3)[[1]]
R <- diag(c(1.5,2.3,1.7)) %*% matrix(c(1,0.5,0.6,0.5,1,0.4,0.6,0.4,1), 3, 3) %*% diag(c(1.5,2.3,1.7)); rownames(R) <- colnames(R) <- paste0("tr", 1:dim(R)[1])
root_state <- rep(0, dim(R)[1])
traits <- mvSIM(tree, model="BM1", param=list(sigma=R, theta=root_state)) 

#simulate thresholded traits at the population level
n_indiv <- 20
threshold <- rep(0, dim(traits)[2])
traits_indiv <- lapply(1:dim(traits)[1], function(tip) rmvnorm(n = n_indiv, mean = traits[tip,], sigma = R)); names(traits_indiv) <- rownames(traits)
traits_indiv_discr <- lapply(1:dim(traits)[1], function(tip) t(sapply(1:n_indiv, function(indiv) as.numeric(traits_indiv[[tip]][indiv,] > threshold)))); names(traits_indiv_discr) <- names(traits_indiv) 

#now try to infer data-generating tree
iter <- 1e3
nTraits <- dim(R)[1]
thin <- 2 # set the thinning rate
MCMCTree <- unroot(rtree(n=numTips)) ##TREE SHOULD BE UNROOTED
MCMCTree$edge.length <- rep(1, (2*numTips-3))
MCMCTree$tip.label <- sample(populations, length(populations), replace = F)
sigmaBM <- sig

trees <- rmtree(N = iterMCMC/thin, n = nTaxa, rooted = F)
LLs <- vector(length = iterMCMC)
est_traits <- array(data = 0, dim = c(nTaxa, nTraits, iter / thin), dimnames = list(tree$tip.label, rownames(R), 1:(iter/thin)))

for (i in 1:iter){
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
  if(i %% thin == 0){
    MCMCTrees[[i/thin]] <- MCMCTree
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

multi_pmvnorm <- function(mean, sigma, binaries, algorithm = c("GenzBretz", "Miwa")[1]){
  upper <- binaries
  upper[] <- 0
  lower <- upper
  upper[binaries == 1] <- Inf
  lower[binaries == 0] <- -Inf
  return(sapply(1:dim(binaries)[1], function(indiv) pmvnorm(lower = lower[indiv,], upper = upper[indiv,], mean = mean, sigma = sigma, algorithm = algorithm)[1]))
}

uniqueSP <- function(indivs){
  SPs <- table(sapply(1:length(indivs[,1]), function(indiv) paste0(indivs[indiv,], collapse = "")))
  counts <- as.numeric(SPs)
  SPs <- t(sapply(1:length(SPs), function(SP) as.numeric(unlist(strsplit(names(SPs)[SP], "")))))
  return(list(SPs = SPs, counts = counts))
}

prune <- function(tree){ 
  ntips <- length(tree$tip.label)
  contrast <- matrix(0, nrow = ntips-1, ncol = ntips); colnames(contrast) <- tree$tip.label
  BLength <- rep(0, ntips - 1)
  
  traitTransforms <- matrix(0, nrow = 2 * ntips - 2, ncol = ntips)
  traitTransforms[1:ntips, 1:ntips] <- diag(ntips)
  rownames(traitTransforms) <- c(tree$tip.label, rep(NA, ntips-2))
  colnames(traitTransforms) <- tree$tip.label
  for(i in 1:(ntips - 1)){
    #Find two sister tips
    #preterminal nodes
    ptn <- tree$edge[sapply(1:length(tree$tip.label), function(t) which(tree$edge[,2] == t)),1]
    uniq <- unique(ptn)
    parent <- uniq[which(sapply(1:length(uniq), function (x) sum(uniq[x] == ptn)) > 1)][1]
    sisters <- tree$edge[tree$edge[,1] == parent, 2]
    sisters <- intersect(sisters, 1:(length(tree$tip.label)))
    
    #calculate contrast between two sister nodes
    contrast[i,] <- traitTransforms[tree$tip.label[sisters[1]],] - traitTransforms[tree$tip.label[sisters[2]],] 
    
    #calculate BL along contrast
    BLength[i] <- tree$edge.length[tree$edge[,2] == sisters[1]] + tree$edge.length[tree$edge[,2] == sisters[2]]
    
    #calculate multivariate normal density of contrast along Branch Length
    
    if(length(tree$tip.label) > 2){   
      tipName = paste(tree$tip.label[sisters], collapse = "+")
      tt_index <- which(is.na(rownames(traitTransforms)))[1]
      rownames(traitTransforms)[tt_index] <- tipName
      #Fix Trait Matrix
      nodeValues <- (traitTransforms[tree$tip.label[sisters[1]],]*tree$edge.length[tree$edge[,2] == sisters[2]] + 
                       traitTransforms[tree$tip.label[sisters[2]],]*tree$edge.length[tree$edge[,2] == sisters[1]]) / BLength[i]
      
      
      traitTransforms[tt_index,] <- nodeValues
      
      #Fix Tree
      newTip <- list(edge=matrix(c(2,1),1,2),
                     tip.label=tipName,
                     edge.length=(prod(tree$edge.length[tree$edge[,2] == sisters[1]], tree$edge.length[tree$edge[,2] == sisters[2]]) / BLength[i]) - tree$edge.length[tree$edge[,2] == sisters[1]],
                     Nnode=1)
      class(newTip)<-"phylo"
      oldw <- getOption("warn")
      options(warn = -1)
      tree <- bind.tree(x = tree, y = newTip, where = sisters[1])
      options(warn = oldw)
      tree <- drop.tip(tree, sisters[2])
    }
  }
  return(list(contrastTransformationMatrix = contrast, contrastBranchLengths = BLength, nodalTraitConditionals = traitTransforms))
}

treeCV <- vcv.phylo(tree)
prunePCV <- function(treeCV){ 
  ntips <- dim(treeCV)[1]
  contrast <- matrix(0, nrow = ntips-1, ncol = ntips); colnames(contrast) <- tree$tip.label
  BLength <- rep(0, ntips - 1)
  
  traitTransforms <- matrix(0, nrow = 2 * ntips - 2, ncol = ntips)
  traitTransforms[1:ntips, 1:ntips] <- diag(ntips)
  rownames(traitTransforms) <- c(rownames(treeCV), rep(NA, ntips-2))
  colnames(traitTransforms) <- rownames(treeCV)
  for(i in 1:(ntips - 1)){
    #Find two sister tips
    #preterminal nodes
    sisters_i <- which(treeCV == max(treeCV - diag(diag(treeCV))), arr.ind = T)[1,]
    sisters <- rownames(treeCV)[sisters_i]
    
    #calculate contrast between two sister nodes
    contrast[i,] <- traitTransforms[sisters[1],] - traitTransforms[sisters[2],] 
    
    #calculate BL along contrast
    sister_bl <- c(treeCV[sisters_i[1], sisters_i[1]], treeCV[sisters_i[2], sisters_i[2]]) - treeCV[sisters_i[1], sisters_i[2]]
    BLength[i] <- sum(sister_bl)
    
    #calculate multivariate normal density of contrast along Branch Length
    
    if(dim(treeCV)[1] > 2){   
      tipName = paste(sisters, collapse = "+")
      tt_index <- which(is.na(rownames(traitTransforms)))[1]
      rownames(traitTransforms)[tt_index] <- tipName
      #Fix Trait Matrix
      nodeValues <- (traitTransforms[sisters[1],] * sister_bl[2] + 
                       traitTransforms[sisters[2],] * sister_bl[1]) / BLength[i]
      
      traitTransforms[tt_index,] <- nodeValues
      
      #Fix Tree
      treeCV[sisters_i[1], sisters_i[1]] <- treeCV[sisters_i[1], sisters_i[1]] - sister_bl[1] + sister_bl[1] * sister_bl[2] / BLength[i] 
      rownames(treeCV)[sisters_i[1]] <- colnames(treeCV)[sisters_i[1]] <- tipName
      treeCV <- treeCV[-sisters_i[2], -sisters_i[2]]
      
    }
  }
  return(list(contrastTransformationMatrix = contrast, contrastBranchLengths = BLength, nodalTraitConditionals = traitTransforms))
}

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

BMpruneLLuvtChol <- function(rawTraits, sig, tree, printContrast = F, printLL = F, printBL = F){ #this function evaluates univariate normal densities
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
    if(printContrast){print(contrast)}
    
    #calculate BL separating two sister nodes
    BLength <- tree$edge.length[tree$edge[,2] == sisters[1]] + tree$edge.length[tree$edge[,2] == sisters[2]]
    if(printBL){print(BLength)}
    
    #calculate multivariate normal density of contrast along Branch Length
    LLs[i] <- sum(dnorm(contrast / BLength^0.5, mean = 0, sd = 1, log = T)) - log(detcov*(BLength^dim))/2
    if(printLL){print(LLs[i])}
    
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
