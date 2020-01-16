#inference under the phylogenetics multivariate brownian threshold model

#load compute libraries
library(gpuR)
library(Rcpp)

#load phylogenetic libraries
library(mvMORPH)
library(TESS)

tree <- tess.sim.taxa(n = 1, nTaxa = 8, lambda = 1, mu = 0, max = 1E3)[[1]]
R <- diag(c(1.5,2.3,1.7)) %*% matrix(c(1,0.5,0.6,0.5,1,0.4,0.6,0.4,1), 3, 3) %*% diag(c(1.5,2.3,1.7))
root_state <- rep(0, dim(R)[1])
traits <- mvSIM(tree, model="BM1", param=list(sigma=R, theta=root_state)) 


prunes <- prune(tree)
contrasts <- prunes[[1]] %*% traits[rownames(prunes[[3]])[1:length(tree$tip.label)],]

sum(sapply(1:(length(tree$tip.label)-1), function(x) dmvnorm(contrasts[x,], sigma = R*prunes[[2]][x], log = T)))
BMpruneLL(traits = traits, sig = R, tree = tree)





tree <- pbtree(n = 8); plot(tree); axisPhylo()
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
                     edge.length=(prod(tree$edge.length[tree$edge[,2] == sisters[1]], tree$edge.length[tree$edge[,2] == sisters[2]])/BLength) - tree$edge.length[tree$edge[,2] == sisters[1]],
                     Nnode=1)
      class(newTip)<-"phylo"
      oldw <- getOption("warn")
      options(warn = -1)
      tree <- bind.tree(x = tree, y = newTip, where = sisters[1])
      options(warn = oldw)
      tree <- drop.tip(tree, sisters[2])
    }
  }
  return(list(contrast, BLength, traitTransforms))
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
    print(contrast)
    
    #calculate BL separating two sister nodes
    BLength <- tree$edge.length[tree$edge[,2] == sisters[1]] + tree$edge.length[tree$edge[,2] == sisters[2]]
    print(paste0(BLength))
    
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
