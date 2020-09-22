library(MCMCpack)
library(microbenchmark)
library(phytools)
library(phangorn)
library(mvMORPH)
library(tictoc)

meansNexus <- function (characterMatrix, outputFilePath, traitsOnColumns = T){
  if (traitsOnColumns) {ntax = dim(characterMatrix)[1]; nchar = dim(characterMatrix)[2]}
  else {ntax = dim(characterMatrix)[2]; nchar = dim(characterMatrix)[1]}
  sink(outputFilePath, append=F)
  cat("#NEXUS\n\nBegin data;\n\tDimensions ntax=", ntax, "nchar=", nchar, ";\n\tFormat datatype=Continuous missing=?;\n\tMatrix\n\n")
  for(i in 1:length(characterMatrix[,1])) {cat(rownames(characterMatrix)[i], characterMatrix[i,], "\n")}
  cat(";\nEnd;")
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

betaMove <- function(simplex, index, conc, offset = 1){
  
  proposed_value <- rbeta(1, shape1 = simplex[index] * length(simplex) * conc + offset, shape2 = (1 - simplex[index]) * length(simplex) * conc + offset)
  proposed_simplex <- simplex
  proposed_simplex[index] <- proposed_value
  proposed_simplex[-index] <- proposed_simplex[-index] / sum(proposed_simplex[-index]) * (1-proposed_value) 
  
  forward_prob <- dbeta(proposed_value, shape1 = simplex[index] * length(simplex) * conc + offset, shape2 = (1 - simplex[index]) * length(simplex) * conc + offset, log = T)
  backward_prob <- dbeta(simplex[index], shape1 = proposed_value * length(simplex) * conc + offset, shape2 = (1 - proposed_value) * length(simplex) * conc + offset, log = T)
  
  return(list(prop = proposed_simplex, log_prop_ratio = backward_prob - forward_prob))
  
}

slideMove <- function(value, window){
  
  prop <- runif(1, min = value - window/2, max = value + window/2)
  return(list(prop = prop, log_prop_ratio = 0))
  
}

scaleMove <- function(value, window){
  sign_value <- sign(value)
  prop <- log(abs(value)) + runif(1, min = -window/2, max = window/2)
  return(list(prop = exp(prop) * sign_value, log_prop_ratio = 0))
   
}

prunePCV <- function(treeCV){ 
  ntips <- dim(treeCV)[1]
  contrast <- matrix(0, nrow = ntips-1, ncol = ntips); colnames(contrast) <- colnames(treeCV)
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

updateLi <- function(curr_Li, sd_index, change_ratio){
  if(length(sd_index) == 1){
    curr_Li[,sd_index] <- curr_Li[,sd_index] / change_ratio
  } else if(length(sd_index) == nrow(curr_Li)){
    curr_Li <- curr_Li %*% diag(sqrt(change_ratio))
  }
  return(curr_Li)
}

updateTransformedTraits <- function(means, Li){
  return(t(Li %*% t(means)))
}

updateTransformedContrasts_rate <- function(means, Li, prunes){
  transformed_contrasts <- prunes[[1]] %*% updateTransformedTraits(means, Li)[rownames(prunes[[3]])[1:nTaxa],]
  return(transformed_contrasts)
}

updateTransformedContrasts_tree <- function(transformed_traits, prunes){
  transformed_contrasts <- prunes[[1]] %*% transformed_traits[rownames(prunes[[3]])[1:nTaxa],]
  return(transformed_contrasts)
}

updateDetcov <- function(detcov, change_ratio){
  return(detcov * prod(1 / change_ratio))
}

nTaxa <- 6
tree <- pbtree(n = nTaxa, nsim = 1)
treeCV <- vcv.phylo(tree)
prunes <- prunePCV(treeCV)
n_traits <- 200
R <- rethinking::rlkjcorr(1, n_traits, 1)
sds <- diag(runif(n = n_traits, 1, 10))
sds <- diag(rep(2, n_traits))
sds2 <- sds; sds2[3,3] <- 5
Rm <- (sds) %*% R %*% (sds)
Rm2 <- (sds2) %*% R %*% (sds2)

upC <- chol(Rm)
lowC <- t(upC)

# solve(t(chol((sds2) %*% R %*% (sds2)))) - updateLi(Li, 3, 2.5)

Li <- solve(lowC)
detcov <- det(Rm)

updateDetcov(det(Rm), 2.5)
det(Rm2)

means <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = n_traits, sigma=Rm, mu= rep(0,n_traits)))

BMpruneLLuvtChol(means, Rm, tree)

transformed_traits <- t(Li %*% t(means)); rownames(transformed_traits) <- rownames(means)

t(Li %*% t(means))

transformed_contrasts <- prunes[[1]] %*% transformed_traits[rownames(prunes[[3]])[1:nTaxa],]
sum(sapply(1:nrow(transformed_contrasts), function(x) sum(dnorm(transformed_contrasts[x,] / prunes[[2]][x]^0.5, mean = 0, sd = 1, log = T)) - log(detcov*(prunes[[2]][x]^n_traits))/2))

BM_LL_PRECOMPUTE <- function(transformed_contrasts, prunes, detcov){
  n_traits <- ncol(transformed_contrasts)
  sum(sapply(1:nrow(transformed_contrasts), function(x) sum(dnorm(transformed_contrasts[x,] / prunes[[2]][x]^0.5, mean = 0, sd = 1, log = T)) - log(detcov*(prunes[[2]][x]^n_traits))/2))
}

BM_LL_PRECOMPUTE(transformed_contrasts, prunes, detcov)

#check for insensitivity to rooting
tree <- rtree(nTaxa, rooted = T)
BM_LL_PRECOMPUTE(updateTransformedContrasts_tree(transformed_traits = transformed_traits, prunes = prunePCV(vcv.phylo(tree))), prunePCV(vcv.phylo(tree)), detcov)
BM_LL_PRECOMPUTE(updateTransformedContrasts_tree(transformed_traits = transformed_traits, prunes = prunePCV(vcv.phylo(midpoint(tree)))), prunePCV(vcv.phylo(midpoint(tree))), detcov)
BM_LL_PRECOMPUTE(updateTransformedContrasts_tree(transformed_traits = transformed_traits, prunes = prunePCV(vcv.phylo(unroot(tree)))), prunePCV(vcv.phylo(unroot(tree))), detcov)


BMpruneLLuvtChol(means, Rm, tree)

############
### MCMC ###
############

#simulate data
nTaxa <- 5
n_traits <- 300
tree <- pbtree(n = nTaxa, nsim = 1)
R <- rethinking::rlkjcorr(1, n_traits, 1)
sds <- diag(runif(n = n_traits, 1, 5))
Rm <- (sds) %*% R %*% (sds)
means <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = n_traits, sigma=Rm, mu= rep(0,n_traits)))

#write test data to file
# meansNexus(means, outputFilePath = "~/data/test_rMCMC.nex")
# ratesTSV(R, "~/data/test_rMCMC.tsv")

#specify priors
trait_rates_conc_params <- rep(1, n_traits)
bls_conc_params <- rep(1, 2*nTaxa-3)
tree_prior <- "uniform"
TL_log_prior <- c(1,1)

#initialize chain
curr_tree <- rtree(n = nTaxa, rooted = F)
curr_relative_rates <- c(rdirichlet(1, alpha = trait_rates_conc_params))
curr_rates <- curr_relative_rates * n_traits
curr_Rm <- diag(sqrt(curr_rates)) %*% R %*% diag(sqrt(curr_rates))
curr_detcov <- det(curr_Rm)
curr_Li <- solve(t(chol(curr_Rm)))
curr_TL_log<- rnorm(1,TL_log_prior[1], TL_log_prior[2]) 
curr_TL <- 10^curr_TL_log
curr_relative_bls <- rdirichlet(1, alpha = bls_conc_params)
curr_bls <- curr_TL * curr_relative_bls
curr_tree$edge.length <- curr_bls
curr_transformed_traits <- t(curr_Li %*% t(means)); rownames(curr_transformed_traits) <- rownames(means)
curr_prunes <- prunePCV(vcv.phylo(curr_tree))
curr_transformed_contrasts <- updateTransformedContrasts_tree(curr_transformed_traits, prunes = curr_prunes)
curr_ll <- BM_LL_PRECOMPUTE(transformed_contrasts = curr_transformed_contrasts, curr_prunes, detcov = curr_detcov)
curr_ll
BMpruneLLuvtChol(means, sig = diag(sqrt(curr_rates)) %*% R %*% diag(sqrt(curr_rates)), tree = curr_tree)

#specify chain params
n_iter <- 2E6
thin <- 1E3
n_out <- round(n_iter/thin)
burnin_prop = 0.4
proposal_probs <- c(0.2, 0.3,0.2,0.3) #c("topology", "rates", "TL", "BLs")

#initialize chain record
trees <- rmtree(N = n_out, n = nTaxa); trees[[1]] <- curr_tree
rates <- matrix(0, nrow = n_out, ncol = n_traits); rates[1,] <- curr_rates
lls <- rep(0, n_out); lls[1] <- curr_ll
n_prop <- rep(0, 4)
n_accept <- rep(0, 4) 

#initialize floating proposals
prop_tree <- curr_tree
prop_relative_rates <- curr_relative_rates
prop_rates <- curr_rates
prop_Rm <- curr_Rm
prop_detcov <- curr_detcov
prop_Li <- curr_Li
prop_TL_log <- curr_TL_log
prop_TL <- curr_TL
prop_relative_bls <- curr_relative_bls
prop_bls <- curr_bls
prop_transformed_traits <- curr_transformed_traits
prop_prunes <- curr_prunes
prop_transformed_contrasts <- curr_transformed_contrasts
prop_ll <- curr_ll


debug_bugs <- T
tic()
for(i in 2:n_iter){
  
  #progress bar
  if(i %% (n_iter / 100) == 0){
    cat(paste0(i / n_iter * 100, "%  : \n"))
    toc()
    tic()
    }
  
  #propose a move and calculate prior density and 
  type_of_move <- sample(c("topology", "rates", "TL", "BLs"), size = 1, prob = proposal_probs)
  
  if(type_of_move == "topology"){
    n_prop[1] <- n_prop[1] + 1
    type_of_topology_move <- sample(c("NNI", "SPR"), 1)
    type_of_topology_move <- "NNI"
    if(type_of_topology_move == "NNI"){
      prop_tree <- rNNI(curr_tree)
    }
    if(type_of_topology_move == "SPR"){
      prop_tree <- rSPR(curr_tree)
    }
    prop_prunes <- prunePCV(vcv.phylo(prop_tree))
    prop_transformed_contrasts <- updateTransformedContrasts_tree(transformed_traits = curr_transformed_traits, prunes = prop_prunes)
    prop_detcov <- curr_detcov
    log_prop_ratio <- 0
    log_prior_ratio <- 0
  }
  
  if(type_of_move == "BLs"){
    n_prop[2] <- n_prop[2] + 1
    prop_tree <- curr_tree
    which_BL <- sample(1:length(curr_tree$edge.length), 1)
    curr_bls <- prop_tree$edge.length
    curr_relative_bls <- curr_bls / curr_TL
    prop_bl_info <- betaMove(curr_relative_bls, index = which_BL, conc = 5)
    prop_relative_bls <- prop_bl_info$prop
    prop_bls <- prop_relative_bls * curr_TL
    prop_tree$edge.length <- prop_bls
    prop_prunes <- prunePCV(vcv.phylo(prop_tree))
    prop_transformed_contrasts <- updateTransformedContrasts_tree(transformed_traits = curr_transformed_traits, prunes = prop_prunes)
    log_prop_ratio <- prop_bl_info$log_prop_ratio
    prop_detcov <- curr_detcov
    log_prior_ratio <- 0
  }
  
  if(type_of_move == "TL"){
    n_prop[3] <- n_prop[3] + 1
    prop_tree <- curr_tree
    curr_TL <- sum(curr_tree$edge.length)
    curr_TL_log <- log10(curr_TL)
    prop_TL_log <- slideMove(curr_TL_log, window = 0.5)$prop
    prop_TL <- 10^prop_TL_log
    prop_tree$edge.length <- prop_tree$edge.length / curr_TL * prop_TL
    prop_prunes <- curr_prunes 
    prop_prunes$contrastBranchLengths <- prop_prunes$contrastBranchLengths / curr_TL * prop_TL
    prop_transformed_contrasts <- curr_transformed_contrasts
    log_prop_ratio <- 0
    log_prior_ratio <- dnorm(prop_TL_log, mean = TL_log_prior[1], TL_log_prior[2], log = T) - dnorm(curr_TL_log, mean = TL_log_prior[1], TL_log_prior[2], log = T) 
    prop_detcov <- curr_detcov
  }
  
  if(type_of_move == "rates"){
    n_prop[4] <- n_prop[4] + 1
    which_Rate <- sample(1:n_traits, 1)
    curr_relative_rates <- curr_rates / n_traits
    prop_relative_rates_info <- betaMove(curr_relative_rates, index = which_Rate, conc = 5)
    log_prop_ratio <- prop_relative_rates_info$log_prop_ratio
    log_prior_ratio <- 0
    prop_rates <- prop_relative_rates_info$prop * n_traits
    
    prop_change_ratio <- curr_rates / prop_rates
    prop_Li <- updateLi(curr_Li = curr_Li, sd_index = 1:n_traits, change_ratio = prop_change_ratio)
    prop_detcov <- updateDetcov(curr_detcov, change_ratio = prop_change_ratio)
    
    prop_transformed_traits <- updateTransformedTraits(means, prop_Li)
    prop_transformed_contrasts <- updateTransformedContrasts_tree(transformed_traits = prop_transformed_traits, prunes = curr_prunes)
    prop_prunes <- curr_prunes
    log_prior_ratio <- 0
    
  }
  
  #compute likelihood ratio
  prop_ll  <- BM_LL_PRECOMPUTE(transformed_contrasts = prop_transformed_contrasts, prop_prunes, detcov = prop_detcov)
  log_ll_ratio <- prop_ll - curr_ll
  
  #accept or reject
  log_accept_prob <- log_prop_ratio + log_ll_ratio + log_prior_ratio
  
  if(debug_bugs){
    # log_accept_prob <- -Inf
    print(type_of_move)
    print(c(BM_LL_PRECOMPUTE(transformed_contrasts = prop_transformed_contrasts, prop_prunes, detcov = prop_detcov),
            BMpruneLLuvtChol(means, sig = diag(sqrt(prop_rates)) %*% R %*% diag(sqrt(prop_rates)), tree = prop_tree)))
  }
  
  if(log(runif(1, 0, 1)) < log_accept_prob){
    n_accept[which(type_of_move == c("topology", "rates", "TL", "BLs"))] <- n_accept[which(type_of_move == c("topology", "rates", "TL", "BLs"))] + 1
    curr_tree <- prop_tree
    curr_relative_rates <- prop_relative_rates
    curr_rates <- prop_rates
    curr_Rm <- prop_Rm
    curr_detcov <- prop_detcov
    curr_Li <- prop_Li
    curr_TL_log <- prop_TL_log
    curr_TL <- prop_TL
    curr_relative_bls <- prop_relative_bls
    curr_bls <- prop_bls
    curr_transformed_traits <- prop_transformed_traits
    curr_prunes <- prop_prunes
    curr_transformed_contrasts <- prop_transformed_contrasts
    curr_ll <- prop_ll
  } else {
    prop_tree <- curr_tree
    prop_relative_rates <- curr_relative_rates
    prop_rates <- curr_rates
    prop_Rm <- curr_Rm
    prop_detcov <- curr_detcov
    prop_Li <- curr_Li
    prop_TL_log <- curr_TL_log
    prop_TL <- curr_TL
    prop_relative_bls <- curr_relative_bls
    prop_bls <- curr_bls
    prop_transformed_traits <- curr_transformed_traits
    prop_prunes <- curr_prunes
    prop_transformed_contrasts <- curr_transformed_contrasts
    prop_ll <- curr_ll
  }
  
  #record given thinning 
  if(i %% thin == 0){
    trees[[i/thin]] <- curr_tree
    rates[i/thin,] <- curr_rates
    lls[i/thin] <- curr_ll
    
  }
      
}

#discard burnin
trees <- trees[(burnin_prop*n_out):n_out]
rates <- rates[(burnin_prop*n_out):n_out,]
lls <- lls[(burnin_prop*n_out):n_out]

n_accept / n_prop 

plot(diag(sds^2), apply(rates, 2, mean))
effectiveSize(rates)
effectiveSize(lls)
plot(sapply(1:length(trees), function(x) sum(trees[[x]]$edge.length)), type = "l")
hist(sapply(1:length(trees), function(x) sum(trees[[x]]$edge.length)))
mean(sapply(1:length(trees), function(x) sum(trees[[x]]$edge.length)))
apply(rates, 2, mean)
plot(maxCladeCred(trees))
plot(tree)
RF.dist(maxCladeCred(trees), tree, normalize = T)

# rb_trees <- read.tree("~/output/test_rMCMC.trees")
# 
# startsWith2 <- function(x, pre){apply(sapply(1:length(pre), function(n) startsWith(x, pre[n])), 1, any)}
# greedyCT <- function(trees){
#   tipNames <- trees[[1]]$tip.label
#   tipNums <- paste0("t", 1:length(tipNames))
#   for(i in 1:length(trees)){
#     trees[[i]]$tip.label <- tipNums[match(trees[[i]]$tip.label, tipNames)]
#   }
#   ape::write.tree(trees, file = "gct_trees.txt")
#   writeLines("gct_trees.txt\nY", con = "consense_script.txt")
#   if(file.exists("outfile")){file.remove("outfile")}
#   if(file.exists("outtree")){file.remove("outtree")}
#   system(paste0("/Applications/phylip-3.695/exe/consense < consense_script.txt"), intern = T)
#   gct <- ape::read.tree("outtree")
#   if(file.exists("outfile")){file.remove("outfile")}
#   if(file.exists("outtree")){file.remove("outtree")}
#   if(file.exists("gct_trees.txt")){file.remove("gct_trees.txt")}
#   if(file.exists("consense_script.txt")){file.remove("consense_script.txt")}
#   gct$tip.label <- tipNames[match(gct$tip.label, tipNums)]
#   return(gct)
# }
# convNum2Str <- function(nums, key){
#   sapply(1:length(nums), function(x) key[nums[x]])
# }
# prop.part.df <- function(trees, cutoff = 0.01, bs = T){
#   if(class(trees) == "multiPhylo"){
#     if(trees[[1]]$Nnode == (length(trees[[1]]$tip.label) - 1)){
#       trees <- unroot(trees) #unroot rooted trees
#     }
#     tipLabs <- trees[[1]]$tip.label
#     numTrees <- length(trees)
#     if(bs) {
#       out <- as.prop.part(bitsplits(trees))
#     } else {
#       out <- prop.part(trees)
#     }
#     outList <- as.list.data.frame(out)
#     pps <- attributes(outList)$number/numTrees
#     props <- data.frame(pps, as.matrix(outList)); colnames(props) <- c("postProbs", "clade")
#     props <- props[order(-pps),]
#     props <- props[props[,1] > cutoff,]
#     rownames(props) <- 1:nrow(props)
#     props$cladeNames <- sapply(1:length(props[,1]), function(x) sort(convNum2Str(props$clade[[x]], attributes(outList)$labels)))
#     props <- props[,c(1,3)]
#     if(!bs) {
#       props <- props[-1,]
#     }
#     allClades <- c(props$cladeNames, lapply(1:length(props$postProbs), function(x) sort(setdiff(tipLabs, props$cladeNames[[x]]))))
#     if(any(duplicated(allClades))){
#       # print("duplicate clades found")
#       dupClades <- allClades[duplicated(allClades)]
#       dupCladesComp <- lapply(1:length(dupClades), function(x) sort(setdiff(tipLabs, dupClades[[x]])))
#       matchedDupes <- cbind(1:length(dupClades), sapply(1:length(dupClades), function(x) which(sapply(1:length(dupCladesComp), function(y) setequal(dupClades[[x]], dupCladesComp[[y]])))))
#       dupes <- matrix(data = 0, nrow = length(matchedDupes[,1])/2, ncol = 2)
#       for(i in 1:length(matchedDupes[,1])){
#         pair <- sort(matchedDupes[i,])
#         if(!any(sapply(1:length(dupes[,1]), function(x) pair == dupes[x,]))){
#           dupes[min(which(apply(dupes, 1, sum) == 0)),] <- pair
#         }
#       }
#       for(i in 1:length(dupes[,1])){
#         clade1 <- dupClades[dupes[i,1]][[1]]
#         clade2 <- dupClades[dupes[i,2]][[1]]
#         propsInd1 <- which(sapply(1:length(props[,1]), function(x) setequal(clade1, props$cladeNames[[x]])))
#         propsInd2 <- which(sapply(1:length(props[,1]), function(x) setequal(clade2, props$cladeNames[[x]])))
#         props[propsInd1,1] <- props[propsInd1,1] + props[propsInd2,1]
#         props <- props[-propsInd2,]
#       }
#       props <- props[order(props[,1], decreasing = T),]
#     }
#     props
#   } else if (class(trees) == "phylo"){
#     if(trees$Nnode == (length(trees$tip.label) - 1)){
#       trees <- unroot(trees) #unroot rooted trees
#     }
#     tipLabs <- trees$tip.label
#     numTrees <- 1
#     if(bs) {
#       out <- as.prop.part(bitsplits(trees))
#     } else {
#       out <- prop.part(trees)
#     }
#     outList <- as.list.data.frame(out)
#     pps <- attributes(outList)$number/numTrees
#     props <- data.frame(pps, as.matrix(outList)); colnames(props) <- c("postProbs", "clade")
#     props <- props[order(-pps),]
#     props <- props[props[,1] > cutoff,]
#     rownames(props) <- 1:nrow(props)
#     props$cladeNames <- sapply(1:length(props[,1]), function(x) sort(convNum2Str(props$clade[[x]], attributes(outList)$labels)))
#     props <- props[,c(1,3)]
#     if(!bs) {
#       props <- props[-1,]
#     }
#     props
#   }
# }
# compareTrees <- function(biparts1, biparts2, tipLabs, returnCladeNames = F){
#   matchProbs <- matrix(0, nrow = sum(length(biparts1[,1]), length(biparts2[,1])), ncol = ifelse(returnCladeNames, 3, 2))
#   counter <- 1
#   biparts2$notSeen <- 1
#   for(clade in 1:length(biparts1[,2])){
#     cladeName <- biparts1[,2][[clade]]
#     isThere <- sapply(1:length(biparts2[,1]), function(x) identical(cladeName, biparts2[x,2][[1]]))
#     altCladeName <- sort(setdiff(tipLabs, cladeName))
#     shortCladeName <- paste0(sort(shorter(cladeName, altCladeName)), collapse = ", ")
#     orIsThere <- sapply(1:length(biparts2[,1]), function(x) identical(altCladeName, biparts2[x,2][[1]]))
#     if(any(isThere)){
#       biparts2$notSeen[isThere] <- 0
#       if(returnCladeNames){
#         matchProbs[counter,] <- c(biparts1[,1][clade], biparts2[,1][isThere], shortCladeName)
#       } else {
#         matchProbs[counter,] <- c(biparts1[,1][clade], biparts2[,1][isThere])
#       }
#     } else if (any(orIsThere)) {
#       biparts2$notSeen[orIsThere] <- 0
#       if(returnCladeNames){
#         matchProbs[counter,] <- c(biparts1[,1][clade], biparts2[,1][orIsThere], shortCladeName)
#       } else {
#         matchProbs[counter,] <- c(biparts1[,1][clade], biparts2[,1][orIsThere])
#       }
#     } else {
#       if(returnCladeNames){
#         matchProbs[counter,] <- c(biparts1[,1][[clade]], 0, shortCladeName)
#       } else {
#         matchProbs[counter,] <- c(biparts1[,1][[clade]], 0)
#       }
#     }
#     counter <- counter + 1
#   }
#   if(sum(biparts2$notSeen) > 0){
#     if(returnCladeNames){
#       cladeNames <- sapply(1:length(biparts2[biparts2$notSeen == 1,2]), function(bip) 
#         paste0(sort(shorter(biparts2[biparts2$notSeen == 1,2][[bip]], sort(setdiff(tipLabs, biparts2[biparts2$notSeen == 1,2][[bip]])))), collapse = ", "))
#       matchProbs[counter:(counter+sum(biparts2$notSeen)-1),] <- cbind(0, cbind(biparts2[biparts2$notSeen == 1,1], cladeNames)) #needs to be an sapply
#     } else {
#       matchProbs[counter:(counter+sum(biparts2$notSeen)-1),] <- cbind(0, biparts2[biparts2$notSeen == 1,1])
#     }
#   }
#   if(returnCladeNames){
#     matchProbs <- matchProbs[rowSums(matchProbs == c(0,0)) != ifelse(returnCladeNames,3,2),]
#   }
#   matchProbs
# }
# clade_prob <- function(tipNames, trees, allTipNames = NA, partfreqs = NA){
#   if(is.na(partfreqs)[1]){
#     clade_monophyly <- prop.part.df(trees, cutoff = 0.00001)
#     if(is.na(allTipNames)[1]){
#       allTipNames <- trees[[1]]$tip.label
#     }
#   } else {
#     clade_monophyly <- partfreqs
#     if(is.na(allTipNames)[1]){
#       allTipNames <- unique(unlist(bipart_probs$cladeNames))
#     }
#   }
#   cladeName <- tipNames
#   clade <- list(allTipNames[startsWith2(allTipNames, cladeName)], allTipNames[!startsWith2(allTipNames, cladeName)])
#   clade_ind <- (sapply(1:length(clade_monophyly$cladeNames), function(x) setequal(clade_monophyly$cladeNames[[x]], clade[[1]]) |
#                          setequal(clade_monophyly$cladeNames[[x]], clade[[2]])))
#   if(all(!clade_ind)){return(0)}
#   # clade_monophyly$cladeNames[clade_ind]
#   return(clade_monophyly$postProbs[clade_ind])
# }
# internal_nodes_probs <- function(tree, trees_to_use){
#   bipart_probs <- prop.part.df(trees_to_use, cutoff = 0.00001)
#   ti <- 1:length(tree$tip.label)
#   tips <- sapply(1:Nnode(tree), function(node) as.numeric(na.omit(ti[getDescendants(tree, length(tree$tip.label) + node)])))
#   tips <- sapply(1:length(tips), function(node) tree$tip.label[tips[[node]]])
#   node_pps <- sapply(1:length(tips), function(node) clade_prob(tipNames = tips[[node]], allTipNames = tree$tip.label, partfreqs = bipart_probs))
#   return(node_pps)
# }
# shorter <- function(x1, x2){
#   if(length(x1) > length(x2)){
#     return(x2)
#   } else{
#     return(x1)
#   }
# }
# plot(compareTrees(prop.part.df(trees), prop.part.df(rb_trees), tipLabs = tree$tip.label)); abline(0,1)
