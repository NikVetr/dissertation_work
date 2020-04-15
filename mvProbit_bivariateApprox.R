setwd("/Volumes/1TB/Bailey/")
library(pbivnorm)
library(Matrix)
library(mvtnorm)
library(phangorn)
library(parallel)
library(doParallel)
library(foreach)
library(data.table)
library(readr)
library(microbenchmark)
library(TESS)
library(mvMORPH)
library(questionr)
library(tictoc)
?pbivnorm

rlkj <- function (K, eta = 1) {
  alpha <- eta + (K - 2)/2
  r12 <- 2 * rbeta(1, alpha, alpha) - 1
  R <- matrix(0, K, K)
  R[1, 1] <- 1
  R[1, 2] <- r12
  R[2, 2] <- sqrt(1 - r12^2)
  if (K > 2) 
    for (m in 2:(K - 1)) {
      alpha <- alpha - 0.5
      y <- rbeta(1, m/2, alpha)
      z <- rnorm(m, 0, 1)
      z <- z/sqrt(crossprod(z)[1])
      R[1:m, m + 1] <- sqrt(y) * z
      R[m + 1, m + 1] <- sqrt(1 - y)
    }
  return(crossprod(R))
}

multi_pmvnorm <- function(mean, sigma, binaries, algorithm = c("GenzBretz", "Miwa")[1]){
  upper <- binaries
  upper[] <- 0
  lower <- upper
  upper[binaries == 1] <- Inf
  lower[binaries == 0] <- -Inf
  return(sapply(1:dim(binaries)[1], function(indiv) pmvnorm(lower = lower[indiv,], upper = upper[indiv,], mean = mean, sigma = sigma, algorithm = algorithm)[1]))
}

multi_pmvnorm_optim <- function(mean, sigma, binaries, algorithm = c("GenzBretz", "Miwa")[1]){
  upper <- binaries
  upper[] <- 0
  lower <- upper
  upper[binaries == 1] <- Inf
  lower[binaries == 0] <- -Inf
  return(sapply(1:dim(binaries)[1], function(indiv) pmvnorm(lower = as.numeric(lower[indiv,]), 
                                                            upper = as.numeric(upper[indiv,]), 
                                                            mean = mean, sigma = sigma, algorithm = algorithm)[1]))
}

# multi_pbivnorm_optim <- function(mean, sigma, binaries){
#   if(dim(sigma)[1] != 2){
#     stop("not a bivariate normal -- did you mean to use multi_pmbnorm_optim()?")
#   }
#   rho <- sigma[1,2]
#   upper <- binaries
#   upper[] <- -mean
#   lower <- upper <- t(upper)
#   upper[binaries == 1] <- Inf
#   lower[binaries == 0] <- -Inf
#   
#   return(sapply(1:dim(binaries)[1], function(indiv) pbivnorm(x = t(cbind(as.numeric(upper[indiv,]), 
#                                                                          (as.numeric(lower[indiv,])))), rho = rho)[1]))
# }

uniqueSP <- function(indivs){
  SPs <- table(sapply(1:length(indivs[,1]), function(indiv) paste0(indivs[indiv,], collapse = "")))
  counts <- as.numeric(SPs)
  SPs <- t(sapply(1:length(SPs), function(SP) as.numeric(unlist(strsplit(names(SPs)[SP], "")))))
  return(list(SPs = SPs, counts = counts))
}

#first let's resimulate data
nTaxa <- 5
d_traits <- 5
tree <- tess.sim.taxa(n = 1, nTaxa = nTaxa, lambda = 1, mu = 0, max = 1E3)[[1]]
target_tree_height <- 2
tree$edge.length <- tree$edge.length * target_tree_height / node.depth.edgelength(tree)[1]

R <- rlkj(d_traits); rownames(R) <- colnames(R) <- paste0("tr", 1:dim(R)[1])
root_state <- rep(0, dim(R)[1])
root_state <- rnorm(d_traits)
traits <- mvSIM(tree, model="BM1", param=list(sigma=R, theta=root_state)) 

#simulate thresholded traits at the population level
n_indiv <- 5000
threshold <- rep(0, dim(traits)[2])
traits_indiv <- lapply(1:dim(traits)[1], function(tip) rmvnorm(n = n_indiv, mean = traits[tip,], sigma = R)); names(traits_indiv) <- rownames(traits)
traits_indiv_discr <- lapply(1:dim(traits)[1], function(tip) t(sapply(1:n_indiv, function(indiv) as.numeric(traits_indiv[[tip]][indiv,] > threshold)))); names(traits_indiv_discr) <- names(traits_indiv) 

#identify unique site patterns
traits_indiv_discr_unique <- lapply(1:dim(traits)[1], function(tip) uniqueSP(traits_indiv_discr[[tip]]))

#find probs of these site patterns and multiply by freq
fast = TRUE
if(!fast){
  logLLs <- sapply(1:length(traits_indiv_discr_unique), function(tip) sum(log(multi_pmvnorm(mean = traits[tip,], sigma = R, binaries = traits_indiv_discr_unique[[tip]]$SPs)) * traits_indiv_discr_unique[[tip]]$counts))
}

#identify *bivariate* unique site patterns
tr1 <- 1
tr2 <- 2
logLLs <- sapply(1:nTaxa, function(tip) sum(log(
  multi_pmvnorm(mean = traits[tip,c(tr1,tr2)], sigma = R[c(tr1,tr2),c(tr1,tr2)], binaries = traits_indiv_discr_unique[[tip]]$SPs)) * traits_indiv_discr_unique[[tip]]$counts))

#now let's try to optimize

logLL_bivProb <- function(means, cor, traits_indiv_discr_unique){
  sum(sapply(1:dim(means)[1], function(tip) sum(log(
    multi_pmvnorm(mean = means[tip,], sigma = matrix(c(1,cor,cor,1),2,2), 
                  binaries = traits_indiv_discr_unique[[tip]]$SPs)) * traits_indiv_discr_unique[[tip]]$counts)))
}

dat = data.frame(tr1 = unlist(lapply(1:nTaxa, function(tip) traits_indiv_discr_unique[[tip]]$SPs[,tr1])),
                 tr2 = unlist(lapply(1:nTaxa, function(tip) traits_indiv_discr_unique[[tip]]$SPs[,tr2])),
                 nSP = unlist(lapply(1:nTaxa, function(tip) traits_indiv_discr_unique[[tip]]$counts)),
                 tip = unlist(lapply(1:nTaxa, function(tip) rep(tip, length(traits_indiv_discr_unique[[tip]]$counts)))))


par <- c(rep(c(0,0), nTaxa), 30, 0)

logLL_bivProb_optim <- function(par, dat){
  howbig <- length(par)
  rho <- par[howbig]
  reg <- par[howbig - 1]
  par <- matrix(par[-c(howbig-1, howbig)], ncol = 2, byrow = T)
  logLL <- -sum(sapply(1:dim(par)[1], function(tip) sum(log(multi_pmvnorm_optim(mean = par[tip,], sigma = matrix(c(1,rho,rho,1),2,2), 
                  binaries = dat[dat$tip == tip,c(1,2)])) * dat[dat$tip == tip,3]))) 
  # regTerm <- - sum(dnorm(x = as.vector(par), mean = 0, sd = sqrt(reg), log = T))
  regTerm <- -sum(sapply(1:dim(par)[1], function(tip) dmvnorm(x = par[tip,], matrix(c(1,rho,rho,1),2,2) * reg, mean = c(0,0), log = T)))
  returnVal <- logLL + regTerm
  if(returnVal == Inf){returnVal <- 1E10} #hacky temporary fix to make optim's generalized "L-BFGS-B" method stop complaining about singularities :/
  return(returnVal)
}

logLL_bivProb_optim(par, dat)

method <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent")[4]
t_s <- Sys.time()
optim_out <- optim(par = par, logLL_bivProb_optim, dat = dat, method = method, control = list(trace = 3, REPORT = 10),
                   upper = c(rep(Inf, length(par)-1), 1), lower = c(rep(-Inf, length(par)-2), 0, -1))
Sys.time() - t_s
# optim_out <- optim(par = par, logLL_bivProb_optim, dat = dat, method = method)
optim_est <- optim_out$par
plot(c(t(traits[,c(tr1,tr2)])), optim_est[-c(length(optim_est)-1, length(optim_est))]); 
abline(a = 0,b = 1, col = "green", lwd = 2)
paste0("true r = ", round(R[tr1, tr2],3), ", optim r = ", round(optim_est[length(optim_est)],3), collapse = "")
paste0("regularizing variance = ", round(optim_est[length(optim_est)-1], 3))

#goal -- alternate optimization under nearPD'd correlations and mean liabs
#can also initialize from random values and take global min of optim_out$value

#########################################################
#what if we try to optimize all the correlations jointly?
#########################################################

dat <- as.data.frame(do.call(rbind, sapply(1:nTaxa, function(tip) traits_indiv_discr_unique[[tip]]$SPs)))
colnames(dat) <- paste0("tr", 1:d_traits)
dat2 = data.frame(nSP = unlist(lapply(1:nTaxa, function(tip) traits_indiv_discr_unique[[tip]]$counts)),
                 tip = unlist(lapply(1:nTaxa, function(tip) rep(tip, length(traits_indiv_discr_unique[[tip]]$counts)))))
dat <- cbind(dat, dat2)

#bad init
par <- c(rep(0, nTaxa * d_traits), 10, rep(0, choose(d_traits, 2)))

#better init
init_cor <- cor(do.call(rbind, lapply(1:nTaxa, function(x) (traits_indiv_discr[[x]]))))
init_means <- t(sapply(1:nTaxa, function(tip) apply(X = (traits_indiv_discr[[tip]]), MARGIN = 2, FUN = base::mean)))
init_means <- qnorm(init_means)
nsamp_ptip <- sapply(1:length(traits_indiv_discr), function(tip) dim(traits_indiv_discr[[tip]])[1])
std_norm_5050_invariance <- sapply(1:length(nsamp_ptip), function(tip) qnorm(p = 1-0.5^(1/nsamp_ptip[tip]))) 
for(i in 1:nrow(init_means)){
  init_means[i, init_means[i,] == -Inf] <- std_norm_5050_invariance[i]
  init_means[i, init_means[i,] == Inf] <- -std_norm_5050_invariance[i]
}
plot(init_means, traits)
plot(init_cor, R)
par <- c(c(t(init_means)), 10, init_cor[upper.tri(init_cor)])


logLL_multivProb_optim_joint <- function(par, dat){
  hb <- length(par)
  nTaxa <- length(unique(dat$tip))
  d_trait <- ncol(dat) - 2
  cor <- diag(d_trait)
  cor[upper.tri(cor)] <- par[(nTaxa*d_trait+2) : hb];
  cor <- cor + t(cor) - diag(d_trait)
  reg <- par[nTaxa * d_trait + 1]
  means <- matrix(par[1:(nTaxa*d_trait)], ncol = d_trait, nrow = nTaxa, byrow = T)
  logLL <- -sum(sapply(1:dim(means)[1], function(tip) sum(log(multi_pmvnorm_optim(mean = means[tip,], sigma = cor, 
                                                                                binaries = dat[dat$tip == tip, c(1:d_trait)])) * dat[dat$tip == tip,d_trait+1]))) 
  # regTerm <- - sum(dnorm(x = as.vector(par), mean = 0, sd = sqrt(reg), log = T))
  regTerm <- -sum(sapply(1:dim(means)[1], function(tip) dmvnorm(x = means[tip,], cor * reg, mean = rep(0, d_trait), log = T)))
  returnVal <- logLL + regTerm
  if(returnVal == Inf){returnVal <- 1E10} #hacky temporary fix to make optim's generalized "L-BFGS-B" method stop complaining about singularities :/
  return(returnVal)
}

logLL_multivProb_optim_joint(par, dat)

logLL_bivProb_optim_joint <- function(par, dat){
  hb <- length(par)
  nTaxa <- length(unique(dat$tip))
  d_trait <- ncol(dat) - 2
  cor <- diag(d_trait)
  cor[upper.tri(cor)] <- par[(nTaxa*d_trait+2) : hb];
  cor <- cor + t(cor) - diag(d_trait)
  reg <- par[nTaxa * d_trait + 1]
  means <- matrix(par[1:(nTaxa*d_trait)], ncol = d_trait, nrow = nTaxa, byrow = T)
  logLL <- sum(sapply(1:(d_trait-1), function(rowi) sum(sapply((rowi+1):d_trait, 
                                 function(colj) -sum(sapply(1:nTaxa, function(tip) sum(log(multi_pmvnorm_optim(mean = means[tip,c(rowi,colj)], sigma = matrix(c(1,cor[rowi,colj],cor[rowi,colj],1),2,2), 
                                                                                  binaries = dat[dat$tip == tip, c(rowi, colj)])) * dat[dat$tip == tip,d_trait+1])))))))
  # regTerm <- - sum(dnorm(x = as.vector(par), mean = 0, sd = sqrt(reg), log = T))
  regTerm <- -sum(sapply(1:dim(means)[1], function(tip) dmvnorm(x = means[tip,], (as.matrix(nearPD(cor, corr = T)$mat)) * reg, mean = rep(0, d_trait), log = T)))
  returnVal <- logLL + regTerm
  if(returnVal == Inf){returnVal <- 1E10} #hacky temporary fix to make optim's generalized "L-BFGS-B" method stop complaining about singularities :/
  return(returnVal)
}

logLL_bivProb_optim_joint(par, dat)


microbenchmark(times = 10,
  logLL_multivProb_optim_joint(par, dat),
  logLL_bivProb_optim_joint(par, dat)
)

method <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent")[4]
t_s <- Sys.time()
par <- c(rep(0, nTaxa * d_traits), 30, rep(0, choose(d_traits, 2)))

optim_out <- optim(par = par, logLL_bivProb_optim_joint, dat = dat, method = method, control = list(trace = 3, REPORT = 1),
                   upper = c(rep(Inf, nTaxa * d_traits + 1), rep(1, choose(d_traits, 2))), 
                   lower = c(rep(-Inf, nTaxa * d_traits), 0, rep(-1, choose(d_traits, 2))))
Sys.time() - t_s
# optim_out <- optim(par = par, logLL_bivProb_optim, dat = dat, method = method)
optim_est <- optim_out$par
R_est <- diag(d_traits)
R_est[upper.tri(R_est)] <- optim_est[(nTaxa*d_traits+2) : length(optim_est)];
R_est <- R_est + t(R_est) - diag(d_traits)
reg <- optim_est[nTaxa * d_traits + 1]
paste0("regularizing variance = ", round(reg, 3))
est_means <- matrix(optim_est[1:(nTaxa*d_traits)], ncol = d_traits, nrow = nTaxa, byrow = T)
par(mfrow = c(1,2))
plot(est_means, traits); abline(0,1)
plot(R_est, R, xlim = c(-1,1), ylim = c(-1,1)); abline(0,1)

#TODO: use more efficient vectorized form of the Genz algorithm
#optimize just one pairwise correlation at a time, then average the n-1 pop means?
#optimize some subset of the pop means at a time, then average?
#iterate over all pairwise correlations, and then within each iterate over subsets of means, optimize subsets
#then condition on means and iterate over individual means at random while conditioning on the rest + correlation
#then condition on means and do one final pass over correlation
#Do 1-d optimization over every parameter of the par vector, keep going till they converge; Maybe spend more time on the correlations; Add conditioned on pars to dat list

##################################################################
#let's try iterating over one parameter of the par vector at a time
##################################################################

#first let's resimulate data
nTaxa <- 9
d_traits <- 150
n_indiv <- 35

tree <- tess.sim.taxa(n = 1, nTaxa = nTaxa, lambda = 1, mu = 0, max = 1E3)[[1]]
target_tree_height <- 2
tree$edge.length <- tree$edge.length * target_tree_height / node.depth.edgelength(tree)[1]

R <- rlkj(d_traits); rownames(R) <- colnames(R) <- paste0("tr", 1:dim(R)[1])
root_state <- rep(0, dim(R)[1])
root_state <- rnorm(d_traits)
traits <- mvSIM(tree, model="BM1", param=list(sigma=R, theta=root_state)) 

#simulate thresholded traits at the population level
threshold <- rep(0, dim(traits)[2])
traits_indiv <- lapply(1:dim(traits)[1], function(tip) rmvnorm(n = n_indiv, mean = traits[tip,], sigma = R)); names(traits_indiv) <- rownames(traits)
traits_indiv_discr <- lapply(1:dim(traits)[1], function(tip) t(sapply(1:n_indiv, function(indiv) as.numeric(traits_indiv[[tip]][indiv,] > threshold)))); names(traits_indiv_discr) <- names(traits_indiv) 

#identify unique site patterns
traits_indiv_discr_unique <- lapply(1:dim(traits)[1], function(tip) uniqueSP(traits_indiv_discr[[tip]]))

uniqueSP_redux <- function(indivs, freqs){
  indivs_str <- sapply(1:length(indivs[,1]), function(indiv) paste0(indivs[indiv,], collapse = "")) 
  uniq_indivs_str <- unique(indivs_str)
  counts <- sapply(1:length(uniq_indivs_str), function(string) sum(freqs[indivs_str == uniq_indivs_str[string]]))
  SPs <- t(apply((do.call(rbind, (strsplit(uniq_indivs_str, "")))), 1, as.numeric))
  return(list(SPs = SPs, counts = counts))
}

dat <- as.data.frame(do.call(rbind, lapply(1:nTaxa, function(tip) traits_indiv_discr_unique[[tip]]$SPs)))
colnames(dat) <- paste0("tr", 1:d_traits)
dat2 = data.frame(nSP = unlist(lapply(1:nTaxa, function(tip) traits_indiv_discr_unique[[tip]]$counts)),
                  tip = unlist(lapply(1:nTaxa, function(tip) rep(tip, length(traits_indiv_discr_unique[[tip]]$counts)))))

dat <- dat_orig <- cbind(dat, dat2)

#bad init
par <- c(rep(0, nTaxa * d_traits), 10, rep(0, choose(d_traits, 2)))

#better init
init_cor <- cor(do.call(rbind, lapply(1:nTaxa, function(x) (traits_indiv_discr[[x]]))))
init_cor[is.na(init_cor)] <- 0
init_means <- t(sapply(1:nTaxa, function(tip) apply(X = (traits_indiv_discr[[tip]]), MARGIN = 2, FUN = base::mean)))
init_means <- qnorm(init_means)
rownames(init_means) <- rownames(traits); colnames(init_means) <- colnames(traits)
nsamp_ptip <- sapply(1:length(traits_indiv_discr), function(tip) dim(traits_indiv_discr[[tip]])[1])
std_norm_5050_invariance <- sapply(1:length(nsamp_ptip), function(tip) qnorm(p = 1-0.5^(1/nsamp_ptip[tip]))) 
for(i in 1:nrow(init_means)){
  init_means[i, init_means[i,] == -Inf] <- std_norm_5050_invariance[i]
  init_means[i, init_means[i,] == Inf] <- -std_norm_5050_invariance[i]
}
# plot(init_means, traits); abline(0,1)
# plot(init_cor, R,xlim = c(-1,1), ylim = c(-1,1)); abline(0,1)
par <- c(c(t(init_means)), 10, init_cor[upper.tri(init_cor)])

#figure out shapw of regularizing term
maha <- function(x1, x2, vcvm, squared = FALSE){
  squaredDist <- t(x1-x2) %*% solve(vcvm) %*% (x1-x2)
  if(!squared){
    squaredDist
  } else {
    squaredDist^0.5
  }
} 

mahaMatrix <- function(data, vcvm, squared = FALSE){
  mat <- matrix (nrow = nrow(data), ncol = nrow(data), 0)
  for(i in 1:nrow(mat)){
    for(j in i:nrow(mat)){
      mat[i,j] <- maha(data[i,], data[j,], vcvm, squared = squared)
    }
  }
  mat <- mat + t(mat)
  rownames(mat) <- colnames(mat) <-rownames(data)
  mat
}

#get mahalonibis distance matrix
mahDistMat_init <- mahaMatrix(init_means, init_cor, squared = F)

#compute tree
njTree_init <- nj(mahDistMat_init) #can use outgroup for rooting
upgmaTree_init <- upgma(mahDistMat_init)

#find tree vcv and factorize
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
treeCV_init <- vcv.phylo(upgmaTree_init)
treeCV_init <- (cov2cor(treeCV_init) + diag(nTaxa)) / 2 #weight with star phylogeny
prunes_init <- prunePCV(treeCV_init)

par_to_opt_ind <- 26
par_to_opt <- par[par_to_opt_ind]
k_par <- par
k_par[par_to_opt_ind] <- NA  
dat <- list(dat_orig, k_par, prunes_init, rownames(traits))
logLL_bivProb_optim_single(par_to_opt, dat)

logLL_bivProb_optim_single <- function(par_to_opt, dat){ #TMB package may help with speedup
  treeReg <- T
  tipNames <- dat[[4]]
  prunes <- dat[[3]]
  k_par <- dat[[2]]
  dat <- dat[[1]]
  par_to_opt_ind <- which(is.na(k_par))
  par <- k_par
  par[par_to_opt_ind] <- par_to_opt
  
  hb <- length(par)
  nTaxa <- length(unique(dat$tip))
  d_trait <- ncol(dat) - 2
  cor <- diag(d_trait)
  cor[upper.tri(cor)] <- k_par[(nTaxa*d_trait+2) : hb];
  reg <- par[nTaxa * d_trait + 1]
  means <- matrix(par[1:(nTaxa*d_trait)], ncol = d_trait, nrow = nTaxa, byrow = T)
  rownames(means) <- tipNames
  
  if(par_to_opt_ind > 0 & par_to_opt_ind <= (nTaxa * d_trait)){
    par_type <- "mean"
    mean_trait <- par_to_opt_ind %% d_trait
    if(mean_trait == 0){
      mean_trait <- d_trait
    }
    mean_tip <- ceiling(par_to_opt_ind / d_trait)
  } else if (par_to_opt_ind == (nTaxa * d_trait + 1)){
    par_type = "reg"  
  } else {
    par_type = "corr"
    trs <- which(is.na(cor), arr.ind = T)
    tr1 = trs[1]
    tr2 = trs[2]
    cor[tr1, tr2] <- par_to_opt
  }
  
  cor <- cor + t(cor) - diag(d_trait)
  # psd_cor <- cov2cor(as.matrix(nearPD(cor)$mat))
  
  if(par_type == "mean"){
    #check biv dens of trait mean with all others traits for that tip
    cont_traits_of_tip <- means[mean_tip, -mean_trait]
    corrs_of_trait <- cor[mean_trait, -mean_trait]
    logLL <- -sum(unlist(sapply(1:(d_trait-1), function(trait) log(multi_pmvnorm_optim(mean = c(par_to_opt, cont_traits_of_tip[trait]),
                                                                           sigma = matrix(c(1,corrs_of_trait[trait],corrs_of_trait[trait],1),2,2),
                                                                           binaries = uniqueSP_redux(dat[dat$tip == mean_tip, c(mean_trait, (1:d_trait)[-mean_trait][trait])], dat[dat$tip == mean_tip,d_trait+1])$SPs)) 
                                                                           * uniqueSP_redux(dat[dat$tip == mean_tip, c(mean_trait, (1:d_trait)[-mean_trait][trait])], dat[dat$tip == mean_tip,d_trait+1])$counts)))
     
    
    if(treeReg){
      contrasts <- prunes[[1]] %*% means[rownames(prunes[[3]])[1:nTaxa],]
      regTerm <- -sum(sapply(1:(nTaxa-1), function(x) dnorm(x = contrasts[x,], mean = rep(0,d_trait), 
                                                              sd = sqrt(reg*prunes[[2]][x]), log = T)))
    } else {
      regTerm <- -sum(sapply(1:(d_trait-1), function(trait) mvtnorm::dmvnorm(x = c(par_to_opt, cont_traits_of_tip[trait]),
                                                                            sigma = matrix(c(1,corrs_of_trait[trait],corrs_of_trait[trait],1),2,2) * reg,
                                                                            mean = c(0,0), log = T)))
    }
    
  } else if(par_type == "reg"){
    logLL <- 0
    
    if(treeReg){
      contrasts <- prunes[[1]] %*% means[rownames(prunes[[3]])[1:nTaxa],]
      regTerm <- -sum(sapply(1:(nTaxa-1), function(x) dnorm(x = contrasts[x,], mean = rep(0,d_trait),
                                                                sd = sqrt(reg*prunes[[2]][x]), log = T)))
    } else {
      regTerm <- -sum(dnorm(x = means, sd = sqrt(reg), mean = 0, log = T)) 
    }
    
  } else if(par_type == "corr"){
    cont_traits <- means[, c(tr1, tr2)]
    logLL <- -sum(unlist(sapply(1:nTaxa, function(tip) log(multi_pmvnorm_optim(mean = cont_traits[tip,],
                                                                                sigma = matrix(c(1,par_to_opt,par_to_opt,1),2,2),
                                                                                binaries = uniqueSP_redux(dat[dat$tip == tip, c(tr1, tr2)], dat[dat$tip == tip,d_trait+1])$SPs)) 
                                                                                * uniqueSP_redux(dat[dat$tip == tip, c(tr1, tr2)], dat[dat$tip == tip,d_trait+1])$counts)))
    
    #FIX REDUNDANT UNIQUESP_REDUX
    regTerm <- 0
    
  }
  
  returnVal <- logLL + regTerm
  if(returnVal == Inf){returnVal <- 1E10} #hacky temporary fix to make optim's generalized "L-BFGS-B" method stop complaining about singularities :/
  return(returnVal)
}

nparam <- length(par)
tr_cors_to_opt <- sort(sample(1:d_traits, replace = F, size = 3))
indMat <- matrix(0, d_traits, d_traits)
indMat[upper.tri(indMat)] <- (nTaxa * d_traits + 1) + 1:choose(d_traits, 2)
pars_to_opt_ind <- c(indMat[tr_cors_to_opt[1], tr_cors_to_opt[2]],
                     indMat[tr_cors_to_opt[1], tr_cors_to_opt[3]],
                     indMat[tr_cors_to_opt[2], tr_cors_to_opt[3]])
pars_to_opt <- par[pars_to_opt_ind]
k_par <- par
k_par[pars_to_opt_ind] <- NA  
dat <- list(dat_orig, k_par, prunes_init, rownames(traits))
logLL_trivProb_optim_double(pars_to_opt, dat)

logLL_trivProb_optim_double <- function(pars_to_opt, dat){ #TMB package may help with speedup
  treeReg <- T
  tipNames <- dat[[4]]
  prunes <- dat[[3]]
  k_par <- dat[[2]]
  dat <- dat[[1]]
  pars_to_opt_ind <- which(is.na(k_par))
  par <- k_par
  par[pars_to_opt_ind] <- pars_to_opt
  
  hb <- length(par)
  nTaxa <- length(unique(dat$tip))
  d_trait <- ncol(dat) - 2
  cor <- diag(d_trait)
  cor[upper.tri(cor)] <- k_par[(nTaxa*d_trait+2) : hb];
  reg <- par[nTaxa * d_trait + 1]
  means <- matrix(par[1:(nTaxa*d_trait)], ncol = d_trait, nrow = nTaxa, byrow = T)
  rownames(means) <- tipNames
  
  par_type = "corr"
  trs <- sort(unique(c(which(is.na(cor), arr.ind = T))))
  
  cor[trs[1], trs[2]] <- pars_to_opt[1]
  cor[trs[1], trs[3]] <- pars_to_opt[2]
  cor[trs[2], trs[3]] <- pars_to_opt[3]

  cor <- cor + t(cor) - diag(d_trait)

  cont_traits <- means[, trs]
  cormat <- cor[trs, trs]
  logLL <- -sum(unlist(sapply(1:nTaxa, function(tip) log(multi_pmvnorm_optim(mean = cont_traits[tip,],
                                                                             sigma = cormat,
                                                                             binaries = uniqueSP_redux(dat[dat$tip == tip, trs], dat[dat$tip == tip,d_trait+1])$SPs)) 
                              * uniqueSP_redux(dat[dat$tip == tip, trs], dat[dat$tip == tip,d_trait+1])$counts)))
  

  regTerm <- 0
  
  returnVal <- logLL + regTerm
  if(returnVal == Inf){returnVal <- 1E10} #hacky temporary fix to make optim's generalized "L-BFGS-B" method stop complaining about singularities :/
  return(returnVal)
}

library(microbenchmark)
fast <- T
if(!fast){
  for(i in 1:length(par)){
    print(i)
    par_to_opt_ind <- i
    par_to_opt <- par[par_to_opt_ind]
    k_par <- par
    k_par[par_to_opt_ind] <- NA
    dat <- list(dat_orig, k_par)
    print(microbenchmark(logLL_bivProb_optim_single(par_to_opt, dat), times = 10))
  }
}

method <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent")[4]

par_to_opt_ind <- 1
par_to_opt <- par[par_to_opt_ind]
k_par <- par
k_par[par_to_opt_ind] <- NA  
dat <- list(dat_orig, k_par, prunes_init)
logLL_bivProb_optim_single(par_to_opt, dat)

t_s <- Sys.time()
if(!fast){
  optim_out <- optim(par = par_to_opt, logLL_bivProb_optim_single, dat = dat, method = method, control = list(trace = 0, REPORT = 10),
                   upper = c(rep(Inf, nTaxa * d_traits + 1), rep(1, choose(d_traits, 2))), 
                   lower = c(rep(-Inf, nTaxa * d_traits), 0, rep(-1, choose(d_traits, 2))))
}
Sys.time() - t_s

#now iterate over the elements of par

nrounds <- 4
nparam <- length(par)
prunes <- prunes_init
for(i in 1:nrounds){
  cat(paste0("\n", "round ", i, ": "))
  param_ord <- sample(c(1:nparam, (nTaxa * d_traits + 2):nparam)) #double emphasis on correlation parameters
  for(par_to_opt_ind in param_ord){
    type_of_param <- ifelse(par_to_opt_ind <= nTaxa * d_traits, "m", ifelse(par_to_opt_ind == (nTaxa * d_traits + 1), "v", "c"))
    cat(paste0(par_to_opt_ind, type_of_param, "-"))
    par_to_opt <- par[par_to_opt_ind]
    k_par <- par
    k_par[par_to_opt_ind] <- NA  
    dat <- list(dat_orig, k_par, prunes, rownames(traits))
    upper = c(rep(Inf, nTaxa * d_traits), Inf, rep(1, choose(d_traits, 2)))[par_to_opt_ind] 
    lower = c(rep(-Inf, nTaxa * d_traits), 0.01, rep(-1, choose(d_traits, 2)))[par_to_opt_ind] 
    optim_out <- optim(par = par_to_opt, logLL_bivProb_optim_single, dat = dat, method = method, control = list(trace = 0, REPORT = 10),
                       upper = upper, lower = lower)
    cat(paste0("(", round(abs((par[par_to_opt_ind] - optim_out$par) / par[par_to_opt_ind] * 100)), ")-"))
    par[par_to_opt_ind] <- optim_out$par
  }
  
  #update tree
  optim_est <- par
  R_est <- diag(d_traits)
  R_est[upper.tri(R_est)] <- optim_est[(nTaxa*d_traits+2) : length(optim_est)];
  R_est <- as.matrix(nearPD(R_est + t(R_est) - diag(d_traits), corr = T)$mat)
  est_means <- matrix(optim_est[1:(nTaxa*d_traits)], ncol = d_traits, nrow = nTaxa, byrow = T, dimnames = list(rownames(traits)))
  
  mahDistMat <- mahaMatrix(est_means, R_est, squared = F)
  upgmaTree <- upgma(mahDistMat)
  treeCV <- vcv.phylo(upgmaTree)
  treeCV <- (cov2cor(treeCV) + diag(nTaxa)) / 2
  prunes <- prunePCV(treeCV)
  
}

nrounds <- 2
nparam <- length(par)
indMat <- matrix(0, d_traits, d_traits)
indMat[upper.tri(indMat)] <- (nTaxa * d_traits + 1) + 1:choose(d_traits, 2)
for(i in 1:nrounds){
  cat(paste0("\n\n", "round ", i, ": "))
  
  trs_to_opt <- sample(1:d_traits, replace = F)
  trs_to_opt <- c(trs_to_opt, trs_to_opt[1:(3 - length(trs_to_opt) %% 3)])
  tr_cors_to_opt <- matrix(trs_to_opt, nrow = 3)
  tr_cors_to_opt <- sapply(1:ncol(tr_cors_to_opt), function(j) sort(tr_cors_to_opt[,j]))
  
  trs_to_opt <- t(which(upper.tri(diag(d_traits)), arr.ind = T))
  trs_to_opt_3 <- sapply(1:choose(d_traits, 2), function(j) sample((1:d_traits)[-trs_to_opt[,j]], size = 1))
  tr_cors_to_opt <- rbind(trs_to_opt, trs_to_opt_3)
  tr_cors_to_opt <- as.matrix(sapply(1:ncol(tr_cors_to_opt), function(j) sort(tr_cors_to_opt[,j])))
  
  pars_to_opt_ind_ord <- sapply(1:ncol(tr_cors_to_opt), function(j) 
                     c(indMat[tr_cors_to_opt[,j][1], tr_cors_to_opt[,j][2]],
                       indMat[tr_cors_to_opt[,j][1], tr_cors_to_opt[,j][3]],
                       indMat[tr_cors_to_opt[,j][2], tr_cors_to_opt[,j][3]]))
  
  
  for(par_set in 1:ncol(pars_to_opt_ind_ord)){
    type_of_param <- "c"
    cat(paste0(paste0(par_set, collapse = "-"), type_of_param, "-"))
    
    pars_to_opt_ind <- pars_to_opt_ind_ord[,par_set]
    pars_to_opt <- par[pars_to_opt_ind]
    k_par <- par
    k_par[pars_to_opt_ind] <- NA  
    dat <- list(dat_orig, k_par, prunes, rownames(traits))
    
    upper = c(rep(Inf, nTaxa * d_traits), Inf, rep(1, choose(d_traits, 2)))[pars_to_opt_ind] 
    lower = c(rep(-Inf, nTaxa * d_traits), 0, rep(-1, choose(d_traits, 2)))[pars_to_opt_ind] 
    
    optim_out <- optim(par = pars_to_opt, logLL_trivProb_optim_double, dat = dat, method = method, control = list(trace = 0, REPORT = 1),
                       upper = upper, lower = lower)
    if(optim_out$value == 1e+10){
      optim_out <- optim(par = c(0,0,0), logLL_trivProb_optim_double, dat = dat, method = method, control = list(trace = 0, REPORT = 1),
                         upper = upper, lower = lower)
    }

    cat(paste0("(", round(mean(abs(par[pars_to_opt_ind] - optim_out$par)), 2), ")-"))
    
    par[pars_to_opt_ind] <- optim_out$par
  }
}

optim_est <- par
R_est <- diag(d_traits)
R_est[upper.tri(R_est)] <- optim_est[(nTaxa*d_traits+2) : length(optim_est)];
R_est <- as.matrix(nearPD(R_est + t(R_est) - diag(d_traits), corr = T)$mat)
reg <- optim_est[nTaxa * d_traits + 1]
paste0("regularizing variance = ", round(reg, 3))
est_means <- matrix(optim_est[1:(nTaxa*d_traits)], ncol = d_traits, nrow = nTaxa, byrow = T)
par(mfrow = c(2,2))
plot(init_means, traits, main = paste0("r = ", round(cor(c(traits), c(init_means)), 3))); abline(0,1)
plot(init_cor[upper.tri(R)], R[upper.tri(R)], xlim = c(-1,1), ylim = c(-1,1), xlab = "Pearson's 'r's", ylab = "true correlations",
     main = paste0("r = ", round(cor(c(R[upper.tri(R)]), c(init_cor[upper.tri(R)])), 3))); abline(0,1)
plot(est_means, traits, main = paste0("r = ", round(cor(c(traits), c(est_means)), 3))); abline(0,1)
plot(R_est[upper.tri(R)], R[upper.tri(R)], xlim = c(-1,1), ylim = c(-1,1), xlab = "estimated correlations", ylab = "true correlations",
     main = paste0("r = ", round(cor(c(R[upper.tri(R)]), c(R_est[upper.tri(R)])), 3))); abline(0,1)

######################################################################
######################################################################
###OKIE-DOKIE NEW ATTEMPT NOW LET'S TRY IT WITH MULTIPLE THRESHOLDS###
######################################################################
######################################################################


#first let's resimulate data
nTaxa <- 10
d_traits <- 20
n_indiv <- sample(15:40, nTaxa, replace = T)
n_thresholds <- sample(1:5, d_traits, T) #for variable numbers of thresholds, just need to buffer with 'Inf's on the right
weighPCV <- F

tree <- tess.sim.taxa(n = 1, nTaxa = nTaxa, lambda = 1, mu = 0, max = 1E3)[[1]]
target_tree_height <- 2
tree$edge.length <- tree$edge.length * target_tree_height / node.depth.edgelength(tree)[1]

#sample thresholds from scaled lognormal distribution

logn_mu <- 0.5
logn_sd <- 0.25
threshold_spacings <- sapply(1:d_traits, function(trait) exp(rnorm(n = n_thresholds[trait]-1, mean = logn_mu, sd = logn_sd)) / (n_thresholds[trait]-1))
threshold_spacings[n_thresholds == 1] <- Inf
max_thresholds_count <- max(n_thresholds)
threshold_spacings <- t(sapply(1:length(threshold_spacings), function(trait) c(threshold_spacings[[trait]], rep(Inf, max_thresholds_count - 1 - length(threshold_spacings[[trait]])))))
# thresholds <- cbind(rep(0, d_traits), matrix(runif(n = (n_thresholds - 1) * d_traits, min = 0.25, max = 1), nrow = d_traits, ncol = n_thresholds-1))
thresholds <- cbind(rep(0, d_traits), threshold_spacings)
thresholds <- t(apply(thresholds, 1, cumsum))
threshold_diffs <- t(sapply(1:d_traits, function(trait) diff(thresholds[trait,])))
threshold_diffs[is.na(threshold_diffs)] <- Inf
R <- rlkj(d_traits); rownames(R) <- colnames(R) <- paste0("tr", 1:dim(R)[1])

root_state <- sapply(1:length(n_thresholds), function(trait) ifelse(n_thresholds[trait] %% 2 == 1, thresholds[trait,(median(1:n_thresholds[trait]))], 
                                                      mean(c(thresholds[trait,(median(1:n_thresholds[trait]) - 0.5)], thresholds[trait,(median(1:n_thresholds[trait]) + 0.5)]))))

#for use with a single threshold
# if(n_thresholds %% 2 == 1){root_state <- thresholds[,(median(1:n_thresholds))]} else {root_state <- (thresholds[,(median(1:n_thresholds) - 0.5)] + thresholds[,(median(1:n_thresholds) + 0.5)]) / 2}

traits <- mvSIM(tree, model="BM1", param=list(sigma=R, theta=root_state)) 

#simulate thresholded traits at the population level
traits_indiv <- lapply(1:nTaxa, function(tip) rmvnorm(n = n_indiv[tip], mean = traits[tip,], sigma = R)); names(traits_indiv) <- rownames(traits)
traits_indiv_discr <- lapply(1:nTaxa, function(tip) t(sapply(1:n_indiv[tip], function(indiv) sapply(1:d_traits, function(trait) sum(traits_indiv[[tip]][indiv,trait] - thresholds[trait,] > 0))))); 
names(traits_indiv_discr) <- names(traits_indiv) 

#identify unique site patterns
traits_indiv_discr_unique <- lapply(1:dim(traits)[1], function(tip) uniqueSP(traits_indiv_discr[[tip]]))

dat <- as.data.frame(do.call(rbind, lapply(1:nTaxa, function(tip) traits_indiv_discr_unique[[tip]]$SPs)))
colnames(dat) <- paste0("tr", 1:d_traits)
dat2 = data.frame(nSP = unlist(lapply(1:nTaxa, function(tip) traits_indiv_discr_unique[[tip]]$counts)),
                  tip = unlist(lapply(1:nTaxa, function(tip) rep(tip, length(traits_indiv_discr_unique[[tip]]$counts)))))

dat <- dat_orig <- cbind(dat, dat2)

#bad init
# par <- c(rep(0, nTaxa * d_traits), 10, rep(0, choose(d_traits, 2)), c(sapply(0:(n_thresholds-1), function(x) rep(x, d_trait))))

#better init
thresholds_init <- matrix(rep(0:(max_thresholds_count-1), d_traits), nrow = d_traits, ncol = max_thresholds_count,  byrow = T) / 2
thresholds_init <- t(sapply(1:d_traits, function(trait) c(thresholds_init[tip,1:n_thresholds[trait]], rep(Inf, max_thresholds_count - n_thresholds[trait]))))
threshold_diffs_init <- t(sapply(1:d_traits, function(trait) diff(thresholds_init[trait,]))); threshold_diffs_init[is.na(threshold_diffs_init)] <- Inf
n_thresh <- sapply(1:length(thresholds_init[,1]), function(trait) length((thresholds_init[trait,])))
locations <- cbind(rep(-0.5, d_traits), thresholds_init, rep(Inf, d_traits))
infs_first_index <- sapply(1:d_traits, function(trait) which(locations[trait,] == Inf)[1])
for(i in 1:d_traits){locations[i,infs_first_index] <- locations[i,infs_first_index-1] + 0.5}
init_cor <- cor(do.call(rbind, lapply(1:nTaxa, function(x) (traits_indiv_discr[[x]]))))
init_cor[is.na(init_cor)] <- 0

init_means <- t(sapply(1:nTaxa, function(tip) 
  sapply(1:length(traits_indiv_discr[[tip]][1,]), function(trait) mean(
    sapply(1:length(traits_indiv_discr[[tip]][,1]), function(indiv) 
      mean(locations[trait, c(traits_indiv_discr[[tip]][indiv,][trait] + 1, traits_indiv_discr[[tip]][indiv,][trait] + 2)]))
      )
    )
  )
)

rownames(init_means) <- rownames(traits); colnames(init_means) <- colnames(traits)

par <- c(c(t(init_means)), 10, init_cor[upper.tri(init_cor)], c(threshold_diffs_init))


#figure out shapw of regularizing term
maha <- function(x1, x2, vcvm, squared = FALSE){
  squaredDist <- t(x1-x2) %*% solve(vcvm) %*% (x1-x2)
  if(!squared){
    squaredDist
  } else {
    squaredDist^0.5
  }
} 

mahaMatrix <- function(data, vcvm, squared = FALSE){
  mat <- matrix (nrow = nrow(data), ncol = nrow(data), 0)
  for(i in 1:nrow(mat)){
    for(j in i:nrow(mat)){
      mat[i,j] <- maha(data[i,], data[j,], vcvm, squared = squared)
    }
  }
  mat <- mat + t(mat)
  rownames(mat) <- colnames(mat) <-rownames(data)
  mat
}

#get mahalonibis distance matrix
mahDistMat_init <- mahaMatrix(init_means, init_cor, squared = F)

#compute tree
njTree_init <- nj(mahDistMat_init) #can use outgroup for rooting
upgmaTree_init <- upgma(mahDistMat_init)

#find tree vcv and factorize
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
treeCV_init <- cov2cor(vcv.phylo(upgmaTree_init))
if(weighPCV){treeCV_init <- (treeCV_init + diag(nTaxa)) / 2} #weight with star phylogeny
prunes_init <- prunePCV(treeCV_init)

multi_pmvnorm_optim_multithresh <- function(mean, sigma, ordinals, thresholds, algorithm = c("GenzBretz", "Miwa")[1]){
    n_traits <- dim(thresholds)[1]
    if(is.null(n_traits)){n_traits <- 1; thresholds <- c(rep(-Inf, n_traits), thresholds, rep(Inf, n_traits))} else{
    thresholds <- cbind(rep(-Inf, n_traits), thresholds, rep(Inf, n_traits))}
    
    if(n_traits == 2){
    upper <- t(sapply(1:length(ordinals[,1]), function(pat) c(thresholds[1,ordinals[pat,1]+2], thresholds[2,ordinals[pat,2]+2])))
    lower <- t(sapply(1:length(ordinals[,1]), function(pat) c(thresholds[1,ordinals[pat,1]+1], thresholds[2,ordinals[pat,2]+1])))
    
    to_return <- sapply(1:dim(ordinals)[1], function(indiv) pmvnorm(lower = as.numeric(lower[indiv,]), 
                                                                    upper = as.numeric(upper[indiv,]), 
                                                                    mean = mean, sigma = sigma, algorithm = algorithm))
    to_return[to_return < 0] <- 0
    return(to_return)
    } else if (n_traits == 1){
      upper <- sapply(1:length(ordinals), function(pat) c(thresholds[ordinals[pat]+2]))
      lower <- sapply(1:length(ordinals), function(pat) c(thresholds[ordinals[pat]+1]))
      to_return <- pnorm(upper, mean = mean, sd = 1) - pnorm(lower, mean = mean, sd = 1)
      to_return[to_return < 0] <- 0
      return(to_return)
    }
}

pbivnorm_infs <- function(upper1, upper2, rhos){
  nInf <- (upper1 == -Inf | upper2 == -Inf)
  Inf1 <- upper1 == Inf
  Inf2 <- upper2 == Inf
  n2_1 <- Inf1 & !Inf2
  n1_2 <- Inf2 & !Inf1
  include <- !(nInf | Inf1 | Inf2)
  
  results <- rep(0, length(include))
  results[Inf1 & Inf2] <- 1
  results[n2_1] <- pnorm(q = upper2[n2_1])
  results[n1_2] <- pnorm(q = upper1[n1_2])
  results[include] <- pbivnorm(upper1[include], upper2[include], rhos[include])
  return(results)
}

uniqueSP_fast <- function(indivs, freqs){
  newFreqs <- as.data.frame(wtd.table(indivs[,1], indivs[,2], weights = freqs))
  newFreqs <- (newFreqs[newFreqs$Freq > 0,])
  newFreqs$Var1 <- as.numeric(levels(newFreqs$Var1))[newFreqs$Var1]
  newFreqs$Var2 <- as.numeric(levels(newFreqs$Var2))[newFreqs$Var2]
  return(list(SPs = as.matrix(newFreqs[,1:2]), counts = newFreqs[,3]))
}


getUniqueSitePatterns <- function(par_to_opt, dat){ #TMB package may help with speedup
  vectorized <- T
  univThresh <- F
  treeReg <- T
  n_thresh <- dat[[5]]
  tipNames <- dat[[4]]
  prunes <- dat[[3]]
  k_par <- dat[[2]]
  dat <- dat[[1]]
  par_to_opt_ind <- which(is.na(k_par))
  par <- k_par
  par[par_to_opt_ind] <- par_to_opt
  
  hb <- length(par)
  nTaxa <- length(unique(dat$tip))
  d_trait <- ncol(dat) - 2
  cor <- diag(d_trait)
  cor[upper.tri(cor)] <- k_par[(nTaxa*d_trait+2) : (nTaxa*d_trait+1+choose(d_trait, 2))];
  reg <- par[nTaxa * d_trait + 1]
  means <- matrix(par[1:(nTaxa*d_trait)], ncol = d_trait, nrow = nTaxa, byrow = T)
  rownames(means) <- tipNames
  threshold_mat <- cbind(rep(0, d_trait), matrix(par[(nTaxa*d_trait+2+choose(d_trait, 2)) : ((nTaxa*d_trait+1+choose(d_trait, 2)) + sum(n_thresh-1))], nrow = d_trait, ncol = n_thresh[1]-1))
  
  if(length(par_to_opt_ind) == "1"){
    if(par_to_opt_ind > 0 & par_to_opt_ind <= (nTaxa * d_trait)){
      par_type <- "mean"
      mean_trait <- par_to_opt_ind %% d_trait
      if(mean_trait == 0){
        mean_trait <- d_trait
      }
      mean_tip <- ceiling(par_to_opt_ind / d_trait)
    } else if (par_to_opt_ind == (nTaxa * d_trait + 1)){
      par_type = "reg"  
    } else if (par_to_opt_ind >= (nTaxa*d_trait+2) & par_to_opt_ind <= (nTaxa*d_trait+1+choose(d_trait, 2))){
      par_type = "corr"
      trs <- which(is.na(cor), arr.ind = T)
      tr1 = trs[1]
      tr2 = trs[2]
      cor[tr1, tr2] <- par_to_opt
    } else if(par_to_opt_ind > (nTaxa*d_trait+1+choose(d_trait, 2)) & par_to_opt_ind <= ((nTaxa*d_trait+1+choose(d_trait, 2)) + sum(n_thresh-1))){
      par_type <- "thresh"
      threshold_mat_unknown <- cbind(rep(0, d_trait), matrix(k_par[(nTaxa*d_trait+2+choose(d_trait, 2)) : ((nTaxa*d_trait+1+choose(d_trait, 2)) + sum(n_thresh-1))], nrow = d_trait, ncol = n_thresh[1]-1))
      thresh_trait <- which(is.na(threshold_mat_unknown), arr.ind = T)
      thresh_ind <- thresh_trait[2]
      thresh_trait <- thresh_trait[1]
      
    }
  } else if(all(par_to_opt_ind > (nTaxa*d_trait+1+choose(d_trait, 2)) & par_to_opt_ind <= ((nTaxa*d_trait+1+choose(d_trait, 2)) + sum(n_thresh-1)))){
    par_type <- "thresh"
    threshold_mat_unknown <- cbind(rep(0, d_trait), matrix(k_par[(nTaxa*d_trait+2+choose(d_trait, 2)) : ((nTaxa*d_trait+1+choose(d_trait, 2)) + sum(n_thresh-1))], nrow = d_trait, ncol = n_thresh[1]-1))
    thresh_trait <- which(is.na(threshold_mat_unknown), arr.ind = T)
    thresh_ind <- thresh_trait[,2]
    thresh_trait <- thresh_trait[1,1]
  }
  
  cor <- cor + t(cor) - diag(d_trait)
  threshold_mat <- t(sapply(1:nrow(threshold_mat), function(trait) cumsum(threshold_mat[trait,])))
  
  if(par_type == "mean"){
        uniqueSitePatterns <- sapply(1:(d_trait-1), function(trait) uniqueSP_fast(dat[dat$tip == mean_tip, c(mean_trait, (1:d_trait)[-mean_trait][trait])], dat[dat$tip == mean_tip,d_trait+1]))
  } else if(par_type == "corr"){
        uniqueSitePatterns <- sapply(1:nTaxa, function(tip)  uniqueSP_fast(dat[dat$tip == tip, c(tr1, tr2)], dat[dat$tip == tip,d_trait+1]))
  }  else if(par_type == "thresh"){
        uniqueSitePatterns <- lapply(1:nTaxa, function(tip) sapply(1:(d_trait-1), function(trait) uniqueSP_fast(dat[dat$tip == tip, c(thresh_trait, (1:d_trait)[-thresh_trait][trait])], dat[dat$tip == tip,d_trait+1])))
  }
  return(uniqueSitePatterns)
}

logLL_bivProb_optim_multithresh <- function(par_to_opt, dat){ #TMB package may help with speedup
  vectorized <- T
  univThresh <- F
  treeReg <- T
  mvBMreg <- F
  enforcePSD <- F
  if(vectorized & length(dat) >= 6){
    uniqueSitePatterns_precomp <- dat[[6]]
  }
  n_thresh <- dat[[5]]
  tipNames <- dat[[4]]
  prunes <- dat[[3]]
  k_par <- dat[[2]]
  dat <- dat[[1]]
  par_to_opt_ind <- which(is.na(k_par))
  par <- k_par
  par[par_to_opt_ind] <- par_to_opt
  
  hb <- length(par)
  nTaxa <- length(unique(dat$tip))
  d_trait <- ncol(dat) - 2
  cor <- diag(d_trait)
  cor[upper.tri(cor)] <- k_par[(nTaxa*d_trait+2) : (nTaxa*d_trait+1+choose(d_trait, 2))];
  reg <- par[nTaxa * d_trait + 1]
  means <- matrix(par[1:(nTaxa*d_trait)], ncol = d_trait, nrow = nTaxa, byrow = T)
  rownames(means) <- tipNames
  threshold_mat <- cbind(rep(0, d_trait), matrix(par[(nTaxa*d_trait+2+choose(d_trait, 2)) : ((nTaxa*d_trait+1+choose(d_trait, 2)) + sum(n_thresh-1))], nrow = d_trait, ncol = n_thresh[1]-1))
  
  if(length(par_to_opt_ind) == "1"){
    if(par_to_opt_ind > 0 & par_to_opt_ind <= (nTaxa * d_trait)){
      par_type <- "mean"
      mean_trait <- par_to_opt_ind %% d_trait
      if(mean_trait == 0){
        mean_trait <- d_trait
      }
      mean_tip <- ceiling(par_to_opt_ind / d_trait)
    } else if (par_to_opt_ind == (nTaxa * d_trait + 1)){
      par_type = "reg"  
    } else if (par_to_opt_ind >= (nTaxa*d_trait+2) & par_to_opt_ind <= (nTaxa*d_trait+1+choose(d_trait, 2))){
      par_type = "corr"
      trs <- which(is.na(cor), arr.ind = T)
      tr1 = trs[1]
      tr2 = trs[2]
      cor[tr1, tr2] <- par_to_opt
    } else if(par_to_opt_ind > (nTaxa*d_trait+1+choose(d_trait, 2)) & par_to_opt_ind <= ((nTaxa*d_trait+1+choose(d_trait, 2)) + sum(n_thresh-1))){
      par_type <- "thresh"
      threshold_mat_unknown <- cbind(rep(0, d_trait), matrix(k_par[(nTaxa*d_trait+2+choose(d_trait, 2)) : ((nTaxa*d_trait+1+choose(d_trait, 2)) + sum(n_thresh-1))], nrow = d_trait, ncol = n_thresh[1]-1))
      thresh_trait <- which(is.na(threshold_mat_unknown), arr.ind = T)
      thresh_ind <- thresh_trait[2]
      thresh_trait <- thresh_trait[1]
      
    }
  } else if(all(par_to_opt_ind > (nTaxa*d_trait+1+choose(d_trait, 2)) & par_to_opt_ind <= ((nTaxa*d_trait+1+choose(d_trait, 2)) + sum(n_thresh-1)))){
    par_type <- "thresh"
    threshold_mat_unknown <- cbind(rep(0, d_trait), matrix(k_par[(nTaxa*d_trait+2+choose(d_trait, 2)) : ((nTaxa*d_trait+1+choose(d_trait, 2)) + sum(n_thresh-1))], nrow = d_trait, ncol = n_thresh[1]-1))
    thresh_trait <- which(is.na(threshold_mat_unknown), arr.ind = T)
    thresh_ind <- thresh_trait[,2]
    thresh_trait <- thresh_trait[1,1]
  }
  
    cor <- cor + t(cor) - diag(d_trait)
    threshold_mat <- t(sapply(1:nrow(threshold_mat), function(trait) cumsum(threshold_mat[trait,])))
  
  # psd_cor <- cov2cor(as.matrix(nearPD(cor)$mat))
    if(enforcePSD & par_type == "corr"){
      if(!is.positive.definite(cor)){
        return(1E100)
      }
    }
    
  if(par_type == "mean"){
    #check biv dens of trait mean with all others traits for that tip
    
    cont_traits_of_tip <- means[mean_tip, -mean_trait]
    corrs_of_trait <- cor[mean_trait, -mean_trait]
    threshs_of_other_traits <- threshold_mat[-mean_trait,]
    
    if(!vectorized){
      
      logLL <- -sum(unlist(sapply(1:(d_trait-1), function(trait) log(multi_pmvnorm_optim_multithresh(mean = c(par_to_opt, cont_traits_of_tip[trait]),
                                                                                         sigma = matrix(c(1,corrs_of_trait[trait],corrs_of_trait[trait],1),2,2),
                                                                                         ordinals = uniqueSP_redux(dat[dat$tip == mean_tip, c(mean_trait, (1:d_trait)[-mean_trait][trait])], dat[dat$tip == mean_tip,d_trait+1])$SPs,
                                                                                         thresholds = rbind(threshold_mat[mean_trait,], threshs_of_other_traits[trait,]))) * 
                                    uniqueSP_redux(dat[dat$tip == mean_tip, c(mean_trait, (1:d_trait)[-mean_trait][trait])], dat[dat$tip == mean_tip,d_trait+1])$counts)))
    } else if(vectorized){ ## alternatively, in vectorized form
      
      if(!exists("uniqueSitePatterns_precomp")){
        uniqueSitePatterns <- sapply(1:(d_trait-1), function(trait) uniqueSP_fast(dat[dat$tip == mean_tip, c(mean_trait, (1:d_trait)[-mean_trait][trait])], dat[dat$tip == mean_tip,d_trait+1]))
      } else {
        uniqueSitePatterns <-uniqueSitePatterns_precomp
      }
      
      #convert to standard normal
      nOrds <- sapply(1:(d_trait-1), function(trait) nrow(uniqueSitePatterns[1,trait][[1]]))
      nOrds_traits <- rep(1:(d_trait-1), nOrds)
      new_std_thresholds <- cbind(rep(-Inf, d_trait), sapply(1:ncol(threshold_mat), function(col) threshold_mat[,col] - means[mean_tip,]), rep(Inf, d_trait))
      SPs <- do.call(rbind, uniqueSitePatterns[1,])
      # uppers <- t(sapply(1:length(nOrds_traits), function(pat) c(new_std_thresholds[mean_trait,SPs[pat,1]+2], new_std_thresholds[nOrds_traits[pat] + ifelse(nOrds_traits[pat] < mean_trait, 0, 1),SPs[pat,2]+2])))
      # lowers <- t(sapply(1:length(nOrds_traits), function(pat) c(new_std_thresholds[mean_trait,SPs[pat,1]+1], new_std_thresholds[nOrds_traits[pat] + ifelse(nOrds_traits[pat] < mean_trait, 0, 1),SPs[pat,2]+1])))
      uppers_lowers <- t(sapply(1:length(nOrds_traits), function(pat) c(new_std_thresholds[mean_trait,SPs[pat,1]+c(2,1)], new_std_thresholds[nOrds_traits[pat] + ifelse(nOrds_traits[pat] < mean_trait, 0, 1), SPs[pat,2]+c(2,1)])))
      uppers <- uppers_lowers[,c(1,3)]
      lowers <- uppers_lowers[,c(2,4)]
      
      rhos = rep(corrs_of_trait, nOrds)
      counts = unlist(uniqueSitePatterns[2,])
      
      
      probs <- pbivnorm_infs(uppers[,1], uppers[,2], rhos) - 
        pbivnorm_infs(lowers[,1], uppers[,2], rhos) - 
        pbivnorm_infs(uppers[,1], lowers[,2], rhos) + 
        pbivnorm_infs(lowers[,1], lowers[,2], rhos)
      probs[probs < 0] <- 0 #get rid of floating point errors
      logLL <- -sum(log(probs) * counts)

    }
      # TODO make regularizing term multivariate normal
    if(treeReg){
      
      contrasts <- prunes[[1]] %*% means[rownames(prunes[[3]])[1:nTaxa],]
      
      if(!mvBMreg){
        regTerm <- -sum(sapply(1:(nTaxa-1), function(x) dnorm(x = contrasts[x,], mean = rep(0,d_trait), #mean 0 is marginalizes over all mean states
                                                            sd = sqrt(reg*prunes[[2]][x]), log = T)))
      } else if(mvBMreg){
        #equiv to above
        # regTerm <- -sum(sapply(1:(nTaxa-1), function(x) mvtnorm::dmvnorm(x = contrasts[x,], mean = rep(0,d_trait), #mean 0 is marginalizes over all mean states
        #                                                                  sigma = (reg*prunes[[2]][x]*diag(d_trait)), log = T)))
        # 
        # npdcor <- cov2cor(as.matrix(nearPD(cor)$mat))
        # regTerm <- -sum(sapply(1:(nTaxa-1), function(x) mvtnorm::dmvnorm(x = contrasts[x,], mean = rep(0,d_trait), #mean 0 is marginalizes over all mean states
        #                                                                  sigma = (reg*prunes[[2]][x]*npdcor), log = T)))
        # 
        # -sum(sapply(1:(nTaxa-1), function(x) mvtnorm::dmvnorm(x = contrasts[x,], mean = rep(0,d_trait), #mean 0 is marginalizes over all mean states
        #                                                       sigma = (reg*prunes[[2]][x]*npdcor), log = T)))
        # 
        # -sum(sapply(1:(nTaxa-1), function(x) mvtnorm::dmvnorm(x = contrasts[x,], mean = rep(0,d_trait), #mean 0 is marginalizes over all mean states
        #                                                       sigma = (reg*prunes[[2]][x]*diag(d_trait)), log = T)))
        # 
        # -sum(sapply(1:(nTaxa-1), function(x) mvtnorm::dmvnorm(x = contrasts[x,], mean = rep(0,d_trait), #mean 0 is marginalizes over all mean states
        #                                                       sigma = (reg*prunes[[2]][x]*(diag(d_trait)+npdcor)/2), log = T)))
        # 
        # 
        # univContrasts <-(contrasts) %*% chol(npdcor)
        # detR <- det(npdcor)
        # LLs[i] <- sum(dnorm(contrast / BLength^0.5, mean = 0, sd = 1, log = T)) 
        # 
        # regTerm <- -sum(sapply(1:(nTaxa-1), function(x) sum(dnorm(x = univContrasts[x,] / sqrt(reg*prunes[[2]][x]), mean = rep(0,d_trait), 
        #                                                       sd = 1, log = T)) - log(detR*(sqrt(reg*prunes[[2]][x])^length(univContrasts[x,])))/2))
        
      }
    } else {
      regTerm <- -sum(sapply(1:(d_trait-1), function(trait) mvtnorm::dmvnorm(x = c(par_to_opt, cont_traits_of_tip[trait]),
                                                                             sigma = matrix(c(1,corrs_of_trait[trait],corrs_of_trait[trait],1),2,2) * reg,
                                                                             mean = c(0,0), log = T)))
    }
    
  } else if(par_type == "reg"){
    logLL <- 0
    
    if(treeReg){
      # TODO make regularizing term multivariate normal
      contrasts <- prunes[[1]] %*% means[rownames(prunes[[3]])[1:nTaxa],]
      regTerm <- -sum(sapply(1:(nTaxa-1), function(x) dnorm(x = contrasts[x,], mean = rep(0,d_trait),
                                                            sd = sqrt(reg*prunes[[2]][x]), log = T)))
    } else {
      regTerm <- -sum(dnorm(x = means, sd = sqrt(reg), mean = 0, log = T)) 
    }
    
  } else if(par_type == "corr"){
    cont_traits <- means[, c(tr1, tr2)]
    corrs_threshes <- thresholds[c(tr1, tr2),]
     if(!vectorized){
       logLL <- -sum(unlist(sapply(1:nTaxa, function(tip) log(multi_pmvnorm_optim_multithresh(mean = cont_traits[tip,],
                                                                                                     sigma = matrix(c(1,par_to_opt,par_to_opt,1),2,2),
                                                                                                     ordinals = uniqueSP_redux(dat[dat$tip == tip, c(tr1, tr2)], dat[dat$tip == tip,d_trait+1])$SPs,
                                                                                                     thresholds = corrs_threshes)) * 
                                                                                                     uniqueSP_redux(dat[dat$tip == tip, c(tr1, tr2)], dat[dat$tip == tip,d_trait+1])$counts)))
      } else if(vectorized){ ## alternatively, in vectorized form
        
        if(!exists("uniqueSitePatterns_precomp")){
          uniqueSitePatterns <- sapply(1:nTaxa, function(tip)  uniqueSP_fast(dat[dat$tip == tip, c(tr1, tr2)], dat[dat$tip == tip,d_trait+1]))
        } else {
          uniqueSitePatterns <-uniqueSitePatterns_precomp
        }
        
        #convert to standard normal
        nOrds <- sapply(1:nTaxa, function(tip) nrow(uniqueSitePatterns[1,tip][[1]])) #number of unique ordinals per tip
        nOrds_tips <- rep(1:nTaxa, nOrds)
        new_std_thresholds <- lapply(1:nTaxa, function(tip) cbind(rep(-Inf, 2), sapply(1:ncol(corrs_threshes), function(col) corrs_threshes[,col] - cont_traits[tip,]), rep(Inf, 2)))
        SPs <- do.call(rbind, uniqueSitePatterns[1,])
        uppers_lowers <- t(sapply(1:length(nOrds_tips), function(pat) c(new_std_thresholds[[nOrds_tips[pat]]][1, SPs[pat,1]+c(2,1)], 
                                                                        new_std_thresholds[[nOrds_tips[pat]]][2, SPs[pat,2]+c(2,1)])))
        uppers <- uppers_lowers[,c(1,3)]
        lowers <- uppers_lowers[,c(2,4)]
        
        rhos = rep(par_to_opt, sum(nOrds))
        counts = unlist(uniqueSitePatterns[2,])
        probs <- pbivnorm_infs(uppers[,1], uppers[,2], rhos) - 
          pbivnorm_infs(lowers[,1], uppers[,2], rhos) - 
          pbivnorm_infs(uppers[,1], lowers[,2], rhos) + 
          pbivnorm_infs(lowers[,1], lowers[,2], rhos)
        probs[probs < 0] <- 0 #get rid of floating point errors
        logLL <- -sum(log(probs) * counts)
        
      }
    
    #FIX REDUNDANT UNIQUESP_REDUX
    regTerm <- 0
    
  } else if(par_type == "thresh"){
    # threshold_mat
    # thresh_ind 
    # thresh_trait
    cont_traits_of_thresh <- means[, thresh_trait]
    
    
    if(univThresh){
      logLL <- -sum(unlist(sapply(1:(nTaxa), function(tip) sum(log(multi_pmvnorm_optim_multithresh(mean = cont_traits_of_thresh[tip],
                                                                                                         sigma = 1,
                                                                                                         ordinals = as.integer(names(table(rep(dat[dat$tip == tip, thresh_trait], dat[dat$tip == tip, d_trait+1])))),
                                                                                                         thresholds = threshold_mat[thresh_trait,]))) * 
                                                                                                         as.vector(table(rep(dat[dat$tip == tip, thresh_trait], dat[dat$tip == tip, d_trait+1]))))))
      

    } else {
      corrs_of_thresh <- cor[thresh_trait, -thresh_trait]
      
      if(!vectorized){
        cont_traits_of_other_traits <- means[, -thresh_trait]
        threshs_of_other_traits <- threshold_mat[-thresh_trait,]
        logLL <- sum(sapply(1:nTaxa, function(tip) -sum(unlist(sapply(1:(d_trait-1), function(trait) log(multi_pmvnorm_optim_multithresh(mean = c(cont_traits_of_thresh[tip], cont_traits_of_other_traits[tip,trait]),
                                                                                                     sigma = matrix(c(1,corrs_of_thresh[trait],corrs_of_thresh[trait],1),2,2),
                                                                                                     ordinals = uniqueSP_redux(dat[dat$tip == tip, c(thresh_trait, (1:d_trait)[-thresh_trait][trait])], dat[dat$tip == tip,d_trait+1])$SPs,
                                                                                                     thresholds = rbind(threshold_mat[thresh_trait,], threshs_of_other_traits[trait,]))) * 
                                                                  uniqueSP_redux(dat[dat$tip == tip, c(thresh_trait, (1:d_trait)[-thresh_trait][trait])], dat[dat$tip == tip,d_trait+1])$counts)))))
      } else if(vectorized){ ## alternatively, in vectorized form
           
        if(!exists("uniqueSitePatterns_precomp")){
          uniqueSitePatterns <- lapply(1:nTaxa, function(tip) sapply(1:(d_trait-1), function(trait) uniqueSP_fast(dat[dat$tip == tip, c(thresh_trait, (1:d_trait)[-thresh_trait][trait])], dat[dat$tip == tip,d_trait+1])))
        } else {
          uniqueSitePatterns <-uniqueSitePatterns_precomp
        }
        
        #convert to standard normal
        nOrds <- lapply(1:nTaxa, function(tip) sapply(1:(d_trait-1), function(trait) nrow(uniqueSitePatterns[[tip]][1,trait][[1]]))) #number of unique ordinals per tip per trait match
        nOrds_traits <- lapply(1:nTaxa, function(tip) rep(1:(d_trait-1), nOrds[[tip]]))
        new_std_thresholds <- lapply(1:nTaxa, function(tip) cbind(rep(-Inf, d_trait), sapply(1:ncol(threshold_mat), function(col) threshold_mat[,col] - means[tip,]), rep(Inf, d_trait)))
        SPs <- lapply(1:nTaxa, function(tip) do.call(rbind, uniqueSitePatterns[[tip]][1,]))
        uppers_lowers <- lapply(1:nTaxa, function(tip) t(sapply(1:length(nOrds_traits[[tip]]), 
                                         function(pat) c(new_std_thresholds[[tip]][thresh_trait, SPs[[tip]][pat,1]+c(2,1)], 
                                                         new_std_thresholds[[tip]][nOrds_traits[[tip]][pat] + ifelse(nOrds_traits[[tip]][pat] < thresh_trait, 0, 1), SPs[[tip]][pat,2]+c(2,1)]))))
        uppers <- do.call(rbind, lapply(1:nTaxa, function(tip) uppers_lowers[[tip]][,c(1,3)]))
        lowers <- do.call(rbind, lapply(1:nTaxa, function(tip) uppers_lowers[[tip]][,c(2,4)]))
        
        rhos = corrs_of_thresh[unlist(nOrds_traits)]
        counts = unlist(lapply(1:nTaxa, function(tip) unlist(uniqueSitePatterns[[tip]][2,])))
        
        probs <- pbivnorm_infs(uppers[,1], uppers[,2], rhos) - 
          pbivnorm_infs(lowers[,1], uppers[,2], rhos) - 
          pbivnorm_infs(uppers[,1], lowers[,2], rhos) + 
          pbivnorm_infs(lowers[,1], lowers[,2], rhos)
        probs[probs < 0] <- 0 #get rid of floating point errors
        logLL <- -sum(log(probs) * counts)
        
      }
    }
    
    regTerm <- 0
  }

  returnVal <- logLL + regTerm
  if(returnVal == Inf){returnVal <- 1E10} #hacky temporary fix to make optim's generalized "L-BFGS-B" method stop complaining about singularities :/
  return(returnVal)
}

thresh_trait_ind <- sample(1:d_traits, 1)
par_to_opt_ind <- nTaxa*d_traits+2+choose(d_traits, 2) + 0:(max_thresholds_count-2) * d_traits + thresh_trait_ind - 1
par_to_opt_ind <- par_to_opt_ind[1:(n_thresholds[thresh_trait_ind]-1)]
# par <- c(c(t(traits)), 10, R[upper.tri(R)], c(threshold_diffs_init)) #try using true par values to check for good behavior
par_to_opt_ind <- 1
par_to_opt <- par[par_to_opt_ind]
type_of_param <- ifelse(par_to_opt_ind[1] <= nTaxa * d_traits, "m", 
                 ifelse(par_to_opt_ind[1] == (nTaxa * d_traits + 1), "v",
                 ifelse(par_to_opt_ind[1] > (nTaxa * d_traits + 1) & par_to_opt_ind[1] <= (nTaxa * d_traits + 1 + choose(d_traits, 2)), "c", "t")))
print(type_of_param)
k_par <- par
k_par[par_to_opt_ind] <- NA  
dat <- list(dat_orig, k_par, prunes_init, rownames(traits), n_thresh)
dat <- list(dat_orig, k_par, prunes_init, rownames(traits), n_thresh, getUniqueSitePatterns(par_to_opt, dat)) 
logLL_bivProb_optim_multithresh(par_to_opt, dat)
logLL_bivProb_optim_multithresh(0.95, dat)

# sapply(1:1000/100, function(asdf) logLL_bivProb_optim_multithresh(asdf, dat))

fast = T
method <- "L-BFGS-B"
if(!fast){
  t_s <- Sys.time()
  optim_out <- optim(par = par_to_opt, logLL_bivProb_optim_multithresh, dat = dat, method = method, control = list(trace = 0, REPORT = 10),
                     upper = c(rep(Inf, nTaxa * d_traits + 1), rep(1, choose(d_traits, 2)), rep(Inf, (max_thresholds_count-1) * d_traits))[par_to_opt_ind], 
                     lower = c(rep(-Inf, nTaxa * d_traits), 0, rep(-1, choose(d_traits, 2)), rep(0, (max_thresholds_count-1) * d_traits))[par_to_opt_ind])
  Sys.time() - t_s
  optim_out$par
  threshold_diffs[,thresh_trait_ind]
}

tic()
nrounds <- 5
ncore <- 1; mclapply_over_foreach <- T;
nparam <- length(par)
prunes <- prunes_init
for(i in 1:nrounds){
  cat(paste0("\n", "round ", i, ": "))
  param_ord <- as.list(sample(c(1:(nTaxa*d_traits+1+choose(d_traits, 2))))) #univ optim, double emphasis on correlation parameters and threshold diffs
  # param_ord <- as.list(sample(c(1:nparam))) #univ optim, double emphasis on correlation parameters and threshold diffs
  thresh_trait_ind <- sample(which(n_thresholds > 1))
  par_to_opt_ind <- lapply(thresh_trait_ind, function(trait) (nTaxa*d_traits+2+choose(d_traits, 2) + 0:(max_thresholds_count-2) * d_traits + trait - 1)[0:(n_thresholds[trait]-1)])
  if(T){
    param_ord <- c(par_to_opt_ind, param_ord)
  } else {
    param_ord <- sample(c(par_to_opt_ind, param_ord))
  }
  num_iters_passed <- 0
  if(ncore == 1){
    for(par_to_opt_ind in param_ord){
      type_of_param <- ifelse(par_to_opt_ind[1] <= nTaxa * d_traits, "m", 
                       ifelse(par_to_opt_ind[1] == (nTaxa * d_traits + 1), "v",
                       ifelse(par_to_opt_ind[1] > (nTaxa * d_traits + 1) & par_to_opt_ind[1] <= (nTaxa * d_traits + 1 + choose(d_traits, 2)), "c", 
                        "t")))
      cat(paste0("round ", i, ", ", round(num_iters_passed / (length(param_ord) / ncore) * 100), "%: ", type_of_param, "-"))
      num_iters_passed = num_iters_passed + 1
      # cat(paste0(par_to_opt_ind, type_of_param, "-"))
      par_to_opt <- par[par_to_opt_ind]
      k_par <- par
      k_par[par_to_opt_ind] <- NA  
      dat <- list(dat_orig, k_par, prunes, rownames(traits), n_thresh)
      dat <- list(dat_orig, k_par, prunes_init, rownames(traits), n_thresh, getUniqueSitePatterns(par_to_opt, dat)) 
      upper = c(rep(Inf, nTaxa * d_traits + 1), rep(0.99, choose(d_traits, 2)), rep(Inf, (max_thresholds_count-1) * d_traits))[par_to_opt_ind]
      lower = c(rep(-Inf, nTaxa * d_traits), 0.01, rep(-0.99, choose(d_traits, 2)), rep(0, (max_thresholds_count-1) * d_traits))[par_to_opt_ind]
      optim_out <- optim(par = par_to_opt, logLL_bivProb_optim_multithresh, dat = dat, method = method, control = list(trace = 0, REPORT = 1),
                         upper = upper, lower = lower)
      cat(paste0("(", round(mean(abs(par[par_to_opt_ind] - optim_out$par)), 2), ")-"))
      par[par_to_opt_ind] <- optim_out$par
    }
  } else {
    if(mclapply_over_foreach){
      param_ord <- c(param_ord, as.list(sample(1:nparam, size = ncore - length(param_ord) %% ncore, replace = F)))
      for(par_to_opt_ind in 1:(length(param_ord) / ncore)){
        mc_par_to_opt_ind <- (par_to_opt_ind-1)*ncore + 1:ncore
        mc_par_to_opt_ind <- param_ord[mc_par_to_opt_ind]
        type_of_param <- paste0(sapply(1:ncore, function(param_type) ifelse(mc_par_to_opt_ind[[param_type]][1] <= nTaxa * d_traits, "m",
                                ifelse(mc_par_to_opt_ind[[param_type]][1] == (nTaxa * d_traits + 1), "v",
                                       ifelse(mc_par_to_opt_ind[[param_type]][1] > (nTaxa * d_traits + 1) & mc_par_to_opt_ind[[param_type]][1] <= (nTaxa * d_traits + 1 + choose(d_traits, 2)), "c",
                                              "t")))), collapse = "-")
        cat(paste0("round ", i, ", ", round(par_to_opt_ind / (length(param_ord) / ncore) * 100), "%: ", type_of_param, "-"))
        par_to_opt <- lapply(1:ncore, function(core) par[mc_par_to_opt_ind[[core]]])
        k_par <- lapply(1:ncore, function(core) par)
        for(core in 1:ncore){
          k_par[[core]][mc_par_to_opt_ind[[core]]] <- NA
        }
        dat <- lapply(1:ncore, function(core) list(dat_orig, k_par[[core]], prunes, rownames(traits), n_thresh))
        upper = c(rep(Inf, nTaxa * d_traits + 1), rep(0.99, choose(d_traits, 2)), rep(Inf, (max_thresholds_count-1) * d_traits))
        upper = lapply(1:ncore, function(core) upper[mc_par_to_opt_ind[[core]]])
        lower = c(rep(-Inf, nTaxa * d_traits), 0.01, rep(-0.99, choose(d_traits, 2)), rep(0, (max_thresholds_count-1) * d_traits))
        lower = lapply(1:ncore, function(core) lower[mc_par_to_opt_ind[[core]]])
        optim_out <- mclapply(1:ncore, function(core) optim(par = par_to_opt[[core]], logLL_bivProb_optim_multithresh, dat = dat[[core]],
                           method = method, control = list(trace = 0, REPORT = 1),
                           upper = upper[[core]], lower = lower[[core]]), mc.preschedule = F, mc.cores = ncore)
        change <- mean(sapply(1:ncore, function(core) round(mean(abs(par[mc_par_to_opt_ind[[core]]] - optim_out[[core]]$par)), 2)))
        cat(paste0("(", round(change,3), ")-"))
        for(core in 1:ncore){
          par[mc_par_to_opt_ind[[core]]] <- optim_out[[core]]$par
        }
      }
    } else {
      fwrite(list(par), file = "probit_params.txt")
  
      cl <- makeCluster(8, outfile="")
      registerDoParallel(cl)
      getDoParWorkers()
      
      foreach(m=1:length(param_ord), .packages = c("pbivnorm", "Matrix", "mvtnorm", "phangorn", "parallel", "data.table")) %dopar% {
        par <- fread("probit_params.txt")$V1
        par_to_opt_ind <- param_ord[[m]]
        type_of_param <- ifelse(par_to_opt_ind[1] <= nTaxa * d_traits, "m", 
                                ifelse(par_to_opt_ind[1] == (nTaxa * d_traits + 1), "v",
                                       ifelse(par_to_opt_ind[1] > (nTaxa * d_traits + 1) & par_to_opt_ind[1] <= (nTaxa * d_traits + 1 + choose(d_traits, 2)), "c", 
                                              "t")))
        cat(paste0(par_to_opt_ind, type_of_param, "-"))
        par_to_opt <- par[par_to_opt_ind]
        k_par <- par
        k_par[par_to_opt_ind] <- NA  
        dat <- list(dat_orig, k_par, prunes, rownames(traits), n_thresh)
        upper = c(rep(Inf, nTaxa * d_traits + 1), rep(0.999, choose(d_traits, 2)), rep(Inf, (max_thresholds_count-1) * d_traits))[par_to_opt_ind]
        lower = c(rep(-Inf, nTaxa * d_traits), 0.01, rep(-0.999, choose(d_traits, 2)), rep(0, (max_thresholds_count-1) * d_traits))[par_to_opt_ind]
        optim_out <- optim(par = par_to_opt, logLL_bivProb_optim_multithresh, dat = dat, method = method, control = list(trace = 0, REPORT = 1),
                           upper = upper, lower = lower)
        cat(paste0("(", round(mean(abs(par[par_to_opt_ind] - optim_out$par)), 2), ")-"))
        
        par <- fread("probit_params.txt")$V1
        par[par_to_opt_ind] <- optim_out$par
        
        fwrite(list(par), file = "probit_params.txt")
      }
    }
    stopCluster(cl)
    
  }
  
  
  if(ncore > 1){
    par <- fread("probit_params.txt")$V1
  }

  #update tree
  optim_est <- par
  R_est <- diag(d_traits)
  R_est[upper.tri(R_est)] <- optim_est[(nTaxa*d_traits+2) : (nTaxa * d_traits + 1 + choose(d_traits, 2))];
  R_est <- as.matrix(nearPD(R_est + t(R_est) - diag(d_traits), corr = T)$mat)
  est_means <- matrix(optim_est[1:(nTaxa*d_traits)], ncol = d_traits, nrow = nTaxa, byrow = T, dimnames = list(rownames(traits)))
  
  mahDistMat <- mahaMatrix(est_means, R_est, squared = F)
  upgmaTree <- upgma(mahDistMat)
  treeCV <- cov2cor(vcv.phylo(upgmaTree))
  if(weighPCV){treeCV <- (treeCV + diag(nTaxa)) / 2} #weight with star phylogeny
  prunes <- prunePCV(treeCV)
  
}
toc()

if(ncore > 1){
  par <- fread("probit_params.txt")$V1
}  
optim_est <- par

R_est <- diag(d_traits)
R_est[upper.tri(R_est)] <- optim_est[(nTaxa*d_traits+2) : (nTaxa * d_traits + 1 + choose(d_traits, 2))];
R_est <- as.matrix(nearPD(R_est + t(R_est) - diag(d_traits), corr = T)$mat)

reg <- optim_est[nTaxa * d_traits + 1]
paste0("regularizing variance = ", round(reg, 3))

est_means <- matrix(optim_est[1:(nTaxa*d_traits)], ncol = d_traits, nrow = nTaxa, byrow = T)

thresholds_est <- cbind(rep(0, d_traits), matrix(optim_est[(nTaxa*d_traits+2+choose(d_traits, 2)) : ((nTaxa*d_traits+1+choose(d_traits, 2)) + sum(d_traits*(max_thresholds_count-1)))], nrow = d_traits, ncol = max_thresholds_count[1]-1))
thresholds_est <- t(sapply(1:nrow(thresholds_est), function(trait) cumsum(thresholds_est[trait,])))

par(mfrow = c(2,3))
plot(init_means, traits, main = paste0("r = ", round(cor(c(traits), c(init_means)), 3))); abline(0,1)
plot(init_cor[upper.tri(R)], R[upper.tri(R)], xlim = c(-1,1), ylim = c(-1,1), xlab = "Pearson's 'r's", ylab = "true correlations",
     main = paste0("r = ", round(cor(c(R[upper.tri(R)]), c(init_cor[upper.tri(R)])), 3))); abline(0,1)
plot.new()
plot(est_means, traits, main = paste0("r = ", round(cor(c(traits), c(est_means)), 3))); abline(0,1)
plot(R_est[upper.tri(R)], R[upper.tri(R)], xlim = c(-1,1), ylim = c(-1,1), xlab = "estimated correlations", ylab = "true correlations",
     main = paste0("r = ", round(cor(c(R[upper.tri(R)]), c(R_est[upper.tri(R)])), 3))); abline(0,1)
thresholds_est_vec <- as.vector(thresholds_est[,-1]); thresholds_est_vec <- thresholds_est_vec[thresholds_est_vec != Inf]
thresholds_vec <- as.vector(thresholds[,-1]); thresholds_vec <- thresholds_vec[thresholds_vec != Inf]
plot(thresholds_est_vec, thresholds_vec, main = paste0("r = ", round(cor(c(thresholds_vec), c(thresholds_est_vec)), 3))); abline(0,1)

#use linear model but have it be on the lognormal scale

par(mfrow = c(1,1))
hist((abs(est_means - traits)))
deviants <- as.matrix(table(as.data.frame(which((abs(est_means - traits)) > 1, arr.ind = T))))
deviants
apply(deviants, 2, sum)
plot(tree)
plot(upgmaTree)

invariants <- as.matrix(table(as.data.frame(which(t(sapply(1:nTaxa, function(tip) apply(traits_indiv_discr_unique[[tip]]$SPs, 2, var) == 0)), arr.ind = T))))
invariants
deviants
par(mfrow = c(1,2)); plot(deviants); plot(invariants)

deviants <- cbind(deviants, matrix(0,nrow(deviants), ncol = length(setdiff(colnames(invariants), colnames(deviants))), dimnames = list(rownames(deviants), setdiff(colnames(invariants), colnames(deviants)))))
deviants <- deviants[,order(as.numeric(colnames(deviants)))]
deviants <- rbind(deviants, matrix(0,ncol(deviants), nrow = length(setdiff(rownames(invariants), rownames(deviants))), dimnames = list(setdiff(rownames(invariants), rownames(deviants)), colnames(deviants))))
deviants <- deviants[order(as.numeric(rownames(deviants))),]
invariants - deviants
sum(invariants)
