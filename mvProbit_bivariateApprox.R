library(pbivnorm)
library(Matrix)
library(mvtnorm)
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

uniqueSP <- function(indivs){
  SPs <- table(sapply(1:length(indivs[,1]), function(indiv) paste0(indivs[indiv,], collapse = "")))
  counts <- as.numeric(SPs)
  SPs <- t(sapply(1:length(SPs), function(SP) as.numeric(unlist(strsplit(names(SPs)[SP], "")))))
  return(list(SPs = SPs, counts = counts))
}

#first let's resimulate data
nTaxa <- 5
d_traits <- 8
tree <- tess.sim.taxa(n = 1, nTaxa = nTaxa, lambda = 1, mu = 0, max = 1E3)[[1]]
target_tree_height <- 2
tree$edge.length <- tree$edge.length * target_tree_height / node.depth.edgelength(tree)[1]

R <- rlkj(d_traits); rownames(R) <- colnames(R) <- paste0("tr", 1:dim(R)[1])
root_state <- rep(0, dim(R)[1])
root_state <- rnorm(d_traits)
traits <- mvSIM(tree, model="BM1", param=list(sigma=R, theta=root_state)) 

#simulate thresholded traits at the population level
n_indiv <- 50
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
plot(init_means, traits); abline(0,1)
plot(init_cor, R); abline(0,1)
par <- c(c(t(init_means)), 10, init_cor[upper.tri(init_cor)])


par_to_opt_ind <- 50
par_to_opt <- par[par_to_opt_ind]
k_par <- par
k_par[par_to_opt_ind] <- NA  
dat <- list(dat, k_par)
  
logLL_bivProb_optim_single <- function(par_to_opt, dat){
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
  
  if(par_to_opt_ind > 0 & par_to_opt_ind <= (nTaxa * d_traits)){
    par_type <- "mean"
    mean_trait <- par_to_opt_ind %% d_trait
    mean_tip <- round(par_to_opt_ind / d_trait)
  } else if (par_to_opt_ind == (nTaxa * d_traits + 1)){
    par_type = "reg"  
  } else {
    par_type = "corr"
    trs <- which(is.na(cor), arr.ind = T)
    tr1 = trs[1]
    tr2 = trs[2]
    cor[tr1, tr2] <- par_to_opt
  }
  
  cor <- cor + t(cor) - diag(d_trait)
  
  if(par_type == "mean"){
    #check biv dens of trait mean with all others traits for that tip
    cont_traits_of_tip <- means[mean_tip, -mean_trait]
    corrs_of_trait <- cor[mean_trait, -mean_trait]
    logLL <- -sum(sapply(1:(d_trait-1), function(trait) log(multi_pmvnorm_optim(mean = c(par_to_opt, cont_traits_of_tip[trait]),
                                                                           sigma = matrix(c(1,corrs_of_trait[trait],corrs_of_trait[trait],1),2,2),
                                                                           binaries = dat[dat$tip == mean_tip, c(tr1, tr2)]))) * dat[dat$tip == mean_tip,d_trait+1])
    regTerm <- -sum(sapply(1:(d_trait-1), function(trait) mvtnorm::dmvnorm(x = c(par_to_opt, cont_traits_of_tip[trait]),
                                                                          sigma = matrix(c(1,corrs_of_trait[trait],corrs_of_trait[trait],1),2,2) * reg,
                                                                          mean = c(0,0), log = T)))
  } else if(par_type == "reg"){
    logLL <- 0
    regTerm <- -sum(dnorm(x = means, sd = sqrt(par_to_opt), mean = 0, log = T))
    
  } else if(par_type == "cor"){
    cont_traits <- means[, c(tr1, tr2)]
    logLL <- -sum(sapply(1:nTaxa, function(tip) log(multi_pmvnorm_optim(mean = cont_traits[tip,],
                                                                                sigma = matrix(c(1,par_to_opt,par_to_opt,1),2,2),
                                                                                binaries = dat[dat$tip == tip, c(tr1, tr2)]))) * dat[dat$tip == tip,d_trait+1])
    ## NEED NEW UNIQUE SITE PATTERNS FOR BINARY TRAIT, SINCE OLD UNIQUE SITE PATTERNS ARE FOR THE WHOLE SHEBANG
    regTerm <- -sum(sapply(1:(d_trait-1), function(trait) mvtnorm::dmvnorm(x = c(par_to_opt, cont_traits_of_tip[trait]),
                                                                           sigma = matrix(c(1,corrs_of_trait[trait],corrs_of_trait[trait],1),2,2) * reg,
                                                                           mean = c(0,0), log = T)))
  }
  
  logLL <- sum(sapply(1:(d_trait-1), function(rowi) sum(sapply((rowi+1):d_trait, 
                                                               function(colj) -sum(sapply(1:nTaxa, function(tip) sum(log(multi_pmvnorm_optim(mean = means[tip,c(rowi,colj)], sigma = matrix(c(1,cor[rowi,colj],cor[rowi,colj],1),2,2), 
                                                                                                                                             binaries = dat[dat$tip == tip, c(rowi, colj)])) * dat[dat$tip == tip,d_trait+1])))))))
  # regTerm <- - sum(dnorm(x = as.vector(par), mean = 0, sd = sqrt(reg), log = T))
  regTerm <- -sum(sapply(1:dim(means)[1], function(tip) dmvnorm(x = means[tip,], (as.matrix(nearPD(cor, corr = T)$mat)) * reg, mean = rep(0, d_trait), log = T)))
  returnVal <- logLL + regTerm
  if(returnVal == Inf){returnVal <- 1E10} #hacky temporary fix to make optim's generalized "L-BFGS-B" method stop complaining about singularities :/
  return(returnVal)
}



method <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent")[4]

t_s <- Sys.time()
optim_out <- optim(par = par, logLL_bivProb_optim_joint_single, dat = dat_a, method = method, control = list(trace = 3, REPORT = 1),
                   upper = c(rep(Inf, nTaxa * d_traits + 1), rep(1, choose(d_traits, 2))), 
                   lower = c(rep(-Inf, nTaxa * d_traits), 0, rep(-1, choose(d_traits, 2))))
Sys.time() - t_s

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
