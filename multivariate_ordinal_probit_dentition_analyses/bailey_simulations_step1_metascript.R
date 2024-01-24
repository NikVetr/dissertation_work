rep_i <- 1
optim_rep_i <- 1

nrounds <- 12
noSmartStart = F
useStarPhylogeny <- T
doEDA <- F
jointThreshAndMeans <- F; jTM_startAtRound <- 6
useUPGMAforNRounds <- 0
impose_constraints <- T 
weighPCV <- F
MCAR <- F #missing data is coded with state '9'
ignoreImputedData <- F
nrounds_to_not_impute <- 4
stochastic_imputation <- T
MNAR_Imputation <- T
updateAllAtOnce <- T #impute each missing value as you iterate through them, or all at once at the end?
ignoreOtherImputedStates <- T
use_conditional_probs_for_observed_data <- F
use_conditional_probs_for_unobserved_data <- T
use_imputed_data_for_corrs <- F
useUPGMA <- F
use_univariate_normals_for_conditionals <- T
use_observed_numerator_ALL_conditionals_denominator <- T
use_conditional_multinomial_approx <- T
add_individuals_to_missing <- F
adjust_conds_at_pop_level <- F
use_conditional_expectation <- F
useNearPDforTree <- F

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
library(gtools)
library(MCMCpack)
library(microbenchmark)
library(phytools)
library(phangorn)
library(mvMORPH)
library(tictoc)
library(mvtnorm)

logLL_bivProb_optim_multithresh <- function(par_to_opt, dat){ #TMB package may help with speedup
  lognormal_prior_on_beta_shape_params_for_corrs <- c(0,1) #let's not infer the shape of the beta distribution,
  # but just fix it to the appropriate beta of the marginal LKJ -- ehhh this is too informative
  disp_lognormal <- 1
  vectorized <- T
  univMeans <- T
  univThresh <- T
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
  
  #for lognormal regularization
  # thresh_reg_mu <- par[(nTaxa*d_trait+1+choose(d_trait, 2)) + sum(n_thresh-1) + 1]
  # thresh_reg_sd <- par[(nTaxa*d_trait+1+choose(d_trait, 2)) + sum(n_thresh-1) + 2]
  
  #for exponential regularization
  lambda <- par[(nTaxa*d_trait+1+choose(d_trait, 2)) + sum(n_thresh-1) + 1]
  
  #for regularizing the correlations
  cor_beta_shape_params <- par[nTaxa * d_trait + 1 + choose(d_trait, 2) + sum(n_thresh-1) + c(2,3)] 
  
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
  } else if(any(par_to_opt_ind > (nTaxa*d_trait+1+choose(d_trait, 2)) & par_to_opt_ind <= ((nTaxa*d_trait+1+choose(d_trait, 2)) + sum(n_thresh-1)))){
    par_type <- "thresh-means"
    threshold_mat_unknown <- cbind(rep(0, d_trait), matrix(k_par[(nTaxa*d_trait+2+choose(d_trait, 2)) : ((nTaxa*d_trait+1+choose(d_trait, 2)) + sum(n_thresh-1))], nrow = d_trait, ncol = n_thresh[1]-1))
    thresh_trait <- which(is.na(threshold_mat_unknown), arr.ind = T)
    thresh_ind <- thresh_trait[,2]
    thresh_trait <- thresh_trait[1,1]
  }
  
  if(length(intersect(par_to_opt_ind, (((nTaxa*d_trait+1+choose(d_trait, 2)) + sum(n_thresh-1)) + c(1,2)))) > 0){
    par_type <- "thresh_reg"
  }
  
  if(length(intersect(par_to_opt_ind, (((nTaxa*d_trait+1+choose(d_trait, 2)) + sum(n_thresh-1)) + 3))) > 0){
    par_type <- "corr_reg"
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
    
    if(univMeans){
      ordinals <- as.data.frame(wtd.table(dat[dat$tip == mean_tip,mean_trait], weights = dat$nSP[dat$tip == mean_tip]))
      stateProbs <- log(computeUnivNormalStateProbs(thresholds = threshold_mat[mean_trait,], mean = means[mean_tip, mean_trait]))
      logLL <- -sum(stateProbs[as.integer(as.character(ordinals[,1])) + 1] * ordinals[,2])
    }
    
    if(!vectorized & !univMeans){
      
      logLL <- -sum(unlist(sapply(1:(d_trait-1), function(trait) log(multi_pmvnorm_optim_multithresh(mean = c(par_to_opt, cont_traits_of_tip[trait]),
                                                                                                     sigma = matrix(c(1,corrs_of_trait[trait],corrs_of_trait[trait],1),2,2),
                                                                                                     ordinals = uniqueSP_redux(dat[dat$tip == mean_tip, c(mean_trait, (1:d_trait)[-mean_trait][trait])], dat[dat$tip == mean_tip,d_trait+1])$SPs,
                                                                                                     thresholds = rbind(threshold_mat[mean_trait,], threshs_of_other_traits[trait,]))) * 
                                    uniqueSP_redux(dat[dat$tip == mean_tip, c(mean_trait, (1:d_trait)[-mean_trait][trait])], dat[dat$tip == mean_tip,d_trait+1])$counts)))
    } else if(vectorized & !univMeans){ ## alternatively, in vectorized form
      
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
        # regTerm <- -sum(sapply(1:(nTaxa-1), function(x) dnorm(x = contrasts[x,], mean = rep(0,d_trait), #mean 0 is marginalizes over all mean states
        #                                                       sd = sqrt(reg*prunes[[2]][x]), log = T)))
        
        regTerm <- -sum(sapply(1:(nTaxa-1), function(x) dnorm(x = contrasts[x,mean_trait], mean = 0, #mean 0 is marginalizes over all mean states
                                                              sd = sqrt(reg*prunes[[2]][x]), log = T)))
        
      } else if(mvBMreg){
        pairs <- cbind(rep(mean_trait, d_trait-1), (1:d_trait)[-mean_trait])
        regTerm <- sum(sapply(1:nrow(pairs), function(pair) -sum(sapply(1:(nTaxa-1), function(x) mvtnorm::dmvnorm(x = contrasts[x,pairs[pair,]], mean = rep(0,2), #mean 0 is marginalizes over all mean states
                                                                                                                  sigma = cor[pairs[pair,], pairs[pair,]] * (reg*prunes[[2]][x]), log = T)))))
        
        
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
      if(!mvBMreg){
        contrasts <- prunes[[1]] %*% means[rownames(prunes[[3]])[1:nTaxa],]
        regTerm <- -sum(sapply(1:(nTaxa-1), function(x) dnorm(x = contrasts[x,], mean = rep(0,d_trait),
                                                              sd = sqrt(reg*prunes[[2]][x]), log = T)))
      } else if(mvBMreg){
        pairs <- as.matrix(which(upper.tri(cor), arr.ind = T))
        regTerm <- sum(sapply(1:nrow(pairs), function(pair) -sum(sapply(1:(nTaxa-1), function(x) mvtnorm::dmvnorm(x = contrasts[x,pairs[pair,]], mean = rep(0,2), 
                                                                                                                  sigma = cor[pairs[pair,], pairs[pair,]] * (reg*prunes[[2]][x]), log = T)))))
      }
    } else {
      regTerm <- -sum(dnorm(x = means, sd = sqrt(reg), mean = 0, log = T)) 
    }
    
  } else if(par_type == "corr"){
    cont_traits <- means[, c(tr1, tr2)]
    corrs_threshes <- threshold_mat[c(tr1, tr2),]
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
    
    regTerm <- -(dbeta(x = (par_to_opt + 1) / 2, shape1 = cor_beta_shape_params[1], shape2 =  cor_beta_shape_params[2], log = TRUE))
    
  } else if(par_type == "thresh" | par_type == "thresh-means"){
    # threshold_mat
    # thresh_ind 
    # thresh_trait
    cont_traits_of_thresh <- means[, thresh_trait]
    
    
    if(univThresh){
      
      ordinals <- lapply(1:nTaxa, function(tip) as.data.frame(wtd.table(dat[dat$tip == tip, thresh_trait], weights = dat$nSP[dat$tip == tip])))
      stateProbs <- lapply(1:nTaxa, function(tip) log(computeUnivNormalStateProbs(thresholds = threshold_mat[thresh_trait,], 
                                                                                  mean = means[tip, thresh_trait])))
      logLL <- -sum(sapply(1:nTaxa, function(tip) sum((stateProbs[[tip]][as.integer(as.character(ordinals[[tip]][,1]))+1]) * ordinals[[tip]][,2])))
      
      
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
    
    # regTerm <- 0
    
    #introduce regularizing term
    #assume joint optimizatiton of thresholds separations
    if(par_type == "thresh"){
      if(abs(sum(par_to_opt) - sum(diff(threshold_mat[thresh_trait,][threshold_mat[thresh_trait,] != Inf]))) > 1E-6){
        if(length(par_to_opt) != (length(threshold_mat[thresh_trait,][threshold_mat[thresh_trait,] != Inf])-1)){
          stop("something has gone terribly wrong")
        } else {
          sum_of_spacings <- sum(diff(threshold_mat[thresh_trait,][threshold_mat[thresh_trait,] != Inf]))
          regTerm <- -sum(dexp(diff(threshold_mat[thresh_trait,][threshold_mat[thresh_trait,] != Inf]), lambda, TRUE))
        } 
      } else {
        sum_of_spacings <- sum(par_to_opt)
        regTerm <- -sum(dexp(par_to_opt, lambda, TRUE))
        
      }
    } else if(par_type == "thresh-means"){
      if(abs(sum(par_to_opt[-(1:nTaxa)]) - sum(diff(threshold_mat[thresh_trait,][threshold_mat[thresh_trait,] != Inf]))) > 1E-6){
        sum_of_spacings <- sum(diff(threshold_mat[thresh_trait,][threshold_mat[thresh_trait,] != Inf]))
        regTerm <- -sum(dexp(diff(threshold_mat[thresh_trait,][threshold_mat[thresh_trait,] != Inf]), lambda, TRUE))
      } else {
        sum_of_spacings <- sum(par_to_opt[-(1:nTaxa)])
        regTerm <- -sum(dexp(par_to_opt[-(1:nTaxa)], lambda, TRUE))
        
      }
    }
    # regTerm <- -dnorm(x = log(sum_of_spacings), mean = thresh_reg_mu, sd = thresh_reg_sd, log = T) #lognormal regularization
    
  } else if(par_type == "thresh_reg"){
    logLL <- 0
    thresh_spacings <- unlist(lapply(1:d_trait, function(trait) diff(threshold_mat[trait,][threshold_mat[trait,] != Inf])))
    thresh_spacings <- thresh_spacings[thresh_spacings > 1E-6]
    regTerm <- -sum(dexp(thresh_spacings, lambda, TRUE))
    
    
    #lognormal regularization
    # thresh_sums_of_spacings <- sapply(1:d_trait, function(trait) sum(diff(threshold_mat[trait,][threshold_mat[trait,] != Inf])))
    # thresh_sums_of_spacings <- thresh_sums_of_spacings[thresh_sums_of_spacings > 1E-6]
    # regTerm <- -sum(dnorm(x = log(thresh_sums_of_spacings), mean = thresh_reg_mu, sd = thresh_reg_sd, log = T))
    
  } else if(par_type == "corr_reg"){
    all_the_cors <- cor[upper.tri(cor)]
    all_the_cors <- (all_the_cors + 1) / 2
    logLL <- -sum(dbeta(x = all_the_cors, shape1 = cor_beta_shape_params[1], shape2 =  cor_beta_shape_params[2], log = TRUE))
    regTerm <- -sum(dnorm(x = log(cor_beta_shape_params-disp_lognormal), mean = lognormal_prior_on_beta_shape_params_for_corrs[1], sd = lognormal_prior_on_beta_shape_params_for_corrs[2], log = T))
    # logLL <- 0
    # regTerm <- 0
  }
  
  returnVal <- logLL + regTerm
  if(returnVal == Inf){returnVal <- 1E10} #hacky temporary fix to make optim's generalized "L-BFGS-B" method stop complaining about singularities :/
  return(returnVal)
}
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
uniqueSP <- function(indivs){
  SPs <- table(sapply(1:length(indivs[,1]), function(indiv) paste0(indivs[indiv,], collapse = "")))
  counts <- as.numeric(SPs)
  SPs <- t(sapply(1:length(SPs), function(SP) as.numeric(unlist(strsplit(names(SPs)[SP], "")))))
  return(list(SPs = SPs, counts = counts))
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
} #checks for buggy edge cases and corrects
uniqueSP_fast <- function(indivs, freqs){
  newFreqs <- as.data.frame(wtd.table(indivs[,1], indivs[,2], weights = freqs))
  newFreqs <- (newFreqs[newFreqs$Freq > 0,])
  newFreqs$Var1 <- as.numeric(levels(newFreqs$Var1))[newFreqs$Var1]
  newFreqs$Var2 <- as.numeric(levels(newFreqs$Var2))[newFreqs$Var2]
  return(list(SPs = as.matrix(newFreqs[,1:2]), counts = newFreqs[,3]))
}
hasNine <- function(x){any(x==9)}
uniqueSP_fast_ignoreImputed <- function(indivs, freqs){
  newFreqs <- as.data.frame(wtd.table(indivs[,1], indivs[,2], weights = freqs))
  newFreqs <- (newFreqs[newFreqs$Freq > 0,])
  newFreqs$Var1 <- as.numeric(levels(newFreqs$Var1))[newFreqs$Var1]
  newFreqs$Var2 <- as.numeric(levels(newFreqs$Var2))[newFreqs$Var2]
  newFreqs <- newFreqs[!apply(newFreqs[,1:2], 1, hasNine),]
  return(list(SPs = as.matrix(newFreqs[,1:2]), counts = newFreqs[,3]))
}
getUniqueSitePatterns <- function(par_to_opt, dat){ #precompute unique site patterns
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
  
  
  uniqueSitePatterns <- NA
  
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
  } else if(any(par_to_opt_ind > (nTaxa*d_trait+1+choose(d_trait, 2)) & par_to_opt_ind <= ((nTaxa*d_trait+1+choose(d_trait, 2)) + sum(n_thresh-1)))){
    par_type <- "thresh-means"
    threshold_mat_unknown <- cbind(rep(0, d_trait), matrix(k_par[(nTaxa*d_trait+2+choose(d_trait, 2)) : ((nTaxa*d_trait+1+choose(d_trait, 2)) + sum(n_thresh-1))], nrow = d_trait, ncol = n_thresh[1]-1))
    thresh_trait <- which(is.na(threshold_mat_unknown), arr.ind = T)
    thresh_ind <- thresh_trait[,2]
    thresh_trait <- thresh_trait[1,1]
  }
  if(length(intersect(par_to_opt_ind, (((nTaxa*d_trait+1+choose(d_trait, 2)) + sum(n_thresh-1)) + c(1,2)))) > 0){
    par_type <- "thresh_reg"
  }
  
  cor <- cor + t(cor) - diag(d_trait)
  threshold_mat <- t(sapply(1:nrow(threshold_mat), function(trait) cumsum(threshold_mat[trait,])))
  
  if(par_type == "mean"){
    uniqueSitePatterns <- sapply(1:(d_trait-1), function(trait) uniqueSP_fast(dat[dat$tip == mean_tip, c(mean_trait, (1:d_trait)[-mean_trait][trait])], dat[dat$tip == mean_tip,d_trait+1]))
  } else if(par_type == "corr"){
    uniqueSitePatterns <- sapply(1:nTaxa, function(tip)  uniqueSP_fast(dat[dat$tip == tip, c(tr1, tr2)], dat[dat$tip == tip,d_trait+1]))
  }  else if(par_type == "thresh" | par_type == "thresh-means"){
    uniqueSitePatterns <- lapply(1:nTaxa, function(tip) sapply(1:(d_trait-1), function(trait) uniqueSP_fast(dat[dat$tip == tip, c(thresh_trait, (1:d_trait)[-thresh_trait][trait])], dat[dat$tip == tip,d_trait+1])))
  }
  return(uniqueSitePatterns)
}
getUniqueSitePatternsIgnoreImputed <- function(par_to_opt, dat){ #precompute unique site patterns
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
  
  uniqueSitePatterns <- NA
  
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
  } else if(any(par_to_opt_ind > (nTaxa*d_trait+1+choose(d_trait, 2)) & par_to_opt_ind <= ((nTaxa*d_trait+1+choose(d_trait, 2)) + sum(n_thresh-1)))){
    par_type <- "thresh-means"
    threshold_mat_unknown <- cbind(rep(0, d_trait), matrix(k_par[(nTaxa*d_trait+2+choose(d_trait, 2)) : ((nTaxa*d_trait+1+choose(d_trait, 2)) + sum(n_thresh-1))], nrow = d_trait, ncol = n_thresh[1]-1))
    thresh_trait <- which(is.na(threshold_mat_unknown), arr.ind = T)
    thresh_ind <- thresh_trait[,2]
    thresh_trait <- thresh_trait[1,1]
  }
  
  if(length(intersect(par_to_opt_ind, (((nTaxa*d_trait+1+choose(d_trait, 2)) + sum(n_thresh-1)) + c(1,2)))) > 0){
    par_type <- "thresh_reg"
  }
  
  cor <- cor + t(cor) - diag(d_trait)
  threshold_mat <- t(sapply(1:nrow(threshold_mat), function(trait) cumsum(threshold_mat[trait,])))
  
  if(par_type == "mean"){
    uniqueSitePatterns <- sapply(1:(d_trait-1), function(trait) uniqueSP_fast_ignoreImputed(dat[dat$tip == mean_tip, c(mean_trait, (1:d_trait)[-mean_trait][trait])], dat[dat$tip == mean_tip,d_trait+1]))
  } else if(par_type == "corr"){
    uniqueSitePatterns <- sapply(1:nTaxa, function(tip)  uniqueSP_fast_ignoreImputed(dat[dat$tip == tip, c(tr1, tr2)], dat[dat$tip == tip,d_trait+1]))
  }  else if(par_type == "thresh" | par_type == "thresh-means"){
    uniqueSitePatterns <- lapply(1:nTaxa, function(tip) sapply(1:(d_trait-1), function(trait) uniqueSP_fast_ignoreImputed(dat[dat$tip == tip, c(thresh_trait, (1:d_trait)[-thresh_trait][trait])], dat[dat$tip == tip,d_trait+1])))
  }
  return(uniqueSitePatterns)
}
recompileMeansCorrsThreshes <- function(par, dat){
  
  n_thresh <- dat[[5]]
  tipNames <- dat[[4]]
  prunes <- dat[[3]]
  k_par <- dat[[2]]
  dat <- dat[[1]]
  
  hb <- length(par)
  nTaxa <- length(unique(dat$tip))
  d_trait <- ncol(dat) - 2
  cor <- diag(d_trait)
  cor[upper.tri(cor)] <- k_par[(nTaxa*d_trait+2) : (nTaxa*d_trait+1+choose(d_trait, 2))];
  reg <- par[nTaxa * d_trait + 1]
  means <- matrix(par[1:(nTaxa*d_trait)], ncol = d_trait, nrow = nTaxa, byrow = T)
  rownames(means) <- tipNames
  threshold_mat <- cbind(rep(0, d_trait), matrix(par[(nTaxa*d_trait+2+choose(d_trait, 2)) : ((nTaxa*d_trait+1+choose(d_trait, 2)) + sum(n_thresh-1))], nrow = d_trait, ncol = n_thresh[1]-1))
  cor <- cor + t(cor) - diag(d_trait)
  threshold_mat <- t(sapply(1:nrow(threshold_mat), function(trait) cumsum(threshold_mat[trait,])))
  
  return(list(means = means, cor = cor, threshold_mat = threshold_mat))
  
}
invlogit <- function(x){exp(x) / (1 + exp(x))}
logit <- function(p){log(p / (1 - p))}
softmax <- function(x) exp(x) / sum(exp(x))
allEqualNinf <- function(nums){
  return(all(nums == -Inf))
}
evaluateStates <- function(missinds, possible_states, tipMeans, otherStates, other_missing_traits_in_individual, threshold_mat, cor, dimAdj = T, ignoreOtherImputedStates = T){
  #1st index is person, 2nd index of missinds is trait
  
  #if ignoring other imputed states, change the other states to be correct
  if(ignoreOtherImputedStates & (length(other_missing_traits_in_individual) > 0)){
    if(length(other_missing_traits_in_individual) + 1 == length(tipMeans)){
      return(rep(log(1 / length(possible_states)), length(possible_states)))
    }
    missinds[2] <- missinds[2] - sum(missinds[2] > other_missing_traits_in_individual)
    tipMeans <- tipMeans[-other_missing_traits_in_individual]
    otherStates <- otherStates[,-other_missing_traits_in_individual]
    threshold_mat <- threshold_mat[-other_missing_traits_in_individual,]
    cor <- cor[-other_missing_traits_in_individual, -other_missing_traits_in_individual]
  }
  
  bivMeans <- cbind(rep(tipMeans[missinds[2]], length(tipMeans)-1), tipMeans[-missinds[2]])
  personStates <- cbind(as.vector(sapply(possible_states, function(state) rep(state, length(tipMeans)-1))),
                        rep(otherStates[missinds[1], -missinds[2]], length(possible_states)))
  cent_threshold_mat <- cbind(-Inf, threshold_mat - tipMeans, Inf)
  
  other_lower_inds <- other_upper_inds <- cbind((1:length(tipMeans))[-missinds[2]], otherStates[missinds[1], -missinds[2]]+2) #first column is trait ind, second is thresh ind
  other_lower_inds[,2] <- other_lower_inds[,2] - 1
  uppers <- cbind(cent_threshold_mat[missinds[2],][possible_states+2][personStates[,1]+1],
                  rep(sapply(1:(length(tipMeans)-1), function(trait) cent_threshold_mat[other_upper_inds[trait,1], other_upper_inds[trait,2]]), length(possible_states)))
  lowers <- cbind(cent_threshold_mat[missinds[2],][possible_states+1][personStates[,1]+1],
                  rep(sapply(1:(length(tipMeans)-1), function(trait) cent_threshold_mat[other_lower_inds[trait,1], other_lower_inds[trait,2]]), length(possible_states)))
  rhos <- rep(cor[missinds[2],-missinds[2]], length(possible_states))
  
  probs <- pbivnorm_infs(uppers[,1], uppers[,2], rhos) - 
    pbivnorm_infs(lowers[,1], uppers[,2], rhos) - 
    pbivnorm_infs(uppers[,1], lowers[,2], rhos) + 
    pbivnorm_infs(lowers[,1], lowers[,2], rhos)
  
  probs[probs < 0] <- 0 #get rid of floating point errors
  probs <- log(probs)
  if(dimAdj){
    #need all other pairwise combis... or do I? differences between them do not change with other states!
    probs <- sapply(possible_states, function(state) (probs[(state*(length(tipMeans)-1)+1:(length(tipMeans)-1))]))
    if(class(probs) == "numeric"){probs <- t(as.matrix(probs))}
    if(length(which(apply(probs, 1, allEqualNinf))) > 0){probs <- probs[-which(apply(probs, 1, allEqualNinf)),]}
    probs <- (apply(probs, 2, sum) / (nrow(probs) + choose(dim(cor)[1] - 1, 2)))
    probs <- probs / 2 * dim(cor)[1]
    
  } else{
    probs <- sapply(possible_states, function(state) mean(probs[(state*(length(tipMeans)-1)+1:(length(tipMeans)-1))]))
  }
  return(probs)
  
}
bayesian_normalize <- function(stateProbs, missing_cond_probs){
  unnorm_probs <- stateProbs * missing_cond_probs
  return(unnorm_probs / sum(unnorm_probs))
}
MNAR_imputation <- function(imputed_probs, maxStates){
  buffNegInfs <- sapply(1:length(imputed_probs), function(tip) t(sapply(1:length(imputed_probs[[tip]][[2]]), 
                                                                        function(ind) c(imputed_probs[[tip]][[2]][[ind]], rep(-Inf, maxStates - length(imputed_probs[[tip]][[2]][[ind]]))))))
  buffNegInfs <- do.call(rbind, buffNegInfs)
  allStateProbs <- t(apply(buffNegInfs, 1, softmax))
  bimp_pars <- c(0.31, 0.75) #base prob missingness & logit-scale effect on base prob
  bimp_dat <- allStateProbs
  bimp_pars <- optim( par = bimp_pars, bimp_optim_func, bimp_dat = bimp_dat, method = method, control = list(trace = 0, REPORT = 1),
                      upper = c(1-exp(-4),10), lower = c(exp(-4),-10))$par
  bimp_par <- 0.5
  bimp_par <- optim( par = bimp_par, bimp_optim_func_just_disp, bimp_dat = bimp_dat, method = method, control = list(trace = 0, REPORT = 1),
                     upper = 5, lower = -5)$par
  return(bimp_pars)
}
bimp_optim_func <- function(bimp_pars, bimp_dat){
  base_state <- 3
  state_specific_missing_probs <- invlogit(logit(bimp_pars[1]) + cumsum(rep(bimp_pars[2], ncol(bimp_dat))) - base_state * bimp_pars[2])
  updated_state_probs <- t(sapply(1:nrow(bimp_dat), function(obs) bayesian_normalize(bimp_dat[obs,], state_specific_missing_probs)))
  probs_of_preferred_state <- sapply(1:nrow(updated_state_probs), function(obs) updated_state_probs[obs,which.max(updated_state_probs[obs,])])
  return(-sum(log(probs_of_preferred_state)))
}
bimp_optim_func_just_disp <- function(bimp_par, bimp_dat){
  base_state <- 3
  base_prob <- 0.3
  state_specific_missing_probs <- invlogit(logit(base_prob) + cumsum(rep(bimp_par, ncol(bimp_dat))) - base_state * bimp_par)
  updated_state_probs <- t(sapply(1:nrow(bimp_dat), function(obs) bayesian_normalize(bimp_dat[obs,], state_specific_missing_probs)))
  probs_of_preferred_state <- sapply(1:nrow(updated_state_probs), function(obs) updated_state_probs[obs,which.max(updated_state_probs[obs,])])
  return(-sum(log(probs_of_preferred_state)))
}
countObsStates <- function(states, maxStates = 7){
  sttab <- (table(states))
  obs_states <- rep(0, maxStates)
  obs_states[as.numeric(attr(sttab, which = "dimnames")$states) + 1] <- as.integer(sttab)
  return(obs_states)
}
countUnobsStates <- function(trait_inds, state_probs, maxStates = 7, ntraits, use_opt = F){
  conv_buff_probs <- t(sapply(1:length(state_probs), function(foo) softmax(c(state_probs[[foo]], rep(-Inf, maxStates - length(state_probs[[foo]]))))))
  if(use_opt){
    maxInds <- apply(conv_buff_probs, 1, which.max)
    conv_buff_probs <- t(sapply(1:length(maxInds), function(st) c(rep(0, maxInds[st] - 1), 1, rep(0, maxStates - maxInds[st]))))
  }
  unique_trait_inds <- sort(unique(trait_inds))
  trait_state_probs <- matrix(0, nrow = ntraits, ncol = maxStates)
  trait_state_probs[unique_trait_inds,] <- t(sapply(unique_trait_inds, function(trait) apply(conv_buff_probs[trait_inds == trait,], 2, sum)))
  trait_state_probs
}
rbeta_flat_update_sample_matrix <- function(matrix_a, matrix_b){
  samp_probs <- sapply(1:nrow(matrix_a), function(rowi) 
    sapply(1:ncol(matrix_a), function(colj) 
      rbeta(1, shape1 = 1 + max(matrix_a[rowi, colj], -0.99), shape2 = 1 + max(matrix_b[rowi, colj], -0.99))
    )
  )
  return(t(samp_probs))
}
find_missing_probs <- function(n_obs_trait_state_per_tip, n_imputed_trait_state_per_tip, per_pop = F, per_trait = T, raw = F){
  if(raw){
    if(per_pop){
      return(lapply(1:length(n_obs_trait_state_per_tip), function(tip) n_imputed_trait_state_per_tip[[tip]] / (n_imputed_trait_state_per_tip[[tip]] + n_obs_trait_state_per_tip[[tip]])))
    } else {
      n_obs_trait_state_per_tip_total <- n_obs_trait_state_per_tip[[1]]
      for(i in 2:length(n_obs_trait_state_per_tip)){
        n_obs_trait_state_per_tip_total <- n_obs_trait_state_per_tip_total + n_obs_trait_state_per_tip[[i]]
      }
      n_imputed_trait_state_per_tip_total <- n_imputed_trait_state_per_tip[[1]]
      for(i in 2:length(n_imputed_trait_state_per_tip)){
        n_imputed_trait_state_per_tip_total <- n_imputed_trait_state_per_tip_total + n_imputed_trait_state_per_tip[[i]]
      }
      if(!per_trait){
        return(sapply(1:7, function(x) sum(n_imputed_trait_state_per_tip_total[,x]) / sum((n_imputed_trait_state_per_tip_total[,x] + n_obs_trait_state_per_tip_total[,x]))))
      } else {
        return(n_imputed_trait_state_per_tip_total / (n_imputed_trait_state_per_tip_total + n_obs_trait_state_per_tip_total)) 
      }
    }
  } else {
    #update beta distribution with binomial count
    if(per_pop){
      return(lapply(1:length(n_obs_trait_state_per_tip), function(tip) 
        rbeta_flat_update_sample_matrix(n_imputed_trait_state_per_tip[[tip]], n_obs_trait_state_per_tip[[tip]])))
    } else {
      n_obs_trait_state_per_tip_total <- n_obs_trait_state_per_tip[[1]]
      for(i in 2:length(n_obs_trait_state_per_tip)){
        n_obs_trait_state_per_tip_total <- n_obs_trait_state_per_tip_total + n_obs_trait_state_per_tip[[i]]
      }
      n_imputed_trait_state_per_tip_total <- n_imputed_trait_state_per_tip[[1]]
      for(i in 2:length(n_imputed_trait_state_per_tip)){
        n_imputed_trait_state_per_tip_total <- n_imputed_trait_state_per_tip_total + n_imputed_trait_state_per_tip[[i]]
      }
      if(!per_trait){
        tr_imp <- sapply(1:ncol(n_imputed_trait_state_per_tip_total), function(x) sum(n_imputed_trait_state_per_tip_total[,x]))
        tr_obs <- sapply(1:ncol(n_obs_trait_state_per_tip_total), function(x) sum(n_obs_trait_state_per_tip_total[,x]))
        return(rbeta_flat_update_sample_matrix(t(as.matrix(tr_imp)), t(as.matrix(tr_obs))))
      } else {
        return(rbeta_flat_update_sample_matrix(n_imputed_trait_state_per_tip_total, n_obs_trait_state_per_tip_total)) 
      }
    }
  }
}
nudgeCorHigher <- function(r){
  invlogit(logit(abs(r)) + 1) * sign(r)
}
computeUnivNormalStateProbs <- function(thresholds, mean, sd = 1, npop = NA){
  thresholds <- c(-Inf, thresholds, Inf)
  lower <- thresholds[1:(length(thresholds) - 1)]
  upper <- thresholds[2:length(thresholds)]
  if(is.na(npop)){
    return(pnorm(upper, mean = mean, sd = sd) - pnorm(lower, mean = mean, sd = sd))
  } else {
    return((pnorm(upper, mean = mean, sd = sd) - pnorm(lower, mean = mean, sd = sd)) * npop)
  }
}
computeConditionalsLessObserved <- function(conditionals, observed, npop = NA, addPop = F){
  rawDiff <- conditionals - observed
  if(min(rawDiff) > 0 | is.na(npop)){
    return(rawDiff)
  } else {
    if(addPop){
      condProbs <- conditionals / npop
      extraMissing <- -min(rawDiff / condProbs, na.rm = T)
      new_npop <- npop + extraMissing
      return(condProbs * new_npop - observed)
    } else {
      rawDiff[rawDiff < 0] <- 0
      return(rawDiff)
    }
  }
}
approxMultinomialCond <- function(probs, obs_counts, total_count, min_sample_size = 10, raw = F, binomial_rescue = NA, mcmc_rescue = 5, mcmc_start_niter = 1E4,
                                  mcmc_thin = 1E2, mcmc_burnin = 0.2, mcmc_swap_num_divisor = 5, troubleshoot = F){
  if(sum(obs_counts) == total_count){
    #hack to return 0s when there is no missing data... may not work with non-default parameter settings
    if(!raw){return(rep(0, length(obs_counts)))} else{return(matrix(0, nrow = min_sample_size, ncol = length(obs_counts)))}
  }
  samps <- t(rmultinom(n = min_sample_size, size = total_count, prob = probs))
  diffs <- t(sapply(1:nrow(samps), function(simsamp) samps[simsamp,] - obs_counts))
  valid <- sapply(1:nrow(samps), function(simsamp) all(diffs[simsamp,] >= 0))
  valid_diffs <- diffs[valid,]
  needed_extra <- min_sample_size / sum(valid) * min_sample_size
  order_mag_needed_extra <- 1
  if(needed_extra == Inf){
    needed_extra = 2^order_mag_needed_extra*min_sample_size
    order_mag_needed_extra <- order_mag_needed_extra + 1
  }
  while(needed_extra > min_sample_size){
    extra_samps <- t(rmultinom(n = needed_extra, size = total_count, prob = probs))
    extra_diffs <- t(sapply(1:nrow(extra_samps), function(simsamp) extra_samps[simsamp,] - obs_counts))
    extra_valid <- sapply(1:nrow(extra_samps), function(simsamp) all(extra_diffs[simsamp,] >= 0))
    valid_diffs <- rbind(valid_diffs, extra_diffs[extra_valid,])
    valid <- sum(extra_valid, valid)
    needed_extra <- min_sample_size / valid * min_sample_size
    if(needed_extra == Inf){
      needed_extra = 2^order_mag_needed_extra*min_sample_size
      order_mag_needed_extra <- order_mag_needed_extra + 1
      if(!is.na(binomial_rescue)){
        if(order_mag_needed_extra > binomial_rescue){
          print("unable to approximate multinomial probability -- resorting to binomial probability")
          ntraits <- length(probs)
          exp_extras <- sapply(1:ntraits, function(marg_ind)
            sum(diff(pbinom(q = obs_counts[marg_ind]:(total_count-sum(obs_counts[-marg_ind])), size = total_count, prob = probs[marg_ind])) * 1:(total_count-sum(obs_counts)))
          )
          if(!raw){return(exp_extras / sum(exp_extras))} else{
            return(t(sapply(1:min_sample_size, function(foo) exp_extras)))
          }
        }
      }
      if(!is.na(mcmc_rescue)){
        if(order_mag_needed_extra > mcmc_rescue){
          print("unable to approximate multinomial probability via monte carlo simulation -- resorting to ~markov chain~ monte carlo simulation")
          exp_extras <- approxMultinomialCondMCMC(probs = probs, obs_counts = obs_counts, total_count = total_count, 
                                                  min_ess = min_sample_size, start_niter = mcmc_start_niter, thin = mcmc_thin, burnin = mcmc_burnin, exp_count = raw,
                                                  swap_num_divisor = mcmc_swap_num_divisor, troubleshoot = troubleshoot)
          if(raw){return(exp_extras)} else{
            return(apply(exp_extras, 2, mean))
          }
        }
      }
      
      
    }
  }
  
  #approximate the probabilities of the remaining, unobserved categories
  if(!raw){
    cond_probs <- apply(valid_diffs, 2, sum)
    cond_probs <- cond_probs / sum(cond_probs) 
    return(cond_probs)
  } else {
    if(class(valid_diffs) == "matrix"){
      return(valid_diffs[1:min_sample_size,])
    } else {
      return(valid_diffs)
    }
  }
}
approxMultinomialCondMCMC <- function(probs, obs_counts, total_count, min_ess = 100, start_niter = 1E5, thin = 1E1, burnin = 0.2, 
                                      exp_count = T, samples = T, swap_num_divisor = 5, thinAuto = T, troubleshoot = F, extraIterProp = 2){
  ntraits <- length(probs)
  size_obs <- sum(obs_counts)
  samp_counts <- matrix(0, nrow = start_niter / thin, ncol = ntraits)
  curr_state <- obs_counts + t(rmultinom(n = 1, size = total_count - size_obs, prob = probs))
  curr_prob <- dmultinom(curr_state, prob = probs, log = T)
  accept_count <- 0
  for(i in 1:start_niter){
    if(troubleshoot){if(i %% thin == 0){cat(paste0(i / thin, " "))}}
    amount_to_swap <- sample(size = 1, 1:floor((total_count - size_obs)/swap_num_divisor))
    inds_to_swap <- sample(size = 2, 1:ntraits, replace = F)
    prop_state <- curr_state
    prop_state[inds_to_swap[1]] <- prop_state[inds_to_swap[1]] + amount_to_swap
    prop_state[inds_to_swap[2]] <- prop_state[inds_to_swap[2]] - amount_to_swap
    if(all(prop_state >= obs_counts)){
      prop_prob <- dmultinom(prop_state, prob = probs, log = T)
      if(exp(prop_prob - curr_prob) > runif(1,0,1)){
        curr_state <- prop_state
        curr_prob <- prop_prob
        accept_count <- accept_count + 1
      }
    }# else {prop_prob <- -Inf} ...implied
    
    if(i %% thin == 0){
      if(troubleshoot){cat(paste0("r", i / thin, " "))}
      samp_counts[i/thin,] <- curr_state
    }
  }
  if(troubleshoot){print("we made it to the first ess check")}
  samp_counts <- samp_counts[ceiling(burnin*nrow(samp_counts)):nrow(samp_counts),]
  if(troubleshoot){cat("samp_counts:", " "); print(samp_counts)}
  diff_counts <- t(sapply(1:nrow(samp_counts), function(count) samp_counts[count,] - obs_counts))
  if(troubleshoot){cat("diff_counts:", " "); print(diff_counts)}
  if(troubleshoot){print(paste(sum(obs_counts), total_count))}
  ess <- coda::effectiveSize(samp_counts[,apply(diff_counts,2,sum) != 0])
  if(troubleshoot){cat("ess:", " "); print(ess)}
  if(troubleshoot){print("we made it past the first ess check")}
  if(all(ess > min_ess)){
    if(exp_count){
      if(samples){
        return(diff_counts)
      } else {
        return(apply(diff_counts, 2, mean))
      }
    } else {
      est_cond_probs <- apply(diff_counts, 2, mean)
      return(est_cond_probs / sum(est_cond_probs))
    }
  } else {
    while(min(ess) < min_ess){
      if(min(ess) == 0){
        ess[ess==0] <- 1
        if(troubleshoot){print("doing the adjustment for ess = 0")}
        
      }
      followup_niter <- floor((min_ess / min(ess)) * start_niter * extraIterProp)
      if(troubleshoot){print(paste0("starting over with ", followup_niter," samples..."))}
      if(thinAuto){
        thin <- floor(followup_niter / 1E4)
      }
      followup_samp_counts <- matrix(0, nrow = floor(followup_niter / thin), ncol = ntraits)
      for(i in 1:followup_niter){
        amount_to_swap <- sample(size = 1, 1:floor((total_count - size_obs)/swap_num_divisor))
        inds_to_swap <- sample(size = 2, 1:ntraits, replace = F)
        prop_state <- curr_state
        prop_state[inds_to_swap[1]] <- prop_state[inds_to_swap[1]] + amount_to_swap
        prop_state[inds_to_swap[2]] <- prop_state[inds_to_swap[2]] - amount_to_swap
        if(all(prop_state >= obs_counts)){
          prop_prob <- dmultinom(prop_state, prob = probs, log = T)
          if(exp(prop_prob - curr_prob) > runif(1,0,1)){
            curr_state <- prop_state
            curr_prob <- prop_prob
          }
        }
        if(i %% thin == 0){
          followup_samp_counts[i/thin,] <- curr_state
        }
      }
      followup_diff_counts <- t(sapply(1:nrow(followup_samp_counts), function(count) followup_samp_counts[count,] - obs_counts))
      diff_counts <- rbind(diff_counts, followup_diff_counts)
      ess <- coda::effectiveSize(diff_counts[,apply(diff_counts,2,sum) != 0])
    }
    if(exp_count){
      if(samples){
        return(diff_counts)
      } else {
        return(apply(diff_counts, 2, mean))
      }
    } else {
      est_cond_probs <- apply(diff_counts, 2, mean)
      return(est_cond_probs / sum(est_cond_probs))
    }
  }
}
MAP_tree <- function(trees, return_posteriorMean_BLs = T, n_trees_to_return = 1){
  all_the_trees <- trees
  unique_topologies <- trees[[1]]
  same_trees <- sapply(1:length(trees), function(tree) RF.dist(trees[[tree]], trees[[1]])) == 0
  if(n_trees_to_return == "all" & return_posteriorMean_BLs){
    unique_topologies <- maxCladeCred(trees[same_trees])
  } else {
    unique_topologies <- trees[[1]]
  }
  n_per_unique_topology <- sum(same_trees)
  trees <- trees[!same_trees]
  # while(max(n_per_unique_topology) < length(trees)){
  while(length(trees) > 0){
    same_trees <- sapply(1:length(trees), function(tree) RF.dist(trees[[tree]], trees[[1]])) == 0
    if(n_trees_to_return == "all" & return_posteriorMean_BLs){
      unique_topologies <- c(unique_topologies, maxCladeCred(trees[same_trees]))
    } else {
      unique_topologies <- c(unique_topologies, trees[[1]])
    }
    n_per_unique_topology <- c(n_per_unique_topology, sum(same_trees))
    trees <- trees[!same_trees]
  }
  map_tree = unique_topologies[[which.max(n_per_unique_topology)]]
  if(return_posteriorMean_BLs & n_trees_to_return != "all"){
    map_tree = maxCladeCred(all_the_trees[sapply(1:length(all_the_trees), function(tree) RF.dist(all_the_trees[[tree]], map_tree) == 0)])
  }
  if(n_trees_to_return == 1){
    return(list(map_tree = map_tree, prob = max(n_per_unique_topology) / sum(c(n_per_unique_topology, length(trees)))))
  }
  if(return_posteriorMean_BLs & n_trees_to_return == "all"){
    unique_topologies <- unique_topologies[order(n_per_unique_topology, decreasing = T)]
    n_per_unique_topology <- sort(n_per_unique_topology, decreasing = T)
    return(lapply(1:length(n_per_unique_topology), function(tree) list(tree = unique_topologies[tree], 
                                                                       prob = n_per_unique_topology[tree] / length(all_the_trees))))
  }
}
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
startsWith2 <- function(x, pre){apply(sapply(1:length(pre), function(n) startsWith(x, pre[n])), 1, any)}
greedyCT <- function(trees){
  tipNames <- trees[[1]]$tip.label
  tipNums <- paste0("t", 1:length(tipNames))
  for(i in 1:length(trees)){
    trees[[i]]$tip.label <- tipNums[match(trees[[i]]$tip.label, tipNames)]
  }
  ape::write.tree(trees, file = "gct_trees.txt")
  writeLines("gct_trees.txt\nY", con = "consense_script.txt")
  if(file.exists("outfile")){file.remove("outfile")}
  if(file.exists("outtree")){file.remove("outtree")}
  system(paste0("/Applications/phylip-3.695/exe/consense < consense_script.txt"), intern = T)
  gct <- ape::read.tree("outtree")
  if(file.exists("outfile")){file.remove("outfile")}
  if(file.exists("outtree")){file.remove("outtree")}
  if(file.exists("gct_trees.txt")){file.remove("gct_trees.txt")}
  if(file.exists("consense_script.txt")){file.remove("consense_script.txt")}
  gct$tip.label <- tipNames[match(gct$tip.label, tipNums)]
  return(gct)
}
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
compareTrees <- function(biparts1, biparts2, tipLabs, returnCladeNames = F){
  matchProbs <- matrix(0, nrow = sum(length(biparts1[,1]), length(biparts2[,1])), ncol = ifelse(returnCladeNames, 3, 2))
  counter <- 1
  biparts2$notSeen <- 1
  for(clade in 1:length(biparts1[,2])){
    cladeName <- biparts1[,2][[clade]]
    isThere <- sapply(1:length(biparts2[,1]), function(x) identical(cladeName, biparts2[x,2][[1]]))
    altCladeName <- sort(setdiff(tipLabs, cladeName))
    shortCladeName <- paste0(sort(shorter(cladeName, altCladeName)), collapse = ", ")
    orIsThere <- sapply(1:length(biparts2[,1]), function(x) identical(altCladeName, biparts2[x,2][[1]]))
    if(any(isThere)){
      biparts2$notSeen[isThere] <- 0
      if(returnCladeNames){
        matchProbs[counter,] <- c(biparts1[,1][clade], biparts2[,1][isThere], shortCladeName)
      } else {
        matchProbs[counter,] <- c(biparts1[,1][clade], biparts2[,1][isThere])
      }
    } else if (any(orIsThere)) {
      biparts2$notSeen[orIsThere] <- 0
      if(returnCladeNames){
        matchProbs[counter,] <- c(biparts1[,1][clade], biparts2[,1][orIsThere], shortCladeName)
      } else {
        matchProbs[counter,] <- c(biparts1[,1][clade], biparts2[,1][orIsThere])
      }
    } else {
      if(returnCladeNames){
        matchProbs[counter,] <- c(biparts1[,1][[clade]], 0, shortCladeName)
      } else {
        matchProbs[counter,] <- c(biparts1[,1][[clade]], 0)
      }
    }
    counter <- counter + 1
  }
  if(sum(biparts2$notSeen) > 0){
    if(returnCladeNames){
      cladeNames <- sapply(1:length(biparts2[biparts2$notSeen == 1,2]), function(bip) 
        paste0(sort(shorter(biparts2[biparts2$notSeen == 1,2][[bip]], sort(setdiff(tipLabs, biparts2[biparts2$notSeen == 1,2][[bip]])))), collapse = ", "))
      matchProbs[counter:(counter+sum(biparts2$notSeen)-1),] <- cbind(0, cbind(biparts2[biparts2$notSeen == 1,1], cladeNames)) #needs to be an sapply
    } else {
      matchProbs[counter:(counter+sum(biparts2$notSeen)-1),] <- cbind(0, biparts2[biparts2$notSeen == 1,1])
    }
  }
  if(returnCladeNames){
    matchProbs <- matchProbs[rowSums(matchProbs == c(0,0)) != ifelse(returnCladeNames,3,2),]
  }
  matchProbs
}
clade_prob <- function(tipNames, trees, allTipNames = NA, partfreqs = NA){
  if(is.na(partfreqs)[1]){
    clade_monophyly <- prop.part.df(trees, cutoff = 0.00001)
    if(is.na(allTipNames)[1]){
      allTipNames <- trees[[1]]$tip.label
    }
  } else {
    clade_monophyly <- partfreqs
    if(is.na(allTipNames)[1]){
      allTipNames <- unique(unlist(bipart_probs$cladeNames))
    }
  }
  cladeName <- tipNames
  clade <- list(allTipNames[startsWith2(allTipNames, cladeName)], allTipNames[!startsWith2(allTipNames, cladeName)])
  clade_ind <- (sapply(1:length(clade_monophyly$cladeNames), function(x) setequal(clade_monophyly$cladeNames[[x]], clade[[1]]) |
                         setequal(clade_monophyly$cladeNames[[x]], clade[[2]])))
  if(all(!clade_ind)){return(0)}
  # clade_monophyly$cladeNames[clade_ind]
  return(clade_monophyly$postProbs[clade_ind])
}
internal_nodes_probs <- function(tree, trees_to_use){
  bipart_probs <- prop.part.df(trees_to_use, cutoff = 0.00001)
  ti <- 1:length(tree$tip.label)
  tips <- sapply(1:Nnode(tree), function(node) as.numeric(na.omit(ti[getDescendants(tree, length(tree$tip.label) + node)])))
  tips <- sapply(1:length(tips), function(node) tree$tip.label[tips[[node]]])
  node_pps <- sapply(1:length(tips), function(node) clade_prob(tipNames = tips[[node]], allTipNames = tree$tip.label, partfreqs = bipart_probs))
  return(node_pps)
}
shorter <- function(x1, x2){
  if(length(x1) > length(x2)){
    return(x2)
  } else{
    return(x1)
  }
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
} #slowest step by far
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
updateLi2 <- function(curr_Li, sd_index, change_ratio){
  curr_Li[,sd_index] <- curr_Li[,sd_index] * sqrt(change_ratio[1])
  curr_Li[,-sd_index] <- curr_Li[,-sd_index] * sqrt(change_ratio[2])
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
BM_LL_PRECOMPUTE <- function(transformed_contrasts, prunes, detcov){
  n_traits <- ncol(transformed_contrasts)
  sum(sapply(1:nrow(transformed_contrasts), function(x) sum(dnorm(transformed_contrasts[x,] / prunes[[2]][x]^0.5, mean = 0, sd = 1, log = T)) - log(detcov*(prunes[[2]][x]^n_traits))/2))
}
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
condMVN <- function(means, cov, obs, inds){
  A <- cov[-inds, -inds]
  C <- cov[inds, inds]
  B <- cov[-inds, inds]
  condCov <- A - B%*%solve(C)%*%t(B)
  condMean <- as.vector(means[-inds] + B%*%solve(C)%*%(obs - means[inds]))
  return(list(means = condMean, cov = condCov))
}
estimate_root_brownian_bridge <- function(tree, means, rate_matrix){
  n_traits <- nrow(rate_matrix)
  artificial_root_tip <- sample(tree$tip.label, 1)
  tree_root <- add.tips(tree, tips = "ROOT", where = length(tree$tip.label) + 1, edge.length = 0)
  tree_root <- reroot(tree_root, node.number = which(tree_root$tip.label == artificial_root_tip))
  tree_cov <- vcv.phylo(tree_root)
  tree_cov <- tree_cov[rownames(tree_cov) != artificial_root_tip, rownames(tree_cov) != artificial_root_tip]
  data_root_cov <- kronecker(tree_cov, rate_matrix)
  artificial_root_mean <- rep(inferred_empirical_means[rownames(inferred_empirical_means) == artificial_root_tip,], nrow(tree_cov))
  all_other_means <- inferred_empirical_means[rownames(inferred_empirical_means) != artificial_root_tip,]
  all_other_means <- all_other_means[colnames(tree_cov)[colnames(tree_cov) != "ROOT"],]
  all_other_means <- c(t(all_other_means))
  root_conditional <- condMVN(means = artificial_root_mean, cov = data_root_cov, obs = all_other_means, inds = 1:(n_traits*(nrow(tree_cov)-1)))
  return(root_conditional)
}
avmode <- function(x) {
  if(all(is.na(x))){return(1)}
  ux <- unique(x[!is.na(x)])
  ux[which.max(tabulate(match(as.vector(x), ux)))]
}
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(c(alpha, beta))
}
marginal_beta_from_LKJ <- function(eta, dim){
  #this function retrieves the marginal rescaled beta distribution of each correlation parameter of an LKJ distribution with specified eta
  return(rep(eta - 1 + dim / 2, 2))
} 

#read in data from previous step
load(file = paste0("traits_indiv_discr/traits_indiv_discr_", rep_i))
load(file = paste0("obs/obs_", rep_i))
load(file = paste0("constraint_arrays/constraintArray_", rep_i))

npops <- nTaxa <- length(obs)
d_traits <- ncol(obs[[1]])
pops <- names(obs)

#get useful variables
trait_names <- paste0("tr", 1:d_traits)
n_indiv_emp <- c(119, 51, 84, 40, 40, 17, 135, 198); n_indiv_emp_popmatch <- c("NEAND", "AUS-NG", "EUR", "WAS", "SAS", "NEAS", "SSAF", "AMER")
n_thresholds <- c(5, 6, 4, 6, 7, 5, 3, 6, 4, 3, 6, 3, 4, 2, 2, 2, 4, 3, 3, 2, 4, 2, 2, 2, 1, 4, 3, 1, 8, 5, 5, 7, 4, 3, 3, 3, 3, 5, 5, 5, 7, 4, 3, 3, 5, 3, 5, 5, 5, 7, 6, 4, 2, 3, 3, 3, 1, 2, 4, 3, 3, 3, 3, 2, 5, 2, 3, 2, 1, 2, 4, 2, 4, 2, 3, 3, 3, 2, 3, 3, 4, 2, 4, 4, 2, 3, 6, 5, 5, 4, 4, 5, 3, 2, 3, 6, 5, 5, 4, 4, 4, 6, 3, 2, 3, 7, 5, 5, 4, 4, 4, 3, 4, 3, 3, 3, 3, 3)
max_thresholds_count <- max(n_thresholds)

########################################
### now we are ready to do inference ###
########################################

#first pass initialization of all missing states to their tip's median value
for(tip in pops){
  for(indiv in 1:length(traits_indiv_discr[[tip]][,1])){
    whichNA <- which(is.na(traits_indiv_discr[[tip]][indiv,]))
    traits_indiv_discr[[tip]][indiv,whichNA] <- apply(as.matrix(traits_indiv_discr[[tip]][,whichNA]), 2, avmode)
  }
}

#identify unique site patterns
traits_indiv_discr_unique <- lapply(pops, function(tip) uniqueSP(traits_indiv_discr[[tip]]))
# names(traits_indiv_discr_unique) <- pops
dat <- as.data.frame(do.call(rbind, lapply(1:nTaxa, function(tip) traits_indiv_discr_unique[[tip]]$SPs)))
colnames(dat) <- trait_names
dat2 = data.frame(nSP = unlist(lapply(1:nTaxa, function(tip) traits_indiv_discr_unique[[tip]]$counts)),
                  tip = unlist(lapply(1:nTaxa, function(tip) rep(tip, length(traits_indiv_discr_unique[[tip]]$counts)))))

dat <- dat_orig <- cbind(dat, dat2)

#better init
thresholds_init <- matrix(rep(0:(max_thresholds_count-1), d_traits), nrow = d_traits, ncol = max_thresholds_count,  byrow = T) / 2
thresholds_init <- t(sapply(1:d_traits, function(trait) c(thresholds_init[trait,1:n_thresholds[trait]], rep(Inf, max_thresholds_count - n_thresholds[trait]))))
threshold_diffs_init <- t(sapply(1:d_traits, function(trait) diff(thresholds_init[trait,]))); threshold_diffs_init[is.na(threshold_diffs_init)] <- Inf
n_thresh <- sapply(1:length(thresholds_init[,1]), function(trait) length((thresholds_init[trait,])))
locations <- cbind(rep(-0.5, d_traits), thresholds_init, rep(Inf, d_traits))
infs_first_index <- sapply(1:d_traits, function(trait) which(locations[trait,] == Inf)[1])
for(i in 1:d_traits){locations[i,infs_first_index] <- locations[i,infs_first_index-1] + 0.5}
init_cor <- cor(do.call(rbind, lapply(pops, function(x) (traits_indiv_discr[[x]]))))
init_cor[is.na(init_cor)] <- 0

#initialize beta regularizing param
init_cor_beta <- init_cor[upper.tri(init_cor)]
init_cor_beta <- (init_cor_beta + 1) / 2
init_cor_beta <- estBetaParams(mean(init_cor_beta), var(init_cor_beta))
init_cor_beta <- marginal_beta_from_LKJ(1, d_traits)
init_cor_beta <- c(10,10)

#initialize means
init_means <- t(sapply(pops, function(tip) 
  sapply(1:length(traits_indiv_discr[[tip]][1,]), function(trait) mean(
    sapply(1:length(traits_indiv_discr[[tip]][,1]), function(indiv) 
      mean(locations[trait, c(traits_indiv_discr[[tip]][indiv,][trait] + 1, traits_indiv_discr[[tip]][indiv,][trait] + 2)]))
  )
  )
)
)

rownames(init_means) <- pops; colnames(init_means) <- trait_names

# thresh_mu_init <- 1
# thresh_sd_init <- 100
lambda_init <- 0.01

par <- c(c(t(init_means)), 10, init_cor[upper.tri(init_cor)], c(threshold_diffs_init), lambda_init, init_cor_beta)

#get the data that ignores imputed values ready
traits_indiv_discr9 <- traits_indiv_discr
for(tip in pops){
  for(indiv in 1:length(traits_indiv_discr[[tip]][,1])){
    whichNA <- which(!(obs[[tip]][indiv,]))
    traits_indiv_discr9[[tip]][indiv,whichNA] <- 9 #here stands for a missing trait to be ignored
  }
}
traits_indiv_discr_unique9 <- lapply(pops, function(tip) uniqueSP(traits_indiv_discr9[[tip]]))
# names(traits_indiv_discr_unique9) <- pops
dat9 <- as.data.frame(do.call(rbind, lapply(1:nTaxa, function(tip) traits_indiv_discr_unique9[[tip]]$SPs)))
colnames(dat9) <- trait_names
dat29 = data.frame(nSP = unlist(lapply(1:nTaxa, function(tip) traits_indiv_discr_unique9[[tip]]$counts)),
                   tip = unlist(lapply(1:nTaxa, function(tip) rep(tip, length(traits_indiv_discr_unique9[[tip]]$counts)))))

dat9 <- dat_orig9 <- cbind(dat9, dat29)

#get mahalonibis distance matrix
mahDistMat_init <- mahaMatrix(init_means, init_cor, squared = F)

#compute tree
njTree_init <- nj(mahDistMat_init) #can use outgroup for rooting
njTree_init$edge.length[njTree_init$edge.length < 0] <- 0
njTree_init <- reroot(tree = njTree_init, node.number = which(njTree_init$tip.label == "NEAND"), position = njTree_init$edge.length[which(njTree_init$edge[,2] == which(njTree_init$tip.label == "NEAND"))] / 2)
plot(njTree_init)
upgmaTree_init <- upgma(mahDistMat_init)
if(useUPGMA){
  treeCV_init <- cov2cor(vcv.phylo(upgmaTree_init))
} else {
  treeCV_init <- cov2cor(vcv.phylo(njTree_init))
}
if(weighPCV){treeCV_init <- (treeCV_init + diag(nTaxa)) / 2} #weight with star phylogeny
if(useStarPhylogeny){
  treeCV_init[upper.tri(treeCV_init)] <- 0
  treeCV_init[lower.tri(treeCV_init)] <- 0
}
prunes_init <- prunePCV(treeCV_init)

#TEST IF ALGORITHM WORKS
thresh_trait_ind <- sample(1:d_traits, 1)
par_to_opt_ind <- nTaxa*d_traits+2+choose(d_traits, 2) + 0:(max_thresholds_count-2) * d_traits + thresh_trait_ind - 1
par_to_opt_ind <- par_to_opt_ind[1:(n_thresholds[thresh_trait_ind]-1)]
#thresholds AND all the means for that threshold?
par_to_opt_ind <- c(thresh_trait_ind + d_traits*0:(nTaxa-1), par_to_opt_ind) #thresholds AND all the means for that threshold?
# par <- c(c(t(traits)), 10, R[upper.tri(R)], c(threshold_diffs_init)) #try using true par values to check for good behavior
# par_to_opt_ind <- length(par) + c(-1,0)
par_to_opt_ind <- 5 #mean
par_to_opt_ind <- (nTaxa * d_traits + 2) #correlation
par_to_opt_ind <- (nTaxa * d_traits + 1 + choose(d_traits, 2) + (max_thresholds_count-1) * d_traits) + c(2,3) #correlation regularizing beta
par_to_opt <- par[par_to_opt_ind]
type_of_param <- ifelse(par_to_opt_ind[1] <= nTaxa * d_traits, ifelse(par_to_opt_ind[length(par_to_opt_ind)] > nTaxa * d_traits, "tr-m", "m"), 
                        ifelse(par_to_opt_ind[1] == (nTaxa * d_traits + 1), "v",
                               ifelse(par_to_opt_ind[1] > (nTaxa * d_traits + 1) & par_to_opt_ind[1] <= (nTaxa * d_traits + 1 + choose(d_traits, 2)), "c", 
                                      ifelse(par_to_opt_ind[1] > (nTaxa * d_traits + 1 + choose(d_traits, 2)) & 
                                               par_to_opt_ind[1] <= (nTaxa * d_traits + 1 + choose(d_traits, 2) + (max_thresholds_count-1) * d_traits), "t", 
                                             ifelse(par_to_opt_ind == (nTaxa * d_traits + 1 + choose(d_traits, 2) + (max_thresholds_count-1) * d_traits) + 1, "t-r",
                                                    ifelse(all.equal(sort(par_to_opt_ind), (nTaxa * d_traits + 1 + choose(d_traits, 2) + (max_thresholds_count-1) * d_traits) + c(2,3)), "c-r", ""))))))
print(type_of_param)
k_par <- par
k_par[par_to_opt_ind] <- NA  
dat <- list(dat_orig, k_par, prunes_init, pops, n_thresh)
dat <- list(dat_orig, k_par, prunes_init, pops, n_thresh, getUniqueSitePatterns(par_to_opt, dat)) 
logLL_bivProb_optim_multithresh(par_to_opt = par_to_opt, dat =  dat)
logLL_bivProb_optim_multithresh(rep(100,2), dat)

method <- "L-BFGS-B"



################################################################################
################################################################################
################################################################################

############################ Start Inference Here ##############################

################################################################################
################################################################################
################################################################################

ncore <- 1; mclapply_over_foreach <- T;
nparam <- length(par)
prunes <- prunes_init

if(noSmartStart){
  starTree <- diag(length(njTree_init$tip.label))
  colnames(starTree) <- rownames(starTree) <- (njTree_init$tip.label)
  prunes <- prunePCV(starTree*10)
  par <- par <- c(rep(0, length(c(t(init_means)))), 10, rep(0, length(init_cor[upper.tri(init_cor)])), 
                  c(threshold_diffs_init), lambda_init, init_cor_beta)
}


for(i in 1:nrounds){
  cat(paste0("\n", "round ", i, ": "))
  param_ord <- as.list(sample(c(1:(nTaxa*d_traits+1+choose(d_traits, 2))))) #univ optim, double emphasis on correlation parameters and threshold diffs
  # param_ord <- as.list(sample(c(1:nparam))) #univ optim, double emphasis on correlation parameters and threshold diffs
  thresh_trait_ind <- sample(which(n_thresholds > 1))
  par_to_opt_ind <- lapply(thresh_trait_ind, function(trait) (nTaxa*d_traits+2+choose(d_traits, 2) + 0:(max_thresholds_count-2) * d_traits + trait - 1)[0:(n_thresholds[trait]-1)])
  if(jointThreshAndMeans & i >= jTM_startAtRound){
    par_to_opt_ind <- lapply(thresh_trait_ind, function(trait) c(trait + d_traits*0:(nTaxa-1),
                                                                 (nTaxa*d_traits+2+choose(d_traits, 2) + 0:(max_thresholds_count-2) * d_traits + trait - 1)[0:(n_thresholds[trait]-1)]))
  }
  if(T){
    param_ord <- c(par_to_opt_ind, param_ord)
  } else {
    param_ord <- sample(c(par_to_opt_ind, param_ord))
  }
  
  #add in the threshold regularizing terms at the end
  param_ord <- (c(param_ord, list(((nTaxa*d_traits+1+choose(d_traits, 2)) + sum(n_thresh-1)) + 1)))
  
  #add in the correlation regularizing terms at the end
  # param_ord <- (c(param_ord, list(((nTaxa*d_traits+1+choose(d_traits, 2)) + sum(n_thresh-1)) + c(2,3))))
  
  num_iters_passed <- 0
  for(par_to_opt_ind in param_ord){
    type_of_param <- ifelse(par_to_opt_ind[1] <= nTaxa * d_traits, ifelse(par_to_opt_ind[length(par_to_opt_ind)] > nTaxa * d_traits, "tr-m", "m"), 
                            ifelse(par_to_opt_ind[1] == (nTaxa * d_traits + 1), "v",
                                   ifelse(par_to_opt_ind[1] > (nTaxa * d_traits + 1) & par_to_opt_ind[1] <= (nTaxa * d_traits + 1 + choose(d_traits, 2)), "c", 
                                          ifelse(par_to_opt_ind[1] > (nTaxa * d_traits + 1 + choose(d_traits, 2)) & 
                                                   par_to_opt_ind[1] <= (nTaxa * d_traits + 1 + choose(d_traits, 2) + (max_thresholds_count-1) * d_traits), "t", 
                                                 ifelse(par_to_opt_ind == (nTaxa * d_traits + 1 + choose(d_traits, 2) + (max_thresholds_count-1) * d_traits) + 1, "t-r",
                                                        ifelse(all.equal(sort(par_to_opt_ind), (nTaxa * d_traits + 1 + choose(d_traits, 2) + (max_thresholds_count-1) * d_traits) + c(2,3)), "c-r", ""))))))
    cat(paste0("round ", i, ", ", round(num_iters_passed / (length(param_ord) / ncore) * 100), "%: ", type_of_param, "-"))
    num_iters_passed = num_iters_passed + 1
    # cat(paste0(par_to_opt_ind, type_of_param, "-"))
    par_to_opt <- par[par_to_opt_ind]
    k_par <- par
    k_par[par_to_opt_ind] <- NA  
    dat <- list(dat_orig, k_par, prunes, pops, n_thresh)
    if(type_of_param != "c" | use_imputed_data_for_corrs){
      corrOK <- T
    } else {
      corrOK <- F
    }
    # uniqueSitePats <- getUniqueSitePatterns(par_to_opt, dat)
    if(!ignoreImputedData & (i > nrounds_to_not_impute) & corrOK){
      uniqueSitePats <- getUniqueSitePatterns(par_to_opt, dat)
    } else if(ignoreImputedData){
      uniqueSitePats <- getUniqueSitePatternsIgnoreImputed(par_to_opt, dat)
    } else if(!ignoreImputedData & (i <= nrounds_to_not_impute) | !corrOK){
      dat9 <- list(dat_orig9, k_par, prunes, pops, n_thresh)
      uniqueSitePats <- getUniqueSitePatternsIgnoreImputed(par_to_opt, dat9)
    }
    dat <- list(dat_orig, k_par, prunes_init, pops, n_thresh, uniqueSitePats) 
    upper = c(rep(Inf, nTaxa * d_traits + 1), rep(1, choose(d_traits, 2)), rep(Inf, (max_thresholds_count-1) * d_traits), c(Inf, Inf, Inf, Inf))[par_to_opt_ind]
    lower = c(rep(-Inf, nTaxa * d_traits), 0.01, rep(-1, choose(d_traits, 2)), rep(0, (max_thresholds_count-1) * d_traits), c(-Inf, 1, 1, 1))[par_to_opt_ind]
    optim_out <- optim(par = par_to_opt, logLL_bivProb_optim_multithresh, dat = dat, method = method, control = list(trace = 0, REPORT = 1),
                       upper = upper, lower = lower)
    if(length(par_to_opt_ind) > 1){
      cat(paste0("(", round(mean(abs(par[par_to_opt_ind] - optim_out$par)), 2), ")-"))
    } else {
      cat(paste0("(", round(par[par_to_opt_ind] - optim_out$par, 2), ")-"))
    }
    par[par_to_opt_ind] <- optim_out$par
  }
  
  
  
  #update tree
  optim_est <- par
  R_est <- diag(d_traits)
  R_est[upper.tri(R_est)] <- optim_est[(nTaxa*d_traits+2) : (nTaxa * d_traits + 1 + choose(d_traits, 2))];
  R_est <- R_est + t(R_est) - diag(d_traits)
  R_est_nearPD <- as.matrix(nearPD(R_est, corr = T, maxit = 1E5)$mat)
  est_means <- matrix(optim_est[1:(nTaxa*d_traits)], ncol = d_traits, nrow = nTaxa, byrow = T, dimnames = list(pops))
  eigensystem <- eigen(R_est)
  V <- eigensystem$vectors
  L <- eigensystem$values
  PC_Scores <- t(t(V) %*% t(est_means))
  est_means_PCs <- sapply(1:which.max(cumsum(L)), function(ev) PC_Scores[,ev] / L[ev]^0.5)
  est_means_PCs <- est_means_PCs[,1:sum(cumsum(eigensystem$values / sum(eigensystem$values)) < 0.99)]
  
  
  if(useNearPDforTree){
    mahDistMat <- mahaMatrix(est_means, R_est_nearPD, squared = F)
  } else {
    mahDistMat <- mahaMatrix(est_means_PCs, diag(ncol(est_means_PCs)), squared = F)
  }
  upgmaTree <- upgma(mahDistMat)
  njTree <- nj(mahDistMat)
  if(weighPCV){treeCV <- (treeCV_init + diag(nTaxa)) / 2} #weight with star phylogeny
  if(useUPGMA | (i <= useUPGMAforNRounds)){
    if(!useStarPhylogeny){
      treeCV <- cov2cor(vcv.phylo(upgmaTree))
    }
  } else {
    njTree$edge.length[njTree$edge.length < 0] <- 0
    njTree <- reroot(tree = njTree, node.number = which(njTree$tip.label == "NEAND"), position = njTree$edge.length[which(njTree$edge[,2] == which(njTree$tip.label == "NEAND"))] / 2)
    plot(njTree)
    if(!useStarPhylogeny){
      treeCV <- cov2cor(vcv.phylo(njTree))
    }
  }
  if(weighPCV){treeCV <- (treeCV + diag(nTaxa)) / 2} #weight with star phylogeny
  if(useStarPhylogeny){
    treeCV <- diag(nTaxa)
    rownames(treeCV) <- colnames(treeCV) <- pops
  }
  prunes <- prunePCV(treeCV)
  
  ############################
  #update imputed missing data
  ############################
  if(!ignoreImputedData){
    #iterate over each missing value and substitute in its likeliest state
    # current_params <- list(means = traits, cor = R, threshold_mat = thresholds) #for debugging
    current_params <- recompileMeansCorrsThreshes(par, dat)
    imputed_probs <- obs
    for(tip in 1:length(obs)){
      
      cat(paste0(ifelse(tip == 1, "\nimputing missing states: ", ""), tip, " "))
      missing <- which(!obs[[tip]], arr.ind = T) #rows are people, cols are traits
      missing_scores <- lapply(1:length(missing[,1]), function(x) NA)
      
      for(missdat in 1:length(missing[,1])){
        missinds <- missing[missdat,]
        other_missing_traits_in_individual <- missing[missing[,1] == missinds[1],2]
        other_missing_traits_in_individual <- setdiff(other_missing_traits_in_individual, missinds[2])
        possible_states <- 0:(n_thresholds[missinds[2]])
        
        #evaluate -logLik for all poss_states for that tip given the tip mean, the thresholds, and the correlation matrix
        scores <- evaluateStates(missinds = as.numeric(missinds),
                                 possible_states = possible_states, 
                                 tipMeans = current_params$means[tip,], 
                                 otherStates = traits_indiv_discr[[tip]], 
                                 threshold_mat = current_params$threshold_mat,
                                 other_missing_traits_in_individual = other_missing_traits_in_individual,
                                 cor = current_params$cor,
                                 ignoreOtherImputedStates = ignoreOtherImputedStates)
        
        
        missing_scores[[missdat]] <- scores
        if(!updateAllAtOnce){
          bestState <- which.max(scores) - 1
          traits_indiv_discr[[tip]][missinds[1], missinds[2]] <- bestState
        }
      }
      
      imputed_probs[[tip]] <- list(missing, missing_scores)
    }
    
    #find state-specific probabilities of missingness
    if(MNAR_Imputation){
      
      maxStates <- max(n_thresholds) + 1
      
      if(!use_conditional_probs_for_observed_data){
        n_obs_trait_state_per_tip <- lapply(1:length(obs), function(tip)
          t(sapply(1:ncol(obs[[tip]]), function(trait) countObsStates(traits_indiv_discr[[tip]][obs[[tip]][,trait],trait], maxStates = maxStates))
          )
        )
        
        #for debugging purposes
        # n_unobs_trait_state_per_tip <- lapply(1:length(obs), function(tip)
        #   t(sapply(1:ncol(obs[[tip]]), function(trait) countObsStates(traits_indiv_discr_true[[tip]][!obs[[tip]][,trait],trait], maxStates = maxStates))
        #   )
        # )
        
      } else if(use_conditional_probs_for_observed_data | (use_observed_numerator_ALL_conditionals_denominator & !use_univariate_normals_for_conditionals)){
        
        obs_probs <- obs
        for(tip in 1:length(obs)){
          
          cat(paste0(ifelse(tip == 1, "\nimputing observed states: ", ""), tip, " "))
          missing <- which(obs[[tip]], arr.ind = T) #rows are people, cols are traits
          observed <- which(!obs[[tip]], arr.ind = T) #rows are people, cols are traits
          observed_scores <- lapply(1:length(observed[,1]), function(x) NA)
          
          for(obsdat in 1:length(observed[,1])){
            obsinds <- observed[obsdat,]
            missing_traits_in_individual <- missing[missing[,1] == obsinds[1],2]
            possible_states <- 0:(n_thresholds[obsinds[2]])
            
            #evaluate -logLik for all poss_states for that tip given the tip mean, the thresholds, and the correlation matrix
            scores <- evaluateStates(missinds = as.numeric(obsinds),
                                     possible_states = possible_states, 
                                     tipMeans = current_params$means[tip,], 
                                     otherStates = traits_indiv_discr[[tip]], 
                                     threshold_mat = current_params$threshold_mat,
                                     other_missing_traits_in_individual = missing_traits_in_individual,
                                     cor = current_params$cor,
                                     ignoreOtherImputedStates = ignoreOtherImputedStates)
            
            
            observed_scores[[missdat]] <- scores
          }
          
          obs_probs[[tip]] <- list(observed, observed_scores)
        }
        
        if(use_conditional_probs_for_observed_data){
          n_obs_trait_state_per_tip <- lapply(1:length(obs_probs), function(tip)
            countUnobsStates(imputed_probs[[tip]][[1]][,2], imputed_probs[[tip]][[2]], maxStates = maxStates, ntraits = d_traits)
          )
        }
      }
      
      if(!use_observed_numerator_ALL_conditionals_denominator){
        
        if(use_conditional_probs_for_unobserved_data){
          
          n_imputed_trait_state_per_tip <- lapply(1:length(imputed_probs), function(tip)
            countUnobsStates(trait_inds = imputed_probs[[tip]][[1]][,2], state_probs = imputed_probs[[tip]][[2]], maxStates = maxStates, ntraits = d_traits)
          )
          
        } else if(!use_conditional_probs_for_unobserved_data){
          
          n_imputed_trait_state_per_tip <- lapply(1:length(imputed_probs), function(tip)
            countUnobsStates(trait_inds = imputed_probs[[tip]][[1]][,2], state_probs = imputed_probs[[tip]][[2]], maxStates = maxStates, ntraits = d_traits, use_opt = T)
          )
          
        }
        
      } else if(use_observed_numerator_ALL_conditionals_denominator){
        if(use_univariate_normals_for_conditionals){
          #TODO perform univariate estimation of all conditional frequencies
          
          conditional_expectation_per_tip <- lapply(1:length(obs), function(tip) t(sapply(1:ncol(current_params$means), function(trait)
            computeUnivNormalStateProbs(mean = current_params$means[tip, trait], thresholds = current_params$threshold_mat[trait,], npop = nrow(obs[[tip]])))
          ))
          
          conditional_probs_per_tip <- lapply(1:length(obs), function(tip) t(sapply(1:ncol(current_params$means), function(trait)
            computeUnivNormalStateProbs(mean = current_params$means[tip, trait], thresholds = current_params$threshold_mat[trait,]))
          ))
          
          if(use_conditional_multinomial_approx){
            
            conditionals_less_observed <-  lapply(1:length(obs), function(tip) t(sapply(1:ncol(current_params$means), function(trait)
              apply(approxMultinomialCond(probs = conditional_probs_per_tip[[tip]][trait,], obs_counts = n_obs_trait_state_per_tip[[tip]][trait,], 
                                          total_count = nrow(obs[[tip]]), min_sample_size = 500, raw = T, mcmc_rescue = 4, troubleshoot = F), 2, mean)
            )))
            
          } else if(use_conditional_expectation){
            
            conditional_distribution_per_tip <- conditional_expectation_per_tip
            
          }
          
        } else{
          
          n_obs_conds_state_per_tip <- lapply(1:length(obs_probs), function(tip)
            countUnobsStates(imputed_probs[[tip]][[1]][,2], imputed_probs[[tip]][[2]], maxStates = maxStates, ntraits = d_traits)
          )
          
          #combine the two sets of conditionals
          conditional_distribution_per_tip <- lapply(1:length(obs), function(tip)
            n_obs_conds_state_per_tip[[tip]] + n_imputed_trait_state_per_tip[[tip]]
          )
        }
      }
      
      if(use_observed_numerator_ALL_conditionals_denominator){
        
        if(!use_conditional_multinomial_approx){
          conditionals_less_observed <- lapply(1:length(obs), function(tip)
            computeConditionalsLessObserved(conditionals = conditional_distribution_per_tip[[tip]], observed = n_obs_trait_state_per_tip[[tip]], 
                                            npop = ifelse(adjust_conds_at_pop_level, nrow(obs[[tip]]), NA), addPop = add_individuals_to_missing)
          )
        }
        trait_specific_missing_probs <- find_missing_probs(n_obs_trait_state_per_tip, conditionals_less_observed, per_pop = F, raw = F, per_trait = T)
        print("estimated missing probs:")
        print(find_missing_probs(n_obs_trait_state_per_tip, conditionals_less_observed, per_pop = F, raw = T, per_trait = F))
        if(i == nrounds){
          write.table(find_missing_probs(n_obs_trait_state_per_tip, conditionals_less_observed, per_pop = F, raw = T, per_trait = F), 
                      file = paste0("estimated_missing_probs_replicate_", rep_i, "_optimrep_", optim_rep_i))
        }
        
        
      } else {
        
        trait_specific_missing_probs <- find_missing_probs(n_obs_trait_state_per_tip, n_imputed_trait_state_per_tip, per_pop = F, raw = F)
        print("estimated missing probs:")
        print(find_missing_probs(n_obs_trait_state_per_tip, n_imputed_trait_state_per_tip, per_pop = F, raw = F, per_trait = F))
        
      }
    }
    
    if(updateAllAtOnce){
      for(tip in 1:length(obs)){
        missing <- imputed_probs[[tip]][[1]]
        missing_scores <- imputed_probs[[tip]][[2]]
        for(missdat in 1:length(missing[,1])){
          missinds <- missing[missdat,]
          possible_states <- 0:(n_thresholds[missinds[2]])
          scores <- missing_scores[[missdat]]
          if(stochastic_imputation){
            stateProbs <- softmax(scores)
            if(MNAR_Imputation){
              missing_cond_probs <- trait_specific_missing_probs[missinds[2],1:length(stateProbs)]
            } else{
              missing_cond_probs <- rep(0.1, length(stateProbs))
            }
            stateProbs <- bayesian_normalize(stateProbs, missing_cond_probs)
            if(impose_constraints){
              if(any(stateProbs * constraintArray[[tip]][missinds[1], missinds[2],][possible_states + 1] != 0)){
                stateProbs <- bayesian_normalize(stateProbs, constraintArray[[tip]][missinds[1], missinds[2],][possible_states + 1])
              }
            }
            bestState <- sample(0:(length(stateProbs)-1), size = 1, prob = stateProbs)
          } else {
            bestState <- which.max(scores) - 1
          }
          traits_indiv_discr[[tip]][missinds[1], missinds[2]] <- bestState
        }
      }
    }
    
    traits_indiv_discr_unique <- lapply(pops, function(tip) uniqueSP(traits_indiv_discr[[tip]]))
    dat <- as.data.frame(do.call(rbind, lapply(1:npops, function(tip) traits_indiv_discr_unique[[tip]]$SPs)))
    colnames(dat) <- trait_names
    dat2 = data.frame(nSP = unlist(lapply(1:npops, function(tip) traits_indiv_discr_unique[[tip]]$counts)),
                      tip = unlist(lapply(1:npops, function(tip) rep(tip, length(traits_indiv_discr_unique[[tip]]$counts)))))
    dat_orig <- cbind(dat, dat2)
  }
  
}

#save current parameter estimates to disk
save(current_params, file = paste0("data_sim/param_estimates_replicate_", rep_i, "_optimrep_", optim_rep_i))
q(save = "no")