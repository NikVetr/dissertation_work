#set random seed 
set.seed(12345)

#set dimension of multinomial distribution
ntraits <- 10

#specify number of monte carlo replicates
n_mc_reps <- 1E5

#specify size of the observed multinomial sample that is being conditioned on 
size_obs <- 50

#specify the total number of draws from the multinomial distribution
size_tot <- 100

#specify the probabilities of each of the categories in the total distribution
probs <- seq(0.2, 0.8, length.out = ntraits)
probs <- probs / sum(probs)

#sample the observed values to be conditioned on that is specifically distinct from
#the overall multinomial distribution
obs <- t(rmultinom(n = 1, size = size_obs, prob = rep(0.5, ntraits)))

#perform monte-carlo simulation of total multinomial distribution
samps <- t(rmultinom(n = n_mc_reps, size = size_tot, prob = probs))

#identify which of these distributions are compatible with the observed values
diffs <- t(sapply(1:n_mc_reps, function(simsamp) samps[simsamp,] - obs))
valid <- sapply(1:n_mc_reps, function(simsamp) all(diffs[simsamp,] >= 0))

#subset the remaining amount of the total distribution compatible with the observed counts
valid_samps <- diffs[valid,]

#approximate the probabilities of the remaining, unobserved categories
cond_probs <- apply(valid_samps, 2, sum)
cond_counts <- cond_probs / nrow(valid_samps)
cond_probs <- cond_probs / sum(cond_probs)

#print these probabilities to screen
print("unconditional probabilities: "); cat(probs)
print("probabilities conditional on obs: "); cat(cond_probs)

#can we estimate the probs marginally?
marg_conds_ests <- rep(0,ntraits)
marg_nmc <- 1E5
for(i in 1:ntraits){
  marg_ind <- i
  marg_samps <- rbinom(n = marg_nmc, size = size_tot, prob = probs[marg_ind])
  marg_samps <- cbind(marg_samps, size_tot - marg_samps)
  marg_diffs <- marg_samps - cbind(rep(obs[marg_ind], marg_nmc), rep(sum(obs[-marg_ind]), marg_nmc))
  # marg_diffs <- marg_samps - cbind(rep(obs[marg_ind], marg_nmc), rep(0, marg_nmc))
  marg_valid_diffs <- marg_diffs[marg_diffs[,1] >= 0,]
  marg_cond_probs <- apply(marg_valid_diffs, 2, sum)
  marg_cond_probs <- marg_cond_probs / marg_nmc #sum(marg_cond_probs)
  marg_conds_ests[i] <- marg_cond_probs[1]
}
marg_conds_ests <- marg_conds_ests / sum(marg_conds_ests)
marg_conds_ests
cond_probs

#is there an exact solution to the binomial question?
exp_extras <- sapply(1:ntraits, function(marg_ind)
  sum(diff(pbinom(q = obs[marg_ind]:size_tot, size = size_tot, prob = probs[marg_ind])) * 1:(size_tot-obs[marg_ind]))
)
exp_extras / sum(exp_extras)

#higher counts are impossible in full multinomial distribution... maybe it's not the total size I need to look at, but the remaining size? 
exp_extras <- sapply(1:ntraits, function(marg_ind)
  sum(diff(pbinom(q = obs[marg_ind]:(size_tot-sum(obs[-marg_ind])), size = size_tot, prob = probs[marg_ind])) * 1:(size_tot-sum(obs)))
)
exp_extras / sum(exp_extras)


#what if we generate our multinomial draws iteratively? Maybe we can randomly choose our starting state? In proportion to their probs??
nmc_it <- 2000
valid_samps <- matrix(0, nrow = nmc_it, ncol = ntraits)
valid_samps[,1] <- size_tot
for(i in 1:(ntraits-1)){
  print(i)
  prob <- probs[i] / sum(probs[i:ntraits])
  samps <- rbinom(n = nmc_it, size = valid_samps[,i],  prob =  prob)
  if(any(samps < obs[i]) | any((valid_samps[,i] - samps) < sum(obs[(i+1):ntraits]))){
    some_invalid <- T 
    reroll <- unique(c(which(samps < obs[i]), which((valid_samps[,i] - samps) <= sum(obs[(i+1):ntraits]))))
    while(some_invalid){
      samps[reroll] <- rbinom(n = length(reroll), size = valid_samps[reroll,i],  prob =  prob)
      reroll <- reroll[unique(c(which(samps[reroll] < obs[i]), which((valid_samps[,i] - samps[reroll]) <= sum(obs[(i+1):ntraits]))))]
      if(length(reroll) == 0){
        some_invalid <- F
      }
    }
  }
  valid_samps[,i+1] <- valid_samps[,i] - samps
  valid_samps[,i] <- samps
}
valid_diffs <- t(sapply(1:nmc_it, function(simsamp) valid_samps[simsamp,] - obs))
cond_probs_it <- apply(valid_diffs, 2, sum)
cond_probs_it <- cond_probs_it / sum(cond_probs_it)
cond_probs_it
cond_probs


#let's write a function to do the brute force approach
approxMultinomialCond <- function(probs, obs_counts, total_count, min_sample_size = 10, raw = F){
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
approxMultinomialCond(probs, obs, total_count = 100, min_sample_size = 100, raw = F)
microbenchmark(approxMultinomialCond(probs, obs, total_count = 100, min_sample_size = 500, raw = T))

##################################################
# let's write a more elegant metropolis sampler? #
##################################################

min_ess = 500
n_iter <- 1E5
thin <- 1E1
thin_print <- n_iter / 100
samp_counts <- matrix(0, nrow = n_iter / thin, ncol = ntraits)
curr_state <- obs + t(rmultinom(n = 1, size = size_tot - size_obs, prob = probs))
curr_prob <- dmultinom(curr_state, prob = probs, log = T)
accept_count <- 0
for(i in 1:n_iter){
  if(i %% thin_print == 0){cat(paste0(i / thin_print, " "))}
  amount_to_swap <- sample(size = 1, 1:floor((size_tot - size_obs)/5))
  inds_to_swap <- sample(size = 2, 1:ntraits, replace = F)
  prop_state <- curr_state
  prop_state[inds_to_swap[1]] <- prop_state[inds_to_swap[1]] + amount_to_swap
  prop_state[inds_to_swap[2]] <- prop_state[inds_to_swap[2]] - amount_to_swap
  if(all(prop_state >= obs)){
    prop_prob <- dmultinom(prop_state, prob = probs, log = T)
  } else {
    prop_prob <- -Inf
  }
  if(exp(prop_prob - curr_prob) > runif(1,0,1)){
    curr_state <- prop_state
    curr_prob <- prop_prob
    accept_count <- accept_count + 1
  }
  if(i %% thin == 0){
    samp_counts[i/thin,] <- curr_state
  }
}
accept_count / n_iter
samp_counts <- samp_counts[ceiling(0.2*nrow(samp_counts)):nrow(samp_counts),]
diff_counts <- t(sapply(1:nrow(samp_counts), function(count) samp_counts[count,] - obs))
exp_extras_mcmc <- apply(diff_counts, 2, mean)
exp_extras_mcmc / sum(exp_extras_mcmc)
cond_probs
exp_extras_mcmc
cond_counts
coda::effectiveSize(diff_counts)


#ok, let's wrap it into a function
approxMultinomialCondMCMC <- function(probs, obs_counts, total_count, min_ess = 100, start_niter = 1E5, thin = 1E1, burnin = 0.2, exp_count = T, samples = T, swap_num_divisor = 5){
  ntraits <- length(probs)
  size_obs <- sum(obs_counts)
  samp_counts <- matrix(0, nrow = start_niter / thin, ncol = ntraits)
  curr_state <- obs_counts + t(rmultinom(n = 1, size = total_count - size_obs, prob = probs))
  curr_prob <- dmultinom(curr_state, prob = probs, log = T)
  accept_count <- 0
  for(i in 1:start_niter){
    # if(i %% thin_print == 0){cat(paste0(i / thin_print, " "))}
    amount_to_swap <- sample(size = 1, 1:floor((total_count - size_obs)/swap_num_divisor))
    inds_to_swap <- sample(size = 2, 1:ntraits, replace = F)
    prop_state <- curr_state
    prop_state[inds_to_swap[1]] <- prop_state[inds_to_swap[1]] + amount_to_swap
    prop_state[inds_to_swap[2]] <- prop_state[inds_to_swap[2]] - amount_to_swap
    if(all(prop_state >= obs)){
      prop_prob <- dmultinom(prop_state, prob = probs, log = T)
      if(exp(prop_prob - curr_prob) > runif(1,0,1)){
        curr_state <- prop_state
        curr_prob <- prop_prob
        accept_count <- accept_count + 1
      }
    }# else {prop_prob <- -Inf} ...implied
    
    if(i %% thin == 0){
      samp_counts[i/thin,] <- curr_state
    }
  }
  samp_counts <- samp_counts[ceiling(burnin*nrow(samp_counts)):nrow(samp_counts),]
  diff_counts <- t(sapply(1:nrow(samp_counts), function(count) samp_counts[count,] - obs))
  ess <- coda::effectiveSize(diff_counts)
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
     followup_niter <- floor((min_ess / min(ess)) * start_niter)
     followup_samp_counts <- matrix(0, nrow = floor(followup_niter / thin), ncol = ntraits)
     for(i in 1:followup_niter){
       amount_to_swap <- sample(size = 1, 1:floor((total_count - size_obs)/swap_num_divisor))
       inds_to_swap <- sample(size = 2, 1:ntraits, replace = F)
       prop_state <- curr_state
       prop_state[inds_to_swap[1]] <- prop_state[inds_to_swap[1]] + amount_to_swap
       prop_state[inds_to_swap[2]] <- prop_state[inds_to_swap[2]] - amount_to_swap
       if(all(prop_state >= obs)){
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
     followup_diff_counts <- t(sapply(1:nrow(followup_samp_counts), function(count) followup_samp_counts[count,] - obs))
     diff_counts <- rbind(diff_counts, followup_diff_counts)
     ess <- coda::effectiveSize(diff_counts)
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

approxMultinomialCondMCMC(probs = probs, obs_counts = obs, total_count = size_tot, min_ess = 100, start_niter = 1000)
