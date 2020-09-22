
betaMove <- function(simplex, index, conc, offset = 1){
  #concentration param needed to increasingly concentrate probability dist on current value, but offset needed bc sometimes these are still <<1
  proposed_value <- rbeta(1, shape1 = simplex[index] * length(simplex) * conc + offset, shape2 = (1 - simplex[index]) * length(simplex) * conc + offset) #sample from beta
  proposed_simplex <- simplex
  proposed_simplex[index] <- proposed_value
  proposed_simplex[-index] <- proposed_simplex[-index] / sum(proposed_simplex[-index]) * (1-proposed_value) #rescale to unit sum
  
  forward_prob <- dbeta(proposed_value, shape1 = simplex[index] * length(simplex) * conc + offset, shape2 = (1 - simplex[index]) * length(simplex) * conc + offset, log = T)
  backward_prob <- dbeta(simplex[index], shape1 = proposed_value * length(simplex) * conc + offset, shape2 = (1 - proposed_value) * length(simplex) * conc + offset, log = T)
  
  return(list(prop = proposed_simplex, log_prop_ratio = (forward_prob - backward_prob)))
  
}

s1 <- 3 #shape parameter #1
s2 <- 5 #shape parameter #2
rbeta(1, s1, s2) #random draw from this distribution

curr_nums <- c(0.5,0.5)
curr_dens <- dbeta(curr_nums[1], s1, s2, log = T)

niter <- 5E6
thin = 1E3
nums <- matrix(0, ncol = 2, nrow = niter / thin)
for(i in 1:niter){
  if(i %% thin == 0){cat(paste0(i/thin, " "))}
  num_to_prop <- sample(1:2, 1) #which of the two indices of the simplex to make a proposal to?
  prop_set <- betaMove(curr_nums, index = num_to_prop, conc = 2) #get proposed value and log proposal ratio
  prop_nums <- prop_set$prop #these are the proposed values
  prop_rat <- prop_set$log_prop_ratio #this is the proposal ratio
  prop_dens <- dbeta(prop_nums[1], s1, s2, log = T) #find the proposed probability in the target distribution
  dens_rat <- prop_dens - curr_dens #find the log ratio of target densities
  log_accept_prob <- prop_rat + dens_rat #find the log acceptance probability
  if(log(runif(1, 0, 1)) < log_accept_prob){ #roll a 1d\infty to see if move is accepted
    curr_nums <- prop_nums
    curr_dens <- prop_dens
  }
  if(i %% thin == 0){
    nums[i/thin,] <- curr_nums
  }
}

hist(nums[,1], col = rgb(1,0,0,0.5), probability = T, ylim = c(0,3))
hist(rbeta(1E4, s1, s2), col = rgb(0,1,0,0.5), probability = T, add = T)
#...sometimes works, sometimes does not

#################################
# let's try a sliding move out? #
#################################

slideMove <- function(value, window, min = -Inf, max = Inf){
  
  low_bound <- max(value - window/2, min)
  high_bound <- min(value + window/2, max)
  prop <- runif(1, min = low_bound, max = high_bound)
  
  forward_log_prob <- log(high_bound - low_bound)
  backward_low_bound <- max(prop - window/2, min)
  backward_high_bound <- min(prop + window/2, max)
  backward_log_prob <- log(backward_high_bound - backward_low_bound)
  
  return(list(prop = prop, log_prop_ratio = forward_log_prob - backward_log_prob ))
  
}

s1 <- 3
s2 <- 5
rbeta(1, s1, s2)

curr_nums <- c(0.5,0.5)
curr_dens <- dbeta(curr_nums[1], s1, s2, log = T)

niter <- 5E6
thin = 5E2
nums <- matrix(0, ncol = 2, nrow = niter / thin)
for(i in 1:niter){
  if(i %% thin == 0){cat(paste0(i/thin, " "))}
  num_to_prop <- sample(1:2, 1)
  prop_set <- slideMove(curr_nums[num_to_prop], window = 0.1, min = 0, max = 1)$prop
  if(num_to_prop == 1){
    prop_nums <- c(prop_set, 1-prop_set)
  }  else {
    prop_nums <- c(1-prop_set, prop_set)
  }
  prop_rat <- 0
  prop_dens <- dbeta(prop_nums[1], s1, s2, log = T) #target density function will auto-reject proposals outside (0,1)
  dens_rat <- prop_dens - curr_dens
  log_accept_prob <- prop_rat + dens_rat
  if(log(runif(1, 0, 1)) < log_accept_prob){
    curr_nums <- prop_nums
    curr_dens <- prop_dens
  }
  if(i %% thin == 0){
    nums[i/thin,] <- curr_nums
  }
}

analytical_distribution <- rbeta(1E5, s1, s2)
hist(nums[,1], col = rgb(1,0,0,0.5), probability = T, ylim = c(0,3), breaks = 10)
hist(analytical_distribution, col = rgb(0,1,0,0.5), probability = T, add = T, breaks = 10)


coda::effectiveSize(c(nums[,1], analytical_distribution[1:nrow(nums)]))
coda::gelman.diag(list(coda::mcmc(nums[,1]), coda::mcmc(analytical_distribution[1:nrow(nums)])))

#yes this does appear to work more consistently

####################################
#OK let's try a higher-dim simplex??
####################################

sps <- c(2,4,1,7,4,2,6,2,11)

curr_nums <- runif(n = length(sps), min = 0, max = 1)
curr_nums <- curr_nums / sum(curr_nums)
curr_dens <- DirichletReg::ddirichlet(x = t(as.matrix(curr_nums)), alpha = sps, log = T)

niter <- 5E6
thin = 1E3
n_out <- niter / thin
burnin_prop <- 0.3
nums <- matrix(0, ncol = length(sps), nrow = n_out)
for(i in 1:niter){
  if(i %% thin == 0){cat(paste0(i/thin, " "))}
  num_to_prop <- sample(1:length(sps), 1) #which of the two indices of the simplex to make a proposal to?
  prop_set <- betaMove(curr_nums, index = num_to_prop, conc = 2) #get proposed value and log proposal ratio
  prop_nums <- prop_set$prop #these are the proposed values
  prop_rat <- prop_set$log_prop_ratio #this is the proposal ratio
  prop_dens <- DirichletReg::ddirichlet(x = t(as.matrix(prop_nums)), alpha = sps, log = T) #find the proposed probability in the target distribution
  dens_rat <- prop_dens - curr_dens #find the log ratio of target densities
  log_accept_prob <- prop_rat + dens_rat #find the log acceptance probability
  if(log(runif(1, 0, 1)) < log_accept_prob){ #roll a 1d\infty to see if move is accepted
    curr_nums <- prop_nums
    curr_dens <- prop_dens
  }
  if(i %% thin == 0){
    nums[i/thin,] <- curr_nums
  }
}
nums <- nums[(n_out * burnin_prop):n_out,]
ess <- coda::effectiveSize(nums)
analytical_distribution <- DirichletReg::rdirichlet(1E5, sps)
par(mfrow = c(ceiling(sqrt(length(sps))), ceiling(sqrt(length(sps)))))
for(i in 1:length(sps)){
  hist(nums[,i], col = rgb(1,0,0,0.5), probability = T, breaks = 10, 
       main = latex2exp::TeX(paste0("index ", i, ", ESS = ", round(ess[i]), ", $\\alpha$ = ", sps[i])))
  hist(analytical_distribution[,i], col = rgb(0,1,0,0.5), probability = T, add = T, breaks = 10)
}
coda::effectiveSize(rbind(nums, analytical_distribution[1:nrow(nums),]))
coda::gelman.diag(list(coda::mcmc(nums), coda::mcmc(analytical_distribution[1:nrow(nums),])))

####################################
#and a sliding move on this one??
####################################

sps <- c(2,4,1,7,4,2,6,2,11)

curr_nums <- runif(n = length(sps), min = 0, max = 1)
curr_nums <- curr_nums / sum(curr_nums)
curr_dens <- DirichletReg::ddirichlet(x = t(as.matrix(curr_nums)), alpha = sps, log = T)

niter <- 5E6
thin = 1E3
n_out <- niter / thin
burnin_prop <- 0.3
nums <- matrix(0, ncol = length(sps), nrow = n_out)
for(i in 1:niter){
  if(i %% thin == 0){cat(paste0(i/thin, " "))}
  num_to_prop <- sample(1:length(sps), 1) #which of the two indices of the simplex to make a proposal to?
  prop_set <- slideMove(curr_nums[num_to_prop], window = 0.05)
  prop_rat <- prop_set$log_prop_ratio
  prop_set <- prop_set$prop
  prop_nums <- curr_nums
  prop_nums[num_to_prop] <- prop_set
  prop_nums[-num_to_prop] <- prop_nums[-num_to_prop] / sum(prop_nums[-num_to_prop]) * (1-prop_set)
  
  if(prop_set > 0 & prop_set < 1){
    prop_dens <- DirichletReg::ddirichlet(x = t(as.matrix(prop_nums)), alpha = sps, log = T) #find the proposed probability in the target distribution
  } else {
    prop_dens = -Inf
  }
  dens_rat <- prop_dens - curr_dens #find the log ratio of target densities
  log_accept_prob <- prop_rat + dens_rat #find the log acceptance probability
  if(log(runif(1, 0, 1)) < log_accept_prob){ #roll a 1d\infty to see if move is accepted
    curr_nums <- prop_nums
    curr_dens <- prop_dens
  }
  if(i %% thin == 0){
    nums[i/thin,] <- curr_nums
  }
}
nums <- nums[(n_out * burnin_prop):n_out,]
ess <- coda::effectiveSize(nums)
analytical_distribution <- DirichletReg::rdirichlet(1E5, sps)
par(mfrow = c(ceiling(sqrt(length(sps))), ceiling(sqrt(length(sps)))))
for(i in 1:length(sps)){
  hist(nums[,i], col = rgb(1,0,0,0.5), probability = T, breaks = 100, 
       main = latex2exp::TeX(paste0("index ", i, ", ESS = ", round(ess[i]), ", $\\alpha$ = ", sps[i])))
  hist(analytical_distribution[,i], col = rgb(0,1,0,0.5), probability = T, add = T, breaks = 100)
}

coda::effectiveSize(rbind(nums, analytical_distribution[1:nrow(nums),]))

coda::gelman.diag(list(coda::mcmc(nums[,1:9]), coda::mcmc(analytical_distribution[1:nrow(nums),1:9])))


#confirm slide move is ok in uniform context


slideMove <- function(value, window, min = -Inf, max = Inf){
  
  low_bound <- max(value - window/2, min)
  high_bound <- min(value + window/2, max)
  prop <- runif(1, min = low_bound, max = high_bound)
  
  forward_log_prob <- log(high_bound - low_bound)
  backward_low_bound <- max(prop - window/2, min)
  backward_high_bound <- min(prop + window/2, max)
  backward_log_prob <- log(backward_high_bound - backward_low_bound)
  
  return(list(prop = prop, log_prop_ratio = - forward_log_prob + backward_log_prob ))
  
}


curr_nums = 0
curr_dens <- dnorm(curr_nums, 0, 1)

niter <- 1E6
thin = 5E2
nums <- matrix(0, ncol = 1, nrow = niter / thin)
for(i in 1:niter){
  if(i %% thin == 0){cat(paste0(i/thin, " "))}
  prop_set <- slideMove(curr_nums, window = 0.5, min = -0.5, max = 0.5)
  prop_nums <- prop_set$prop
  prop_rat = prop_set$log_prop_ratio
  prop_dens <- dnorm(prop_nums, 0, 1) #find the proposed probability in the target distribution
  dens_rat <- prop_dens - curr_dens #find the log ratio of target densities
  log_accept_prob <- prop_rat + dens_rat #find the log acceptance probability
  if(log(runif(1, 0, 1)) < log_accept_prob){ #roll a 1d\infty to see if move is accepted
    curr_nums <- prop_nums
    curr_dens <- prop_dens
  }
  if(i %% thin == 0){
    nums[i/thin,] <- curr_nums
  }
}
hist(nums, breaks = 10, probability = T)
mean(nums)
sd(nums)
x <- rnorm(1E6, 0, 1)
x <- x[x<0.5]
x <- x[x>-0.5]

hist(x, breaks = 10)
sd(x)
sd(nums)

#rejection sample sliding window move

curr_nums = 0
curr_dens <- dnorm(curr_nums, 0, 1)

niter <- 5E6
thin = 5E2
nums <- matrix(0, ncol = 1, nrow = niter / thin)
for(i in 1:niter){
  if(i %% thin == 0){cat(paste0(i/thin, " "))}
  prop_set <- slideMove(curr_nums, window = 0.5, min = -Inf, max = Inf)
  prop_nums <- prop_set$prop
  prop_rat = prop_set$log_prop_ratio
  if(prop_nums > -1 & prop_nums < 1){
    prop_dens <- dnorm(prop_nums, 0, 1) #find the proposed probability in the target distribution
  } else {
    prop_dens <- -Inf
  }
  dens_rat <- prop_dens - curr_dens #find the log ratio of target densities
  log_accept_prob <- prop_rat + dens_rat #find the log acceptance probability
  if(log(runif(1, 0, 1)) < log_accept_prob){ #roll a 1d\infty to see if move is accepted
    curr_nums <- prop_nums
    curr_dens <- prop_dens
  }
  if(i %% thin == 0){
    nums[i/thin,] <- curr_nums
  }
}
hist(nums, breaks = 10, probability = T)
mean(nums)
sd(nums)
x <- rnorm(1E6, 0, 1)
x <- x[x<1]
x <- x[x>-1]

hist(x, breaks = 10)
sd(x)
sd(nums)

#######################################################################
#######################################################################
#let's try to sample from a uniform (0,1) with either distribution eh?
#######################################################################
#######################################################################

curr_nums <- c(0.5, 0.5)
curr_dens <- 1

niter <- 5E5
thin = 1E2
nums <- matrix(0, ncol = 2, nrow = niter / thin)
for(i in 1:niter){
  if(i %% thin == 0){cat(paste0(i/thin, " "))}
  num_to_prop <- sample(1:2, 1) #which of the two indices of the simplex to make a proposal to?
  prop_set <- betaMove(curr_nums, index = num_to_prop, conc = 2) #get proposed value and log proposal ratio
  prop_nums <- prop_set$prop #these are the proposed values
  prop_rat <- prop_set$log_prop_ratio #this is the proposal ratio
  prop_dens <- 1
  dens_rat <- prop_dens - curr_dens #find the log ratio of target densities
  log_accept_prob <- prop_rat + dens_rat #find the log acceptance probability
  if(log(runif(1, 0, 1)) < log_accept_prob){ #roll a 1d\infty to see if move is accepted
    curr_nums <- prop_nums
    curr_dens <- prop_dens
  }
  if(i %% thin == 0){
    nums[i/thin,] <- curr_nums
  }
}

hist(nums[,1], col = rgb(1,0,0,0.5), probability = T, ylim = c(0,3))
#...sometimes works, sometimes does not

#and now a sliding move?

slideMove <- function(value, window, min = -Inf, max = Inf){
  
  low_bound <- max(value - window/2, min)
  high_bound <- min(value + window/2, max)
  prop <- runif(1, min = low_bound, max = high_bound)
  
  forward_log_prob <- log(1 / (high_bound - low_bound))
  backward_low_bound <- max(prop - window/2, min)
  backward_high_bound <- min(prop + window/2, max)
  backward_log_prob <- log(1 / (backward_high_bound - backward_low_bound))
  
  return(list(prop = prop, log_prop_ratio = backward_log_prob - forward_log_prob))
  
}

curr_nums <- c(runif(1,0,1),0.5)
curr_nums[2] <- 1-curr_nums[1]
curr_dens <- 1

niter <- 5E6
thin = 1E2
nums <- matrix(0, ncol = 2, nrow = niter / thin)
for(i in 1:niter){
  if(i %% thin == 0){cat(paste0(i/thin, " "))}
  num_to_prop <- sample(1:2, 1)
  prop_set <- slideMove(curr_nums[num_to_prop], window = 0.1, min = 0, max = 1)
  if(num_to_prop == 1){
    prop_nums <- c(prop_set$prop, 1-prop_set$prop)
  }  else {
    prop_nums <- c(1-prop_set$prop, prop_set$prop)
  }
  prop_rat <- prop_set$log_prop_ratio
  prop_dens <- 1
  dens_rat <- prop_dens - curr_dens
  log_accept_prob <- prop_rat + dens_rat
  if(log(runif(1, 0, 1)) < log_accept_prob){
    curr_nums <- prop_nums
    curr_dens <- prop_dens
  }
  if(i %% thin == 0){
    nums[i/thin,] <- curr_nums
  }
}

hist(nums[,1], col = rgb(1,0,0,0.5), probability = T, ylim = c(0,3), breaks = 10)
