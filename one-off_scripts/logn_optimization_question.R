
#initialize simulation parameters
d_traits <- 1E2 #number of characters to iptimize
logn_mu <- 0.5 #mean of normal to be exponentiated
logn_sd <- 0.5 #sd of normal to be exponentiated
n_samps <- 1E3 #number of replicates
resid_sd <- rep(0, n_samps) #initialize vector to store residual 'sd's
resid_mu <- rep(0, n_samps) #initialize vector to store residual 'mu's
method <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent")[4]

for(i in 1:n_samps){
  
  if(i %% (n_samps/100) == 0){cat(paste0(i/(n_samps/100), "%... "))}
  
  #sample number of thresholds for simulated traits
  n_thresholds <- sample(1:5, d_traits, T) 
  
  #sample sum of threshold spacings from lognormal distribution, and multiply by dirichlet proportions
  threshold_spacings <- sapply(1:d_traits, function(trait) as.vector(exp(rnorm(n = 1, mean = logn_mu, sd = logn_sd)) * rdirichlet(1, alpha = rep(1, (n_thresholds[trait]-1)))))
  
  #store sums for later confirmation
  thresh_sums <- sapply(1:length(threshold_spacings), function(trait) sum(threshold_spacings[[trait]]))
  
  #traits with just one threshold have infinitely far apart thresholds
  threshold_spacings[n_thresholds == 1] <- Inf
  
  #to identify dimensionality of matrix for storage
  max_thresholds_count <- max(n_thresholds)
  
  #pad traits with submaximal thresholds counts with infinities on the right
  threshold_spacings <- t(sapply(1:length(threshold_spacings), function(trait) c(threshold_spacings[[trait]], rep(Inf, max_thresholds_count - 1 - length(threshold_spacings[[trait]])))))
  
  #set first threshold at 0
  threshold_spacings <- cbind(rep(0, d_traits), threshold_spacings)
  
  #convert threshold spacings to thresholds
  threshold_mat <- t(apply(threshold_spacings, 1, cumsum))
  
  #go back and deconvert into spacings and then sum
  thresh_sums_of_spacings <- sapply(1:d_traits, function(trait) sum(diff(threshold_mat[trait,][threshold_mat[trait,] != Inf])))
  
  #confirm these are the same as previously calculated thresh_sums
  all((thresh_sums_of_spacings - thresh_sums) < 1E-6)
  
  #get rid of '0's while accounting for floating point errors
  dat <- thresh_sums_of_spacings[thresh_sums_of_spacings > 1E-6]
  
  #initialize the two parameters
  pars_to_opt <- c(1,1)
  
  #specify function to be optimized
  dens_in_logn <- function(pars_to_opt, dat){-sum(dnorm(x = log(dat), mean = pars_to_opt[1], sd = pars_to_opt[2], log = T))}
  
  #confirm function is operational
  dens_in_logn(pars_to_opt, dat)
  
  #perform optimization
  pars <- optim(par = pars_to_opt, dens_in_logn, dat = dat, method = "L-BFGS-B", upper = c(Inf, Inf), lower = c(-Inf, 1E-3), control = list(trace = 0, REPORT = 1))$par
  
  #store residual sd from true sd 
  resid_mu[i] <- pars[1] - logn_mu
  resid_sd[i] <- pars[2] - logn_sd 
  
}

par(mfrow = c(1,2))
hist(resid_sd)
sum(resid_sd > 0) / n_samps

hist(resid_mu)
sum(resid_mu > 0) / n_samps
