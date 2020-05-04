#check met hastings operation on correlation matrix with wishart move
library(MCMCpack)
shape_1 <- 15
shape_2 <- 86
par(mfrow = c(2,2))
hist(rbeta(1E6, shape1 = shape_1, shape2 = shape_2) * 2 - 1, breaks = 1E2, main = "target distribution", xlim = c(-1,1), xlab = "")

corr_init <- 0.001
n_iter <- 2E4
corrs <- rep(0, n_iter)
corrs[1] <- corr_init
df_tune <- 20
useCorrMatDens <- F
for(i in 2:n_iter){
  if(i %% (n_iter / 100) == 0){cat(paste0(i / n_iter * 100, "%  "))}
  cov_prop <- (rwish(v = df_tune, S = matrix(c(1,corrs[i-1], corrs[i-1], 1), 2, 2)) / df_tune)
  corr_prop <- cov2cor(cov_prop)[1,2]
  if(useCorrMatDens){
    log_prop_ratio <-  log(dwish(W = matrix(c(1,corrs[i-1], corrs[i-1], 1), 2, 2) * df_tune, v = df_tune, S = cov2cor(cov_prop))) - 
      log(dwish(W = cov2cor(cov_prop) * df_tune, v = df_tune, S = matrix(c(1,corrs[i-1], corrs[i-1], 1), 2, 2)))
  } else {
    log_prop_ratio <-  log(dwish(W = matrix(c(1,corrs[i-1], corrs[i-1], 1), 2, 2) * df_tune, v = df_tune, S = cov_prop)) - 
      log(dwish(W = cov_prop * df_tune, v = df_tune, S = matrix(c(1,corrs[i-1], corrs[i-1], 1), 2, 2)))
  }
  log_dens_ratio <- dbeta(x = (corr_prop + 1) / 2, shape1 = shape_1, shape2 = shape_2, log = T) - 
    dbeta(x = (corrs[i-1] + 1) / 2, shape1 = shape_1, shape2 = shape_2, log = T)
  log_accept_prob <- log_prop_ratio + log_dens_ratio
  if(log(runif(1, 0, 1)) < log_accept_prob){
    corrs[i] <- corr_prop
  } else {
    corrs[i] <- corrs[i-1]
  }
}

hist(corrs[(n_iter / 5) : n_iter], breaks = 1E2, main = "sampled distribution", xlim = c(-1,1))
plot(corrs, type = "l", main = "trace-plot"); lines(corrs[1:(n_iter/5)], col = 2, ylim = c(-1,1)); abline(v = n_iter / 5, lty = 2)
quantile(corrs[(n_iter / 5) : n_iter], probs = c(1:9)/10)
qbeta(p = c(1:9)/10, shape1 = shape_1, shape2 = shape_2) * 2 - 1
plot(x = quantile(corrs[(n_iter / 5) : n_iter], probs = c(1:99)/100),
     y = qbeta(p = c(1:99)/100, shape1 = shape_1, shape2 = shape_2) * 2 - 1,
     type = "l", main = "Q-Q Plot", xlab = "samples", ylab = "target"); abline(a = 0, b = 1, lty = 2, col = 3)


#how well does a blooming onion proposal work?

blooming_onion <- function(cor = 1, hmean = 0, eta = 1, d){
  r <- diag(rep(1,d))
  if(!identical(cor, 1)){d_s <- dim(cor)[1]} else {d_s <- 1}
  B <- eta + (d-1)/2
  if(d_s == 1){
    u <- rbeta(1, shape1 = B,shape2 = B)
    rij <- 2*u-1
    r[1,2] <- r[2,1] <- rij
  } else {
    r[1:d_s,1:d_s] <- cor
    alpha <- alpha - (d_s - 2) / 2
  }
  for(m in (d_s):(d-1)){
    B <- B - 1/2
    y <- rbeta(1, m/2, B)
    u <- rnorm(m, mean = hmean, sd = 1) 
    u <- u / sqrt(sum(u^2)) 
    w <- sqrt(y)*u
    q <- as.vector(t(chol(r[1:m,1:m])) %*% w)
    r[(m+1),1:m] <- r[1:m, (m+1)] <- q
  }
  return(r)
}

r <- rlkj(10)
lc <- t(chol(r))
lc_i <- solve(lc)
target_w = lc_i[1:9, 1:9] %*% r[10, 1:9]
target_u = target_w / crossprod(target_w)[1,1]
target_y = crossprod(target_w)[1,1]
y <- target_y
u <- as.vector(target_u)
u <- u / sqrt(sum(u^2)) 
w <- sqrt(y)*u
q <- as.vector(lc[1:9, 1:9] %*% w)
q - r[10, 1:9]

y_range <- 1:99/100
for(y in y_range){
  u <- as.vector(target_u)
  u <- u / sqrt(sum(u^2)) 
  w <- sqrt(y)*u
  q <- as.vector(lc[1:9, 1:9] %*% w)
  print(paste0("true_y = ", round(target_y, 2), ", current_y = ", y, " diff in y = ", round(target_y - y, 2), 
               ", diff. in corrs. = ", paste0(round(q - r[10, 1:9], digits = 2), collapse = " ")))
}

blooming_onion_2 <- function (cor = 1, d, eta = 1) {
  if(identical(cor, 1)){d_s <- 1}else{d_s <- dim(cor)[1]}
  R <- matrix(0, d, d)
  alpha <- eta + (d - 2) / 2
  if(d_s == 1){
    r12 <- 2 * rbeta(1, alpha, alpha) - 1
    R[1, 1] <- 1
    R[1, 2] <- r12
    R[2, 2] <- sqrt(1 - r12^2)
  } else {
    R[1:d_s, 1:d_s] <- chol(cor)
    alpha <- alpha - (d_s - 1) / 2
  }
  if (d > 2) {
    for (m in (d_s):(d - 1)) {
      alpha <- alpha - 0.5
      y <- rbeta(1, m/2, alpha)
      z <- rnorm(m, 0, 1)
      z <- z/sqrt(crossprod(z)[1])
      R[1:m, m + 1] <- sqrt(y) * z
      R[m + 1, m + 1] <- sqrt(1 - y)
    }
  }
  return(crossprod(R))
}


blooming_onion_3 <- function (cor, varB, varN) {
  d <- dim(cor)[1]
  R <- matrix(0, d, d)
  m <- d - 1
  chol_cor <- chol(cor)
  R[1:m, 1:m] <- chol_cor[1:m,1:m]
  target_chol <- chol_cor[1:d,d]
  target_y <- 1 - (target_chol[d]^2)
  target_z <- target_chol[1:m] / sqrt(target_y)
  s1 <- ((1-target_y) / varB - 1 / target_y) * target_y^2
  s2 <- s1 * (1 / target_y - 1)
  y <- rbeta(1, shape1 = s1, shape2 = s2)
  z <- rnorm(m, target_z, sqrt(varN))
  z <- z/sqrt(crossprod(z)[1])
  R[1:m, m + 1] <- sqrt(y) * z
  R[m + 1, m + 1] <- sqrt(1 - y)
  return(crossprod(R))
}

blooming_onion_3(cor, varB = 0.001, varN = 0.001)

min(eigen(blooming_onion_2(cor = matrix(c(1,0.6,0.6,1), 2, 2), d = 100, hmean = 1))$values)
dim <- 20; hist(blooming_onion_2(cor = matrix(c(1,0.6,0.6,1), 2, 2), d = dim, hmean = 0)[upper.tri(diag(dim))])
hist(replicate(1E4, blooming_onion_2(d = 2, eta = 1)[1,2]))


dim <- 5
cor_orig <- blooming_onion_2(d = dim)
cor_sub <- cor_orig[1:(dim-1), 1:(dim-1)]
cor_new <- blooming_onion_2(d = dim, cor = cor_sub, hmean = chol(cor_orig)[1:(dim-1),dim], hsd = 0.01)
max(abs(cor_orig - cor_new))


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


#dirichlet proposal mechanism
library(rethinking)
library(MCMCpack)

tp <- 100 #set the tuning parameter
dim <- 10 # set the correlation matrix dimension
cor <- rlkjcorr(1, dim) #sample correlation matrix from a flat lkj
cc <- chol(cor) #compute upper right cholesky factor of sampled matrix
apply(cc, 2, crossprod) #verify that columns are of unit length
new_squares <- sapply(1:dim, function(x) rdirichlet(1, tp * cc[,x]^2)) #sample new columns from tuned dirichlet
prop_cc <- sqrt(new_squares) #compute the element-wise square root
prop_cc[cc < 0] <- -prop_cc[cc < 0] # coerce back into original sign
prop_cor <- crossprod(prop_cc) #compute proposed correlation matrix
all(eigen(prop_cor)$values > 0) #confirm positive semidefiniteness
all(abs(diag(prop_cor) - 1) < 1E-6) #confirm unit diagonal
cor - prop_cor #compute deviations in correlation from original matrix
max(abs(cor - prop_cor)) #compute size of maximum deviation

#now let's simulate to find typical behavior
tp <- 10^(0:16/2)
nrep <- 300
dim <- 10

max_abs_diff <- matrix(0, length(tp), nrep)

for(i in 1:length(tp)){
  cor <- replicate(nrep, rlkjcorr(1, dim))
  cc <- lapply(1:nrep, function(x) chol(cor[,,x]))
  new_squares <- lapply(1:nrep, function(y)
    sapply(1:dim, function(x) rdirichlet(1, tp[i] * cc[[y]][,x]^2)))
  prop_cc <- lapply(1:nrep, function(x) sqrt(new_squares[[x]]))
  for(j in 1:nrep){
    prop_cc[[j]][cc[[j]] < 0] <- -prop_cc[[j]][cc[[j]] < 0]
  }
  prop_cor <- lapply(1:nrep, function(x) crossprod(prop_cc[[x]]))
  max_abs_diff[i,] <- sapply(1:nrep, function(x) max(abs(cor[,,x] - prop_cor[[x]])))
}
avg_max_abs_diff <- apply(max_abs_diff, 1, mean)
plot(tp, avg_max_abs_diff, type = "l", log = "x", 
     xlab = "dirichlet tuning parameter", ylab = "average maximum absolute difference in correlation",
     main = "comparing original and proposed correlation matrices")

#now let's try to run the mcmc sampler again
shape_1 <- 15
shape_2 <- 10
par(mfrow = c(2,2))
hist(rbeta(1E6, shape1 = shape_1, shape2 = shape_2) * 2 - 1, breaks = 1E2, main = "target distribution", xlim = c(-1,1), xlab = "")

corr_init <- 0.001
n_iter <- 1E4
corrs <- rep(0, n_iter)
corrs[1] <- corr_init
dirichlet_tp <- 20

dirichlet_prop <- function(cor, tp, i_weight = 2, offset_1 = 0, offset_2 = 0){
  dim <- dim(cor)[1]
  cor <- (cor + i_weight) / (i_weight + 1)
  cc <- chol(cor)
  cc2 <- cc^2
  new_squares <- sapply(1:dim, function(x) rdirichlet(1, tp * (cc2[,x] + ifelse(cc2[,x] == 0, 0, offset_1)) + ifelse(cc2[,x] == 0, 0, offset_2))) 
  prop_cc <- sqrt(new_squares) 
  prop_cc[cc < 0] <- -prop_cc[cc < 0]
  #prop_cc <- prop_cc * matrix(sample(c(-1, 1), dim^2, replace = T), dim, dim) #does randomizing the signs help?
  # inds <- sample(1:dim, 2)
  # prop_cc[inds] <- -prop_cc[inds]
  prop_cor <- crossprod(prop_cc) 
  prop_cor <- prop_cor * (i_weight + 1) - i_weight
  log_prop_ratio <- sapply(1:dim, function(x) log(ddirichlet(x = cc2[1:x,x], alpha = ((new_squares[,x] + ifelse(new_squares[,x] == 0, 0, offset_1)) * tp + ifelse(new_squares[,x] == 0, 0, offset_2))[1:x])) - 
                             log(ddirichlet(x = new_squares[1:x,x], alpha = ((cc2[,x] + ifelse(cc2[,x] == 0, 0, offset_1)) * tp + ifelse(cc2[,x] == 0, 0, offset_2))[1:x])))
  return(list(sample = prop_cor, log_prop_ratio = log_prop_ratio))
}

for(i in 2:n_iter){
  if(i %% (n_iter / 10) == 0){cat(paste0(i / n_iter * 100, "%  "))}
  corr_prop_rat <- dirichlet_prop(matrix(c(1,corrs[i-1], corrs[i-1], 1), 2, 2), tp = dirichlet_tp, i_weight = 10)
  corr_prop <- corr_prop_rat$sample[1,2]
  log_prop_ratio <-  sum(corr_prop_rat$log_prop_ratio)
  log_dens_ratio <- dbeta(x = (corr_prop + 1) / 2, shape1 = shape_1, shape2 = shape_2, log = T) - 
    dbeta(x = (corrs[i-1] + 1) / 2, shape1 = shape_1, shape2 = shape_2, log = T)
  log_accept_prob <- log_prop_ratio + log_dens_ratio
  if(log(runif(1, 0, 1)) < log_accept_prob){
    corrs[i] <- corr_prop
  } else {
    corrs[i] <- corrs[i-1]
  }
}

hist(corrs[(n_iter / 5) : n_iter], breaks = 50, main = "sampled distribution", xlim = c(-1,1))
plot(corrs, type = "l", main = "trace-plot"); lines(corrs[1:(n_iter/5)], col = 2, ylim = c(-1,1)); abline(v = n_iter / 5, lty = 2)
quantile(corrs[(n_iter / 5) : n_iter], probs = c(1:9)/10)
qbeta(p = c(1:9)/10, shape1 = shape_1, shape2 = shape_2) * 2 - 1
plot(x = quantile(corrs[(n_iter / 5) : n_iter], probs = c(1:99)/100),
     y = qbeta(p = c(1:99)/100, shape1 = shape_1, shape2 = shape_2) * 2 - 1,
     type = "l", main = "Q-Q Plot", xlab = "samples", ylab = "target"); abline(a = 0, b = 1, lty = 2, col = 3)
legend(x = "bottomright", legend = c("quantiles", "1-to-1 line"), col = c(1,3), lty = c(1,2))

## ok so that works ok, what about the blooming onion? i.e. blooming_onion_2
shape_1 <- 15
shape_2 <- 10
par(mfrow = c(2,2))
hist(rbeta(1E6, shape1 = shape_1, shape2 = shape_2) * 2 - 1, breaks = 1E2, main = "target distribution", xlim = c(-1,1), xlab = "")

corr_init <- 0.001
n_iter <- 4E4
corrs <- rep(0, n_iter)
corrs[1] <- corr_init
dirichlet_tp <- 20

for(i in 2:n_iter){
  if(i %% (n_iter / 100) == 0){cat(paste0(i / n_iter * 100, "%  "))}
  corr_prop <- blooming_onion_2(cor = 1, d = 2, eta = 1)[1,2]
  log_prop_ratio <-  sum(corr_prop_rat$log_prop_ratio)
  log_dens_ratio <- dbeta(x = (corr_prop + 1) / 2, shape1 = shape_1, shape2 = shape_2, log = T) - 
    dbeta(x = (corrs[i-1] + 1) / 2, shape1 = shape_1, shape2 = shape_2, log = T)
  log_accept_prob <- log_prop_ratio + log_dens_ratio
  if(log(runif(1, 0, 1)) < log_accept_prob){
    corrs[i] <- corr_prop
  } else {
    corrs[i] <- corrs[i-1]
  }
}

hist(corrs[(n_iter / 5) : n_iter], breaks = 1E2, main = "sampled distribution", xlim = c(-1,1))
plot(corrs, type = "l", main = "trace-plot"); lines(corrs[1:(n_iter/5)], col = 2, ylim = c(-1,1)); abline(v = n_iter / 5, lty = 2)
quantile(corrs[(n_iter / 5) : n_iter], probs = c(1:9)/10)
qbeta(p = c(1:9)/10, shape1 = shape_1, shape2 = shape_2) * 2 - 1
plot(x = quantile(corrs[(n_iter / 5) : n_iter], probs = c(1:99)/100),
     y = qbeta(p = c(1:99)/100, shape1 = shape_1, shape2 = shape_2) * 2 - 1,
     type = "l", main = "Q-Q Plot", xlab = "samples", ylab = "target"); abline(a = 0, b = 1, lty = 2, col = 3)

##########################################
# does it work in more than one dimension?
##########################################

eta <- 25
dim <- 3
par(mfrow = c(2,2))
hist(as.vector(replicate(4E4, rlkj(3, eta)[upper.tri(diag(dim))])), breaks = 1E2, main = "target distribution", xlim = c(-1,1), xlab = "")

corr_init <- diag(dim)
n_iter <- 1E4
corr_mats <- replicate(n_iter, corr_init)

for(i in 2:n_iter){
  if(i %% (n_iter / 10) == 0){cat(paste0(i / n_iter * 100, "%  "))}
  trait_to_prop <- sample(1:dim, 1)
  #permute mat
  corr_prop <- blooming_onion_ind(cor = corr_mats[,,i-1], indices = trait_to_prop, eta = 1)
  log_prop_ratio <-  0
  log_dens_ratio <- dlkjcorr(x = corr_prop, eta = eta, log = T) - 
    dlkjcorr(x = corr_mats[,,i-1], eta = eta, log = T)
  log_accept_prob <- log_prop_ratio + log_dens_ratio
  if(log(runif(1, 0, 1)) < log_accept_prob){
    corr_mats[,,i] <- corr_prop
  } else {
    corr_mats[,,i] <- corr_mats[,,i-1]
  }
}

corrs <- as.vector(sapply(1:n_iter, function(x) corr_mats[,,x][upper.tri(diag(dim))]))
hist(corrs[(n_iter / 5) : n_iter], breaks = 1E2, main = "sampled distribution", xlim = c(-1,1))
plot(corrs, type = "l", main = "trace-plot"); lines(corrs[1:(n_iter/5)], col = 2, ylim = c(-1,1)); abline(v = n_iter / 5, lty = 2)
quantile(corrs[(n_iter / 5) : n_iter], probs = c(1:9)/10)
qbeta(p = c(1:9)/10, shape1 = shape_1, shape2 = shape_2) * 2 - 1
plot(x = quantile(corrs[(n_iter / 5) : n_iter], probs = c(1:99)/100),
     y = qbeta(p = c(1:99)/100, shape1 = shape_1, shape2 = shape_2) * 2 - 1,
     type = "l", main = "Q-Q Plot", xlab = "samples", ylab = "target"); abline(a = 0, b = 1, lty = 2, col = 3)

# what if we don't use the lkj?

blooming_onion_ind <- function (cor = 1, indices = NA, eta = 1) {
  d <- dim(cor)[1]
  if(all(is.na(indices))){indices <- d}
  if(identical(cor, 1)){d_s <- 1}else{d_s <- d - length(indices)}
  R <- matrix(0, d, d)
  R_seed <- cor[-indices, -indices]
  alpha <- eta + (d - 2) / 2
  if(d_s == 1){
    r12 <- 2 * rbeta(1, alpha, alpha) - 1
    R[1, 1] <- 1
    R[1, 2] <- r12
    R[2, 2] <- sqrt(1 - r12^2)
  } else {
    R[1:d_s, 1:d_s] <- chol(R_seed)
    alpha <- alpha - (d_s - 2) / 2
  }
  if (d > 2) {
    for (m in (d_s):(d - 1)) {
      alpha <- alpha - 0.5
      y <- rbeta(1, m/2, alpha)
      z <- rnorm(m, 0, 1)
      z <- z/sqrt(crossprod(z)[1])
      R[1:m, m + 1] <- sqrt(y) * z
      R[m + 1, m + 1] <- sqrt(1 - y)
    }
  }
  R_new <- crossprod(R)
  to_insert <- (d - length(indices) + 1) : d
  new_indices <- order(c((1:d)[-indices], indices))
  R_new <- R_new[new_indices, new_indices]
  return(R_new)
}

dim <- 10
target_corr <- rlkj(dim)
true_corrs <- target_corr[upper.tri(target_corr)]
n_obs <- 100
obs <- rmvnorm(n = n_obs, mean = rep(0, dim), sigma = target_corr)
emp_corrs <- cor(obs)[upper.tri(diag(dim))]
target_corr
par(mfrow = c(3, 2))

corr_init <- diag(dim)
n_iter <- 1E5
corr_mats <- replicate(n_iter, corr_init)

for(i in 2:n_iter){
  if(i %% (n_iter / 10) == 0){cat(paste0(i / n_iter * 100, "%  "))}
  trait_to_prop <- sample(1:dim, 1)
  #permute mat
  corr_prop <- blooming_onion_ind(cor = corr_mats[,,i-1], indices = trait_to_prop, eta = 1)
  log_prop_ratio <-  0
  log_dens_ratio <- sum(dmvnorm(x = obs, mean = rep(0, dim), sigma = corr_prop, log = T)) - 
    sum(dmvnorm(x = obs, mean = rep(0, dim), sigma = corr_mats[,,i-1], log = T))
  log_dens_ratio <- 0 #check if sampling from uniform
  log_accept_prob <- log_prop_ratio + log_dens_ratio
  if(log(runif(1, 0, 1)) < log_accept_prob){
    corr_mats[,,i] <- corr_prop
  } else {
    corr_mats[,,i] <- corr_mats[,,i-1]
  }
}

corrs <- sapply(1:n_iter, function(x) corr_mats[,,x][upper.tri(diag(dim))])
ind_mat <- matrix(1:dim^2, dim, dim)
target_ind_mat <- ind_mat[upper.tri(ind_mat)]
ind_mat <- t(sapply(1:choose(dim, 2), function(x) which(ind_mat == target_ind_mat[x], arr.ind = T)))
for(i in 1:choose(dim, 2)){
  hist(corrs[i,][(n_iter / 5) : n_iter], breaks = 1E2, main = paste0("correlation betw. trait ", ind_mat[i,1], " and trait ", ind_mat[i,2]));
  if(i == 1){legend(x = "topright", lwd = 2, legend = c("true value", "observed value"), col = c("red", "purple"))}
  abline(v = true_corrs[i], col = 2, lwd = 2); abline(v = emp_corrs[i], col = "purple", lwd = 2)
  plot(corrs[i,], type = "l"); lines(corrs[i,][1:(n_iter/5)], col = "grey"); abline(h = true_corrs[i], col = 2, lwd = 2); abline(h = emp_corrs[i], col = "purple", lwd = 2)
}
for(i in 1:choose(dim, 2)){
  hist(corrs[i,][(n_iter / 5) : n_iter], breaks = 1E2, main = paste0("correlation betw. trait ", ind_mat[i,1], " and trait ", ind_mat[i,2]), xlim = c(-1,1));
  if(i == 1){legend(x = "topright", lwd = 2, legend = c("true value", "observed value"), col = c("red", "purple"))}
  abline(v = true_corrs[i], col = 2, lwd = 2); abline(v = emp_corrs[i], col = "purple", lwd = 2)
  plot(corrs[i,], type = "l", ylim = c(-1,1)); lines(corrs[i,][1:(n_iter/5)], col = "grey"); abline(h = true_corrs[i], col = 2, lwd = 2); abline(h = emp_corrs[i], col = "purple", lwd = 2)
}

#looks to, but proposals are too aggressive
#let's try using the dirichlet proposal thing?

normal_chol_prop <- function(cor, chol = F, varN, inds = NA){ ## does not work
  if(!chol){
    cc <- chol(cor)
  } else {
    cc <- cor
  }
  
  if(is.na(inds)){
    inds <- dim(cor)[1]
  }
  
  for(i in length(inds)){
    cc[1:inds[i],inds[i]] <- cc[1:inds[i],inds[i]] + rnorm(inds[i], 0, sqrt(varN))
    cc[1:inds[i],inds[i]] <- cc[1:inds[i],inds[i]] / crossprod(cc[1:inds[i],inds[i]])
  }
  
  if(chol){
    return(list(sample = cc, log_prop_ratio = 0))
  } else {
    return(list(sample = crossprod(cc), log_prop_ratio = 0))
  }
}

dirichlet_prop <- function(cor, tp, i_weight = 2, offset_1 = 0, offset_2 = 0){ #also does not work
  cor_orig <- cor
  dim <- dim(cor)[1]
  plus_mat <- matrix(1.0,dim,dim); diag(plus_mat) <- 1
  pos_def_mat_obtained <- F
  while(!pos_def_mat_obtained){
    cor <- cor_orig
    permute_cols <- sample(1:dim)
    cor <- cor[permute_cols, permute_cols]
    cor <- (cor + plus_mat * i_weight) / (i_weight + 1)
    cc <- chol(cor)
    cc2 <- cc^2
    new_squares <- sapply(1:dim, function(x) rdirichlet(1, tp * (cc2[,x] + ifelse(cc2[,x] == 0, 0, offset_1)) + ifelse(cc2[,x] == 0, 0, offset_2))) 
    prop_cc <- sqrt(new_squares) 
    prop_cc[cc < 0] <- -prop_cc[cc < 0] #need to implement reflecting boundary for this?
    #prop_cc <- prop_cc * matrix(sample(c(-1, 1), dim^2, replace = T), dim, dim) #does randomizing the signs help?
    # inds <- sample(1:dim, 2)
    # prop_cc[inds] <- -prop_cc[inds]
    prop_cor <- crossprod(prop_cc) 
    prop_cor <- prop_cor * (i_weight + 1) - plus_mat * i_weight
    pos_def_mat_obtained <- is.positive.definite(prop_cor)
  }
  prop_cor <- prop_cor[order(permute_cols), order(permute_cols)]
  log_prop_ratio <- sapply(1:dim, function(x) log(ddirichlet(x = cc2[1:x,x], alpha = ((new_squares[,x] + ifelse(new_squares[,x] == 0, 0, offset_1)) * tp + ifelse(new_squares[,x] == 0, 0, offset_2))[1:x])) - 
                             log(ddirichlet(x = new_squares[1:x,x], alpha = ((cc2[,x] + ifelse(cc2[,x] == 0, 0, offset_1)) * tp + ifelse(cc2[,x] == 0, 0, offset_2))[1:x])))
  return(list(sample = prop_cor, log_prop_ratio = log_prop_ratio))
}

dim <- 4
target_corr <- rlkj(dim)
dirichlet_tp <- 20
true_corrs <- target_corr[upper.tri(target_corr)]
n_obs <- 10
obs <- rmvnorm(n = n_obs, mean = rep(0, dim), sigma = target_corr)
emp_corrs <- cor(obs)[upper.tri(diag(dim))]
target_corr
par(mfrow = c(ifelse(dim < 5, choose(dim, 2), 4), 2))

corr_init <- diag(dim)
n_iter <- 1E4
corr_mats <- replicate(n_iter, corr_init)
log_prop_rats <- rep(0, n_iter-1)

for(i in 2:n_iter){
  if(i %% (n_iter / 10) == 0){cat(paste0(i / n_iter * 100, "%  "))}
  corr_prop_rat <- dirichlet_prop(cor = corr_mats[,,i-1], tp = dirichlet_tp, offset_1 = 0, offset_2 = 0, i_weight = 1)
  corr_prop <- corr_prop_rat$sample
  log_prop_ratio <- sum(corr_prop_rat$log_prop_ratio)
  log_prop_rats[i] <- log_prop_ratio
  # log_dens_ratio <- sum(dmvnorm(x = obs, mean = rep(0, dim), sigma = corr_prop, log = T)) - 
  #   sum(dmvnorm(x = obs, mean = rep(0, dim), sigma = corr_mats[,,i-1], log = T))
  log_dens_ratio <- 0 #proposal ratio not being computed correctly? compare to corr_mats <- aperm((rlkjcorr(1E5, K = dim, eta = 1)), c(2,3,1))
  #### or HMMM, maybe invalid matrices are being returned? YESS that is it! Something about the proposal distribution is faulty!
  if(!is.positive.definite(corr_prop))
  {
    cat(".")
    log_dens_ratio <- -Inf
  }
  log_accept_prob <- log_prop_ratio + log_dens_ratio
  if(log(runif(1, 0, 1)) < log_accept_prob){
    corr_mats[,,i] <- corr_prop
  } else {
    corr_mats[,,i] <- corr_mats[,,i-1]
  }
}

corrs <- sapply(1:n_iter, function(x) corr_mats[,,x][upper.tri(diag(dim))])
ind_mat <- matrix(1:dim^2, dim, dim)
target_ind_mat <- ind_mat[upper.tri(ind_mat)]
ind_mat <- t(sapply(1:choose(dim, 2), function(x) which(ind_mat == target_ind_mat[x], arr.ind = T)))
for(i in 1:choose(dim, 2)){
  hist(corrs[i,][(n_iter / 5) : n_iter], breaks = 1E2, main = paste0("correlation betw. trait ", ind_mat[i,1], " and trait ", ind_mat[i,2]));
  if(i == 1){legend(x = "topright", lwd = 2, legend = c("true value", "observed value"), col = c("red", "purple"))}
  abline(v = true_corrs[i], col = 2, lwd = 2); abline(v = emp_corrs[i], col = "purple", lwd = 2)
  plot(corrs[i,], type = "l"); lines(corrs[i,][1:(n_iter/5)], col = "grey"); abline(h = true_corrs[i], col = 2, lwd = 2); abline(h = emp_corrs[i], col = "purple", lwd = 2)
}
for(i in 1:choose(dim, 2)){
  hist(corrs[i,][(n_iter / 5) : n_iter], breaks = 1E2, main = paste0("correlation betw. trait ", ind_mat[i,1], " and trait ", ind_mat[i,2]), xlim = c(-1,1));
  if(i == 1){legend(x = "topright", lwd = 2, legend = c("true value", "observed value"), col = c("red", "purple"))}
  abline(v = true_corrs[i], col = 2, lwd = 2); abline(v = emp_corrs[i], col = "purple", lwd = 2)
  plot(corrs[i,], type = "l", ylim = c(-1,1)); lines(corrs[i,][1:(n_iter/5)], col = "grey"); abline(h = true_corrs[i], col = 2, lwd = 2); abline(h = emp_corrs[i], col = "purple", lwd = 2)
}

#bah so this doesn't seem to work, proposal ratio misspecified maybe?
#let's try to go back to the tuned lkj proposal distribution


blooming_onion_ind_tune <- function (cor = 1, indices = NA, eta = 1) {
  d <- dim(cor)[1]
  chol_cor <- chol(cor)
  if(all(is.na(indices))){indices <- d}
  if(identical(cor, 1)){return(1)}else{d_s <- d - length(indices)}
  R <- matrix(0, d, d)
  R_seed <- cor[-indices, -indices]
  alpha <- eta + (d - 2) / 2
  if(d_s == 1){
    r12 <- 2 * rbeta(1, alpha, alpha) - 1
    R[1, 1] <- 1
    R[1, 2] <- r12
    R[2, 2] <- sqrt(1 - r12^2)
  } else {
    R[1:d_s, 1:d_s] <- chol(R_seed)
    alpha <- alpha - (d_s - 2) / 2
  }
  if (d > 2) {
    for (m in (d_s):(d - 1)) {
      alpha <- alpha - 0.5
      y <- rbeta(1, m/2, alpha)
      z <- rnorm(m, 0, 1)
      z <- z/sqrt(crossprod(z)[1])
      R[1:m, m + 1] <- sqrt(y) * z
      R[m + 1, m + 1] <- sqrt(1 - y)
    }
  }
  R_new <- crossprod(R)
  to_insert <- (d - length(indices) + 1) : d
  new_indices <- order(c((1:d)[-indices], indices))
  R_new <- R_new[new_indices, new_indices]
  return(R_new)
}


blooming_onion_tune <- function (cor, varB = 0.01, varN = 0.1, redrawBeta = T, betaWindow = NA) {
  d <- dim(cor)[1]
  R <- matrix(0, d, d)
  permute_cols <- sample(1:d)
  cor <- cor[permute_cols, permute_cols]
  m <- d - 1
  chol_cor <- chol(cor)
  R[1:m, 1:m] <- chol_cor[1:m,1:m]
  target_chol <- chol_cor[1:d,d]
  target_y <- 1 - (target_chol[d]^2)
  target_z <- target_chol[1:m] / sqrt(target_y)
  if(redrawBeta){
    alpha <- 1 #when eta = 1
    if(is.na(betaWindow)){
      y <- rbeta(1, m/2, alpha)
    } else {
      betaBounds <- c(target_y + betaWindow, target_y - betaWindow)
      betaBounds[betaBounds > 1] <- 1
      betaBounds[betaBounds < -1] <- -1
      betaBounds_prob <- pbeta(q = betaBounds, shape1 = m/2, shape2 = alpha)
      unif_draw <- runif(1)
      beta_subprob <- betaBounds_prob[1] - betaBounds_prob[2]
      betaProb <- betaBounds_prob[2] + beta_subprob * unif_draw
      y <- qbeta(p = betaProb, shape1 = (d-1)/2, shape2 = alpha)
    }
  } else {
    y <- target_y 
  }
  z <- rnorm(m, target_z, sqrt(varN))
  z <- z/sqrt(crossprod(z)[1])
  R[1:m, m + 1] <- sqrt(y) * z
  R[m + 1, m + 1] <- sqrt(1 - y)
  #log orig in new center - log new in orig prop dist
  if(redrawBeta){
    if(is.na(betaWindow)){
      
      # log_prop_ratio <- dbeta(x = target_y, (d-1)/2, alpha) - dbeta(x = y, (d-1)/2, alpha) +
      #   sum(dnorm(x = target_z, mean = z, sd = sqrt(varN))) - sum(dnorm(x = z, mean =target_z, sd = sqrt(varN)))
      
      log_prop_ratio <- 0
      
    } else {
      
      betaBounds_rev <- c(y + betaWindow, y - betaWindow)
      betaBounds_rev[betaBounds_rev > 1] <- 1
      betaBounds_rev[betaBounds_rev < -1] <- -1
      
      #############################################################################
      # is it the truncated beta densities that are needed for the proposal ratio?
      #############################################################################
    
      # betaBounds_prob_rev <- pbeta(q = betaBounds_rev, shape1 = (d-1)/2, shape2 = alpha)
      # beta_subprob_rev <- betaBounds_prob_rev[1] - betaBounds_prob_rev[2]
      #   
      # log_prop_ratio <- dbeta(x = target_y, (d-1)/2, alpha) / beta_subprob_rev - dbeta(x = y, (d-1)/2, alpha) / beta_subprob +
      #   sum(dnorm(x = target_z, mean = z, sd = sqrt(varN))) - sum(dnorm(x = z, mean =target_z, sd = sqrt(varN)))
      
      ####################################
      # what about the uniform densities?
      ####################################
      
      # log_prop_ratio <- log(abs(1 / diff(betaBounds_rev))) - log(abs(1 / diff(betaBounds))) +
      #   sum(dnorm(x = target_z, mean = z, sd = sqrt(varN))) - sum(dnorm(x = z, mean =target_z, sd = sqrt(varN)))

      ###########################################
      # ... what about the beta cumulative probs?
      ###########################################
      
      betaBounds_prob_rev <- pbeta(q = betaBounds_rev, shape1 = (d-1)/2, shape2 = alpha)
      beta_subprob_rev <- betaBounds_prob_rev[1] - betaBounds_prob_rev[2]

      log_prop_ratio <- log(abs(1 / beta_subprob_rev)) - log(abs(1 / beta_subprob)) +
        sum(dnorm(x = target_z, mean = z, sd = sqrt(varN))) - sum(dnorm(x = z, mean =target_z, sd = sqrt(varN)))
      
    }
  } else {
    log_prop_ratio <- sum(dnorm(x = target_z, mean = z, sd = sqrt(varN))) - sum(dnorm(x = z, mean =target_z, sd = sqrt(varN)))
  }
  R <- crossprod(R)
  R <- R[order(permute_cols), order(permute_cols)]
  return(list(sample = R, log_prop_ratio = log_prop_ratio))
}


blooming_onion_tune_ind <- function (cor, varN = 0.1, betaWindow = NA, chol = F, inds = NA) { #also does not work
  d <- dim(cor)[1]
  R <- matrix(0, d, d)
  if(any(is.na(inds))){
    inds <- d
  }
  chol_cor <- chol(cor)
  R <- chol_cor
  log_prop_ratio <- 0
  for(i in inds){
    m <- i - 1
    target_chol <- chol_cor[1:i,i]
    target_y <- 1 - (target_chol[i]^2)
    target_z <- target_chol[1:m] / sqrt(target_y)
    if(!is.na(betaWindow)){
      alpha <- 1 + (d-i)*0.5
      if(betaWindow >= 1){
        y <- rbeta(1, m/2, alpha)
      } else {
        betaBounds <- c(target_y + betaWindow, target_y - betaWindow)
        betaBounds[betaBounds > 1] <- 1
        betaBounds[betaBounds < -1] <- -1
        betaBounds_prob <- pbeta(q = betaBounds, shape1 = m/2, shape2 = alpha)
        unif_draw <- runif(1)
        beta_subprob <- betaBounds_prob[1] - betaBounds_prob[2]
        betaProb <- betaBounds_prob[2] + beta_subprob * unif_draw
        y <- qbeta(p = betaProb, shape1 = (d-1)/2, shape2 = alpha)
      }
    } else {
      y <- target_y 
    }
    
    z <- rnorm(m, target_z, sqrt(varN))
    z <- z/sqrt(crossprod(z)[1])
    R[1:m, m + 1] <- sqrt(y) * z
    R[m + 1, m + 1] <- sqrt(1 - y)
  
      if(!is.na(betaWindow)){
        if(betaWindow >= 1){
          
          log_prop_ratio <- log_prop_ratio + 0
          
        } else {
          
          betaBounds_rev <- c(y + betaWindow, y - betaWindow)
          betaBounds_rev[betaBounds_rev > 1] <- 1
          betaBounds_rev[betaBounds_rev < -1] <- -1
          
          betaBounds_prob_rev <- pbeta(q = betaBounds_rev, shape1 = (d-1)/2, shape2 = alpha)
          beta_subprob_rev <- betaBounds_prob_rev[1] - betaBounds_prob_rev[2]
          
          log_prop_ratio <- log_prop_ratio + log(abs(1 / beta_subprob_rev)) - log(abs(1 / beta_subprob)) 
        }
        
    } else {
      log_prop_ratio <- log_prop_ratio + 0
    }
  }
  R <- crossprod(R)
  
  return(list(sample = R, log_prop_ratio = log_prop_ratio))
}

blooming_onion_tune_clean <- function (cor, varN = 0.1, betaWindow = NA) {
  d <- dim(cor)[1]
  R <- matrix(0, d, d)
  permute_cols <- sample(1:d)
  cor <- cor[permute_cols, permute_cols]
  m <- d - 1
  chol_cor <- chol(cor)
  R[1:m, 1:m] <- chol_cor[1:m,1:m]
  target_chol <- chol_cor[1:d,d]
  target_y <- 1 - (target_chol[d]^2)
  target_z <- target_chol[1:m] / sqrt(target_y)
  if(!is.na(betaWindow)){
    alpha <- 1 
    if(betaWindow >= 1){
      y <- rbeta(1, m/2, alpha)
    } else {
      betaBounds <- c(target_y + betaWindow, target_y - betaWindow)
      betaBounds[betaBounds > 1] <- 1
      betaBounds[betaBounds < -1] <- -1
      betaBounds_prob <- pbeta(q = betaBounds, shape1 = m/2, shape2 = alpha)
      unif_draw <- runif(1)
      beta_subprob <- betaBounds_prob[1] - betaBounds_prob[2]
      betaProb <- betaBounds_prob[2] + beta_subprob * unif_draw
      y <- qbeta(p = betaProb, shape1 = (d-1)/2, shape2 = alpha)
    }
  } else {
    y <- target_y 
  }
  z <- rnorm(m, target_z, sqrt(varN))
  z <- z/sqrt(crossprod(z)[1])
  R[1:m, m + 1] <- sqrt(y) * z
  R[m + 1, m + 1] <- sqrt(1 - y)

  if(redrawBeta){
    if(is.na(betaWindow)){
        
      log_prop_ratio <- 0
      
    } else {
      
      betaBounds_rev <- c(y + betaWindow, y - betaWindow)
      betaBounds_rev[betaBounds_rev > 1] <- 1
      betaBounds_rev[betaBounds_rev < -1] <- -1
 
      betaBounds_prob_rev <- pbeta(q = betaBounds_rev, shape1 = (d-1)/2, shape2 = alpha)
      beta_subprob_rev <- betaBounds_prob_rev[1] - betaBounds_prob_rev[2]
      
      log_prop_ratio <- log(abs(1 / beta_subprob_rev)) - log(abs(1 / beta_subprob)) 
      
    }
  } else {
    log_prop_ratio <- 0
  }
  
  R <- crossprod(R)
  R <- R[order(permute_cols), order(permute_cols)]
  
  return(list(sample = R, log_prop_ratio = log_prop_ratio))
}


dim <- 5
target_corr <- rlkj(dim)
true_corrs <- target_corr[upper.tri(target_corr)]
n_obs <- 40
varN <- 0.5
betaWindow <- 0.2
sampleUniform <- T
obs <- rmvnorm(n = n_obs, mean = rep(0, dim), sigma = target_corr)
emp_corrs <- cor(obs)[upper.tri(diag(dim))]
par(mfrow = c(ifelse(dim < 5, choose(dim, 2), 4), 2))

corr_init <- rlkj(dim) 
n_iter <- 1E6
thin <- 1E2
n_out <- round(n_iter/thin)

corr_mats <- replicate(n_out, diag(dim))
corr_mats[,,1] <- corr_curr <-  corr_init

n_accept <- 0
for(i in 2:n_iter){
  
  #progress bar
  if(i %% (n_iter / 10) == 0){cat(paste0(i / n_iter * 100, "%  "))}
  
  #propose a move
  corr_prop_rat <- blooming_onion_tune_clean(corr_curr, varN = varN, betaWindow = betaWindow)
  corr_prop <- corr_prop_rat$sample
  log_prop_ratio <-  corr_prop_rat$log_prop_ratio
  
  #compute ratio of target densities
  if(sampleUniform){
    log_dens_ratio <- 0 
  } else {
    log_dens_ratio <- sum(dmvnorm(x = obs, mean = rep(0, dim), sigma = corr_prop, log = T)) - 
      sum(dmvnorm(x = obs, mean = rep(0, dim), sigma = corr_mats[,,i-1], log = T))
  }
  
  #accept or reject
  log_accept_prob <- log_prop_ratio + log_dens_ratio
  if(log(runif(1, 0, 1)) < log_accept_prob){
    corr_curr <- corr_prop
    n_accept <- n_accept + 1
  } 
  
  #record given thinning 
  if(i %% thin == 0){
    corr_mats[,,i/thin] <- corr_curr
  }

}

corrs <- sapply(1:n_out, function(x) corr_mats[,,x][upper.tri(diag(dim))])
ind_mat <- matrix(1:dim^2, dim, dim)
target_ind_mat <- ind_mat[upper.tri(ind_mat)]
ind_mat <- t(sapply(1:choose(dim, 2), function(x) which(ind_mat == target_ind_mat[x], arr.ind = T)))
for(i in 1:choose(dim, 2)){
  hist(corrs[i,][(n_out / 5) : n_out], breaks = 1E2, main = paste0("correlation betw. trait ", ind_mat[i,1], " and trait ", ind_mat[i,2]));
  if(i == 1){legend(x = "topright", lwd = 2, legend = c("true value", "observed value"), col = c("red", "purple"))}
  abline(v = true_corrs[i], col = 2, lwd = 2); abline(v = emp_corrs[i], col = "purple", lwd = 2)
  plot(corrs[i,], type = "l"); lines(corrs[i,][1:(n_out/5)], col = "grey"); abline(h = true_corrs[i], col = 2, lwd = 2); abline(h = emp_corrs[i], col = "purple", lwd = 2)
}
for(i in 1:choose(dim, 2)){
  hist(corrs[i,][(n_out / 5) : n_out], breaks = 1E2, main = paste0("correlation betw. trait ", ind_mat[i,1], " and trait ", ind_mat[i,2]), xlim = c(-1,1));
  if(i == 1){legend(x = "topright", lwd = 2, legend = c("true value", "observed value"), col = c("red", "purple"))}
  abline(v = true_corrs[i], col = 2, lwd = 2); abline(v = emp_corrs[i], col = "purple", lwd = 2)
  plot(corrs[i,], type = "l", ylim = c(-1,1)); lines(corrs[i,][1:(n_out/5)], col = "grey"); abline(h = true_corrs[i], col = 2, lwd = 2); abline(h = emp_corrs[i], col = "purple", lwd = 2)
}

if(sampleUniform){
  eta <- 1

  par(mfrow = c(3, 2))
  dets_samp <- (sapply(1:n_out, function(x) det(corr_mats[,,x])))[(n_out / 5) : n_out]
  hist(dets_samp, main = "Sampled Determinants")
  title(main = paste0("dimension = ", dim), outer = T, line = -2, cex.main = 2)
  
  known_lkj_distribution <- aperm((rlkjcorr(n_out * 0.8, K = dim, eta = eta)), c(2,3,1))
  corrs_targ <- sapply(1:(n_out * 0.8), function(x) known_lkj_distribution[,,x][upper.tri(diag(dim))])
  dets_targ <- (sapply(1:(n_out * 0.8), function(x) det(known_lkj_distribution[,,x])))
  hist(dets_targ, main = "Target Determinants")
  
  plot(x = quantile(dets_samp, probs = c(1:99)/100),
       y = quantile(dets_targ, probs = c(1:99)/100),
       type = "l", main = "Q-Q Plot of Determinants", xlab = "samples", ylab = "target"); abline(a = 0, b = 1, lty = 2, col = 3)
  legend(x = "bottomright", legend = c("quantiles", "1-to-1 line"), col = c(1,3), lty = c(1,2))
  
  plot(x = quantile(as.vector(corrs), probs = c(1:99)/100),
       y = quantile(as.vector(corrs_targ), probs = c(1:99)/100),
       type = "l", main = "Q-Q Plot of Marginal Correlations", xlab = "samples", ylab = "target"); abline(a = 0, b = 1, lty = 2, col = 3)
  legend(x = "bottomright", legend = c("quantiles", "1-to-1 line"), col = c(1,3), lty = c(1,2))
  legend(x = "bottomright", legend = c("quantiles", "1-to-1 line"), col = c(1,3), lty = c(1,2))
  
  hist(as.vector(corrs), main = "Marginal Sampled Correlations")
  hist(as.vector(corrs_targ), main = "Marginal Target Correlations")
}

paste0("acceptance ratio = ", n_accept / n_iter)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### let's see if we can poke other indicesas well??????????????
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#nope, does not appear to!?

### ### ### ### ### ### ### ### ### ### ### ### 
### How about in the more efficient formulation?
### ### ### ### ### ### ### ### ### ### ### ### 

choleskalator <- function(mat, ind){
  dim <- dim(mat)[1]
  if(ind == dim){
    return(mat)
  } else {
    ind_removed <- mgcv::choldrop(mat, ind)
    Rn_target <- c((t(mat) %*% mat[,ind])[-ind],1)  
    URn_target_sle <- backsolve(ind_removed, (Rn_target[1:(dim-1)]), transpose = T) 
    URn_target_sle <- c(URn_target_sle, sqrt(1 - crossprod(URn_target_sle))) 
    return(cbind(rbind(ind_removed, rep(0,dim-1)), URn_target_sle))
  }
}
blooming_onion_tune_chol <- function (upper_cholesky_factor, ind, varN = 0.1, betaWindow = NA) {
  redrawBeta <- ifelse(!is.na(betaWindow), T, F)
  d <- dim(upper_cholesky_factor)[1]
  m <- d-1
  R <- matrix(0, d, d)
  chol_cor <- choleskalator(upper_cholesky_factor, ind)
  R[1:m, 1:m] <- chol_cor[1:m,1:m]
  target_chol <- chol_cor[1:d,d]
  target_y <- 1 - (target_chol[d]^2)
  target_z <- target_chol[1:m] / sqrt(target_y)
  if(!is.na(betaWindow)){
    alpha <- 1 
    if(betaWindow >= 1){
      y <- rbeta(1, m/2, alpha)
    } else {
      betaBounds <- c(target_y + betaWindow, target_y - betaWindow)
      betaBounds[betaBounds > 1] <- 1
      betaBounds[betaBounds < -1] <- -1
      betaBounds_prob <- pbeta(q = betaBounds, shape1 = m/2, shape2 = alpha)
      unif_draw <- runif(1)
      beta_subprob <- betaBounds_prob[1] - betaBounds_prob[2]
      betaProb <- betaBounds_prob[2] + beta_subprob * unif_draw
      y <- qbeta(p = betaProb, shape1 = (d-1)/2, shape2 = alpha)
    }
  } else {
    y <- target_y 
  }
  z <- rnorm(m, target_z, sqrt(varN))
  z <- z/sqrt(crossprod(z)[1])
  R[1:m, m + 1] <- sqrt(y) * z
  R[m + 1, m + 1] <- sqrt(1 - y)
  
  if(redrawBeta){
    if(is.na(betaWindow)){
      
      log_prop_ratio <- 0
      
    } else {
      
      betaBounds_rev <- c(y + betaWindow, y - betaWindow)
      betaBounds_rev[betaBounds_rev > 1] <- 1
      betaBounds_rev[betaBounds_rev < -1] <- -1
      
      betaBounds_prob_rev <- pbeta(q = betaBounds_rev, shape1 = (d-1)/2, shape2 = alpha)
      beta_subprob_rev <- betaBounds_prob_rev[1] - betaBounds_prob_rev[2]
      
      log_prop_ratio <- log(abs(1 / beta_subprob_rev)) - log(abs(1 / beta_subprob)) 
      
    }
  } else {
    log_prop_ratio <- 0
  }
  
  return(list(sample = R, log_prop_ratio = log_prop_ratio))
}

dim <- 10
cormat <- rlkj(dim)
cholfac <- chol(cormat)
microbenchmark(blooming_onion_tune_chol(upper_cholesky_factor = cholfac, ind = sample(1:dim, 1), varN = 0.1, betaWindow = 0.1), 
               blooming_onion_tune_clean(cor = cormat, varN = 0.1, betaWindow = 0.1))


dim <- 10
target_corr <- rlkj(dim)
true_corrs <- target_corr[upper.tri(target_corr)]
n_obs <- 40
varN <- 0.5
betaWindow <- 0.2
sampleUniform <- T
obs <- rmvnorm(n = n_obs, mean = rep(0, dim), sigma = target_corr)
emp_corrs <- cor(obs)[upper.tri(diag(dim))]
par(mfrow = c(ifelse(dim < 5, choose(dim, 2), 4), 2))

corr_init <- chol(rlkj(dim))
n_iter <- 2.5E6
thin <- 1E2
n_out <- round(n_iter/thin)

corr_mats <- replicate(n_out, diag(dim))
corr_mats[,,1] <- corr_curr <-  corr_init

inds_init <- 1:dim
inds <- replicate(n_out, rep(0,dim))
inds[,1] <- inds_curr <-  inds_init


n_accept <- 0
for(i in 2:n_iter){
  
  #progress bar
  if(i %% (n_iter / 10) == 0){cat(paste0(i / n_iter * 100, "%  "))}
  
  #propose a move
  ind_to_prop <- sample(1:dim, 1)
  inds_prop <- c(inds_curr[-ind_to_prop], inds_curr[ind_to_prop])
  corr_prop_rat <- blooming_onion_tune_chol(corr_curr, ind = ind_to_prop, varN = varN, betaWindow = betaWindow)
  corr_prop <- corr_prop_rat$sample
  log_prop_ratio <-  corr_prop_rat$log_prop_ratio
  
  #compute ratio of target densities
  if(sampleUniform){
    log_dens_ratio <- 0 
  } else {
    log_dens_ratio <- sum(dmvnorm(x = obs, mean = rep(0, dim), sigma = corr_prop, log = T)) - 
      sum(dmvnorm(x = obs, mean = rep(0, dim), sigma = corr_mats[,,i-1], log = T))
  }
  
  #accept or reject
  log_accept_prob <- log_prop_ratio + log_dens_ratio
  if(log(runif(1, 0, 1)) < log_accept_prob){
    corr_curr <- corr_prop
    inds_curr <- inds_prop
    n_accept <- n_accept + 1
  } 
  
  #record given thinning 
  if(i %% thin == 0){
    corr_mats[,,i/thin] <- corr_curr
    inds[,i/thin] <- inds_curr
  }
  
}

corr_mats <- sapply(1:n_out, function(x) (t(corr_mats[,,x]) %*% corr_mats[,,x])[inds[,x],inds[,x]], simplify = "array")

corrs <- sapply(1:n_out, function(x) corr_mats[,,x][upper.tri(diag(dim))])
ind_mat <- matrix(1:dim^2, dim, dim)
target_ind_mat <- ind_mat[upper.tri(ind_mat)]
ind_mat <- t(sapply(1:choose(dim, 2), function(x) which(ind_mat == target_ind_mat[x], arr.ind = T)))
for(i in 1:choose(dim, 2)){
  hist(corrs[i,][(n_out / 5) : n_out], breaks = 1E2, main = paste0("correlation betw. trait ", ind_mat[i,1], " and trait ", ind_mat[i,2]));
  if(i == 1){legend(x = "topright", lwd = 2, legend = c("true value", "observed value"), col = c("red", "purple"))}
  abline(v = true_corrs[i], col = 2, lwd = 2); abline(v = emp_corrs[i], col = "purple", lwd = 2)
  plot(corrs[i,], type = "l"); lines(corrs[i,][1:(n_out/5)], col = "grey"); abline(h = true_corrs[i], col = 2, lwd = 2); abline(h = emp_corrs[i], col = "purple", lwd = 2)
}
for(i in 1:choose(dim, 2)){
  hist(corrs[i,][(n_out / 5) : n_out], breaks = 1E2, main = paste0("correlation betw. trait ", ind_mat[i,1], " and trait ", ind_mat[i,2]), xlim = c(-1,1));
  if(i == 1){legend(x = "topright", lwd = 2, legend = c("true value", "observed value"), col = c("red", "purple"))}
  abline(v = true_corrs[i], col = 2, lwd = 2); abline(v = emp_corrs[i], col = "purple", lwd = 2)
  plot(corrs[i,], type = "l", ylim = c(-1,1)); lines(corrs[i,][1:(n_out/5)], col = "grey"); abline(h = true_corrs[i], col = 2, lwd = 2); abline(h = emp_corrs[i], col = "purple", lwd = 2)
}

if(sampleUniform){
  eta <- 1
  
  par(mfrow = c(3, 2))
  dets_samp <- (sapply(1:n_out, function(x) det(corr_mats[,,x])))[(n_out / 5) : n_out]
  hist(dets_samp, main = "Sampled Determinants")
  title(main = paste0("dimension = ", dim), outer = T, line = -2, cex.main = 2)
  
  known_lkj_distribution <- aperm((rlkjcorr(n_out * 0.8, K = dim, eta = eta)), c(2,3,1))
  corrs_targ <- sapply(1:(n_out * 0.8), function(x) known_lkj_distribution[,,x][upper.tri(diag(dim))])
  dets_targ <- (sapply(1:(n_out * 0.8), function(x) det(known_lkj_distribution[,,x])))
  hist(dets_targ, main = "Target Determinants")
  
  plot(x = quantile(dets_samp, probs = c(1:99)/100),
       y = quantile(dets_targ, probs = c(1:99)/100),
       type = "l", main = "Q-Q Plot of Determinants", xlab = "samples", ylab = "target"); abline(a = 0, b = 1, lty = 2, col = 3)
  legend(x = "bottomright", legend = c("quantiles", "1-to-1 line"), col = c(1,3), lty = c(1,2))
  
  plot(x = quantile(as.vector(corrs), probs = c(1:99)/100),
       y = quantile(as.vector(corrs_targ), probs = c(1:99)/100),
       type = "l", main = "Q-Q Plot of Marginal Correlations", xlab = "samples", ylab = "target"); abline(a = 0, b = 1, lty = 2, col = 3)
  legend(x = "bottomright", legend = c("quantiles", "1-to-1 line"), col = c(1,3), lty = c(1,2))
  legend(x = "bottomright", legend = c("quantiles", "1-to-1 line"), col = c(1,3), lty = c(1,2))
  
  hist(as.vector(corrs), main = "Marginal Sampled Correlations")
  hist(as.vector(corrs_targ), main = "Marginal Target Correlations")
}

paste0("acceptance ratio = ", n_accept / n_iter)
