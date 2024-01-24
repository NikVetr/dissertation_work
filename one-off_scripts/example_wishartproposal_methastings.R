#check met hastings operation on correlation matrix with wishart move
library(MCMCpack)
shape_1 <- 72
shape_2 <- 15
par(mfrow = c(2,2))
hist(rbeta(1E6, shape1 = shape_1, shape2 = shape_2) * 2 - 1, breaks = 1E2, main = "target distribution", xlim = c(-1,1), xlab = "")

corr_init <- 0.001
n_iter <- 5E4
corrs <- rep(0, n_iter)
corrs[1] <- corr_init
df_tune <- 20
useCorrMatDens <- T
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
