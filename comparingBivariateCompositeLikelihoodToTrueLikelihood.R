library(mvtnorm)
library(pbivnorm)

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

compareCompVsFull <- function(dimr, cor = NA, meanl = -1, meanu = 1, lower = NA, upper = NA){
  if(all(is.na(cor))){
    x <- rlkj(dimr)
  } else {
    x <- cor
  }
  if(all(is.na(lower))){
    lower <- rnorm(dimr, mean = meanl)
  }
  if(all(is.na(upper))){
    upper <- rnorm(dimr, mean = meanu)
  }
  swapi <- lower > upper
  swap <- upper[swapi]
  upper[swapi] <- lower[swapi]
  lower[swapi] <- swap
  trueprob = log(pmvnorm(lower = lower, upper = upper, corr = x)[1])
  
  uppers <- matrix(0, choose(dimr, 2), 2)
  lowers <- matrix(0, choose(dimr, 2), 2)
  rhos <- rep(0, choose(dimr, 2))
  counter <- 0
  for(i in 1:(dimr-1)){
    for(j in (i+1):dimr){
      counter <- counter + 1
      uppers[counter,] <- c(upper[i], upper[j])
      lowers[counter,] <- c(lower[i], lower[j])
      rhos[counter] <- x[i,j]
    }
  }
  
  probs <- pbivnorm_infs(uppers[,1], uppers[,2], rhos) - 
    pbivnorm_infs(lowers[,1], uppers[,2], rhos) - 
    pbivnorm_infs(uppers[,1], lowers[,2], rhos) + 
    pbivnorm_infs(lowers[,1], lowers[,2], rhos)
  
  probs[probs < 0] <- 0 #get rid of floating point errors
  prob <- mean(log(probs))# * dimr * 0.6 + dimr / 8
  return(c(trueprob, prob))
}

dimr <- 50
nrep <- 200
r <- diag(dimr) * 0.5 + 0.5 
comparingTrueVsComp <- t(replicate(nrep, compareCompVsFull(dimr, cor = NA)))
# comparingTrueVsComp[,2] <- comparingTrueVsComp[,2] / dimr
lm(comparingTrueVsComp[,1] ~ comparingTrueVsComp[,2])
# lm(exp(comparingTrueVsComp[,1]) ~ exp(comparingTrueVsComp[,2]))


dev.off()
par(mar = c(5,5,3,3))
plot(comparingTrueVsComp, xlab = "\"full\" multivariate normal integral", ylab = "composite integral obtained from \nall pairwise bivariate integrals")
title(paste0("dimension = ", dimr, ", r = ", round(cor(comparingTrueVsComp)[1,2], 3)))
abline(0, 1, lwd = 2, col = 2)
legend(x = "topleft", legend = "1-to-1 line", col = 2, lwd = 2, lty = 1)

plot(comparingTrueVsComp, xlab = "\"full\" multivariate normal integral", ylab = "composite integral obtained from \nall pairwise bivariate integrals")
# compareCompVsFull(dimr, lower = rep(-3, dimr), upper = rep(3, dimr))

dimrange <- 1:10
intercepts <- rep(0, length(dimrange))
slopes <- rep(0, length(dimrange))
for(dimr in dimrange){
  print(dimr)
  nrep <- 10000
  comparingTrueVsComp <- t(replicate(nrep, compareCompVsFull(dimr*10, cor = NA)))
  lmx <- lm(comparingTrueVsComp[,1] ~ comparingTrueVsComp[,2])
  intercepts[dimr] <- lmx[1]$coefficients[1]
  slopes[dimr] <- lmx[1]$coefficients[2]
}
par(mfrow = c(1,2))
plot(dimrange*10, slopes, xlab = "dimension", ylab = "least squares estimate of slope for \n\"full mvn integral ~ mean pairwise bivn integral\"", main = "slope")
tdmr <- dimrange*10
lm(slopes ~ tdmr)
mean(slopes / (dimrange*10))
plot(dimrange*10, intercepts, xlab = "dimension", ylab = "least squares estimate of intercept for \n\"full mvn integral ~ mean pairwise bivn integral\"", main = "intercept")
lm(intercepts ~ tdmr)
(intercepts / (dimrange*2))
diff(intercepts)
