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

detbivcor <- function(rhos){
  sapply(rhos, function(rho) det(diag(2)*(1-rho) + rho))
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

compareCompVsFull <- function(dimr, cor = NA, meanl = -1, meanu = 1, lower = NA, upper = NA, uniformWindow = F, unifL = c(-1,0), unifU = c(0,1), repError = F){
  if(all(is.na(cor))){
    x <- rlkj(dimr)
  } else {
    x <- cor
    dimr <- dim(x)[1]
  }
  
  if(!uniformWindow){
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
  }
  
  if(uniformWindow){
    lower <- runif(n = dimr, unifL[1], unifL[2])
    upper <- lower + runif(n = dimr, unifU[1], unifU[2])
  }
  
  trueprob <- (pmvnorm(lower = lower, upper = upper, corr = x, algorithm = GenzBretz()))
  errorMag <- attr(trueprob, "error") / trueprob[1]
  trueprob <- log(trueprob[1])
  
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
  prob1 <- mean(log(probs)) #geometric mean
  # prob2 <- log(1 / (mean(1 / (probs)))) #harmonic mean
  # prob3 <- log(mean(probs)) #arithmetic mean
  prob <- prob1
  prob <- prob * dimr / 2
  if(!repError){
    return(c(trueprob, prob))
  } else if(repError){
    return(c(trueprob, prob, errorMag))
  }
}

r2_1to1 <- function(var1, var2, uplow = F){
  if(!uplow){return(1 - (var(var1 - var2) / var(var1)))}
  if(uplow){
    tot <- 1 - (var(var1 - var2) / var(var1))
    low <- var1 > var2
    up <- 1 - (var(var1[!low] - var2[!low]) / var(var1[!low]))
    low <- 1 - (var(var1[low] - var2[low]) / var(var1[low]))
    return(c(tot, up, low))
  }
}

dimr <- 80
nrep <- 500
comparingTrueVsComp <- t(replicate(nrep, compareCompVsFull(dimr, cor = NA)))
# comparingTrueVsComp[,2] <- comparingTrueVsComp[,2] / dimr
# lm(comparingTrueVsComp[,1] ~ comparingTrueVsComp[,2])
cor(comparingTrueVsComp[,1], comparingTrueVsComp[,2])^2
r2_1to1(comparingTrueVsComp[,1], comparingTrueVsComp[,2], uplow = T)
# lm(exp(comparingTrueVsComp[,1]) ~ exp(comparingTrueVsComp[,2]))

par(mar = c(5,5,3,3))
plot(comparingTrueVsComp[,1:2], xlab = "\"full\" multivariate normal integral", ylab = "composite integral obtained from \nall pairwise bivariate integrals")
title(paste0("dimension = ", dimr, ", r = ", round(cor(comparingTrueVsComp)[1,2], 3)))
abline(0, 1, lwd = 2, col = 2)
legend(x = "topleft", legend = "1-to-1 line", col = 2, lwd = 2, lty = 1)


##############################
#### making a nice figure ####
##############################

#approximation will underestimate how tiny very low full mvn probs are
#but this might also be running up against the limits of GenzBretz' abilities

dims <- 2:150
r2s <- matrix(0, nrow = c(length(dims)), ncol = 4)
nrep <- 2000
error_magnitude_threshold <- 0.25

for(dimr in dims){
  cat(paste0(dimr, " "))
  comparingTrueVsComp <- t(replicate(nrep, compareCompVsFull(dimr, cor = NA, uniformWindow = T, unifL = c(-1,0), unifU = c(0,1), repError = T)))
  if(length(which(comparingTrueVsComp == -Inf, arr.ind = T)[,1]) > 0){
    comparingTrueVsComp <- comparingTrueVsComp[-which(comparingTrueVsComp == -Inf, arr.ind = T)[,1],] #remove -Inf errors
  }
  if(length(which(comparingTrueVsComp[,3] > error_magnitude_threshold)) > 0){
    comparingTrueVsComp <- comparingTrueVsComp[-which(comparingTrueVsComp[,3] > error_magnitude_threshold),] #remove high error estimates for full multiv
  }
  if(length(which(comparingTrueVsComp[,3] == error_magnitude_threshold)) > 0){
    comparingTrueVsComp <- comparingTrueVsComp[-which(comparingTrueVsComp[,3] == error_magnitude_threshold),] #remove high error estimates for full multiv
  }
  
  r2s[dimr-1,] <- c(cor(comparingTrueVsComp[,1], comparingTrueVsComp[,2])^2, r2_1to1(comparingTrueVsComp[,1], comparingTrueVsComp[,2], uplow = T))
}

plot(r2s[,3])

# bigdiff <- F
# while(!bigdiff){
#   r <- rlkj(dimr)
#   lower <- runif(n = dimr, unifL[1], unifL[2])
#   upper <- lower + runif(n = dimr, unifU[1], unifU[2])
#   vals <- compareCompVsFull(cor = r, lower = lower, upper = upper) #oh wow the full one is super unstable!
#   if(diff(vals) > 50){
#     bigdiff = T
#   }
# }

# dev.off()


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

# 
# par(mfrow = c(1,1))
# dimr <- 2:20
# corcof <- 0.2
# vals <- sapply(dimr, function(x) compareCompVsFull(dimr = x, cor = (diag(x)*(1-corcof) + corcof), lower = rep(-0.5, x), upper = rep(0.5, x)))
# plot(dimr, vals[1,]); abline()
# dimr2 <- dimr - 2
# lm(vals[1,] ~ (dimr2))
# 
# 
# plot(log(sapply(2:50, function(x) det(diag(x)*(1-corcof) + corcof))))
# lm(log(sapply(dimr, function(x) det(diag(x)*(1-corcof) + corcof)))* 0.5 ~ dimr)
# 
# plot(t(sapply(1:9/10, function(x) compareCompVsFull(dimr = 20, cor = (diag(20)*(1-x) + x), lower = rep(-0.5, 20), upper = rep(0.5, 20)))))
# 
# f <- function(){
#   R <- rlkj(20)
#   c(log(det(R)), sum(log(detbivcor(R[upper.tri(R)]))))
# }
# 
# plot(t(replicate(100, f())))

