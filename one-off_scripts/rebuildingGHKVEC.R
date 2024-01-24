library(bayesm)
library(mvtnorm)

ntraits <- 10
cov <- matrix(data = .2, nrow = ntraits, ncol = ntraits)
diag(cov) <- 1
ltm <- t(chol(cov))
truncation <- rep(0, ntraits)
aboveOrBelow <- rep(1, ntraits)
nreps <- 1e4

ghkvec(ltm, truncation, aboveOrBelow, nreps, HALTON = F)
ghkvecR(ltm, truncation, aboveOrBelow, nreps)
pmvnorm(upper = truncation, mean = rep(0, ntraits), corr = cov, algorithm = GenzBretz())[1]


#code this bad boy up in R

L <- ltm
r <- nreps
trunpt <- truncation
above <- aboveOrBelow

ghk_oneside <- function(L, trunpt, above, r){
  dim <- length(trunpt)
  z <- rep(0, dim)
  res <- 0
  udraw <- runif(n = r*dim, min = 0, max = 1)
  i <- 0
  while(i < r){
    # print(paste0("i is ", i))
    prod <- 1
    j <- 0
    while(j < dim){
      # print(paste0("j is ", j))
      mu <- 0
      k <- 0
      while(k < j){
        # print(paste0("k is ", k))
        mu <- mu + L[k*dim+j+1] * z[k+1] #index the lower cholesky matrix by a single number??
        k <- k + 1
      }
      tpz <- (trunpt[j+1]-mu)/L[j*dim+j+1] #do it here too lol
      if(above[j+1] > 0){
        pa <- 1
        pb <- pnorm(tpz,0,1,1,0)
      } else {
        pb <- 1
        pa <- pnorm(tpz,0,1,1,0) 
      }
      prod <- prod * (pb - pa)
      u <- udraw[i*dim+j+1]
      arg <- u*pb + (1-u)*pa
      if (arg > .999999999){
        arg <- .999999999
      }
      if (arg < .0000000001){
        arg <- .0000000001
      }
      z[j+1] <- qnorm(arg,0,1,1,0)
      j <- j + 1 
    }
    res <- res + prod
    i <- i + 1
  }
  res <- res / r
  res
}

ghkvecR <- function(L, trunpt, above, r){
  dim <- length(above)
  n <- length(trunpt)/dim
  res <- rep(0, n)
  for(i in 1:n){
    res[i] <- ghk_oneside(as.vector(L), trunpt, above, r)
  }
  res
}
