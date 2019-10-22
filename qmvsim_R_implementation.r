library(sfsmisc) #for primes function
library(mvtnorm) 

r <- matrix(c(1, 0.2, 0.3, 0.2, 1, 0.4, 0.3, 0.4, 1), ncol = 3)
a <- c(1,2,3)
b <- c(4,5,6)
qsimvnv(25000, r, a, b)
pmvnorm(lower = a, upper = b, mean = 0, sigma = r)

dim <- 150
r <- matrix(0.5, ncol = dim, nrow = dim); diag(r) <- 1; 
r <- cov2cor(rWishart(n = 1, df = dim*2, Sigma = r)[1:dim, 1:dim,])
a <- rep(0, dim)
b <- sample(2:5, dim, replace = T)
t1 <- Sys.time()
qsimvnv(5e3, r, a, b)
t2 <- Sys.time()
pmvnorm(lower = a, upper = b, mean = 0, sigma = r)
t3 <- Sys.time()
paste0("me = ", t2-t1, ", R = ", t3-t2)


# R = [1 0.2 0.3; 0.2 1 0.4; 0.3 0.4 1]
# a = [1,2,3]
# b = [4,5,6]

r <- matrix(c(5,0.8,0.8,1), ncol = 2)
a <- c(-100,-100)
b <- c(0,0)
qsimvnv(1000, r, a, b)
pmvnorm(lower = a, upper = b, mean = 0, sigma = r)

dim <- 11
x <- matrix(rnorm(100), dim)
r <- cov2cor(x%*%t(x))
a <- rep(0, dim)
b <- sample(2:10, size = dim, replace = T)
qsimvnv(100000, r, a, b)
pmvnorm(lower = a, upper = b, mean = 0, sigma = r, algorithm = GenzBretz())



#functions

#using built in tools lol
erfc <- function (x) { 2 * pnorm(-sqrt(2) * x) }

erfcInv <- function (x) { qnorm(x * 0.5) / -sqrt(2)}

phi <- function(z) {erfc( -z / sqrt(2)) / 2}

phiInv <- function(p) {-sqrt(2)*erfcInv(2*p)}

permchl <- function(corr, upper, lower){
  dc <- dim(corr)[1]
  permCorr <- corr[dc:1, dc:1]
  return(list(t(chol(permCorr)), upper[dc:1], lower[dc:1]))
}

permchl(r,a,b)

qsimvnv <- function(m, r, a, b){
  foo <- permchl(r, a, b); ch <- foo[[1]]; as <- foo[[2]]; bs <- foo[[3]]
  ct <- ch[1,1]
  ai <- as[1]
  bi <- bs[1]
  if(ai > -9*ct){
    if(ai < 9*ct){
      c <- Phi(ai/ct)
    } else {
      c = 1
    }} else {
      c = 0
    }
  
  if(bi > -9*ct){
    if(bi < 9*ct){
      d <- Phi(bi/ct)
    } else {
      d = 1
    }} else {
      d = 0
    }
  
  n <- dim(r)[1]
  ci <- c
  dci <- d - ci
  p = 0
  e = 0
  ps = sqrt(primes(5*n*log(n+1)/4)); q <- t(ps[1:(n-1)])
  ns <- 12
  nv <- trunc(max(m/ns, 1))
  on = rep(1, nv)
  y = matrix(rep(0, (n-1)*nv), nrow = n-1, ncol = nv)
  
  for(j in 1:ns){
    c <- ci*on
    dc <- dci*on
    pv <- dc
    for(i in 2:n){
      x <- abs(2*((q[i-1]%*%1:nv) + runif(n = 1, min = 0, max = 1)) %% 1 - 1)
      y[i-1,] = phiInv(c + x*dc)
      s <- ch[i,(1:(i-1))]%*%y[1:(i-1),]
      ct <- ch[i,i]
      ai <- as[i] - s
      bi <- bs[i] - s
      c = on
      d = c
      c[which(ai < -9*ct)] <- 0
      d[which(bi < -9*ct)] <- 0
      tstl <- which(abs(ai) < 9*ct)
      c[tstl] <- Phi(ai[tstl]/ct)
      tstl <- which(abs(bi) < 9*ct)
      d[tstl] <- Phi(bi[tstl]/ct)
      dc <- d - c
      pv <- pv * dc
    }
    d <- (mean(pv) - p) / j
    p <- p + d
    e <- (j - 2) * e / j + d^2
  }
  e <- 3 * sqrt(e)
  return(list(p, e))
}

qsimvnv(1000, r, a, b)

#implementing functions from scratch 
erfc <- function (x) { 2 * pnorm(-sqrt(2) * x) }

Phi <- function(z) {erfc( -z / sqrt(2)) / 2}

seqcheck <- function(x){
  if(x[1]<=x[length(x)]) {return(x)} else {return(0)}
}

chlrdr <- function(r, a, b){
  ep <- 1E-10
  n <- dim(r)[1]
  c = r
  ap = a
  bp = b
  d = sqrt(diag(c))
  for(i in 1:n){
    if(d[i] > 0){
      c[,i] <- c[,i] / d[i]
      c[i,] <- c[i,] / d[i]
      ap[i] = ap[i] / d[i]
      bp[i] = bp[i] / d[i]
    }
  }
  y <- rep(0, n)
  sqtp <- sqrt(2*pi)
  for(k in 1:n){
    im = k
    print(paste(k, im))
    ckk = 0
    dem = 1
    s = 0
    for(i in k:n){
      if(c[i,i] > ep){
        cii <- sqrt(max(c(c[i,i],0)))
        if(i > 1){
          s <- c[i,(1:(k-1))] %*% y[1:(k-1)] 
        }
        ai <- (ap[i]-s) / cii 
        bi = (bp[i]-s) / cii
        de = Phi(bi) - Phi(ai)
        if(de <= dem){
          ckk <- cii
          dem <- de
          am <- ai
          bm <- bi
          im <- i
        }
      }
    }
    print(paste(k, im))
    if(im > k){
      c[im,im] <- c[k,k]
      ap[c(im,k)] <- ap[c(k,im)]
      bp[c(im,k)] <- bp[c(k,im)]
      c[c(im,k),1:(k-1)] <- c[c(k,im),1:(k-1)] 
      c[((im+1):n),c(im,k)] <- c[((im+1):n),c(k,im)]
      t <- c[(k+1):(im-1),k]
      c[(k+1):(im-1),k] <- t(c[im, (k+1):(im-1)])
      c[im,(k+1):(im-1)] <- t(t)
    }
    c[k,(k+1):n] <- 0
    if(ckk > ep*k){
      c[k,k] <- ckk
      for(i in (k+1):n){
        c[i,k] <- c[i,k] / ckk
        c[i,(k+1):i] <- c[i,(k+1):i] - c[i,k]*t(c[(k+1):i,k])
      }
    if(abs(dem) > ep){
      y[k] <- (exp(-am^2/2) - exp(-bm^2/2))/(sqtp*dem)
      } else {
      if(am < -10){
        y[k] <- bm
      } else if(bm > 10){
        y[k] <- am
      } else {
        y[k] <- (am + bm) / 2
        }
      }
    } else {
      c[k:n,k] <- 0
      y[k] <- 0
    }
  }
  return(c(c, ap, bp))
}

chlrdr(r,a,b)

chlrdr2 <- function(r, a, b){
  ep <- 1E-10
  n <- dim(r)[1]
  c = r
  ap = a
  bp = b
  d = sqrt(diag(c))
  for(i in 1:n){
    if(d[i] > 0){
      c[,i] <- c[,i] / d[i]
      c[i,] <- c[i,] / d[i]
      ap[i] = ap[i] / d[i]
      bp[i] = bp[i] / d[i]
    }
  }
  y <- rep(0, n)
  sqtp <- sqrt(2*pi)
  for(k in 1:n){
    im = k
    print(paste(k, im))
    ckk = 0
    dem = 1
    s = 0
    for(i in k:n){
      if(c[i,i] > ep){
        cii <- sqrt(max(c(c[i,i],0)))
        if(i > 1){
          s <- c[i,(1:(k-1))] %*% y[1:(k-1)] 
        }
        ai <- (ap[i]-s) / cii 
        bi = (bp[i]-s) / cii
        de = Phi(bi) - Phi(ai)
        if(de <= dem){
          ckk <- cii
          dem <- de
          am <- ai
          bm <- bi
          im <- i
        }
      }
    }
    print(paste(k, im))
    if(im > k){
      c[im,im] <- c[k,k]
      ap[c(im,k)] <- ap[c(k,im)]
      bp[c(im,k)] <- bp[c(k,im)]
      c[c(im,k),1:(k-1)] <- c[c(k,im),1:(k-1)] 
      c[seqcheck((im+1):n),c(im,k)] <- c[seqcheck((im+1):n),c(k,im)]
      t <- c[(k+1):(im-1),k]
      c[seqcheck((k+1):(im-1)),k] <- t(c[im, seqcheck((k+1):(im-1))])
      c[im,(k+1):(im-1)] <- t(t)
    }
    c[k,seqcheck((k+1):n)] <- 0
    if(ckk > ep*k){
      c[k,k] <- ckk
      for(i in (k+1):n){
        c[i,k] <- c[i,k] / ckk
        c[i,(k+1):i] <- c[i,seqcheck((k+1):i)] - c[i,k]*t(c[seqcheck((k+1):i),k])
      }
      if(abs(dem) > ep){
        y[k] <- (exp(-am^2/2) - exp(-bm^2/2))/(sqtp*dem)
      } else {
        if(am < -10){
          y[k] <- bm
        } else if(bm > 10){
          y[k] <- am
        } else {
          y[k] <- (am + bm) / 2
        }
      }
    } else {
      c[k:n,k] <- 0
      y[k] <- 0
    }
  }
  return(list(c, ap, bp))
}

chlrdr2(r,a,b)




# primes(n) in matlab returns a row vector of the prime numbers less than or equal to n
# log is base e like in r
# fix(x) returns a vector where each element is rounded towards the nearest integer in the directin of zero, e.g. (-2.1, 4.6) becomes (-2, 4) 
# max returns a scalar in the meaningful case
# ones() and zeros() return matrices of the supplied dimensions filled with 1s and 0s, respectively
# rand returns runif(min = 0, max = 1)
# find() in matlab is which() in R