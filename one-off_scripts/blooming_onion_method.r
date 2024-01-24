dim <- 20
x <- matrix(rnorm(dim^2), dim, dim)
x <- x %*% t(x)
eigen(x)
x <- cov2cor(x)
z <- -solve(x)
diag(z) <- -diag(z)
p <- cov2cor(z)
q <- -p; diag(q) <- -diag(q)
x <- cov2cor(solve(q))
x <- matrix(0.9,10,10); diag(x) <- 1

blooming_onion <- function(cor = 1, hmean = 0, eta = 1, d){
  r <- diag(rep(1,d))
  if(!identical(cor, 1)){d_s <- dim(cor)[1]} else {d_s <- 1}
  B <- eta + (d-2)/2
  if(d_s == 1){
    u <- rbeta(1, shape1 = B,shape2 = B)
    rij <- 2*u-1
    r[1,2] <- r[2,1] <- rij
  } else {
    r[1:d_s,1:d_s] <- cor
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
load("rateMatrix"); cov <- as.matrix(cov)

d <- 500
r <- blooming_onion(cov2cor(cov),0.5,eta=1,d=300)
r <- blooming_onion(1,hmean = 0,eta=1,d=d)
par(mfrow = c(2,1))
hist(r[upper.tri(r)])
# plot(eigen(r)$values, type = "l"); abline(h=0, col = 2, lty = 2)
