cov <- matrix(rnorm(100), 10, 10)
cov <- cov %*% t(cov)
X <- rnorm(10)
U <- chol(cov)
L <- t(U)
Y <- solve(L) %*% X
dmvnorm(x = X, mean = rep(0,10), sigma = cov, log = T)
sum(dnorm(x = Y, log = T, mean = 0, sd = 1)) - log(det(cov))/2
log(prod(dnorm(x = Y, log = F, mean = 0, sd = 1)) / (det(cov))^0.5)

BL <- .5
Yc <- Y/BL^0.5
dmvnorm(x = X, mean = rep(0,10), sigma = BL*cov, log = T)
sum(dnorm(x = Yc, log = T, mean = 0, sd = 1)) - log(det(cov)*BL^10)/2

library(MCMCpack)
cov
rwish(v = 500, cov/500)
solve(riwish(v = 500, solve(cov/500)))

solve(riwish(v = 500, (cov/500)))
solve(rwish(v = 500, 3/500))

sum(3 > sapply(1:1e4, function(x) solve(riwish(v = 500, solve(3/500)))))

sum(3 > sapply(1:1e4, function(x) rwish(v = 500, 3/500)))
solve(riwish(1000, cov))/1000

c <- matrix(c(1,.99,.99,1),2)
sum(0.99 > sapply(1:1e5, function(x) cov2cor(rwish(1e2, c/100))[1,2])) #sampled correlation most often further from 0
mean(sapply(1:1e6, function(x) cov2cor(rwish(1e2, c/100))[1,2])) #but still has expectation sigma
sum(1 > sapply(1:1e5, function(x) (rwish(1e2, c/100))[1,1])) #sampled variance smaller

rwish(v = 500, cov)
cov/(500-dim(cov)[1]-1)

c[1,1]
mean(sapply(1:100000, function(x) rwish(v = 500, S = c/500)[1,1]))


cov <- matrix(rnorm(100), 10)
cov <- cov %*% t(cov)
vecs <- eigen(cov)$vectors
vals <- eigen(cov)$values


round(vecs %*% diag(vals) %*% t(vecs) - cov, digits = 8) 

reords <- sample(1:dim(vecs)[2])
vecs  <- vecs[,reords]
vals <- vals[reords]

round(vecs %*% diag(vals) %*% t(vecs) - cov, digits = 8) 
