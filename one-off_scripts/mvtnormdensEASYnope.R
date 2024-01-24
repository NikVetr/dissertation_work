library(mvtnorm)
#sample a covariance matrix
dim <- 2
C <- matrix(rnorm(n = dim^2), nrow = 2)
R <- C%*%t(C)
R <- cov2cor(R)

#simulate node 1's traits
mea <- rep(0, dim)
vals <- rmvnorm(n=1, mean = mea, sigma = R)
vals <- rep(rmvnorm(n=1, mean = mea, sigma = R)[1],2)

#calculate multivariate normal density
mvdens <- dmvnorm(rep(0, dim), vals, sigma = R, log = T)

#eigentransform the points to PC scores
valsEig <- vals %*% eigen(R)$vectors
meanEig <- mea %*% eigen(R)$vectors

#calculate the product of univariate normal densities of transformed chars
uvdens <- dmvnorm(meanEig, valsEig, sigma = diag(eigen(R)$values), log = T) #done with mvtnorm density where sigma is a diag matrix
vars <- eigen(R)$values
uvdens <- sum(sapply(1:dim, function(x) dnorm(meanEig[x], valsEig[x], sd = (vars[x])^0.5, log = T))) #same as above uvdens
mvdens
uvdens

#calculate the cdf of the multivariate normal distribution
pmvnorm(lower = c(-Inf, -Inf), upper = as.vector(vals), mean = mea, corr = R)
pmvnorm(lower = c(-Inf, -Inf), upper = as.vector(valsEig) / eigen(R)$values^0.5, mean = as.vector(meanEig), corr = diag(1, dim))
prod(pnorm(q = valsEig, mean = meanEig, sd = eigen(R)$values^0.5))/8
#nope
phi <- -sum(vals^2)^0.5
omeg <- seq(-1,0,0.001)
y <- -omeg + phi
norms <- dnorm(x = omeg, mean = 0, sd = eigen(R)$values[1]^0.5)
omega <- omeg[which.min(abs(y-norms))]
yc <- y[which.min(abs(y-norms))]
p1 <- (abs(omega - phi)*yc/2 + pnorm(q = omega, mean = meanEig[1], sd = eigen(R)$values[1]^0.5))

qs <- seq(-1,0,0.0001)
y <- -qs -phi
norms <- dnorm(x = qs, mean = 0, sd = eigen(R)$values[2]^0.5)
Q <- qs[which.min(abs(y-norms))]
yc <- y[which.min(abs(y-norms))]
p2 <- pnorm(-Q, sd = eigen(R)$values[2]^0.5) - pnorm(Q, sd = eigen(R)$values[2]^0.5) - ((yc+phi)*-Q/2 + Q*phi)*2
p1*p2
