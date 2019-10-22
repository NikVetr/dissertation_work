
dim <- 10
x <- matrix(rnorm(dim^2), dim, dim)
cov <- x%*%t(x)
cor <- cov2cor(cov)
cor

L <- t(chol(cor))
sqrt(apply(L^2, 1, sum)) #confirm salient property that rows have unit length (also that diagonal elements all +)
sigma <- 0.01 #tuning parameter
prop <- matrix(0, dim, dim)
prop[7,4] <- 0.5
prop[!upper.tri(prop)] <- rnorm(dim*(dim+1)/2, 0, sigma)
#prop[10,] <- 0

Lp <- L + prop
Lp <- Lp / sqrt(apply(Lp^2, 1, sum))
sqrt(apply(Lp^2, 1, sum))
round(cor - Lp%*%t(Lp), 3)
#proposal ratio tricky -- uncomputable? many ways to get from one cholesky factor to another

#what about proposals to the cholesky factor of the covariance matrix?
L <- t(chol(cov))
sigma <- 0.1 #tuning parameter
prop <- matrix(0, dim, dim)
prop[!upper.tri(prop)] <- rnorm(dim*(dim+1)/2, 0, sigma)

Lp <- L + prop
round(cov - Lp%*%t(Lp), 3)
#just need to ensure that the diagonals are +

#perturb just the diagonals
Lp <- L; diag(Lp) <- 10
round(cov - Lp%*%t(Lp), 3)
round(cov2cor(cov) - cov2cor(Lp%*%t(Lp)), 3)

#perturb everything but diagonals
prop <- matrix(0, dim, dim)
prop[lower.tri(prop)] <- rnorm(dim*(dim-1)/2, 0, sigma)
Lp <- L + prop
round(cov - Lp%*%t(Lp), 3)
round(cov2cor(cov) - cov2cor(Lp%*%t(Lp)), 3)

#perturb last row
prop <- matrix(0, dim, dim)
prop[dim,] <- prop[dim,] + rnorm(dim, 0, sigma)
Lp <- L + prop
round(cov - Lp%*%t(Lp), 3)
round(cov2cor(cov) - cov2cor(Lp%*%t(Lp)), 3)

#perturb last 2 rows
prop <- matrix(0, dim, dim)
prop[(dim-1):dim,] <- prop[(dim-1):dim,] + rnorm(dim*2, 0, sigma)
Lp <- L + prop
round(cov - Lp%*%t(Lp), 3)
round(cov2cor(cov) - cov2cor(Lp%*%t(Lp)), 3)

#perturb all but first row
prop <- matrix(0, dim, dim)
prop[(2):dim,] <- prop[(2):dim,] + rnorm(dim*(dim-1), 0, sigma)
Lp <- L + prop
round(cov - Lp%*%t(Lp), 3)
round(cov2cor(cov) - cov2cor(Lp%*%t(Lp)), 3)

L[1,1]^2 == cov[1,1]