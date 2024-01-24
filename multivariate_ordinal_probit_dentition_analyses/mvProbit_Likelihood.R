library(mvProbit)
## generate a simulated data set
set.seed( 123 )
# number of observations
nObs <- 10
# generate explanatory variables
xMat <- cbind(
  const = rep( 1, nObs ),
  x1 = as.numeric( rnorm( nObs ) > 0 ),
  x2 = as.numeric( rnorm( nObs ) > 0 ),
  x3 = rnorm( nObs ),
  x4 = rnorm( nObs ) )
# model coefficients
beta <- cbind( c( 0.8, 1.2, -1.0, 1.4, -0.8 ),
               c( -0.6, 1.0, 0.6, -1.2, -1.6 ),
               c( 0.5, -0.6, -0.7, 1.1, 1.2 ) )
# covariance matrix of error terms
library( miscTools )
sigma <- symMatrix( c( 1, 0.2, 0.4, 1, -0.1, 1 ) )
# generate dependent variables
yMatLin <- xMat %*% beta
yMat <- ( yMatLin + rmvnorm( nObs, sigma = sigma ) ) > 0
colnames( yMat ) <- paste( "y", 1:3, sep = "" )
# log likelihood values
myData <- as.data.frame( cbind( xMat, yMat ) )
logLikVal <- mvProbitLogLik( cbind( y1, y2, y3 ) ~ x1 + x2 + x3 + x4,
                             coef = c( beta ), sigma = sigma, data = myData, algorithm = TVPACK() )
print( logLikVal )







library(mvProbit)
## generate a simulated data set
set.seed( 123 )
# number of observations
nObs <- 10
# generate explanatory variables
xMat <- cbind(
  const = rep( 1, nObs ),
  x1 = as.numeric( rnorm( nObs ) > 0 ),
  x2 = as.numeric( rnorm( nObs ) > 0 ),
  x3 = rnorm( nObs ),
  x4 = rnorm( nObs ) )

xMat <- cbind(
  const = rep( 1, nObs ))

xMat <- cbind(
  const = rep( 1, nObs ),
  x1 = as.numeric( rnorm( nObs ) > 0 ))

# model coefficients
beta <- cbind( c( 0.8, 1.2, -1.0, 1.4, -0.8 ),
               c( -0.6, 1.0, 0.6, -1.2, -1.6 ),
               c( 0.5, -0.6, -0.7, 1.1, 1.2 ) )

beta <- runif(n = D, min = -2, max = 2)
beta <- rbind(runif(n = D, min = -2, max = 2), runif(n = D, min = -2, max = 2))

# covariance matrix of error terms
library( miscTools )
sigma <- symMatrix( c( 1, 0.2, 0.4, 1, -0.1, 1 ) )
sigma <- Omega
# generate dependent variables
yMatLin <- xMat %*% beta
yMat <- ( yMatLin + rmvnorm( nObs, sigma = sigma ) ) > 0
colnames( yMat ) <- paste( "y", 1:D, sep = "" )
# log likelihood values
myData <- as.data.frame( cbind( xMat, yMat ) )
cat(paste( "y", 1:D, ", ", sep = "" ))
logLikVal <- mvProbitLogLik( cbind( y1,  y2,  y3,  y4,  y5,  y6,  y7,  y8,  y9,  y10,  y11,  y12,  y13,  y14,  y15,  y16,  y17,  y18,  y19,  y20,  y21,  y22,  y23,  y24,  y25,  y26,  y27,  y28,  y29,  y30,  y31,  y32,  y33,  y34,  y35,  y36,  y37,  y38,  y39,  y40 ) ~ x1,
                             coef = c( beta ), sigma = sigma, data = myData, algorithm = GenzBretz() )
print( logLikVal )


library(cubature)

mu1 <- matrix(c(3,3), nrow=2)
sigma1 <- rbind(c(4,-1), c(-1,6))

quadratic <- function(a,b) {
  X <- matrix(c(a,b),nrow=2)
  Q <- (-1/2)*t(X-mu1)%*%solve(sigma1)%*%(X-mu1)
}

NormalPDF <- function(x) {
  f <- (1/(2*pi))*(1/sqrt(det(sigma1)))*exp(quadratic(x[1],x[2]))
}
# Solving for P(1 < X1 < 3, 1 < X2 < 3)
P <- adaptIntegrate( NormalPDF, lowerLimit= sample(1:30, size = 20),  upperLimit=sample(40:60, size = 20))
P$integral

library(mvtnorm)
ntraits <- 100
x <- matrix(rnorm(ntraits^2), ntraits); x <- x %*% t(x)
c <- cov2cor(x); #c <- diag(2)
pmvnorm(lower = rep(0,ntraits), upper = rep(Inf, ntraits), mean = rep(0,ntraits), corr = c)
