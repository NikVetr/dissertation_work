

x_node <- rnorm(n = 4e5, mean = 0, sd = 1) 
x1 <- x_node + rnorm(n = 4e5, mean = 0, sd = 1.5^0.5)
x2 <- x_node + rnorm(n = 4e5, mean = 0, sd = 0.5^0.5)
library(MASS)
bivn <- kde2d(x1, x2, n = 100)
image(bivn, xlim = c(-4,4), ylim = c(-3,3), 
      main = "bivariate density of x1, x2", xlab = "x1", ylab = "x2")
contour(bivn, add = T)


library(mvtnorm)
x1x2 <- rmvnorm(n = 4e5, mean = c(0,0), 
                sigma = matrix(c(2.5, 1, 1, 1.5), nrow = 2))
bivn <- kde2d(x1x2[,1], x1x2[,2], n = 100)
image(bivn, xlim = c(-4,4), ylim = c(-3,3), 
      main = "bivariate density of x1, x2", xlab = "x1", ylab = "x2")
contour(bivn, add = T)
