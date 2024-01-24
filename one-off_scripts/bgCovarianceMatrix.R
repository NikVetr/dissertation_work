library(mvMORPH)
cor <- runif(1,0,1); cor <- rbind(c(1,cor), c(cor,1))
sds <- c(runif(1,0,100), runif(1,0,100))
cov <- diag(sds) %*% cor %*% diag(sds)
ntips <- 5
tree <- rtree(n = ntips, rooted = T)
mean <- c(0,0)
nsamp <- 1e4
traits <- mvSIM(tree, model = "BM1", nsim = nsamp, param = list(theta = mean, sigma = cov))
hist(sapply(1:nsamp, function(x) cor(traits[[x]])[2,1])); abline(v = cor[2,1], col = 2, lwd = 2)
#sample between-group correlations in traits are similar to rate matrix correlation

#concatenating draws returns rate matrix correlation
allTraits <- traits[[1]]; for(i in 2:nsamp) { allTraits <- rbind(allTraits, traits[[i]])}
cor
cor(allTraits) #same thing

#looks multivariate normal
plot(allTraits[1:1e3,1], allTraits[1:1e3,2])

#what is the relation of cov(allTraits) to cov?
cov(allTraits)/cov

#also, the tips come in sets of ntips, but concatenating ignores that nonindependece.
#is it relevant for the sample covariance matrix, though?

#it generally passes e.g. S-W test for multivariate normality
library(mvnormtest)
mshapiro.test(t(allTraits[1:1000,]))
