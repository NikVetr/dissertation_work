library(ape)
library(phangorn)
library(mvtnorm)
library(mvMORPH)
library(adephylo)
library(phytools)

#compute mahalanobis distance between two points in multivariate space given a within-group covariance matrix
maha <- function(x1, x2, vcvm, squared = FALSE){
  squaredDist <- t(x1-x2) %*% solve(vcvm) %*% (x1-x2)
  if(!squared){
    squaredDist
  } else {
    squaredDist^0.5
  }
} 

#compute matrix of mahalanobis distances
mahaMatrix <- function(data, vcvm){
  mat <- matrix (nrow = nrow(data), ncol = nrow(data), 0)
  for(i in 1:nrow(mat)){
    for(j in i:nrow(mat)){
      mat[i,j] <- maha(data[i,], data[j,], vcvm)
    }
  }
  mat <- mat + t(mat)
  rownames(mat) <- colnames(mat) <-rownames(data)
  mat
}

ntraits <- 300
#generate an (arbitrary) covariance matrix
cors <- matrix(nrow = ntraits, ncol = ntraits, data = 0.2)
diag(cors) <- 1
sds <- sample(1:20, ntraits, replace = T)
cov <- diag(sds) %*% cors %*% diag(sds)
#generate a random tree
tree <- rtree(n=8)
treeMat <- as.matrix(distTips(tree))
#simulate traits
traits <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = dim(cov)[1], sigma=round(cov, 4), mu= rep(0, times = dim(cov)[1])))
#compute mahalanobis distance matrix
mm <- mahaMatrix(traits, cov)
# estimate tree via neighbor joining
njTree <- nj(mm)

#off diagonal elements approach 1 as ntraits increases, demonstrating linear proportionality
mm/ntraits/treeMat

#RF distance likewise approaches 0 as ntraits increases
RF.dist(tree, njTree)
