library(ape)
library(phangorn)
library(mvtnorm)
library(adephylo)

maha <- function(x1, x2, vcvm, squared = FALSE){
  squaredDist <- t(x1-x2) %*% solve(vcvm) %*% (x1-x2)
  if(!squared){
    squaredDist
  } else {
    squaredDist^0.5
  }
} 

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

ntraits <- 400
cors <- matrix(nrow = ntraits, ncol = ntraits, data = 0.2)
diag(cors) <- 1
cov <- cors
tree <- pbtree(n = 30)
treeMat <- as.matrix(distTips(tree))
traits <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = dim(cov)[1], sigma=round(cov, 4), mu= rep(0, times = dim(cov)[1])))
mm <- mahaMatrix(traits, cov)
njTree <- nj(mm)
mm/ntraits/treeMat
RF.dist(tree, njTree)
