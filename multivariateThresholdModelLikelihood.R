#threshold model likelihood computation
library(mvtnorm)
library(data.table)

#simulate character evolution on a single branch
ntraits <- 20
startLiab <- rep(0, ntraits) #starting state
bl <- 2 #branch length for brownian motion
R <- matrix(rnorm(ntraits^2), ntraits); R <- R %*% t(R) #generate a brownian motion rate matrix
endLiab <- as.vector(rmvnorm(n = 1, mean = startLiab, sigma = R*bl)) #simulate brownian motion
nindiv <- 30 #number of individuals to sample at the tip
indivLiab <- rmvnorm(n = nindiv, mean = endLiab, sigma = R) #sample liabilities for each individual
threshold <- rep(0, ntraits) #specify the threshold of trait P/A
thresholdMatrix <- matrix(rep(threshold, nindiv), nrow = nindiv, byrow = TRUE)
indivBin <- (indivLiab > thresholdMatrix) #evaluate trait P/A for all traits -- individuals have the trait (TRUE) or do not (FALSE)

#compute likelihood of realized observations on that branch using *true* parameter values and multivariate threshold/BM
indivUpper <- matrix(0, dim(indivBin)[1], dim(indivBin)[2]); indivUpper[indivBin] <- Inf #specify matrix of upper bounds of integration for each individual's character data
indivLower <- matrix(0, dim(indivBin)[1], dim(indivBin)[2]); indivLower[!indivBin] <- -Inf #specify matrix of lower bounds
t1 <- Sys.time()
llIndiv <- sum(sapply(1:nindiv, function(x) log(pmvnorm(lower = indivLower[x,], upper = indivUpper[x,], mean = endLiab, sigma = R)[1]))) #sum over individuals
llBranch <- dmvnorm(endLiab, startLiab, R*bl, log = T) #find likelihood of brownian motion on the branch
print(paste0("log-likelihood = ", sum(llIndiv, llBranch)))
t2 <- Sys.time()
t2-t1

#compute likelihood of using *true* parameter values and univariate threshold/multivariate BM
probPres <- pnorm(q = threshold, mean = endLiab, sd = sqrt(diag(R)), log = F, lower.tail = F) #compute probabilities of trait presence under univariate threshold model
probAbs <- 1-probPres #and of trait absence
indivPres <- apply(indivBin, 2, sum) #compute number of individuals with the trait
indivAbs <- nindiv - indivPres #and number of individuals without
t1 <- Sys.time()
llIndiv.Univ <- sum(indivPres * log(probPres)) + sum(indivAbs * log(probAbs)) #sum over individuals and indepednent traits using element-wise vector multiplication
llBranch <- dmvnorm(endLiab, startLiab, R*bl, log = T)
print(paste0("log-likelihood = ", sum(llIndiv.Univ, llBranch)))
t2 <- Sys.time()
t2-t1

#compare the two calculations
niter <- 100
llcomp <- matrix(0, nrow = niter, ncol = 2)
for(i in 1:niter){
  cat(paste0(i/niter, " "))

  #sample new ending-state liabilities and rate matrices  
  R <- matrix(rnorm(ntraits^2), ntraits); R <- R %*% t(R)
  endLiab <- as.vector(rmvnorm(n = 1, mean = startLiab, sigma = R*bl))
  
  #compute likelihood of realized observations on that branch using *true* parameter values and multivariate threshold/BM
  llIndiv <- sum(sapply(1:nindiv, function(x) log(pmvnorm(lower = indivLower[x,], upper = indivUpper[x,], mean = endLiab, sigma = R)[1])))
  
  #compute likelihood of using *true* parameter values and univariate threshold/multivariate BM
  probPres <- pnorm(q = threshold, mean = endLiab, sd = sqrt(diag(R)), log = F, lower.tail = F)
  probAbs <- 1-probPres
  llIndiv.Univ <- sum(indivPres * log(probPres)) + sum(indivAbs * log(probAbs))
  
  llcomp[i,] <- c(llIndiv, llIndiv.Univ)
}

plot(llcomp, xlab = "true (but still approximated) multivariate threshold log-likelihood", ylab = "log-likelihood approximated with univariate threshold model")
title(paste0("threshold model likelihood under different parameter values: one tip, ", ntraits, " traits, ", nindiv, " individuals per tip"))
abline(0,1,lwd=2,col=2)

library(polycor); library(psych)
tr1 <- indivBin[,7]
tr2 <- indivBin[,6]
polychor(tr1, tr2)
tetrachoric(cbind(tr1,tr2))


#AB <-- this is the contingency matrix shape
#CD

deg2rad <- function(deg) return(deg*pi/180)
rtet <- function(d1, d2) {
  ctable <- table(factor(d1, levels = c(TRUE,FALSE)),factor(d2, levels = c(TRUE,FALSE)))
  
  A <- ctable[1,1]
  B <- ctable[1,2]
  C <- ctable[2,1]
  D <- ctable[2,2]
  
  return(cos(deg2rad(180/(1 + sqrt((A*D)/(B*C))))))
}
rtet(tr1, tr2)

tetmat <- function(indivBin){
  ntraits <- dim(indivBin)[2]
  mat <- matrix(0, ntraits, ntraits)
  for(i in 1:(ntraits-1)){
    for(j in (i+1):ntraits){
      mat[i,j] <- rtet(indivBin[,i], indivBin[,j])
    }
  }
  mat <- mat + t(mat); diag(mat) <- 1
  return(mat)
}

tetmat(indivBin)

#plot relationship to true correlation matrix

#simulate character evolution on a single branch
ntraits <- 50
startLiab <- rep(0, ntraits) #starting state
bl <- 3 #branch length for brownian motion
R <- matrix(rnorm(ntraits^2), ntraits); R <- R %*% t(R) #generate a brownian motion rate matrix
endLiab <- as.vector(rmvnorm(n = 1, mean = startLiab, sigma = R*bl)) #simulate brownian motion
nindiv <- 100 #number of individuals to sample at the tip
indivLiab <- rmvnorm(n = nindiv, mean = endLiab, sigma = R) #sample liabilities for each individual
threshold <- rep(0, ntraits) #specify the threshold of trait P/A
thresholdMatrix <- matrix(rep(threshold, nindiv), nrow = nindiv, byrow = TRUE)
indivBin <- (indivLiab > thresholdMatrix) #evaluate trait P/A for all traits -- individuals have the trait (TRUE) or do not (FALSE)

plot(cov2cor(R)[upper.tri(cov2cor(R))], tetmat(indivBin)[upper.tri(tetmat(indivBin))], xlim = c(-1,1), ylim = c(-1,1), xlab = "true correlations", ylab = "pairwise forumula correlations"); abline(0,1,lwd=2,col=2)

#simulate character evolution on a star phylogeny to find pooled matrix
nbranches <- 30
ntraits <- 30
startLiab <- rep(0, ntraits) #starting state
bl <- 1 #branch length for brownian motion
nindiv <- 100 #number of individuals to sample at the tip
R <- matrix(rnorm(ntraits^2), ntraits); R <- R %*% t(R) #generate a brownian motion rate matrix
indivBin <- matrix(TRUE, nrow = nindiv*nbranches , ncol = ntraits)
for(i in 1:nbranches){
  endLiab <- as.vector(rmvnorm(n = 1, mean = startLiab, sigma = R*bl)) #simulate brownian motion
  indivLiab <- rmvnorm(n = nindiv, mean = endLiab, sigma = R) #sample liabilities for each individual
  threshold <- rep(0, ntraits) #specify the threshold of trait P/A
  thresholdMatrix <- matrix(rep(threshold, nindiv), nrow = nindiv, byrow = TRUE)
  indivBin[(1+(i-1)*nindiv):(i*nindiv),] <- (indivLiab > thresholdMatrix) #evaluate trait P/A for all traits -- individuals have the trait (TRUE) or do not (FALSE)
}
plot(cov2cor(R)[upper.tri(cov2cor(R))], tetmat(indivBin)[upper.tri(tetmat(indivBin))], xlim = c(-1,1), ylim = c(-1,1), xlab = "true correlations", ylab = "pairwise forumula correlations"); abline(0,1,lwd=2,col=2)

#simulate character evolution on a pure-birth phylogeny to find pooled matrix
library(mvMORPH)
ntips <- 30
ntraits <- 30
lambda <- 0.00001
startLiab <- rep(0, ntraits) #starting state
nindiv <- 100 #number of individuals to sample at each tip
R <- matrix(rnorm(ntraits^2), ntraits); R <- R %*% t(R) #generate a brownian motion rate matrix
tree <- pbtree(n = ntips, b = lambda)
endliabs <- mvSIM(tree, param = list(sigma = R))
indivBin <- matrix(TRUE, nrow = nindiv*nbranches , ncol = ntraits)
for(i in 1:ntips){
  endLiab <- as.vector(endliabs[i,]) #simulate brownian motion
  indivLiab <- rmvnorm(n = nindiv, mean = endLiab, sigma = R) #sample liabilities for each individual
  threshold <- rep(0, ntraits) #specify the threshold of trait P/A
  thresholdMatrix <- matrix(rep(threshold, nindiv), nrow = nindiv, byrow = TRUE)
  indivBin[(1+(i-1)*nindiv):(i*nindiv),] <- (indivLiab > thresholdMatrix) #evaluate trait P/A for all traits -- individuals have the trait (TRUE) or do not (FALSE)
}
plot(cov2cor(R)[upper.tri(cov2cor(R))], tetmat(indivBin)[upper.tri(tetmat(indivBin))], xlim = c(-1,1), ylim = c(-1,1), 
     xlab = "true correlations", ylab = "pairwise forumula correlations"); abline(0,1,lwd=2,col=2)
title(paste0(ntips, " tips, ", ntraits, " traits, ", nindiv, " individuals per tip, lambda = ", lambda))
