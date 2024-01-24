library(rethinking)
library(expm)
dimR <- 10
R <- rlkjcorr(1, dimR)
eigensystem <- eigen(R)
V <- eigensystem$vectors
L <- eigensystem$values

rotation = function(x,y){
  u=x/sqrt(sum(x^2))
  v=y-sum(u*y)*u
  v=v/sqrt(sum(v^2))
  cost=sum(x*y)/sqrt(sum(x^2))/sqrt(sum(y^2))
  sint=sqrt(1-cost^2);
  diag(length(x)) - u %*% t(u) - v %*% t(v) + 
    cbind(u,v) %*% matrix(c(cost,-sint,sint,cost), 2) %*% t(cbind(u,v))
}

tuneVars <- 0.01
V1P <-  rnorm(dimR, mean = 0, sd = sqrt(tuneVars))
V1P <- V1P/ sqrt(sum(V1P^2))
rot <- rotation(as.vector(V[,1]), V1P)

#can rejection sample random uniform vectors until you get one in some arbitrary bounds?
distance <- 1
runitvec <- function(dim){v <- rnorm(dim); return(v / sqrt(sum(v^2)))}
min(as.vector(replicate(1E4, as.matrix(dist(rbind(runitvec(dimR), V)))[-1,1])))
#hmmmm... space is just too large, fails for curse of dimensionality reasons

#fractional rotation?
frac <- 100
rot_frac <- expm(logm(rot) / frac)
V[,1] %*% rot - V1P
rotvec <- V[,1]
for(rots in 1:frac){rotvec <- rotvec %*% rot_frac}
rotvec - V1P

#seems to work, and also has desired symmetry
#let's try to rotate eigenvectors then
Vr <- t(t(V) %*% rot_frac)
Vr %*% t(Vr)
V %*% diag(L) %*% t(V) - R
(Vr %*% diag(L) %*% t(Vr)) #ack, not a correlation matrix :'[
