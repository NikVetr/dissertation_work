library(ramcmc)
library(mgcv)

dim <- 200
R <- rlkj(1, dim) #here is our starting correlation matrix
ind <- sample(1:(dim-1), 1) #here is our dimension to tack on to the end
Rr <- R[-ind,-ind]
U <- chol(R) #so this matrix U is what we'd be starting with IRL
Ur <- matrix(0, dim-1, dim-1)
Ur[1:(ind-1), 1:(ind-1)] <- U[1:(ind-1), 1:(ind-1)]
Ur[1:(ind-1),ind:(dim-1)] <- U[1:(ind-1),(ind+1):(dim)]
Ur[ind:(dim-1),ind:(dim-1)] <- chol(t(U[(ind+1):(dim),(ind+1):(dim)]) %*% U[(ind+1):(dim),(ind+1):(dim)] + (t(t(U[ind,(ind+1):dim])) %*% t(U[ind,(ind+1):dim])))
#implementing https://en.wikipedia.org/wiki/Cholesky_decomposition#Rank-one_update will avoid above re-decomposition
all(round(t(Ur) %*% Ur - Rr, 6) == 0) #yep, looks like we are getting the right reduced matric back

#and here it is already implemented: equivalently
Ur[ind:(dim-1),ind:(dim-1)] <- t(chol_update(t(U[(ind+1):(dim),(ind+1):(dim)]), U[ind,(ind+1):dim]))
all(round(t(Ur) %*% Ur - Rr, 6) == 0)

#or just use
Ur <- choldrop(U, ind)
all(round(t(Ur) %*% Ur - Rr, 6) == 0)

#so that's our row/column-deleted cholesky factor, now to add it back in at the end
Rn <- R[c((1:dim)[-ind],ind), c((1:dim)[-ind],ind)] #here is the new matrix we're targetting
URn <- chol(Rn) #and here is its cholesky factor
URn_target <- URn[,dim] #and we're specifically trying to get this column of it, since we've already gotten the reduced form above
Rn_target <- Rn[,dim] #these are the correlations we need, but if we haven't recomposed the entire correlation matrix, we can get them by
Rn_target <- c((t(U) %*% U[,ind])[-ind],1) #I mean, why would we ever want to look at actual correlations?? *scoff*
URn_target_sle <- t(solve(Ur)) %*% (Rn_target[1:(dim-1)]) #to find our final column, we solve a system of linear equations
URn_target_sle <- backsolve(Ur, (Rn_target[1:(dim-1)]), transpose = T) #already in beautifully convenient triangular form, so we can solve in O(n^2)
URn_target_sle <- c(URn_target_sle, sqrt(1 - crossprod(URn_target_sle))) #still need that last element, which we find by enforcing unit length
all(round(URn_target - URn_target_sle, 6) == 0) #huzzah, they are all the same
Un <- cbind(rbind(Ur, rep(0,dim-1)), URn_target_sle) #now let's slap it on to the cholesky factor from the removal step
all(round(t(Un) %*% Un - Rn, 6) == 0) #and check for equivalence!

library(microbenchmark)
microbenchmark(times = 500,
               chol(Rr), 
               chol_update(t(U[(ind+1):(dim),(ind+1):(dim)]), U[ind,(ind+1):dim]), 
               choldrop(U, ind), 
               backsolve(Ur, (Rn_target[1:(dim-1)]), transpose = T))

choleskalator <- function(mat, ind){
  dim <- dim(mat)[1]
  ind_removed <- mgcv::choldrop(mat, ind)
  Rn_target <- c((t(mat) %*% mat[,ind])[-ind],1)  
  URn_target_sle <- backsolve(ind_removed, (Rn_target[1:(dim-1)]), transpose = T) 
  URn_target_sle <- c(URn_target_sle, sqrt(1 - crossprod(URn_target_sle))) 
  return(cbind(rbind(ind_removed, rep(0,dim-1)), URn_target_sle))
}

dims <- round(2^((1:22)/2))+1
times <- matrix(data = 0, nrow = length(dims), ncol = 3, dimnames = list(1:length(dims), 1:3))
times[,1] <- dims; colnames(times)[1] <- "dim"

for(i in 1:length(dims)){
  cat(paste0(i, " "))
  dim <- dims[i]
  ind <- round(dim/2)
  cormat <- rlkj(dim)
  cholfac <- chol(cormat)
  mbm <- microbenchmark(times = 5,
               chol(cormat), 
               choleskalator(cholfac, ind))
  times[i,2] <- mean(mbm$time[mbm$expr == levels(mbm$expr)[1]])
  times[i,3] <- mean(mbm$time[mbm$expr == levels(mbm$expr)[2]])
}

colnames(times)[2] <- levels(mbm$expr)[1]
colnames(times)[3] <- levels(mbm$expr)[2]

plot(times[,1], times[,2], type = "l", col = "red", xlab = "dimensionality", ylab = "time")
lines(times[,1], times[,3], type = "l", col = "green")

times[,2] / times[,3] 
