library(ramcmc)
library(mgcv)
library(microbenchmark)

dim <- 10
rlkj <- function (K, eta = 1) {
  alpha <- eta + (K - 2)/2
  r12 <- 2 * rbeta(1, alpha, alpha) - 1
  R <- matrix(0, K, K)
  R[1, 1] <- 1
  R[1, 2] <- r12
  R[2, 2] <- sqrt(1 - r12^2)
  if (K > 2) 
    for (m in 2:(K - 1)) {
      alpha <- alpha - 0.5
      y <- rbeta(1, m/2, alpha)
      z <- rnorm(m, 0, 1)
      z <- z/sqrt(crossprod(z)[1])
      R[1:m, m + 1] <- sqrt(y) * z
      R[m + 1, m + 1] <- sqrt(1 - y)
    }
  return(crossprod(R))
}
plotMatrix <- function(mobject, size, location, lwd = 2, grid = T, font = 1, cex = 1, rownames = T, colnames = T, title = T, title.label = "Matrix Object"){
  lines(rbind(location, location + c(0,size[2])), lwd = lwd)
  lines(rbind(location, location + c(size[1]/8,0)), lwd = lwd)
  lines(rbind(location + c(0, size[2]), location + c(size[1]/8,size[2])), lwd = lwd)
  lines(rbind(location + c(size[1],0), location + size), lwd = lwd)
  lines(rbind(location + size, location + size - c(size[1]/8,0)), lwd = lwd)
  lines(rbind(location + c(size[1],0), location + c(size[1],0) - c(size[1]/8,0)), lwd = lwd)
  if(grid == T){
    for(i in 1:(dim(mobject)[1]-1)){
      lines(rbind(location + c(0,i*size[2]/dim(mobject)[1]), location + c(size[1], i*size[2]/dim(mobject)[1])))
    }
    for(j in 1:(dim(mobject)[2]-1)){
      lines(rbind(location + c(j*size[1]/dim(mobject)[2],0), location + c(j*size[1]/dim(mobject)[2], size[2])))
    }
  }
  if(class(mobject[1,1]) != "expression" & class(mobject[1,1]) != "character"){mobject <- matrix(as.character(mobject), nrow = dim(mobject)[1], ncol = dim(mobject)[2])}
  for(i in 1:(dim(mobject)[1])){
    for(j in 1:dim(mobject)[2]){
      text(labels = mobject[i,j], x = location[1] + (j-1/2)*size[1]/dim(mobject)[2], y = location[2] + size[2] - (i-1/2)*size[2]/dim(mobject)[1], font = font, cex = cex)
    }
  }
  if(title){
    text(title.label, x = location[1] + size[1]/2, y = location[2] + size[2] + strheight(title.label, font = 2, cex = 1.5)/1.5, cex = 1.5, font = 2)
  }
  if(rownames){
    for(i in 1:dim(mobject)[1]){
      text(rownames(mobject)[i], x = location[1] - strwidth(rownames(mobject)[i])/2 - size[1]/(ncol(mobject)*6), y = location[2] + size[2] - (i-1/2)*size[2]/dim(mobject)[2])
    }
  }
  if(colnames){
    for(i in 1:dim(mobject)[1]){
      text(colnames(mobject)[i], x = location[1] + (i-1/2)*size[1]/dim(mobject)[1], y = location[2] - strheight(colnames(mobject)[i])/2- size[2]/(nrow(mobject)*6))
    }
  }
}
R <- rlkj(K = dim) #here is our starting correlation matrix
plot(1:100, 1:100, col = "white", axes = F, xlab = "", ylab = ""); plotMatrix(mobject = round(R, digits = 2), size = c(75,75), location = c(10,10), title.label = "Correlation Matrix")
ind <- sample(1:(dim-1), 1) #here is our dimension to tack on to the end
Rr <- R[-ind,-ind]
U <- chol(R) #so this matrix U is what we'd be starting with IRL
plot(1:100, 1:100, col = "white", axes = F, xlab = "", ylab = ""); plotMatrix(mobject = round(t(U), digits = 2), size = c(75,75), location = c(10,10), title.label = "Lower Cholesky Factor of Correlation Matrix")
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

dims <- round(2^((1:15)/2))+1
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
