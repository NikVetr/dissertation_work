library(pracma)
d <- 10
A <- rlkj(d)
L <- t(chol(A))
P_inds <- sample(1:d, d, F); 
# P_ind <- sample(1:(d-1), 1, F); 
# P_inds <- 1:10
# P_inds[d] <- P_ind; P_inds[P_ind] <- d
P <- matrix(0, d, d); for(i in 1:d){P[i,P_inds[i]] <- 1}
B <- P %*% A %*% t(P)
Q <- t(chol(B))
Qp <- t(solve(P) %*% Q)
Lp <- t(givens(Qp)$R)
Lp <- t(householder(Qp)$R) #alternatively
all(abs(Lp - L) < 1E-6) #check that cholesky factors identical
all(sapply(1:d, function(x) ((L[,x] - Lp[,x]) < 1E-6) | ((L[,x] + Lp[,x]) < 1E-6))) #or up to a sign, anyway
all(abs(A - Lp %*% t(Lp)) < 1E-6) #check that resultant correlation matrices are identical

data <- matrix(rnorm(d*50, d, 50), d, 50)
library(microbenchmark)
microbenchmark(t(chol(A)) %*% data, L %*% data, givens(Qp), householder(Qp), qr(A))
