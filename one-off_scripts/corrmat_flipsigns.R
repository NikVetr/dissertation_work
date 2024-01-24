n = 200
nrow = 10
ncol = 20
x <- matrix(rnorm(n), nrow, ncol)
R <- cov2cor(x %*% t(x))

Rt <- R
I <- diag(nrow)
while(any(apply(Rt^2 * sign(Rt) - diag(nrow), 1, sum) < 0)){
  i <- which.min(apply(Rt^2 * sign(Rt) - diag(nrow), 1, sum))
  print(sum(apply(Rt^2 * sign(Rt) - diag(nrow), 1, sum)))
  I[i,i] <- I[i,i]*-1
  Rt <- t(R %*% I) %*% I
}
I

load("~/data/smontgom/est_gcor_mat.RData")
gcor <- gcor_mat
gcor[is.na(gcor)] <- 0
gcor <- as.matrix(Matrix::nearPD(gcor, corr = T)$mat)
gcor[abs(gcor) < 0.2] <- 0
Rt <- gcor
I <- diag(nrow(gcor))
while(any(apply(Rt^2 * sign(Rt) - diag(nrow(gcor)), 1, sum) < 0)){
  i <- which.min(apply(Rt^2 * sign(Rt) - diag(nrow(gcor)), 1, sum))
  print(paste0(round(sum(apply(Rt^2 * sign(Rt) - diag(nrow(gcor)), 1, sum)), 2), " ", sum(Rt < 0)))
  I[i,i] <- I[i,i]*-1
  Rt <- t(gcor %*% I) %*% I
}

rownames(gcor)[which(diag(I) < 0)]

# Rt <- sign(R)
# while(any(apply(Rt^2 * sign(Rt) - diag(nrow), 1, sum) < 0)){
#   i <- which.min(apply(Rt^2 * sign(Rt) - diag(nrow), 1, sum))
#   print(sum(Rt < 0))
#   I <- diag(nrow)
#   I[i,i] <- -1
#   Rt <- t(Rt %*% I) %*% I
# }
# Rt <- sign(Rt) * abs(R)
# 


chol(R * 0.5 + diag(nrow) * 0.5) - chol(R)
