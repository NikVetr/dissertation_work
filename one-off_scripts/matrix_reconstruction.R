#if you project data such that the rate matrix is the identity, and then infer a rate matrix, you infer a bit about the rates, too, but not much
dim = 1000
eta = 10
a <- rlkj(dim)
b <- rlkj(dim, eta = eta)
plot(a, cov2cor(t(chol(a)) %*% b %*% (chol(a))))
hist(diag(t(chol(a)) %*% b %*% (chol(a))))
