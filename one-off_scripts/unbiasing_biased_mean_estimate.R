library(truncnorm)

loc_trunc <- 1
data <- rtruncnorm(n = 1E5, a = loc_trunc)
normal_density <- function(par, data){
  -sum(dnorm(x = data, mean = par, log = T))  
}
optim(par = 2, fn = normal_density, data = data, method = "BFGS")
etruncnorm(a = loc_trunc)
#expected value is the optimal location

range_of_integration <- -10000:10000/1000
normal_pds <- dnorm(x = range_of_integration, mean = 0)
total_prob <- sum(normal_pds[-1] * diff(range_of_integration))
total_prob
mean_of_pdf <- sum(normal_pds * range_of_integration) / length(range_of_integration)
mean_of_pdf
where_to_truncate <- 0
trunc_norm_mean_normal_pds <- dnorm(x = range_of_integration, mean = etruncnorm(a = where_to_truncate))
prob_missing <- pnorm(where_to_truncate)
trunc_norm_pds <- normal_pds * as.numeric(range_of_integration > where_to_truncate)
both_pds <- cbind(trunc_norm_pds, trunc_norm_mean_normal_pds)

diff_pds <- apply(both_pds, 1, diff)
diff_pds[diff_pds < 0] <- 0
diff_pds <- diff_pds / sum(diff_pds[-1] * diff(range_of_integration)) * prob_missing


plot(range_of_integration, trunc_norm_pds, type = "l")
lines(range_of_integration, trunc_norm_mean_normal_pds, type = "l", col = 3)
lines(range_of_integration, diff_pds, col = 2)

new_pds <- diff_pds + trunc_norm_pds
plot(range_of_integration, new_pds, type = "l")
total_prob <- sum(new_pds[-1] * diff(range_of_integration))
total_prob
mean_of_pdf <- sum(new_pds * range_of_integration) / sum(new_pds)
mean_of_pdf
etruncnorm(a = where_to_truncate)
etruncnorm(a = where_to_truncate) - mean_of_pdf

### let's make some figures ###

range_of_integration <- -1000:1000/100
normal_pds <- dnorm(x = range_of_integration, mean = 0)

range_of_truncation <- -15:15/5
diff_in_means <- rep(0, length(range_of_truncation))
diff_in_means_true <- rep(0, length(range_of_truncation))
means_of_new_pdfs <-  rep(0, length(range_of_truncation))
  
for(i in 1:length(range_of_truncation)){
  cat(paste0(i, " "))
  where_to_truncate <- range_of_truncation[i]
  mean_of_trunc <- etruncnorm(a = where_to_truncate)
  trunc_norm_mean_normal_pds <- dnorm(x = range_of_integration, mean = mean_of_trunc)
  prob_missing <- pnorm(where_to_truncate)
  trunc_norm_pds <- normal_pds * as.numeric(range_of_integration > where_to_truncate)
  both_pds <- cbind(trunc_norm_pds, trunc_norm_mean_normal_pds)
  diff_pds <- apply(both_pds, 1, diff)
  diff_pds[diff_pds < 0] <- 0
  diff_pds <- diff_pds / sum(diff_pds[-1] * diff(range_of_integration)) * prob_missing
  new_pds <- diff_pds + trunc_norm_pds
  total_prob <- sum(new_pds[-1] * diff(range_of_integration))
  mean_of_pdf <- sum(new_pds * range_of_integration) / sum(new_pds)
  diff_in_means[i] <- (mean_of_pdf - mean_of_trunc)
  diff_in_means_true[i] <- (0 - mean_of_trunc)
  means_of_new_pdfs[i] <- mean_of_pdf
  # print(mean_of_trunc)
  # print(mean_of_pdf)
}
plot(range_of_truncation, diff_in_means_true, type = "l")
lines(range_of_truncation, diff_in_means, col = 2)



nrounds <- 12
diff_in_means_acrounds <-  matrix(0, nrounds, length(range_of_truncation))
diff_in_means_acrounds[1,] <- diff_in_means
means_of_new_pdfs_acrounds <-  matrix(0, nrounds, length(range_of_truncation))
means_of_new_pdfs_acrounds[1,] <- means_of_new_pdfs
# means_of_new_pdfs_acrounds[1,] <- rep(-1, 31)

for(j in 2:nrounds){
  cat(paste0("\n", j, "\n"))
  for(i in 1:length(range_of_truncation)){
    cat(paste0(i, " "))
    where_to_truncate <- range_of_truncation[i]
    mean_of_trunc <- etruncnorm(a = where_to_truncate)
    previous_mean <- means_of_new_pdfs_acrounds[j-1,i]
    trunc_norm_mean_normal_pds <- dnorm(x = range_of_integration, mean = previous_mean)
    prob_missing <- pnorm(where_to_truncate)
    trunc_norm_pds <- normal_pds * as.numeric(range_of_integration > where_to_truncate)
    both_pds <- cbind(trunc_norm_pds, trunc_norm_mean_normal_pds)
    diff_pds <- apply(both_pds, 1, diff)
    diff_pds[diff_pds < 0] <- 0
    diff_pds <- diff_pds / sum(diff_pds[-1] * diff(range_of_integration)) * prob_missing
    new_pds <- diff_pds + trunc_norm_pds
    total_prob <- sum(new_pds[-1] * diff(range_of_integration))
    mean_of_pdf <- sum(new_pds * range_of_integration) / sum(new_pds)
    diff_in_means[i] <- (mean_of_pdf - mean_of_trunc)
    means_of_new_pdfs[i] <- mean_of_pdf
  }
  diff_in_means_acrounds[j,] <- diff_in_means
  means_of_new_pdfs_acrounds[j,] <- means_of_new_pdfs
}

means_of_new_pdfs_acrounds[,11]
for(i in 2:nrounds){
  lines(range_of_truncation, diff_in_means_acrounds[i,], col = 2)
}

# # let's try it with a sigmoidal missingness curve?
# invlogit <- function(x){exp(x)/(1+exp(x))}
# range_of_integration <- -1000:1000/100
# probs_present <- invlogit(range_of_integration)
# normal_pds <- dnorm(x = range_of_integration, mean = 0)
# missing_pdf <- normal_pds * probs_present
# mean_of_missing_pdf <- sum(missing_pdf * range_of_integration) / sum(missing_pdf)
# missing_norm_mean_normal_pds <- dnorm(x = range_of_integration, mean = mean_of_missing_pdf)
# prob_missing <- 1 - sum(missing_pdf[-1] * diff(range_of_integration))
# both_pds <- cbind(missing_pdf, missing_norm_mean_normal_pds)
# diff_pds <- apply(both_pds, 1, diff)
# diff_pds[diff_pds < 0] <- 0
# diff_pds <- diff_pds / sum(diff_pds[-1] * diff(range_of_integration)) * prob_missing
# new_pds <- diff_pds + missing_pdf
# total_prob <- sum(new_pds[-1] * diff(range_of_integration))
# mean_of_pdf <- sum(new_pds * range_of_integration) / sum(new_pds)
# (mean_of_pdf - mean_of_missing_pdf)
# 
# 
# 
# 
# 
# nrounds <- 50
# diff_in_means_acrounds <-  matrix(0, nrounds, length(range_of_truncation))
# diff_in_means_acrounds[1,] <- diff_in_means
# means_of_new_pdfs_acrounds <-  matrix(0, nrounds, length(range_of_truncation))
# means_of_new_pdfs_acrounds[1,] <- means_of_new_pdfs
# 
# for(j in 2:nrounds){
#   cat(paste0("\n", j, "\n"))
#   for(i in 1:length(range_of_truncation)){
#     cat(paste0(i, " "))
#     where_to_truncate <- range_of_truncation[i]
#     mean_of_trunc <- etruncnorm(a = where_to_truncate)
#     previous_mean <- means_of_new_pdfs_acrounds[j-1,i]
#     trunc_norm_mean_normal_pds <- dnorm(x = range_of_integration, mean = previous_mean)
#     prob_missing <- pnorm(where_to_truncate)
#     trunc_norm_pds <- normal_pds * as.numeric(range_of_integration > where_to_truncate)
#     both_pds <- cbind(trunc_norm_pds, trunc_norm_mean_normal_pds)
#     diff_pds <- apply(both_pds, 1, diff)
#     diff_pds[diff_pds < 0] <- 0
#     diff_pds <- diff_pds / sum(diff_pds[-1] * diff(range_of_integration)) * prob_missing
#     new_pds <- diff_pds + trunc_norm_pds
#     total_prob <- sum(new_pds[-1] * diff(range_of_integration))
#     mean_of_pdf <- sum(new_pds * range_of_integration) / sum(new_pds)
#     diff_in_means[i] <- (mean_of_pdf - mean_of_trunc)
#     means_of_new_pdfs[i] <- mean_of_pdf
#   }
#   diff_in_means_acrounds[j,] <- diff_in_means
#   means_of_new_pdfs_acrounds[j,] <- means_of_new_pdfs
# }
# 
# means_of_new_pdfs_acrounds[,11]
# for(i in 2:nrounds){
#   lines(range_of_truncation, diff_in_means_acrounds[i,], col = 2)
# }
