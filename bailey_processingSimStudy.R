setwd("/Volumes/1TB/Bailey")
nreps <- 1:500

n_monomorphic <- nreps
n_ok_traits <- nreps
for(i in nreps){
  load(file = paste0("traits_indiv_discr/traits_indiv_discr_", i))
  all_traits <- do.call(rbind, traits_indiv_discr)
  n_monomorphic[i] <- sum(sapply(1:118, function(colj) length(table(all_traits[,colj])) == 1))
  ldis <- which(sapply(1:118, function(colj) sum(table(all_traits[,colj]) > 10)) > 1)
  n_ok_traits[i] <- length(ldis)
  print(paste0(i, ": ", n_monomorphic[i], ", ", n_ok_traits[i]))
}
quantile(n_monomorphic, 0:20/20)
quantile(118 - n_ok_traits - n_monomorphic, 0:20/20)
for(i in nreps){
  means <- read.table(file = paste0("means/true_means_", i))
  load(file = paste0("data_sim/param_estimates_replicate_", i, "_optimrep_1"))
  print(cor(c(current_params$means), c(as.matrix(means))))
  plot(c(current_params$means), c(as.matrix(means)))
  abline(0, 1)
}

extraFilename <- "_noPanAsian"
corrMat <- as.matrix(read.table(paste0("output/nearPD_Corr", extraFilename, ".txt")))
n_traits <- dim(corrMat)[1]
weight <- 0.02
corrMat <- (corrMat + weight*diag(n_traits)) / (1 + weight)
for(i in nreps){
  load(file = paste0("traits_indiv_discr/traits_indiv_discr_", i))
  all_traits <- do.call(rbind, traits_indiv_discr)
  sapply(1:118, function(colj) sum(table(all_traits[,colj]) > 10)) > 3
  
  +load(file = paste0("data_sim/param_estimates_replicate_", i, "_optimrep_1"))
  cor1 <- current_params$cor
  load(file = paste0("data_sim/param_estimates_replicate_", i, "_optimrep_2"))
  cor2 <- current_params$cor
  load(file = paste0("data_sim/param_estimates_replicate_", i, "_optimrep_3"))
  cor3 <- current_params$cor
  load(file = paste0("data_sim/param_estimates_replicate_", i, "_optimrep_4"))
  cor4 <- current_params$cor
  
  # print(cor(cor1[upper.tri(cor1)], cor2[upper.tri(cor2)]))
  
  cor_est <- as.matrix(Matrix::nearPD((cor1 + cor2 + cor3 + cor4)/4, corr = T)$mat)
  cor_est <- (cor_est + weight*diag(n_traits)) / (1 + weight)
  print(cor(cor_est[upper.tri(cor_est)], corrMat[upper.tri(corrMat)]))
  
  #lotsa data inds
  ldis <- which(sapply(1:118, function(colj) sum(table(all_traits[,colj]) > 5)) > 3)
  cor_est_ld <- cor_est[ldis, ldis]
  corrMat_ld <- corrMat[ldis, ldis]
  print(cor(cor_est_ld[upper.tri(cor_est_ld)], corrMat_ld[upper.tri(corrMat_ld)]))
  plot(cor_est_ld[upper.tri(cor_est_ld)], corrMat_ld[upper.tri(corrMat_ld)])
  
  plot(cor_est[upper.tri(cor_est)], corrMat[upper.tri(corrMat)])
  abline(0, 1)
}


thresholds <- as.matrix(read.table(paste0("output/mean_of_thresholds", extraFilename, ".txt")))
for(i in nreps){
  load(file = paste0("data_sim/param_estimates_replicate_", i, "_optimrep_1"))
  
  print(cor(current_params$threshold_mat[,-1][current_params$threshold_mat[,-1] != Inf], thresholds[,-1][thresholds[,-1] != Inf]))
  plot(current_params$threshold_mat[,-1][current_params$threshold_mat[,-1] != Inf], thresholds[,-1][thresholds[,-1] != Inf])
  abline(0, 1)
}
