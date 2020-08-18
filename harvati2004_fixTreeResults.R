setwd("/Volumes/1TB/Harvati_Empirical/")

library(ape)
library(coda)
library(rwty)
library(phangorn)
library(lattice)
library(latticeExtra)
library(MCMCpack)
library(mvMORPH)
library(phytools)
library(RevGadgets)

meansNexus <- function (characterMatrix, outputFilePath, traitsOnColumns = T){
  if (traitsOnColumns) {ntax = dim(characterMatrix)[1]; nchar = dim(characterMatrix)[2]}
  else {ntax = dim(characterMatrix)[2]; nchar = dim(characterMatrix)[1]}
  sink(outputFilePath, append=F)
  cat("#NEXUS\n\nBegin data;\n\tDimensions ntax=", ntax, "nchar=", nchar, ";\n\tFormat datatype=Continuous missing=?;\n\tMatrix\n\n")
  for(i in 1:length(characterMatrix[,1])) {cat(rownames(characterMatrix)[i], characterMatrix[i,], "\n")}
  cat(";\nEnd;")
  sink()
}

trees_m <- read.tree("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_males_noHomopops_justSpecies_uninf_mvBM_fixTreeInfRM_aggrProps_c1.trees")
trees_f <- read.tree("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_females_noHomopops_justSpecies_uninf_mvBM_fixTreeInfRM_aggrProps_c1.trees")

params_m_1 <- read.table("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_males_noHomopops_justSpecies_uninf_mvBM_fixTreeInfRM_aggrProps_c1.log", header = T)
params_m_2 <- read.table("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_males_noHomopops_justSpecies_uninf_mvBM_fixTreeInfRM_aggrProps_c2.log", header = T)
params_m_3 <- read.table("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_males_noHomopops_justSpecies_uninf_mvBM_mixTreeInfRM_c1.log", header = T)
params_m_4 <- read.table("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_males_noHomopops_justSpecies_uninf_mvBM_mixTreeInfRM_c2.log", header = T)
params_m_all <- rbind(params_m_1, params_m_2, params_m_3, params_m_4)
# sort(coda::effectiveSize(params_m_all))

params_f_1 <- read.table("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_females_noHomopops_justSpecies_uninf_mvBM_fixTreeInfRM_aggrProps_c1.log", header = T)
params_f_2 <- read.table("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_females_noHomopops_justSpecies_uninf_mvBM_fixTreeInfRM_aggrProps_c2.log", header = T)
params_f_3 <- read.table("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_females_noHomopops_justSpecies_uninf_mvBM_fixTreeInfRM_c1.log", header = T)
params_f_4 <- read.table("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_females_noHomopops_justSpecies_uninf_mvBM_fixTreeInfRM_c2.log", header = T)
params_f_all <- rbind(params_f_1, params_f_2, params_f_3, params_f_4)
# sort(coda::effectiveSize(params_f_all))

recomposeRateMatrix <- function(params, sds_keyword = "relative_sd", corrs_keyword = "correlation", ntraits = NA){
  #get the standard deviations
  sds <- params[grep(sds_keyword, colnames(params))]
  if(!all(diff(as.numeric(sapply(1:length(sds), function(x) strsplit(colnames(sds)[[x]], split = "\\.")[[1]][2]))) == 1)){
    stop("error something has gone terrible wrong")
  }
  sds <- diag(as.numeric(sds))
  
  if(is.na(ntraits)){
    ntraits <- nrow(sds)
  }
  
  #get the correlations
  corrmat <- diag(ntraits)
  upper_tri_corrs <- params[grep(corrs_keyword, colnames(params))]
  if(!all(diff(as.numeric(sapply(1:length(upper_tri_corrs), function(x) strsplit(colnames(upper_tri_corrs)[[x]], split = "\\.")[[1]][2]))) == 1)){
    stop("error something done goofed")
  }
  corrmat[lower.tri(corrmat)] <- as.numeric(upper_tri_corrs)
  corrmat <- corrmat + t(corrmat) - diag(ntraits)
  
  #recompose and return
  mvBM_R <- sds %*% corrmat %*% sds
  return(mvBM_R)
}

matchNodes = function(phy){
  
  # get some useful info
  num_tips = length(phy$tip.label)
  num_nodes = phy$Nnode
  tip_indexes = 1:num_tips
  node_indexes = num_tips + num_nodes:1
  
  node_map = data.frame(R=1:(num_tips + num_nodes), Rev=NA, visits=0)
  current_node = phy$Nnode + 2
  k = 1
  t = 1
  
  while(TRUE) {
    
    if ( current_node <= num_tips ) {
      node_map$Rev[node_map$R == current_node] = t
      current_node = phy$edge[phy$edge[,2] == current_node,1]
      t = t + 1
    } else {
      
      if ( node_map$visits[node_map$R == current_node] == 0 ) {
        node_map$Rev[node_map$R == current_node] = node_indexes[k]
        k = k + 1
      }
      node_map$visits[node_map$R == current_node] = node_map$visits[node_map$R == current_node] + 1
      
      if ( node_map$visits[node_map$R == current_node] == 1 ) {
        # go right
        current_node = phy$edge[phy$edge[,1] == current_node,2][2]
      } else if ( node_map$visits[node_map$R == current_node] == 2 ) {
        # go left
        current_node = phy$edge[phy$edge[,1] == current_node,2][1]
      } else if ( node_map$visits[node_map$R == current_node] == 3 ) {
        # go down
        if (current_node == num_tips + 1) {
          break
        } else {
          current_node = phy$edge[phy$edge[,2] == current_node,1]
        }
      }
    }
    
  }
  
  return(node_map[,1:2])
  
}
# trees_f[[1]]$edge.length[order(trees_f[[1]]$edge[,2])][order(matchNodes(trees_f[[1]])[,2])]

getBranchRates <- function(params, tree, keyword = "branch_rates_"){
  br <- params[grep(keyword, colnames(params))]
  if(!all(diff(as.numeric(sapply(1:length(br), function(x) strsplit(colnames(br)[[x]], split = "\\.")[[1]][2]))) == 1)){
    stop("error something has gone terrible wrong")
  }
  br <- as.numeric(br)
  
  #get table of tree edges
  tr_edge <- cbind(edge_ind = 1:length(tree$edge.length), tree$edge)
  tr_edge <- rbind(tr_edge, c(NA, NA, length(tree$tip.label)+1))
  tr_edge <- tr_edge[order(tr_edge[,3]),] #reorder by tip index
  edge_indices_of_desc_tip <- tr_edge[,1]
  
  #construct table of relations between indices
  matcher <- cbind(matchNodes(tree), R_edge_ind = edge_indices_of_desc_tip)
  #reorder to privelege rb
  matcher <- matcher[order(matcher$Rev),]
  #get br's in there
  matcher <- cbind(matcher, rb_br = c(br, NA))
  #reorder by phylo's edge lengths
  matcher <- matcher[order(matcher$R_edge_ind),]
  reordered_brs <- matcher$rb_br
  if(is.na(reordered_brs[length(reordered_brs)])){
    reordered_brs <- reordered_brs[-length(reordered_brs)]
  } else {
    stop("something's afoot with the indexing here")
  }
  return(reordered_brs)
}

ntraits <- 46
n_samps_m <- nrow(params_m_all)
n_samps_f <- nrow(params_f_all) 
rateMatrices_m <- array(data = 0, dim = c(ntraits, ntraits, n_samps_m))
rateMatrices_f <- array(data = 0, dim = c(ntraits, ntraits, n_samps_f))
corrMatrices_m <- array(data = 0, dim = c(ntraits, ntraits, n_samps_m))
corrMatrices_f <- array(data = 0, dim = c(ntraits, ntraits, n_samps_f))
traitRates_m <- matrix(0, ntraits, n_samps_m)
traitRates_f <- matrix(0, ntraits, n_samps_f)
branchRates_m <- matrix(0, 24, n_samps_m)
branchRates_f <- matrix(0, 22, n_samps_f)

for(i in 1:n_samps_m){
  if(i%%10 == 0){cat(paste0(i, " "))}
  rateMatrices_m[,,i] <- recomposeRateMatrix(params_m_all[i,])
  corrMatrices_m[,,i] <- cov2cor(rateMatrices_m[,,i])
  traitRates_m[,i] <- diag(rateMatrices_m[,,i])
  branchRates_m[,i] <- getBranchRates(params_m_all[i,], tree = trees_m[[1]])
}

for(i in 1:n_samps_f){
  if(i%%10 == 0){cat(paste0(i, " "))}
  rateMatrices_f[,,i] <- recomposeRateMatrix(params_f_all[i,])
  corrMatrices_f[,,i] <- cov2cor(rateMatrices_f[,,i])
  traitRates_f[,i] <- diag(rateMatrices_f[,,i])
  branchRates_f[,i] <- getBranchRates(params_f_all[i,], tree = trees_f[[1]])
}

branchRates_m_mean <- apply(branchRates_m, 1, mean)
branchRates_f_mean <- apply(branchRates_f, 1, mean)
traitRates_m_mean <- apply(traitRates_m, 1, mean)
traitRates_f_mean <- apply(traitRates_f, 1, mean)

P_m <- as.matrix(read.table("data2_neanF/Harvati_noHomoPops_males_pooledCov_P.txt"))
P_f <- as.matrix(read.table("data2_neanF/Harvati_noHomoPops_females_pooledCov_P.txt"))
P_m_vars <- diag(P_m)
P_m_vars <- P_m_vars / mean(P_m_vars)
P_f_vars <- diag(P_f)
P_f_vars <- P_f_vars / mean(P_f_vars)
P_m_corrs <- cov2cor(P_m)
P_f_corrs <- cov2cor(P_f)

plot(log(traitRates_m_mean), log(P_m_vars)); abline(0,1)
plot(log(traitRates_f_mean), log(P_f_vars)); abline(0,1)

corrMatrices_m_mean <- corrMatrices_m[,,1]
for(i in 2:n_samps_m){corrMatrices_m_mean <- corrMatrices_m_mean + corrMatrices_m[,,i]}
corrMatrices_m_mean <- corrMatrices_m_mean / n_samps_m

corrMatrices_f_mean <- corrMatrices_f[,,1]
for(i in 2:n_samps_f){corrMatrices_f_mean <- corrMatrices_f_mean + corrMatrices_f[,,i]}
corrMatrices_f_mean <- corrMatrices_f_mean / n_samps_f

plot(corrMatrices_f_mean, corrMatrices_m_mean); abline(0,1)
plot(corrMatrices_f_mean[upper.tri(corrMatrices_f_mean)], P_f_corrs[upper.tri(corrMatrices_f_mean)]); abline(0, 1)
plot(corrMatrices_m_mean[upper.tri(corrMatrices_m_mean)], P_m_corrs[upper.tri(corrMatrices_m_mean)]); abline(0, 1)
