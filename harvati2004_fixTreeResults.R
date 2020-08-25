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

#making figure 4 for the manuscript

tree_m <- trees_m[[1]]
tree_f <- trees_f[[1]]
branchRates_m_mean_rescaled <- branchRates_m_mean / min(branchRates_m_mean)
branchRates_f_mean_rescaled  <- branchRates_f_mean / min(branchRates_f_mean)

tree_m$tip.label[getDescendants(tree_m, tree_m$edge[which(matched_m_branches == 1),2])[getDescendants(tree_m, tree_m$edge[which(matched_m_branches == 1),2]) <= 13]]

coph_plot_rates <- cophylo(tree_m, tree_f)

library(viridis)
cols <- plasma(150)
edge.col<-list(
  left=cols[round(branchRates_m_mean_rescaled * 10)],
  right=cols[round(branchRates_f_mean_rescaled * 10)])
plot(coph_plot_rates, edge.col = edge.col, lwd = 3)

#match branches in one tree to those in the other
corresponding_f_branches <- rep(NA, length(branchRates_m_mean_rescaled))
for(i in 1:length(branchRates_m_mean_rescaled)){
  #these are the male tips descended from the branch subtending their MRCA
  desc <- tree_m$tip.label[getDescendants(tree_m, tree_m$edge[i,2])[getDescendants(tree_m, tree_m$edge[i,2]) <= 13]]
  if(all(desc == "Homo_sapiens")){
    corresponding_f_branches[i] <- NA #for the only sapiens tip
  } else{
    if(any(desc == "Homo_neanderthalensis")){
      desc <- desc[-which(desc == "Homo_neanderthalensis")]
    }
    if(length(desc) > 1){
      corresponding_f_branches[i] <- which(tree_f$edge[,2] == getMRCA(phy = tree_f, tip = desc))
    } else if(length(desc == 1)) {
      corresponding_f_branches[i] <- which(tree_f$edge[,2] == which(tree_f$tip.label == desc))
    } else{
      corresponding_f_branches[i] <- NA #for the neandertal tip
    }
  }
}
plot(branchRates_m_mean_rescaled, branchRates_f_mean_rescaled[corresponding_f_branches])
matched_m_branches <- branchRates_m_mean_rescaled
matched_f_branches <- branchRates_f_mean_rescaled[corresponding_f_branches]
matched_m_branches <- matched_m_branches[!is.na(matched_f_branches)]
matched_f_branches <- matched_f_branches[!is.na(matched_f_branches)]
plot(matched_m_branches, matched_f_branches)
abline(0,1)
# tree$edge.length <- tree$edge.length * getBranchRates(params_m, tree)



#making the actual figure
png(filename = "~/Documents/Harvati_Reanalysis_Manuscript/figures/figure4_final.png", width = 2200, height = 750)
layout(matrix(c(1,1,2), nrow = 1, ncol = 3, byrow = TRUE))
par(xpd=TRUE)

#cophylo plot
plot.cophylo(coph_plot_rates, mar=c(9,1,4,1), lwd = 6, fsize = 3.5, link.lwd=4, edge.col = edge.col)

title("     Males                           Females", cex.main = 4.5, line = 0)

box(which = "figure", lty = 5)
text(labels = "a)", cex = textsiz, x = 0.53, y = 1.065, xpd = T, font = 4, col = "darkred")


xl <- -0.3; yb <- -0.125; xr <- 0.3; yt <- -0.075;
rect(
  head(seq(xl,xr,(xr-xl)/150),-1),
  yb,
  tail(seq(xl,xr,(xr-xl)/150),-1),
  yt,
  col=cols
)
text(labels = c("1x", 2:14, "15x"), y = yb-(yt-yb)/2.5, x = seq(xl,xr,length.out = 15), las=2, cex=2.5)
text(labels = "rate variation", x = (xl+xr)/2, y = yt+(yt-yb)/2.5, font = 2, cex = 3)

#branch rates plot
par(mar=c(8,9,8,2))
plot(matched_m_branches, matched_f_branches, main = "Expected Branch Rates", xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
     cex.axis = 3, cex.main = 4.5, pch = 16, cex = 5, col = rgb(0,0,0,0.5), xlim = c(0,15), ylim = c(0,15))
box(lwd=3)
title(ylab = "females", xlab = "males", cex.lab = 4, line = 5.75)
axis(1, at = 0:15, labels = rep("", 16), lwd = 4, cex.axis = 3, tck = -0.015)
mtext(text = 0:15, side = 1, at = 0:15, cex = 2, line = 2.25)
axis(2, at = 0:15, labels =  rep("", 16), lwd = 4, cex.axis = 3, tck = -0.015)
mtext(text = 0:15, side = 2, at = 0:15, cex = 2, line = 2)
abline(0, 1, lty = 2, lwd = 3, xpd = F)
box(which = "figure", lty = 5)
text(labels = "b)", cex = textsiz, x = 15.5, y = 17, xpd = T, font = 4, col = "darkred")

dev.off()

#making figure 5, comparing the inferred mvBM rate matrix with the P matrix
traitRates_m_mean <- apply(traitRates_m, 1, mean)
traitRates_f_mean <- apply(traitRates_f, 1, mean)
traitRates_m_rescaled <- sapply(1:ncol(traitRates_m), function(samp) traitRates_m[,samp] / (traitRates_m[21,samp]))
traitRates_f_rescaled <- sapply(1:ncol(traitRates_f), function(samp) traitRates_f[,samp] / (traitRates_f[21,samp]))
traitRates_m_89 <- apply(traitRates_m_rescaled, 1, rethinking::HPDI)
traitRates_f_89 <- apply(traitRates_f_rescaled, 1, rethinking::HPDI)


traitRates_m_mean_rescaled <- traitRates_m_mean / min(traitRates_m_mean)
P_m_vars_rescaled <- P_m_vars / min(P_m_vars)
traitRates_f_mean_rescaled <- traitRates_f_mean / min(traitRates_f_mean)
P_f_vars_rescaled <- P_f_vars / min(P_f_vars)

traitRates_m_rescaled_perc <- sapply(1:nrow(traitRates_m_rescaled), function(tr) sum(traitRates_m_rescaled[tr,] > P_m_vars_rescaled[tr]) / ncol(traitRates_m_rescaled))
traitRates_f_rescaled_perc <- sapply(1:nrow(traitRates_f_rescaled), function(tr) sum(traitRates_f_rescaled[tr,] > P_f_vars_rescaled[tr]) / ncol(traitRates_f_rescaled))


corrMatrices_f_vec <- t(sapply(1:n_samps_f, function(cr) corrMatrices_f[,,cr][upper.tri(corrMatrices_f[,,cr])]))
corrMatrices_m_vec <- t(sapply(1:n_samps_m, function(cr) corrMatrices_m[,,cr][upper.tri(corrMatrices_m[,,cr])]))

corrMatrices_f_mean[upper.tri(corrMatrices_f_mean)]
P_f_corrs[upper.tri(corrMatrices_f_mean)]

corrMatrices_f_perc <- sapply(1:length(P_f_corrs[upper.tri(corrMatrices_f_mean)]), 
                              function(cr) sum(corrMatrices_f_vec[,cr] > P_f_corrs[upper.tri(corrMatrices_f_mean)][cr]) / nrow(corrMatrices_f_vec))
corrMatrices_m_perc <- sapply(1:length(P_m_corrs[upper.tri(corrMatrices_m_mean)]), 
                              function(cr) sum(corrMatrices_m_vec[,cr] > P_m_corrs[upper.tri(corrMatrices_m_mean)][cr]) / nrow(corrMatrices_m_vec))







png(filename = "~/Documents/Harvati_Reanalysis_Manuscript/figures/figure5_final.png", width = 2000, height = 500)
layout(matrix(c(1,1,2,2,4,4,5,5,1,1,3,3,4,4,6,6), nrow = 2, ncol = 8, byrow = TRUE))

par(mar=c(8,8,8,3))
plot(log10(traitRates_m_mean_rescaled), log10(P_m_vars_rescaled), pch = 25, xlim = c(0,2.5), ylim = c(0,2.5), 
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", cex.axis = 3, cex.main = 4.5, cex = 2.5, bg = rgb(1,1,1,0.25))
abline(0,1, lty = 2, lwd = 4, col = "grey30")
points(log10(traitRates_f_mean_rescaled), log10(P_f_vars_rescaled), lty = 2, pch = 19, cex = 2.5, col = rgb(0,0,0,0.65)); 
for(i in 1:46){
  segments(x0 = log10(traitRates_m_mean_rescaled)[i], y0 = log10(P_m_vars_rescaled)[i], 
           x1 = log10(traitRates_f_mean_rescaled)[i], y1 = log10(P_f_vars_rescaled)[i], lwd = 2)
}
axis(1, at = 0:5/2, labels = rep("", 6), lwd = 4, cex.axis = 3, tck = -0.015)
mtext(text = sapply(0:5/2, function(i) as.expression(bquote(10^ .(i)))), side = 1, at = 0:5/2 + c(0,0.05,0,0.05,0,0.05), cex = 2, line = 2.25)
axis(2, at = 0:5/2, labels = rep("", 6), lwd = 4, cex.axis = 3, tck = -0.015)
mtext(text = sapply(0:5/2, function(i) as.expression(bquote(10^ .(i)))), side = 2, at = 0:5/2 + c(0,0.05,0,0.05,0,0.05), cex = 2, line = 1)

legend(x = "bottomright", legend = c("females", "males"), pch = c(19, 25), cex = 2.5, box.lty = 3, box.lwd = 2, col = c(rgb(0,0,0,0.65), 1), bg = c(NA, 1))
text(labels = "a)", cex = textsiz, x = 2.665, y = 3.035, xpd = T, font = 4, col = "darkred")

title(xlab = "mvBM relative rates", cex.lab = 4, line = 5.75)
title(ylab = "P-matrix relative variances", cex.lab = 4, line = 4.75)
title(main = "Rates and Variances", cex.main = 5)
box(lwd=4)
box(which = "figure", lty = 5)

par(mar=c(6,6.5,6,1))
par(lwd = 2)
hist(traitRates_m_rescaled_perc, breaks = 5, xlab = "", ylab = "", main = "", cex.axis = 2.5, lwd = 3, ylim = c(0,20), col = "grey75")
text(labels = "b)", cex = textsiz, x = 1.015, y = 26.5, xpd = T, font = 4, col = "darkred")
title(xlab = "percentile", cex.lab = 3, line = 3.75)
title(ylab = "frequency", cex.lab = 3, line = 3.5)
title(main = "Male Rate Percentiles", cex.main = 4)
par(lwd = 1)
box(which = "figure", lty = 5)
par(lwd = 2)
hist(traitRates_f_rescaled_perc, breaks = 5, xlab = "", ylab = "", main = "", cex.axis = 2.5, lwd = 3, ylim = c(0,20), col = "grey75")
title(xlab = "percentile", cex.lab = 3, line = 3.75)
title(ylab = "frequency", cex.lab = 3, line = 3.5)
title(main = "Female Rate Percentiles", cex.main = 4)
par(lwd = 1)
box(which = "figure", lty = 5)

####

par(mar=c(8,8,8,3))
plot(corrMatrices_f_mean[upper.tri(corrMatrices_f_mean)], P_f_corrs[upper.tri(corrMatrices_f_mean)], pch = 19, ylim = c(-1,1), 
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", cex.axis = 3, cex.main = 4.5, cex = 2.5, col = rgb(0,0,0,0.65))
abline(0,1, lty = 2, lwd = 4, col = "grey30")
points(corrMatrices_m_mean[upper.tri(corrMatrices_m_mean)], P_m_corrs[upper.tri(corrMatrices_m_mean)], lty = 2, pch = 25, cex = 2.5, bg = rgb(1,1,1,0.25)); 
# for(i in 1:length(corrMatrices_m_mean[upper.tri(corrMatrices_m_mean)])){
#   segments(x0 = corrMatrices_m_mean[upper.tri(corrMatrices_m_mean)][i], y0 = P_m_corrs[upper.tri(corrMatrices_m_mean)][i], 
#            x1 = corrMatrices_f_mean[upper.tri(corrMatrices_f_mean)][i], y1 = P_f_corrs[upper.tri(corrMatrices_f_mean)][i], lwd = 2)
# }
axis(1, at = -2:2/10, labels = rep("", 6), lwd = 4, cex.axis = 3, tck = -0.015)
mtext(text = -2:2/10, side = 1, at = -2:2/10, cex = 2, line = 2.25)
axis(2, at = -2:2/2, labels = rep("", 5), lwd = 4, cex.axis = 3, tck = -0.015)
mtext(text = -2:2/2, side = 2, at = -2:2/2, cex = 2, line = 1)

legend(x = "bottomright", legend = c("females", "males"), pch = c(19, 25), cex = 2.5, box.lty = 3, box.lwd = 2, col = c(rgb(0,0,0,0.65), 1), bg = c(NA, 1))
text(labels = "c)", cex = textsiz, x = 0.2565, y = 1.4225, xpd = T, font = 4, col = "darkred")

title(xlab = "mvBM correlations", cex.lab = 4, line = 5.75)
title(ylab = "P-matrix correlations", cex.lab = 4, line = 4.75)
title(main = "Correlations", cex.main = 5)
box(lwd=4)
box(which = "figure", lty = 5)

###

par(mar=c(6,6.5,6,1))
par(lwd = 2)
hist(corrMatrices_m_perc, breaks = 5, xlab = "", ylab = "", main = "", cex.axis = 2.5, lwd = 3, ylim = c(0,400), col = "grey75")
text(labels = "d)", cex = textsiz, x = 1.015, y = 525, xpd = T, font = 4, col = "darkred")
title(xlab = "percentile", cex.lab = 3, line = 3.75)
title(ylab = "frequency", cex.lab = 3, line = 3.5)
title(main = "Male Correlation\nPercentiles", cex.main = 4, line = -2.25)
par(lwd = 1)
box(which = "figure", lty = 5)
par(lwd = 2)
hist(corrMatrices_f_perc, breaks = 5, xlab = "", ylab = "", main = "", cex.axis = 2.5, lwd = 3, ylim = c(0,400), col = "grey75")
title(xlab = "percentile", cex.lab = 3, line = 3.75)
title(ylab = "frequency", cex.lab = 3, line = 3.5)
title(main = "Female Correlation\nPercentiles", cex.main = 4, line = -2.25)
par(lwd = 1)
box(which = "figure", lty = 5)

dev.off()





#comparetrees
par(mar=c(8,9,8,2))
plot(comparison_uvBM[,1], comparison_uvBM[,2] + 1/3, main = "Compare-Trees Plot and Histogram", xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
     cex.axis = 3, cex.main = 4.5, pch = 16, cex = 5, col = rgb(0,0,0,0.5), xlim = c(0,1), ylim = c(0,1+1/3))
box(lwd=3)
title(ylab = "molecular", xlab = "morphological", cex.lab = 4, line = 5.75)
axis(1, at = 0:5/5, labels = rep("", 6), lwd = 4, cex.axis = 3, tck = -0.015)
mtext(text = 0:5/5, side = 1, at = 0:5/5, cex = 2, line = 2.25)
axis(2, at = 0:5/5 + 1/3, labels =  rep("", 6), lwd = 4, cex.axis = 3, tck = -0.015)
mtext(text = 0:5/5, side = 2, at = 0:5/5 + 1/3, cex = 2, line = 2)
segments(x0 = 0, x1 = 1, y0 = 1/3, y1 = 4/3, lwd  = 3, lty = 2)

#hists
par(xpd=F)
abline(h = 0.25, lty = 3, lwd = 3)
axis(2, at = 0:4/5/4, labels =  rep("", 5), lwd = 4, cex.axis = 3, tck = -0.0075)
mtext(text = 0:4/4, side = 2, at = 0:4/5/4, cex = 1.25, line = 1, las=2)

highprobs <- discretizeContData(comparison_uvBM[comparison_uvBM[,2] > 0.5,1])
lowprobs <- discretizeContData(comparison_uvBM[comparison_uvBM[,2] < 0.5,1])
width = 0.1
for(i in 1:ncol(highprobs)){
  rect(xleft = highprobs[2,i] - width/2, xright = highprobs[2,i] + width/2, ybottom = 0, ytop = highprobs[1,i] / 5,
       col = rgb(0,0,0,0.5))
}
for(i in 1:ncol(lowprobs)){
  rect(xleft = lowprobs[2,i] - width/2, xright = lowprobs[2,i] + width/2, ybottom = 0, ytop = lowprobs[1,i] / 5,
       col = rgb(red = 1,1,1,0.5))
}


legend(x = c(0.325, 0.675), y = c(0.24, 0.12), 
       legend = c("molecular probability < 0.5", "molecular probability > 0.5"), 
       fill = c(rgb(red = 1,1,1,0.5), rgb(0,0,0,0.5)), cex = 1.8)

box(which = "figure", lty = 5)







corrMatrices_m_mean <- corrMatrices_m[,,1]
for(i in 2:n_samps_m){corrMatrices_m_mean <- corrMatrices_m_mean + corrMatrices_m[,,i]}
corrMatrices_m_mean <- corrMatrices_m_mean / n_samps_m

corrMatrices_f_mean <- corrMatrices_f[,,1]
for(i in 2:n_samps_f){corrMatrices_f_mean <- corrMatrices_f_mean + corrMatrices_f[,,i]}
corrMatrices_f_mean <- corrMatrices_f_mean / n_samps_f

plot(corrMatrices_f_mean, corrMatrices_m_mean); abline(0,1)
plot(corrMatrices_f_mean[upper.tri(corrMatrices_f_mean)], P_f_corrs[upper.tri(corrMatrices_f_mean)]); abline(0, 1)
plot(corrMatrices_m_mean[upper.tri(corrMatrices_m_mean)], P_m_corrs[upper.tri(corrMatrices_m_mean)]); abline(0, 1)
