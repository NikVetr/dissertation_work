setwd("/Volumes/1TB/Harvati")

library(ape)
library(coda)
library(rwty)
library(phangorn)
library(lattice)
library(latticeExtra)
library(MCMCpack)
library(mvMORPH)
library(phytools)

traitNumIter <- (15*(2^(0:4)))*3+1
# traitNumIter <- c(46)
nreps <- 100

rfDists <- matrix(0, nreps, 7); 
colnames(rfDists) <- c("mcc", "upgma", "nj", "pars", "pars-pca", "mk-mvc", "mk-mvc-pca")
rfDists <- lapply(1:length(traitNumIter), function(x) rfDists)

for (j in 1:length(traitNumIter)) { #iterates over trait number
  
  for(i in 1:nreps){ #iterates over replicates
    
    print(i)
    
    replicateNum <- i
    ntraits <- traitNumIter[j]
    fileName <- paste0("simulations_", ntraits, "traits_", "replicate_", replicateNum)
    
    trueTree <- read.nexus(file = paste0("trees/", fileName, "_tree.nex"))

    trees1 <- read.tree(paste0("output/", fileName, "_trees_run_1.trees") )
    trees2 <- read.tree(paste0("output/", fileName, "_trees_run_2.trees") )
    trees <- c(trees1, trees2)

    mccTree <- maxCladeCred(trees)
    njTree <- read.tree(file = paste0("output_distance/", fileName, "_neighborJoining.txt"))
    upgmaTree <- read.tree(file = paste0("output_distance/", fileName, "_upgma.txt"))
    parsTree <- read.tree(file = paste0("output_parsimony/", fileName, "_divergenceCodingParsimony.txt"))
    parsTree_pca <- read.tree(file = paste0("output_parsimony_PCA/", fileName, "_divergenceCodingParsimony.txt"))
    
    mktrees1 <- read.tree(paste0("output_discrete/", fileName, "_MVC_raw_trees_run_1.trees") )
    mktrees2 <- read.tree(paste0("output_discrete/", fileName, "_MVC_raw_trees_run_2.trees") )
    mktrees <- c(mktrees1, mktrees2)
    mccTree_mkRAW <- maxCladeCred(mktrees)
    mkPCAtrees1 <- read.tree(paste0("output_discrete/", fileName, "_MVC_PCs_trees_run_1.trees") )
    mkPCAtrees2 <- read.tree(paste0("output_discrete/", fileName, "_MVC_PCs_trees_run_2.trees") )
    mkPCAtrees <- c(mkPCAtrees1, mkPCAtrees2)
    mccTree_mkPCA <- maxCladeCred(mkPCAtrees)
    
    rfDists[[j]][i,1] <- RF.dist(trueTree, mccTree)
    rfDists[[j]][i,2] <- RF.dist(trueTree, upgmaTree)
    rfDists[[j]][i,3] <- RF.dist(trueTree, njTree)
    if(class(parsTree) == "multiPhylo"){
      rfDists[[j]][i,4] <- mean(sapply(1:length(parsTree), function(x) RF.dist(trueTree, parsTree[[x]])))
    } else if(class(parsTree) == "phylo"){
      rfDists[[j]][i,4] <- RF.dist(trueTree, parsTree)
    }    
    if(class(parsTree_pca) == "multiPhylo"){
      rfDists[[j]][i,5] <- mean(sapply(1:length(parsTree_pca), function(x) RF.dist(trueTree, parsTree_pca[[x]])))
    } else if(class(parsTree_pca) == "phylo"){
      rfDists[[j]][i,5] <- RF.dist(trueTree, parsTree_pca)
    }
    rfDists[[j]][i,6] <- RF.dist(trueTree, mccTree_mkRAW)
    rfDists[[j]][i,7] <- RF.dist(trueTree, mccTree_mkPCA)
  }
}

# save(rfDists, file =  "rfDists.txt")
load(file =  "rfDists.txt")

cols_noalpha <- brewer.pal(7, "Dark2")
alpha = 0.5 * 255
cols <- sapply(1:length(cols_noalpha), function(color) rgb(col2rgb(cols_noalpha[color])[1], 
                                                           col2rgb(cols_noalpha[color])[2], 
                                                           col2rgb(cols_noalpha[color])[3], 
                                                           maxColorValue = 255, alpha = alpha))
meths <- c("mccTree-mvBM", "upgma", "neighbor joining", "parsimony-DC_RAW", "parsimony-DC_PCA", "mccTree-mkRAW", "mccTree-mkPCA")
trimZeros <- function(x, leaveR = 1){ #returns indices of non-0 values from ends trimmed off
  inds <- min( which ( x != 0 )) : max( which( x != 0 )) 
  if(leaveR > 0){
    inds <- c(inds, inds[length(inds)] + 1:leaveR)
  }
  return(inds[inds > 0])
}
png(filename = "~/Documents/Harvati_Reanalysis_Manuscript/figures/figure3_MCC_hists.png", width = 1600, height = 2000)
par(mfrow = c(length(traitNumIter),1), mar = c(4,5,4,3))
for(i in 1:length(traitNumIter)){
  hists <- lapply(1:length(meths), function(meth) hist(rfDists[[i]][,meth], breaks = 2*0:25, plot = F))
  
  plot(hists[[1]], lty="blank", col = cols[1], xlim = c(0,50), ylim = c(0,50), xlab = "Robinson-Foulds Distance", 
       cex.lab = 2, cex.axis = 2, main = "")
  box(which = "plot")
  title(paste0(traitNumIter[i], " Traits"), cex.main = 3)
  abline(v = mean(rfDists[[i]][,1]), col = cols[[1]], lwd = 4)
  if(i==1){
    legend(x = "topright", legend = meths, fill = cols, cex = 2)
    legend(x = "topleft", legend = "Mean RF-Distance", lwd = 4, col = "darkgrey", cex = 2)
  }
  for(meth in 1:length(meths)){
    plot(hists[[meth]], lty="blank", col = cols[meth], add = T); 
    abline(v = mean(rfDists[[i]][,meth]), col = cols[[meth]], lwd = 4)
  }
  
  lines(hists[[1]]$breaks[trimZeros(hists[[1]]$counts)], hists[[1]]$counts[trimZeros(hists[[1]]$counts)], type="s",col=cols[1], add = T, lwd = 2)
  for(meth in 1:length(meths)){
    lines(hists[[meth]]$breaks[trimZeros(hists[[meth]]$counts)], hists[[meth]]$counts[trimZeros(hists[[meth]]$counts)], type="s",col=cols[meth], add = T, lwd = 2)
  }  

}
dev.off()

par(mfrow = c(1,1))
library(RColorBrewer)
cols <- brewer.pal(5, "Dark2")
bw <- 1.5
alph <- 0.25
plot(density(rfDists[[1]][,1], bw = bw), xlim = c(0,50), ylim = c(0,0.1), col = "white", main = "rf-dists to true, data-generating tree", xlab = "Robinson-Foulds Distance")
for(i in 1:5){
  polygon(density(rfDists[[1]][,i], bw = bw), col = add.alpha(cols[i], alpha = alph))
  lines(density(rfDists[[1]][,i], bw = bw), col = cols[i]) 
}
legend(legend = colnames(rfDists), fill = cols, x = "topright")


par(mfrow = c(1,1))
plot(rfDists[,3] - rfDists[,1], type = "l", ylim = c(-10, 48))
lines(rfDists[,2] - rfDists[,1], type = "l", col = 2)
lines(rfDists[,4] - rfDists[,1], type = "l", col = "blue")
lines(rfDists[,5] - rfDists[,1], type = "l", col = "green")
abline(h = 0, lty = 2)
legend(legend = c("neighbor-joining", "upgma", "parsimony", "parsimony-pca"), fill = c(1,2,"blue","green"), x = "topright"); title("difference in RF-dists from MCC Tree RF-dist")
