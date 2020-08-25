library(ape)
library(phytools)
library(phangorn)

colfunc <- colorRampPalette(c("white", "black"))
discretizeContData <- function(data, range = c(0,1), bins = 10){
  bintervals <- seq(from = range[1], to = range[2], length.out = bins+1)
  binMeans <- bintervals[-1] + diff(bintervals) / 2 - bintervals[2]
  discrData <- sapply(data, function(num) sum(num - bintervals > 0))
  discrFreqs <- sapply(1:bins, function(bin) sum(discrData == bin))
  discrFreqs <- discrFreqs / sum(discrFreqs)
  return(rbind(freqs = discrFreqs, binLocations = binMeans))
}
startsWith2 <- function(x, pre){apply(sapply(1:length(pre), function(n) startsWith(x, pre[n])), 1, any)}
greedyCT <- function(trees){
  tipNames <- trees[[1]]$tip.label
  tipNums <- paste0("t", 1:length(tipNames))
  for(i in 1:length(trees)){
    trees[[i]]$tip.label <- tipNums[match(trees[[i]]$tip.label, tipNames)]
  }
  ape::write.tree(trees, file = "gct_trees.txt")
  writeLines("gct_trees.txt\nY", con = "consense_script.txt")
  if(file.exists("outfile")){file.remove("outfile")}
  if(file.exists("outtree")){file.remove("outtree")}
  system(paste0("/Applications/phylip-3.695/exe/consense < consense_script.txt"), intern = T)
  gct <- ape::read.tree("outtree")
  if(file.exists("outfile")){file.remove("outfile")}
  if(file.exists("outtree")){file.remove("outtree")}
  if(file.exists("gct_trees.txt")){file.remove("gct_trees.txt")}
  if(file.exists("consense_script.txt")){file.remove("consense_script.txt")}
  gct$tip.label <- tipNames[match(gct$tip.label, tipNums)]
  return(gct)
}
convNum2Str <- function(nums, key){
  sapply(1:length(nums), function(x) key[nums[x]])
}
prop.part.df <- function(trees, cutoff = 0.01, bs = T){
  if(class(trees) == "multiPhylo"){
    if(trees[[1]]$Nnode == (length(trees[[1]]$tip.label) - 1)){
      trees <- unroot(trees) #unroot rooted trees
    }
    tipLabs <- trees[[1]]$tip.label
    numTrees <- length(trees)
    if(bs) {
      out <- as.prop.part(bitsplits(trees))
    } else {
      out <- prop.part(trees)
    }
    outList <- as.list.data.frame(out)
    pps <- attributes(outList)$number/numTrees
    props <- data.frame(pps, as.matrix(outList)); colnames(props) <- c("postProbs", "clade")
    props <- props[order(-pps),]
    props <- props[props[,1] > cutoff,]
    rownames(props) <- 1:nrow(props)
    props$cladeNames <- sapply(1:length(props[,1]), function(x) sort(convNum2Str(props$clade[[x]], attributes(outList)$labels)))
    props <- props[,c(1,3)]
    if(!bs) {
      props <- props[-1,]
    }
    allClades <- c(props$cladeNames, lapply(1:length(props$postProbs), function(x) sort(setdiff(tipLabs, props$cladeNames[[x]]))))
    if(any(duplicated(allClades))){
      # print("duplicate clades found")
      dupClades <- allClades[duplicated(allClades)]
      dupCladesComp <- lapply(1:length(dupClades), function(x) sort(setdiff(tipLabs, dupClades[[x]])))
      matchedDupes <- cbind(1:length(dupClades), sapply(1:length(dupClades), function(x) which(sapply(1:length(dupCladesComp), function(y) setequal(dupClades[[x]], dupCladesComp[[y]])))))
      dupes <- matrix(data = 0, nrow = length(matchedDupes[,1])/2, ncol = 2)
      for(i in 1:length(matchedDupes[,1])){
        pair <- sort(matchedDupes[i,])
        if(!any(sapply(1:length(dupes[,1]), function(x) pair == dupes[x,]))){
          dupes[min(which(apply(dupes, 1, sum) == 0)),] <- pair
        }
      }
      for(i in 1:length(dupes[,1])){
        clade1 <- dupClades[dupes[i,1]][[1]]
        clade2 <- dupClades[dupes[i,2]][[1]]
        propsInd1 <- which(sapply(1:length(props[,1]), function(x) setequal(clade1, props$cladeNames[[x]])))
        propsInd2 <- which(sapply(1:length(props[,1]), function(x) setequal(clade2, props$cladeNames[[x]])))
        props[propsInd1,1] <- props[propsInd1,1] + props[propsInd2,1]
        props <- props[-propsInd2,]
      }
      props <- props[order(props[,1], decreasing = T),]
    }
    props
  } else if (class(trees) == "phylo"){
    if(trees$Nnode == (length(trees$tip.label) - 1)){
      trees <- unroot(trees) #unroot rooted trees
    }
    tipLabs <- trees$tip.label
    numTrees <- 1
    if(bs) {
      out <- as.prop.part(bitsplits(trees))
    } else {
      out <- prop.part(trees)
    }
    outList <- as.list.data.frame(out)
    pps <- attributes(outList)$number/numTrees
    props <- data.frame(pps, as.matrix(outList)); colnames(props) <- c("postProbs", "clade")
    props <- props[order(-pps),]
    props <- props[props[,1] > cutoff,]
    rownames(props) <- 1:nrow(props)
    props$cladeNames <- sapply(1:length(props[,1]), function(x) sort(convNum2Str(props$clade[[x]], attributes(outList)$labels)))
    props <- props[,c(1,3)]
    if(!bs) {
      props <- props[-1,]
    }
    props
  }
}
compareTrees <- function(biparts1, biparts2, tipLabs){
  matchProbs <- matrix(0, nrow = sum(length(biparts1[,1]), length(biparts2[,1])), ncol = 2)
  counter <- 1
  biparts2$notSeen <- 1
  for(clade in 1:length(biparts1[,2])){
    cladeName <- biparts1[,2][[clade]]
    isThere <- sapply(1:length(biparts2[,1]), function(x) identical(cladeName, biparts2[x,2][[1]]))
    altCladeName <- sort(setdiff(tipLabs, cladeName))
    orIsThere <- sapply(1:length(biparts2[,1]), function(x) identical(altCladeName, biparts2[x,2][[1]]))
    if(any(isThere)){
      biparts2$notSeen[isThere] <- 0
      matchProbs[counter,] <- c(biparts1[,1][clade], biparts2[,1][isThere])
    } else if (any(orIsThere)) {
      biparts2$notSeen[orIsThere] <- 0
      matchProbs[counter,] <- c(biparts1[,1][clade], biparts2[,1][orIsThere])
    } else {
      matchProbs[counter,] <- c(biparts1[,1][[clade]], 0)
    }
    counter <- counter + 1
  }
  if(sum(biparts2$notSeen) > 0){
    matchProbs[counter:(counter+sum(biparts2$notSeen)-1),] <- cbind(0, biparts2[biparts2$notSeen == 1,1])
  }
  matchProbs <- matchProbs[rowSums(matchProbs == c(0,0)) != 2,]
  matchProbs
}
clade_prob <- function(tipNames, trees, allTipNames = NA, partfreqs = NA){
  if(is.na(partfreqs)[1]){
    clade_monophyly <- prop.part.df(trees, cutoff = 0.00001)
    if(is.na(allTipNames)[1]){
      allTipNames <- trees[[1]]$tip.label
    }
  } else {
    clade_monophyly <- partfreqs
    if(is.na(allTipNames)[1]){
      allTipNames <- unique(unlist(bipart_probs$cladeNames))
    }
  }
  cladeName <- tipNames
  clade <- list(allTipNames[startsWith2(allTipNames, cladeName)], allTipNames[!startsWith2(allTipNames, cladeName)])
  clade_ind <- (sapply(1:length(clade_monophyly$cladeNames), function(x) setequal(clade_monophyly$cladeNames[[x]], clade[[1]]) |
                         setequal(clade_monophyly$cladeNames[[x]], clade[[2]])))
  if(all(!clade_ind)){return(0)}
  # clade_monophyly$cladeNames[clade_ind]
  return(clade_monophyly$postProbs[clade_ind])
}
internal_nodes_probs <- function(tree, trees_to_use){
  bipart_probs <- prop.part.df(trees_to_use, cutoff = 0.00001)
  ti <- 1:length(tree$tip.label)
  tips <- sapply(1:Nnode(tree), function(node) as.numeric(na.omit(ti[getDescendants(tree, length(tree$tip.label) + node)])))
  tips <- sapply(1:length(tips), function(node) tree$tip.label[tips[[node]]])
  node_pps <- sapply(1:length(tips), function(node) clade_prob(tipNames = tips[[node]], allTipNames = tree$tip.label, partfreqs = bipart_probs))
  return(node_pps)
}

mol_tree <- read.tree("/Volumes/1TB/Harvati_Empirical/data2_neanF/mol_tree_mcc_spOnly.nex")
tiplabs <- mol_tree$tip.label

#read in and process data for mvBM analysis
dists <- rep(0, 100)
compTrees <- lapply(1:100, function(na) NA)
postDists <- matrix(NA, nrow = 100, ncol = 10002)
filename <- "bothSexes_justSpecies/harvati_PCA_bothSexes_justSpecies_replicate_"
for(i in 1:100){
  cat(paste0(i, " "))
  trees1 <- read.tree(paste0("/Volumes/1TB/Harvati/output/", filename, i, "_c1.trees"))
  trees2 <- read.tree(paste0("/Volumes/1TB/Harvati/output/", filename, i, "_c2.trees"))
  # print(RF.dist(maxCladeCred(trees1), maxCladeCred(trees2)))
  trees <- c(trees1, trees2)
  morph_tree <- maxCladeCred(trees)
  dists[i] <- RF.dist(mol_tree, morph_tree, normalize = T)
  compTrees[[i]] <- compareTrees(prop.part.df(trees), prop.part.df(mol_tree), tiplabs)
  postDists[i,] <- sapply(1:length(trees), function(tree) round(RF.dist(trees[tree], mol_tree, normalize = T), 2))
}

sum(dists >= 0.5) / 100
dists_table <- table(dists)



mvBM_trees1 <- read.tree("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_bothSexes_noHomopops_justSpecies_inf_mvBM_PCA99_c1.trees")
mvBM_trees2 <- read.tree("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_bothSexes_noHomopops_justSpecies_inf_mvBM_PCA99_c2.trees")
mvBM_trees <- c(mvBM_trees1, mvBM_trees2)
comparison_mvBM <- compareTrees(prop.part.df(mvBM_trees), prop.part.df(mol_tree), tiplabs)
mvBM_RFdists <- sapply(1:length(mvBM_trees), function(tree) round(RF.dist(mvBM_trees[tree], mol_tree_mcc_spOnly, normalize = T), 2))
mvBM_table <- table(mvBM_RFdists) / length(mvBM_trees) * 100
compTrees_comb <- do.call(rbind, compTrees)
postDists_tables <- lapply(1:100, function(pd) table(postDists[pd,])  / sum(table(postDists[pd,])) * 100)
















png(filename = "~/Documents/Harvati_Reanalysis_Manuscript/figures/figure3_final.png", width = 2400, height = 800)
par(mfrow = c(1, 3))

#histogram of mcc rf-dists
par(mar=c(10,11,8,2))
plot(dists_table, xlim = c(0,1), xlab = "", ylab = "", ylim = c(0,35), 
     cex.lab = 3.5, cex.axis = 3, lwd = 4.5, xaxt = "n", yaxt = "n", main = "MCC RF-Distances from True Tree", cex.main = 4.25)
lines(dists_table[as.numeric(attr(dists_table, "dimnames")$dists) >= 0.5], col = "darkred", lwd = 4.5)

title(ylab = "percent", cex.lab = 5, line = 7)
title(xlab = "normalized RF-Distance", cex.lab = 5, line = 7.5)
axis(1, at = 1:7/10, labels = rep("", 7), lwd = 4, cex.axis = 4, tck = -0.015)
mtext(text = 1:7/10, side = 1, at = 1:7/10, cex = 3, line = 3.25)
axis(2, at = 0:3*10, labels =  rep("", 4), lwd = 4, cex.axis = 4, tck = -0.015)
mtext(text = 0:3*10, side = 2, at = 0:3*10, cex = 3, line = 2)
box(lwd=4)
legend(x = "topright", legend = c(paste0("µ = ", round(mean(dists), 3), 
                                       ", σ = ", round(sd(dists), 3)), "       (% ≥ 0.5) = 23%"), cex = 3, bty = "n", text.col=c(1, "darkred"))
text(labels = "a)", cex = textsiz, x = 1.0365, y = 39.5, xpd = T, font = 4, col = "darkred")
box(which = "figure", lty = 5)

#comparetrees and histogram
plot(compTrees[[1]][,1], compTrees[[1]][,2] + 1/3, main = "Compare-Trees Plot and Histogram", xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
     cex.axis = 3, cex.main = 4.5, pch = 16, cex = 5, col = rgb(0, 0, 0, 0.05), xlim = c(0,1), ylim = c(0,1+1/3))
jitterPts <- T
for(i in 2:100){
  points(compTrees[[i]][,1], compTrees[[i]][,2] + 1/3 + 
           rexp(length(compTrees[[i]][,2]), rate = 50) * (1-compTrees[[i]][,2]) -
                    rexp(length(compTrees[[i]][,2]), rate = 50) * compTrees[[i]][,2], 
         col = rgb(0, 0, 0, 0.05), cex = 5, pch = 16)
}
points(comparison_mvBM[,1], comparison_mvBM[,2] + 1/3, col = rgb(139/255, 0, 0, 0.75), cex = 5, pch = 1, lwd = 3)
box(lwd=3.25)
title(xlab = "mvBM Bayesian output", cex.lab = 5, line = 7.5)
title(ylab = "true, data-generating tree", cex.lab = 5, line = 7)
axis(1, at = 0:5/5, labels = rep("", 6), lwd = 4, cex.axis = 3, tck = -0.015)
mtext(text = 0:5/5, side = 1, at = 0:5/5, cex = 3, line = 3)
axis(2, at = 0:5/5 + 1/3, labels =  rep("", 6), lwd = 4, cex.axis = 3, tck = -0.015)
mtext(text = 0:5/5, side = 2, at = 0:5/5 + 1/3, cex = 3, line = 1.875)
segments(x0 = 0, x1 = 1, y0 = 1/3, y1 = 4/3, lwd  = 3, lty = 2)

#hists
par(xpd=F)
abline(h = 0.25, lty = 3, lwd = 3)
axis(2, at = 0:4/5/4, labels =  rep("", 5), lwd = 4, cex.axis = 3, tck = -0.0075)
mtext(text = 0:4/4, side = 2, at = 0:4/5/4, cex = 1.75, line = 1, las=2)

highprobs <- discretizeContData(compTrees_comb[compTrees_comb[,2] > 0.5,1])
lowprobs <- discretizeContData(compTrees_comb[compTrees_comb[,2] < 0.5,1])
width = 0.1
for(i in 1:ncol(highprobs)){
  rect(xleft = highprobs[2,i] - width/2, xright = highprobs[2,i] + width/2, ybottom = 0, ytop = highprobs[1,i] / 5,
       col = rgb(0,0,0,0.5))
}
for(i in 1:ncol(lowprobs)){
  rect(xleft = lowprobs[2,i] - width/2, xright = lowprobs[2,i] + width/2, ybottom = 0, ytop = lowprobs[1,i] / 5,
       col = rgb(red = 1,1,1,0.5))
}


legend(x = c(0.325, 0.675), y = c(0.24, 0.1), 
       legend = c("bipartition absent", "bipartition present"), 
       fill = c(rgb(red = 1,1,1,0.5), rgb(0,0,0,0.5)), cex = 2.5)

text(labels = "b)", cex = textsiz, x = 1.0365, y = 1.51, xpd = T, font = 4, col = "darkred")
box(lwd = 3.25)
box(which = "figure", lty = 5)


highprobs <- discretizeContData(comparison_mvBM[comparison_mvBM[,2] > 0.5,1])
lowprobs <- discretizeContData(comparison_mvBM[comparison_mvBM[,2] < 0.5,1])
width = 0.05
for(i in 1:ncol(highprobs)){
  rect(xleft = highprobs[2,i] - width/2, xright = highprobs[2,i] + width/2, ybottom = 0, ytop = highprobs[1,i] / 5,
       col = rgb(0.8,0,0,0.65))
}
for(i in 1:ncol(lowprobs)){
  rect(xleft = lowprobs[2,i] - width/2, xright = lowprobs[2,i] + width/2, ybottom = 0, ytop = lowprobs[1,i] / 5,
       col = rgb(red = 0.4,0,0,0.25))
}


legend(x = c(0.325, 0.675), y = c(0.24, 0.085), 
       legend = c("bipartition absent", "bipartition present"), 
       fill = c(rgb(red = 1,1,1,0.5), rgb(0,0,0,0.5)), cex = 2.5)

text(labels = "b)", cex = textsiz, x = 1.0365, y = 1.51, xpd = T, font = 4, col = "darkred")
box(which = "figure", lty = 5)


#RF-dist of entire posterior




plot(postDists_tables[[1]], xlim = c(0,1), xlab = "", ylab = "", ylim = c(0,60), col = rgb(0,0,0,0.05),
     cex.lab = 4, cex.axis = 3, lwd = 4.5, xaxt = "n", yaxt = "n", main = "RF-Distances from True Tree", cex.main = 4.5)
for(i in 2:100){
  lines(postDists_tables[[i]], col = rgb(0,0,0,0.05), lwd = 4.5)
}
title(xlab = "normalized RF-Distance", cex.lab = 5, line = 7.5)
title(ylab = "percent", cex.lab = 5, line = 7)
axis(1, at = 0:10/10, labels = rep("", 11), lwd = 4, cex.axis = 3, tck = -0.015)
mtext(text = 0:10/10, side = 1, at = 0:10/10, cex = 3, line = 3)
axis(2, at = 0:6*10, labels =  rep("", 7), lwd = 4, cex.axis = 3, tck = -0.015)
mtext(text = 0:6*10, side = 2, at = 0:6*10, cex = 3, line = 1.875)

lines(table(mvBM_RFdists + 0.01) / length(mvBM_trees) * 100, lwd = 4.5, col = rgb(0.75,0,0,0.75))
legend(x = "topright", legend = "empirical result", lwd = 4.5, col = rgb(0.75,0,0,0.75),  cex = 3, box.lty = 3, box.lwd = 2)
box(lwd=3.25)
box(which = "figure", lty = 5)
text(labels = "c)", cex = textsiz, x = 1.0365, y = 67.75, xpd = T, font = 4, col = "darkred")



dev.off()






