library(geomorph)
library(Matrix)
library(phytools)
library(phangorn)

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
  if(all(is.na(partfreqs))){
    clade_monophyly <- prop.part.df(trees, cutoff = 0.00001)
    if(is.na(allTipNames)){
      allTipNames <- trees[[1]]$tip.label
    }
  } else {
    clade_monophyly <- partfreqs
    if(all(is.na(allTipNames))){
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
startsWith2 <- function(x, pre){apply(sapply(1:length(pre), function(n) startsWith(x, pre[n])), 1, any)}

#using the eigenvalue pseudocorrelation correction & a NJ tree
trees1 <- read.tree("/Volumes/macOS/Users/nikolai/output/infPriors_fixRM_Bailey_PCA99_c1.trees")
trees2 <- read.tree("/Volumes/macOS/Users/nikolai/output/infPriors_fixRM_Bailey_PCA99_c2.trees")

trees1 <- read.tree("/Volumes/macOS/Users/nikolai/output/infPriors_fixCorrs_Bailey_PCA99_c1.trees")
trees2 <- read.tree("/Volumes/macOS/Users/nikolai/output/infPriors_fixCorrs_Bailey_PCA99_c2.trees")

trees1 <- read.tree("/Volumes/macOS/Users/nikolai/output/infPriors_Bailey_PCA99_c1.trees")
trees2 <- read.tree("/Volumes/macOS/Users/nikolai/output/infPriors_Bailey_PCA99_c2.trees")

#using the nearPD correction & a star tree
trees1 <- read.tree("/Volumes/macOS/Users/nikolai/output/infPriors_fixRM_Bailey_nearPD_STAR_PCA95_c1.trees")
trees2 <- read.tree("/Volumes/macOS/Users/nikolai/output/infPriors_fixRM_Bailey_nearPD_STAR_PCA95_c2.trees")

trees1 <- read.tree("/Volumes/macOS/Users/nikolai/output/infPriors_fixCorrs_Bailey_nearPD_STAR_PCA95_c1.trees")
trees2 <- read.tree("/Volumes/macOS/Users/nikolai/output/infPriors_fixCorrs_Bailey_nearPD_STAR_PCA95_c2.trees")

trees1 <- read.tree("/Volumes/macOS/Users/nikolai/output/infPriors_Bailey_nearPD_STAR_PCA95_c1.trees")
trees2 <- read.tree("/Volumes/macOS/Users/nikolai/output/infPriors_Bailey_nearPD_STAR_PCA95_c2.trees")

#using the nearPD correction & a star tree & the corr = F option w/ 99% of the PC mass
trees1 <- read.tree("/Volumes/1TB/Bailey/output/infPriors_fixRM_Bailey_nearPD_STAR_nearPDnc_PCA99_c1.trees")
trees2 <- read.tree("/Volumes/1TB/Bailey/output/infPriors_fixRM_Bailey_nearPD_STAR_nearPDnc_PCA99_c2.trees")

#using the nearPD correction, star tree, & corr = T... but only using 50 eigenvectors for computational convenience

trees1 <- read.tree("/Volumes/1TB/Bailey/output/infPriors_fixRM_Bailey_nearPD_STAR_nearPDnc_50PCs_c1.trees")
trees2 <- read.tree("/Volumes/1TB/Bailey/output/infPriors_fixRM_Bailey_nearPD_STAR_nearPDnc_50PCs_c2.trees")

trees <- c(trees1, trees2)

# for(i in 1:length(trees)){
#   trees[i] <- midpoint(trees[i])
# }

par(mfrow = c(1,3))
tiplabs <- trees[[1]]$tip.label
comparison <- compareTrees(prop.part.df(trees1), prop.part.df(trees2), tiplabs)
plot(comparison, main = "same analysis comparetrees"); abline(0,1)
RF.dist(maxCladeCred(trees1), maxCladeCred(trees2))

mcc_tree <- maxCladeCred(trees) 
# gct_tree <- greedyCT(trees)

# mcc_tree <- reroot(tree = mcc_tree, node.number = which(mcc_tree$tip.label == "NEAND"), 
#                    position = mcc_tree$edge.length[which(mcc_tree$edge[,2] == which(mcc_tree$tip.label == "NEAND"))] / 2)
mcc_tree <- midpoint(mcc_tree)

# plot(gct_tree)
plot(mcc_tree, cex = 1.5)

nodelabels(round(internal_nodes_probs(mcc_tree, unroot(trees)) * 100, 0), frame = "circle", bg = "white", cex = 1.25)
