library(geomorph)
library(Matrix)
library(phytools)
library(phangorn)

normalize <- T
noPCA <- T
PCA_covBG <- F
collapseHomo <- T
collapseSubsp <- F


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
  if(is.na(partfreqs)){
    clade_monophyly <- prop.part.df(trees, cutoff = 0.00001)
    if(is.na(allTipNames)){
      allTipNames <- trees[[1]]$tip.label
    }
  } else {
    clade_monophyly <- partfreqs
    if(is.na(allTipNames)){
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
  node_pps <- sapply(1:length(tips), function(node) clade_prob(tips[[node]], allTipNames = tree$tip.label, partfreqs = bipart_probs))
  return(node_pps)
}


data <- readland.nts("/Users/nikolai/data/allPD4.nts")
class(data)
str(data) #landmarks on the rows, dimensions (xyz) on the columns, species on the slices

#mean center all landmarks? landmark by landmark or vs the overall mean?

codes <- attributes(data)$dimnames[[3]]
codes_linnaen <- substr(codes, 1, 5); codes_linnaen <- sort(unique(codes_linnaen))
codes_linnaen <- cbind(codes_linnaen, sapply(1:length(codes_linnaen), function(x) sum(startsWith(codes, paste0(codes_linnaen[x], "M")))),
                       sapply(1:length(codes_linnaen), function(x) sum(startsWith(codes, paste0(codes_linnaen[x], "F")))))
colnames(codes_linnaen) <- c("code", "males", "females")
taxa <- read.csv("/Users/nikolai/data/Harvati2004_Taxa.csv", header = T)
sexMatch <- sapply(1:length(codes_linnaen[,1]), function(y) 
  which(sapply(1:length(taxa[,1]), function(x) 
    all.equal(as.numeric(codes_linnaen[y,2:3]), as.numeric(taxa[x,c("Males", "Females")])) == T))[1])
sexMatch <- unlist(replace(sexMatch, !sapply(sexMatch,length),NA))
codes_linnaen <- cbind(codes_linnaen, as.character(taxa[sexMatch,1]))
codes_linnaen[startsWith(codes_linnaen[,1], "PPPXX"),4] <- "Papio hamadryas papio" #process of elimination!
codes_linnaen[startsWith(codes_linnaen[,1], "HMN"),4] <- "Homo neanderthalensis"
codes_linnaen[startsWith(codes_linnaen[,1], "HMSAM"),4] <- "Homo sapiens (UP)"
codes_linnaen[startsWith(codes_linnaen[,1], "HMSUP"),4] <- "Homo sapiens (UP)"
codes_linnaen[startsWith(codes_linnaen[,1], "HMSAU"),4] <- "Homo sapiens (Australasian)" #right
codes_linnaen[startsWith(codes_linnaen[,1], "HMSTL"),4] <- "Homo sapiens (Australasian)" #right
codes_linnaen[startsWith(codes_linnaen[,1], "HMSEU"),4] <- "Homo sapiens (Eurasian)" #right
codes_linnaen[startsWith(codes_linnaen[,1], "HMSBG"),4] <- "Homo sapiens (Eurasian)" #right
codes_linnaen[startsWith(codes_linnaen[,1], "HMSDG"),4] <- "Homo sapiens (African)"
codes_linnaen[startsWith(codes_linnaen[,1], "HMSES"),4] <- "Homo sapiens (Greenland)"
codes_linnaen[startsWith(codes_linnaen[,1], "HMSSN"),4] <- "Homo sapiens (African)"
codes_linnaen[startsWith(codes_linnaen[,1], "HMS") & is.na(codes_linnaen[,4]),4] <- paste0(
  "Homo sapiens (Unknown ", 1:length(codes_linnaen[startsWith(codes_linnaen[,1], "HMS") & is.na(codes_linnaen[,4]),4]), ")")
if(collapseHomo){
  codes_linnaen[startsWith(codes_linnaen[,1], "HMS"), 4] <- "Homo sapiens"
}
subsp <- unique(codes_linnaen[,4])

# springer_tree <- read.newick("/Users/nikolai/Documents/Springer2012_Primate_Tree.txt")
# 
# library(rethinking)
# data("Primates301_distance_matrix")
# rethinking_tree <- nj(Primates301_distance_matrix)
# root(rethinking_tree, outgroup = "Galeopterus variegatus")
# plot(rethinking_tree)
tenk_trees <- read.nexus("/Users/nikolai/Documents/TreeBlock_10kTrees_Primates_Version3.nex")
tenk_tree <- maxCladeCred(tenk_trees)

# Text S2. Timetrees and RAxML phylograms in Newick format.
# 1. Newick timetree with autocorrelated rates and hard-bounded constraints. One unit = 100 million years. 
# 2. Newick timetree with autocorrelated rates and soft-bounded constraints. One unit = 100 million years. 
# 3. Newick timetree with independent rates and hard-bounded constraints. One unit = 100 million years. 
# 4. Newick timetree with independent rates and soft-bounded constraints. One unit = 100 million years. 
# 5. RAxML phylogram based on the 61199 bp concatenation of 69 nuclear and ten mitochondrial genes. 
# 6. RAxML phylogram based on ten mitochondrial genes.
# 7. RAxML phylogram based on 69 nuclear genes. 

startsWith2 <- function(x, pre){apply(sapply(1:length(pre), function(n) startsWith(x, pre[n])), 1, any)}
# mol_tree <- springer_tree[[5]]
mol_tree <- tenk_tree
# kinda baboon "Papio_hamadryas_kindae" is sister to "Papio_hamadryas_cynocephalus": https://www.sciencedirect.com/science/article/abs/pii/S004724841830232X?via%3Dihub
#so let's add it in

mol_tree <- add.tips(tree = mol_tree, tips = "Papio_kindae", where = "Papio_cynocephalus")
plot(mol_tree)

# trees1 <- read.tree("/Users/nikolai/output/totalCovPCA_fixRM_normalized_collapseHomo_c2.trees")
# trees1 <- read.tree("/Users/nikolai/output/totalCovPCA_uninfRM_c1.trees")

#no PCA analyses
trees1 <- read.tree("output/harvati_noPCA_c1.trees")
trees2 <- read.tree("output/harvati_noPCA_c2.trees")
trees1hp <- read.tree("output/harvati_noPCA_homopops_c1.trees")
trees2hp <- read.tree("output/harvati_noPCA_homopops_c2.trees")
njTree <- read.tree(file = paste0("output/empirical_neighborJoining.txt"))
upgmaTree <- read.tree(file = paste0("output/empirical_upgma.txt"))
mp_tree <- read.tree(file = "maximumParsimonyTree_empirical.txt")
mp_tree_PCA <- read.tree(file = "maximumParsimonyTree_empirical_PCA.txt")
mkTrees_RAW_nhp <- c(read.tree("output/harvati2004_mkModel_nohomopops_raw_c1.trees"), read.tree("output/harvati2004_mkModel_nohomopops_raw_c2.trees"))
mkTrees_PCA_nhp <- c(read.tree("output/harvati2004_mkModel_nohomopops_PCA_c1.trees"), read.tree("output/harvati2004_mkModel_nohomopops_PCA_c2.trees"))


trees <- c(trees1, trees2)
trees_hp <- c(trees1hp, trees2hp)
hptree <- maxCladeCred(trees_hp)
spOnly <- maxCladeCred(trees)
outg <- as.character(unlist(sapply(c("Macaca", "Mandrillus", "Papio"), function(genus) spOnly$tip.label[grepl(genus, spOnly$tip.label)])))
spOnly <- reroot(tree = spOnly, node.number = getMRCA(spOnly, outg), position = 
                   spOnly$edge.length[which(spOnly$edge[,2] ==  getMRCA(spOnly, outg))] / 2); plot(spOnly)
hptree <- reroot(tree = hptree, node.number = getMRCA(hptree, outg), position = 
                   hptree$edge.length[which(hptree$edge[,2] ==  getMRCA(hptree, outg))] / 2); plot(hptree)
newtiplabels_moltree <- setdiff(spOnly$tip.label, mol_tree$tip.label)
oldtiplabels_moltree <- c(sapply(1:8, function(x) paste0(strsplit(newtiplabels_moltree[x], "_")[[1]][c(1,3)], collapse = "_")))
oldtiplabels_moltree[oldtiplabels_moltree == "Homo_NA"] <- "Homo_sapiens_neanderthalensis"
oldtipindices_moltree <- (sapply(1:length(oldtiplabels_moltree), function(tiplab) which(mol_tree$tip.label == oldtiplabels_moltree[tiplab])))
mol_tree$tip.label <- (replace(mol_tree$tip.label, list = oldtipindices_moltree, values = newtiplabels_moltree))

intersect(gsub(x = mol_tree$tip.label, pattern = "_", " "), subsp)
setdiff(subsp, gsub(x = mol_tree$tip.label, pattern = "_", " "))

mol_tree <- keep.tip(phy = mol_tree, tip = spOnly$tip.label)
plot(mol_tree)

dev.off()

###############################
### FIGURE 1: Cophylo Plots ###
###############################

#let's get all the 5k trees in the right format
tenk_trees_filtered <- rep(mol_tree, length(tenk_trees)+2)[-1]
for(tree in 1:length(tenk_trees)){
  if(tree %% 100 == 0){print(tree)}
  tr <- tenk_trees[[tree]]  
  tr <- add.tips(tree = tr, tips = "Papio_kindae", where = "Papio_cynocephalus")
  newtiplabels_moltree <- setdiff(spOnly$tip.label, tr$tip.label)
  oldtiplabels_moltree <- c(sapply(1:8, function(x) paste0(strsplit(newtiplabels_moltree[x], "_")[[1]][c(1,3)], collapse = "_")), "Homo_sapiens_neanderthalensis")
  oldtiplabels_moltree <- oldtiplabels_moltree[oldtiplabels_moltree != "Homo_NA"]
  oldtipindices_moltree <- (sapply(1:length(oldtiplabels_moltree), function(tiplab) which(tr$tip.label == oldtiplabels_moltree[tiplab])))
  tr$tip.label <- (replace(tr$tip.label, list = oldtipindices_moltree, values = newtiplabels_moltree))
  tenk_trees_filtered[[tree]] <- unroot(keep.tip(phy = tr, tip = spOnly$tip.label))
}

# 
# ppMorph <- prop.clades(spOnly, trees, rooted = F)/length(trees)
# # spOnly$node.label <- c(rep(1, Ntip(spOnly) + 1), round(ppMorph,2))
# spOnly$node.label <- round(ppMorph,2)
# plot(spOnly); nodelabels()
# 
# spOnly$node.label <- sapply(1:spOnly$Nnode, function(node) clade_prob(spOnly$tip.label[
#   getDescendants(spOnly, node + length(spOnly$tip.label))[getDescendants(spOnly, node + length(spOnly$tip.label)) <= length(spOnly$tip.label)]], trees = trees))

# ppMol <- prop.clades(mol_tree, trees, rooted = F)/length(trees)
# ppMol[is.na(ppMol)] <- 0
# mol_tree$node.label <- c(rep(1, Ntip(mol_tree) + 1), round(ppMol,2))
# mol_tree$node.label <-round(ppMol,2)
# 
# 
# mol_tree$node.label <- sapply(1:mol_tree$Nnode, function(node) clade_prob(mol_tree$tip.label[
#   getDescendants(mol_tree, node + length(mol_tree$tip.label))[getDescendants(mol_tree, node + length(mol_tree$tip.label)) <= length(mol_tree$tip.label)]], 
#   trees = trees))

clade_prob(c("Pan", "Homo", "Gorilla"), trees)

colfunc <- colorRampPalette(c("white", "black"))

# spOnly$node.label <- 1:20/100
# mol_tree$node.label <- 1:20/100

coph_plot <- cophylo(mol_tree, spOnly, rotate = T)
assoc <- cbind(mol_tree$tip.label[mol_tree$tip.label != "Homo_sapiens"], mol_tree$tip.label[mol_tree$tip.label != "Homo_sapiens"])
assoc <- rbind(assoc, cbind(rep("Homo_sapiens", length(hptree$tip.label[grepl("Homo_sap", hptree$tip.label)])), hptree$tip.label[grepl("Homo_sap", hptree$tip.label)]))
coph_plot_hp <- cophylo(mol_tree, hptree, assoc = assoc, rotate = T)


#vs the mvBM tree
plot_pop_tree <- T
for(i in 1:1){
png(filename = "Documents/Harvati_Reanalysis_Manuscript/figures/figure1.png", width = 1600, height = 800)
  if(plot_pop_tree){
    plot.cophylo(coph_plot_hp, mar=c(5,0,2,0), lwd = 4,  fsize = c(2,1.75), link.lwd=4)
  } else {
    plot.cophylo(coph_plot, mar=c(5,0,2,0), lwd = 4,  fsize = 2, link.lwd=4)
  }
  # df(); plot.cophylo(coph_plot)
  par(xpd=TRUE)
  
  if(!plot_pop_tree){
    col_left <- as.integer(internal_nodes_probs(coph_plot[[1]][[1]], trees)*100)
  } else {
    col_left <- as.integer(internal_nodes_probs(coph_plot_hp[[1]][[1]], trees)*100)
  }
  
  col_left[1] <- 100
  col_left[col_left == 0] <- 1
  col_left <- colfunc(100)[col_left]
  col_left[is.na(col_left)] <- colfunc(100)[1]
  nodelabels.cophylo(pch=21, frame="none", bg=col_left, cex=2.5, which = "left")

  if(!plot_pop_tree){
    col_right <- as.integer(internal_nodes_probs(coph_plot[[1]][[2]], trees)*100)
  } else {
    col_right <- as.integer(internal_nodes_probs(coph_plot_hp[[1]][[2]], trees_hp)*100)
  }
  col_right[1] <- 100
  col_right[col_right == 0] <- 1
  col_right <- colfunc(100)[col_right]
  col_right[is.na(col_right)] <- colfunc(100)[1]
  nodelabels.cophylo(pch=21, frame="none", bg=col_right, cex=2.5, which = "right")
  
  title("       Molecular                                          Morphological", cex.main = 3, line = -0.1)
  title(paste0("RF-Distance = ", RF.dist(mol_tree, spOnly)), cex.main = 1.75, line = -1.2)
  
  xl <- -0.3; yb <- -0.115; xr <- 0.3; yt <- -0.065; 
  rect(
    head(seq(xl,xr,(xr-xl)/100),-1),
    yb,
    tail(seq(xl,xr,(xr-xl)/100),-1),
    yt,
    col=colfunc(100)
  )
  text(labels = seq(from = 0, to = 1, length.out = 11), y = yb-(yt-yb)/2.5, x = seq(xl,xr,length.out = 11), las=2, cex=1.5)
  text(labels = "posterior probability", x = (xl+xr)/2, y = yt+(yt-yb)/2.5, font = 2, cex = 2)

dev.off()
}
clade_prob(c("Macaca_fa", "Macaca_mu"), trees)
clade_prob(c("Homo", "Pan"), trees)
clade_prob(c("Pan"), trees)

RF.dist(mol_tree, spOnly)

#vs the mk tree
spOnly <- maxCladeCred(mkTrees_RAW_nhp)
outg <- as.character(unlist(sapply(c("Macaca", "Mandrillus", "Papio"), function(genus) spOnly$tip.label[grepl(genus, spOnly$tip.label)])))
spOnly <- reroot(tree = spOnly, node.number = getMRCA(spOnly, outg), position = 
                   spOnly$edge.length[which(spOnly$edge[,2] ==  getMRCA(spOnly, outg))] / 2); plot(spOnly)
hptree <- reroot(tree = mkTrees_RAW_hptree, node.number = getMRCA(hptree, outg), position = 
                   hptree$edge.length[which(hptree$edge[,2] ==  getMRCA(hptree, outg))] / 2); plot(hptree)
coph_plot <- cophylo(mol_tree, spOnly, rotate = T)
assoc <- cbind(mol_tree$tip.label[mol_tree$tip.label != "Homo_sapiens"], mol_tree$tip.label[mol_tree$tip.label != "Homo_sapiens"])
assoc <- rbind(assoc, cbind(rep("Homo_sapiens", length(hptree$tip.label[grepl("Homo_sap", hptree$tip.label)])), hptree$tip.label[grepl("Homo_sap", hptree$tip.label)]))
coph_plot_hp <- cophylo(mol_tree, hptree, assoc = assoc, rotate = T)
plot_pop_tree <- F
for(i in 1:1){
  png(filename = "Documents/Harvati_Reanalysis_Manuscript/figures/figure1_mk_raw.png", width = 1600, height = 800)
  if(plot_pop_tree){
    plot.cophylo(coph_plot_hp, mar=c(5,0,2,0), lwd = 4,  fsize = c(2,1.75), link.lwd=4)
  } else {
    plot.cophylo(coph_plot, mar=c(5,0,2,0), lwd = 4,  fsize = 2, link.lwd=4)
  }
  # df(); plot.cophylo(coph_plot)
  par(xpd=TRUE)
  
  if(!plot_pop_tree){
    col_left <- as.integer(internal_nodes_probs(coph_plot[[1]][[1]], mkTrees_PCA_nhp)*100)
  } else {
    col_left <- as.integer(internal_nodes_probs(coph_plot_hp[[1]][[1]], mkTrees_PCA_hp)*100)
  }
  
  col_left[1] <- 100
  col_left[col_left == 0] <- 1
  col_left <- colfunc(100)[col_left]
  col_left[is.na(col_left)] <- colfunc(100)[1]
  nodelabels.cophylo(pch=21, frame="none", bg=col_left, cex=2.5, which = "left")
  
  if(!plot_pop_tree){
    col_right <- as.integer(internal_nodes_probs(coph_plot[[1]][[2]], mkTrees_PCA_nhp)*100)
  } else {
    col_right <- as.integer(internal_nodes_probs(coph_plot_hp[[1]][[2]], mkTrees_PCA_hp)*100)
  }
  col_right[1] <- 100
  col_right[col_right == 0] <- 1
  col_right <- colfunc(100)[col_right]
  col_right[is.na(col_right)] <- colfunc(100)[1]
  nodelabels.cophylo(pch=21, frame="none", bg=col_right, cex=2.5, which = "right")
  
  title("       Molecular                                          Morphological", cex.main = 3, line = -0.1)
  title(paste0("RF-Distance = ", RF.dist(mol_tree, spOnly)), cex.main = 1.75, line = -1.2)
  
  xl <- -0.3; yb <- -0.115; xr <- 0.3; yt <- -0.065; 
  rect(
    head(seq(xl,xr,(xr-xl)/100),-1),
    yb,
    tail(seq(xl,xr,(xr-xl)/100),-1),
    yt,
    col=colfunc(100)
  )
  text(labels = seq(from = 0, to = 1, length.out = 11), y = yb-(yt-yb)/2.5, x = seq(xl,xr,length.out = 11), las=2, cex=1.5)
  text(labels = "posterior probability", x = (xl+xr)/2, y = yt+(yt-yb)/2.5, font = 2, cex = 2)
  
  dev.off()
}


#######################################
## cophylo plots for manuscript figure
#######################################

#compute mcc trees

dev.off()





comparison <- compareTrees(prop.part.df(tenk_trees_filtered), prop.part.df(trees), spOnly$tip.label)
plot(comparison, xlab = "molecular posterior probs", ylab = "morphological posterior probs")
hist(comparison[comparison[,1] > 0.5,2], breaks = 100, col = rgb(0,green = 1,0,0.5), ylim = c(0,25), 
     main = "posterior probability comparison", xlab = "morphological prob")
hist(comparison[comparison[,1] < 0.5,2], breaks = 100, add = T, col = rgb(red = 1,0,0,0.5))
legend(x = "topright", legend = c("molecular probability < 0.5", "molecular probability > 0.5"), fill = c("red", "green"))
