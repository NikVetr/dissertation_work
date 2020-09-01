setwd("~")

library(geomorph)
library(Matrix)
library(phytools)
library(phangorn)

normalize <- T
noPCA <- T
PCA_covBG <- F
collapseHomo <- T
collapseSubsp <- F

shorter <- function(x1, x2){
  if(length(x1) > length(x2)){
    return(x2)
  } else{
    return(x1)
  }
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
compareTrees <- function(biparts1, biparts2, tipLabs, returnCladeNames = F){
  matchProbs <- matrix(0, nrow = sum(length(biparts1[,1]), length(biparts2[,1])), ncol = ifelse(returnCladeNames, 3, 2))
  counter <- 1
  biparts2$notSeen <- 1
  for(clade in 1:length(biparts1[,2])){
    cladeName <- biparts1[,2][[clade]]
    isThere <- sapply(1:length(biparts2[,1]), function(x) identical(cladeName, biparts2[x,2][[1]]))
    altCladeName <- sort(setdiff(tipLabs, cladeName))
    shortCladeName <- paste0(sort(shorter(cladeName, altCladeName)), collapse = ", ")
    orIsThere <- sapply(1:length(biparts2[,1]), function(x) identical(altCladeName, biparts2[x,2][[1]]))
    if(any(isThere)){
      biparts2$notSeen[isThere] <- 0
      if(returnCladeNames){
        matchProbs[counter,] <- c(biparts1[,1][clade], biparts2[,1][isThere], shortCladeName)
      } else {
        matchProbs[counter,] <- c(biparts1[,1][clade], biparts2[,1][isThere])
      }
    } else if (any(orIsThere)) {
      biparts2$notSeen[orIsThere] <- 0
      if(returnCladeNames){
        matchProbs[counter,] <- c(biparts1[,1][clade], biparts2[,1][orIsThere], shortCladeName)
      } else {
        matchProbs[counter,] <- c(biparts1[,1][clade], biparts2[,1][orIsThere])
      }
    } else {
      if(returnCladeNames){
        matchProbs[counter,] <- c(biparts1[,1][[clade]], 0, shortCladeName)
      } else {
        matchProbs[counter,] <- c(biparts1[,1][[clade]], 0)
      }
    }
    counter <- counter + 1
  }
  if(sum(biparts2$notSeen) > 0){
    if(returnCladeNames){
      cladeNames <- sapply(1:length(biparts2[biparts2$notSeen == 1,2]), function(bip) 
        paste0(sort(shorter(biparts2[biparts2$notSeen == 1,2][[bip]], sort(setdiff(tipLabs, biparts2[biparts2$notSeen == 1,2][[bip]])))), collapse = ", "))
      matchProbs[counter:(counter+sum(biparts2$notSeen)-1),] <- cbind(0, cbind(biparts2[biparts2$notSeen == 1,1], cladeNames)) #needs to be an sapply
    } else {
      matchProbs[counter:(counter+sum(biparts2$notSeen)-1),] <- cbind(0, biparts2[biparts2$notSeen == 1,1])
    }
  }
  if(returnCladeNames){
    matchProbs <- matchProbs[rowSums(matchProbs == c(0,0)) != ifelse(returnCladeNames,3,2),]
  }
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
mkTrees_RAW_hp <- c(read.tree("output/harvati2004_mkModel_homopops_raw_c1.trees"), read.tree("output/harvati2004_mkModel_homopops_raw_c2.trees"))
mkTrees_PCA_hp <- c(read.tree("output/harvati2004_mkModel_homopops_PCA_c1.trees"), read.tree("output/harvati2004_mkModel_homopops_PCA_c2.trees"))

jenks6Trees_PCA_nhp <- c(read.tree("output/empirical_46traits_jenks_nohomopops_6_PCs_trees_run_1.trees"), read.tree("output/empirical_46traits_jenks_nohomopops_6_PCs_trees_run_2.trees"))
jenks6Trees_PCA_hp <- c(read.tree("output/empirical_46traits_jenks_6_PCs_trees_run_1.trees"), read.tree("output/empirical_46traits_jenks_6_PCs_trees_run_2.trees"))

jenks2Trees_PCA_nhp <- c(read.tree("output/empirical_46traits_jenks_nohomopops_2_PCs_trees_run_1.trees"), read.tree("output/empirical_46traits_jenks_nohomopops_2_PCs_trees_run_2.trees"))
jenks2Trees_PCA_hp <- c(read.tree("output/empirical_46traits_jenks_2_PCs_trees_run_1.trees"), read.tree("output/empirical_46traits_jenks_2_PCs_trees_run_2.trees"))

jenks4Trees_PCA_nhp <- c(read.tree("output/empirical_46traits_jenks_nohomopops_4_PCs_trees_run_1.trees"), read.tree("output/empirical_46traits_jenks_nohomopops_4_PCs_trees_run_2.trees"))
jenks4Trees_PCA_hp <- c(read.tree("output/empirical_46traits_jenks_4_PCs_trees_run_1.trees"), read.tree("output/empirical_46traits_jenks_4_PCs_trees_run_2.trees"))

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
age_cynocephalus_kindae_split <- 0.6
tenk_trees_filtered <- rep(spOnly, length(tenk_trees)+2)[-1]
for(tree in 1:length(tenk_trees)){
  if(tree %% 100 == 0){print(tree)}
  tr <- tenk_trees[[tree]]  
  tr <- add.tips(tree = tr, tips = "Papio_kindae", where = "Papio_cynocephalus", edge.length = 0)
  newtiplabels_moltree <- setdiff(spOnly$tip.label, tr$tip.label)
  newtiplabels_moltree <- newtiplabels_moltree[-which(newtiplabels_moltree == "Homo_neanderthalensis")]
  newtiplabels_moltree <- c(newtiplabels_moltree, "Homo_neanderthalensis")
  oldtiplabels_moltree <- c(sapply(1:8, function(x) paste0(strsplit(newtiplabels_moltree[x], "_")[[1]][c(1,3)], collapse = "_")), "Homo_sapiens_neanderthalensis")
  oldtiplabels_moltree <- oldtiplabels_moltree[oldtiplabels_moltree != "Homo_NA"]
  oldtipindices_moltree <- (sapply(1:length(oldtiplabels_moltree), function(tiplab) which(tr$tip.label == oldtiplabels_moltree[tiplab])))
  tr$tip.label <- (replace(tr$tip.label, list = oldtipindices_moltree, values = newtiplabels_moltree))
  tr <- (keep.tip(phy = tr, tip = spOnly$tip.label))
  tr$edge.length[which(tr$edge[,2] == which(tr$tip.label == "Papio_hamadryas_kindae"))] <- age_cynocephalus_kindae_split
  tr$edge.length[which(tr$edge[,2] == which(tr$tip.label == "Papio_hamadryas_cynocephalus"))] <- age_cynocephalus_kindae_split
  tr$edge.length[which(tr$edge[,2] == tr$edge[which(tr$edge[,2] == which(tr$tip.label == "Papio_hamadryas_kindae")),1])] <- tr$edge.length[which(tr$edge[,2] == tr$edge[which(tr$edge[,2] == which(tr$tip.label == "Papio_hamadryas_kindae")),1])] - age_cynocephalus_kindae_split
  tenk_trees_filtered[[tree]] <- tr
}
save(tenk_trees_filtered, file = "tenk_trees_filtered")

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
  nodelabels.cophylo(pch=21, frame="none", bg=col_left, cex=3.5, which = "left")

  if(!plot_pop_tree){
    col_right <- as.integer(internal_nodes_probs(coph_plot[[1]][[2]], trees)*100)
  } else {
    col_right <- as.integer(internal_nodes_probs(coph_plot_hp[[1]][[2]], trees_hp)*100)
  }
  col_right[1] <- 100
  col_right[col_right == 0] <- 1
  col_right <- colfunc(100)[col_right]
  col_right[is.na(col_right)] <- colfunc(100)[1]
  nodelabels.cophylo(pch=21, frame="none", bg=col_right, cex=3.5, which = "right")
  
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

#vs the mk tree -- no data transform
spOnly <- maxCladeCred(mkTrees_RAW_nhp)
mkTrees_RAW_hptree <- maxCladeCred(mkTrees_RAW_hp)
outg <- as.character(unlist(sapply(c("Macaca", "Mandrillus", "Papio"), function(genus) spOnly$tip.label[grepl(genus, spOnly$tip.label)])))
spOnly <- reroot(tree = spOnly, node.number = getMRCA(spOnly, outg), position = 
                   spOnly$edge.length[which(spOnly$edge[,2] ==  getMRCA(spOnly, outg))] / 2); plot(spOnly)
hptree <- reroot(tree = mkTrees_RAW_hptree, node.number = getMRCA(mkTrees_RAW_hptree, outg), position = 
                   mkTrees_RAW_hptree$edge.length[which(hptree$edge[,2] ==  getMRCA(mkTrees_RAW_hptree, outg))] / 2); plot(hptree)
coph_plot <- cophylo(mol_tree, spOnly, rotate = T)
assoc <- cbind(mol_tree$tip.label[mol_tree$tip.label != "Homo_sapiens"], mol_tree$tip.label[mol_tree$tip.label != "Homo_sapiens"])
assoc <- rbind(assoc, cbind(rep("Homo_sapiens", length(hptree$tip.label[grepl("Homo_sap", hptree$tip.label)])), hptree$tip.label[grepl("Homo_sap", hptree$tip.label)]))
coph_plot_hp <- cophylo(mol_tree, hptree, assoc = assoc, rotate = T)
plot_pop_tree <- T
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
    col_left <- as.integer(internal_nodes_probs(coph_plot[[1]][[1]], mkTrees_RAW_nhp)*100)
  } else {
    col_left <- as.integer(internal_nodes_probs(coph_plot_hp[[1]][[1]], mkTrees_RAW_nhp)*100)
  }
  
  col_left[1] <- 100
  col_left[col_left == 0] <- 1
  col_left <- colfunc(100)[col_left]
  col_left[is.na(col_left)] <- colfunc(100)[1]
  nodelabels.cophylo(pch=21, frame="none", bg=col_left, cex=3.5, which = "left")
  
  if(!plot_pop_tree){
    col_right <- as.integer(internal_nodes_probs(coph_plot[[1]][[2]], mkTrees_RAW_nhp)*100)
  } else {
    col_right <- as.integer(internal_nodes_probs(coph_plot_hp[[1]][[2]], mkTrees_RAW_hp)*100)
  }
  col_right[1] <- 100
  col_right[col_right == 0] <- 1
  col_right <- colfunc(100)[col_right]
  col_right[is.na(col_right)] <- colfunc(100)[1]
  nodelabels.cophylo(pch=21, frame="none", bg=col_right, cex=3.5, which = "right")
  
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

clade_prob(c("Macaca_fa", "Macaca_mu"), mkTrees_RAW_nhp)
clade_prob(c("Homo", "Pan", "Gorilla"), mkTrees_RAW_nhp)
clade_prob(c("Man", "Pap", "Mac"), mkTrees_RAW_nhp)
clade_prob(c("Gorilla"), mkTrees_RAW_nhp)

#vs the mk tree -- PCA transform
spOnly <- maxCladeCred(mkTrees_PCA_nhp)
mkTrees_PCA_hptree <- maxCladeCred(mkTrees_PCA_hp)
outg <- as.character(unlist(sapply(c("Macaca", "Mandrillus", "Papio"), function(genus) spOnly$tip.label[grepl(genus, spOnly$tip.label)])))
spOnly <- reroot(tree = spOnly, node.number = getMRCA(spOnly, outg), position = 
                   spOnly$edge.length[which(spOnly$edge[,2] ==  getMRCA(spOnly, outg))] / 2); plot(spOnly)
hptree <- reroot(tree = mkTrees_PCA_hptree, node.number = getMRCA(mkTrees_PCA_hptree, outg), position = 
                   mkTrees_PCA_hptree$edge.length[which(hptree$edge[,2] ==  getMRCA(mkTrees_PCA_hptree, outg))] / 2); plot(hptree)
coph_plot <- cophylo(mol_tree, spOnly, rotate = T)
assoc <- cbind(mol_tree$tip.label[mol_tree$tip.label != "Homo_sapiens"], mol_tree$tip.label[mol_tree$tip.label != "Homo_sapiens"])
assoc <- rbind(assoc, cbind(rep("Homo_sapiens", length(hptree$tip.label[grepl("Homo_sap", hptree$tip.label)])), hptree$tip.label[grepl("Homo_sap", hptree$tip.label)]))
coph_plot_hp <- cophylo(mol_tree, hptree, assoc = assoc, rotate = T)
plot_pop_tree <- T
for(i in 1:1){
  png(filename = "Documents/Harvati_Reanalysis_Manuscript/figures/figure1_mk_PCA.png", width = 1600, height = 800)
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
    col_left <- as.integer(internal_nodes_probs(coph_plot_hp[[1]][[1]], mkTrees_PCA_nhp)*100)
  }
  
  col_left[1] <- 100
  col_left[col_left == 0] <- 1
  col_left <- colfunc(100)[col_left]
  col_left[is.na(col_left)] <- colfunc(100)[1]
  nodelabels.cophylo(pch=21, frame="none", bg=col_left, cex=3.5, which = "left")
  
  if(!plot_pop_tree){
    col_right <- as.integer(internal_nodes_probs(coph_plot[[1]][[2]], mkTrees_PCA_nhp)*100)
  } else {
    col_right <- as.integer(internal_nodes_probs(coph_plot_hp[[1]][[2]], mkTrees_PCA_hp)*100)
  }
  col_right[1] <- 100
  col_right[col_right == 0] <- 1
  col_right <- colfunc(100)[col_right]
  col_right[is.na(col_right)] <- colfunc(100)[1]
  nodelabels.cophylo(pch=21, frame="none", bg=col_right, cex=3.5, which = "right")
  
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

clade_prob(c("Macaca_fa", "Macaca_mu"), mkTrees_PCA_nhp)
clade_prob(c("Homo", "Pan", "Gorilla"), mkTrees_PCA_nhp)
clade_prob(c("Man", "Pap", "Mac"), mkTrees_PCA_nhp)
clade_prob(c("Gorilla"), mkTrees_PCA_nhp)
clade_prob(c("Homo"), mkTrees_PCA_nhp)

#vs the JENKS-6 tree -- PCA transform
spOnly <- maxCladeCred(jenks6Trees_PCA_nhp)
jenks6Trees_PCA_hptree <- maxCladeCred(jenks6Trees_PCA_hp)
outg <- as.character(unlist(sapply(c("Macaca", "Mandrillus", "Papio"), function(genus) spOnly$tip.label[grepl(genus, spOnly$tip.label)])))
spOnly <- reroot(tree = spOnly, node.number = getMRCA(spOnly, outg), position = 
                   spOnly$edge.length[which(spOnly$edge[,2] ==  getMRCA(spOnly, outg))] / 2); plot(spOnly)
hptree <- reroot(tree = jenks6Trees_PCA_hptree, node.number = getMRCA(jenks6Trees_PCA_hptree, outg), position = 
                   jenks6Trees_PCA_hptree$edge.length[which(hptree$edge[,2] ==  getMRCA(jenks6Trees_PCA_hptree, outg))] / 2); plot(hptree)
coph_plot <- cophylo(mol_tree, spOnly, rotate = T)
assoc <- cbind(mol_tree$tip.label[mol_tree$tip.label != "Homo_sapiens"], mol_tree$tip.label[mol_tree$tip.label != "Homo_sapiens"])
assoc <- rbind(assoc, cbind(rep("Homo_sapiens", length(hptree$tip.label[grepl("Homo_sap", hptree$tip.label)])), hptree$tip.label[grepl("Homo_sap", hptree$tip.label)]))
coph_plot_hp <- cophylo(mol_tree, hptree, assoc = assoc, rotate = T)
plot_pop_tree <- T
for(i in 1:1){
  png(filename = "Documents/Harvati_Reanalysis_Manuscript/figures/figure1_jenks6_PCA.png", width = 1600, height = 800)
  if(plot_pop_tree){
    plot.cophylo(coph_plot_hp, mar=c(5,0,2,0), lwd = 4,  fsize = c(2,1.75), link.lwd=4)
  } else {
    plot.cophylo(coph_plot, mar=c(5,0,2,0), lwd = 4,  fsize = 2, link.lwd=4)
  }
  # df(); plot.cophylo(coph_plot)
  par(xpd=TRUE)
  
  if(!plot_pop_tree){
    col_left <- as.integer(internal_nodes_probs(coph_plot[[1]][[1]], jenks6Trees_PCA_nhp)*100)
  } else {
    col_left <- as.integer(internal_nodes_probs(coph_plot_hp[[1]][[1]], jenks6Trees_PCA_nhp)*100)
  }
  
  col_left[1] <- 100
  col_left[col_left == 0] <- 1
  col_left <- colfunc(100)[col_left]
  col_left[is.na(col_left)] <- colfunc(100)[1]
  nodelabels.cophylo(pch=21, frame="none", bg=col_left, cex=3.5, which = "left")
  
  if(!plot_pop_tree){
    col_right <- as.integer(internal_nodes_probs(coph_plot[[1]][[2]], jenks6Trees_PCA_nhp)*100)
  } else {
    col_right <- as.integer(internal_nodes_probs(coph_plot_hp[[1]][[2]], jenks6Trees_PCA_hp)*100)
  }
  col_right[1] <- 100
  col_right[col_right == 0] <- 1
  col_right <- colfunc(100)[col_right]
  col_right[is.na(col_right)] <- colfunc(100)[1]
  nodelabels.cophylo(pch=21, frame="none", bg=col_right, cex=3.5, which = "right")
  
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

clade_prob(c("Macaca_fa", "Macaca_mu"), jenks6Trees_PCA_nhp)
clade_prob(c("Homo", "Pan", "Gorilla"), jenks6Trees_PCA_nhp)
clade_prob(c("Man", "Pap", "Mac"), jenks6Trees_PCA_nhp)
clade_prob(c("Gorilla"), jenks6Trees_PCA_nhp)
clade_prob(c("Homo"), mkTrees_PCA_nhp)


#vs the JENKS-2 tree -- PCA transform
spOnly <- maxCladeCred(jenks2Trees_PCA_nhp)
jenks2Trees_PCA_hptree <- maxCladeCred(jenks2Trees_PCA_hp)
outg <- as.character(unlist(sapply(c("Macaca", "Mandrillus", "Papio"), function(genus) spOnly$tip.label[grepl(genus, spOnly$tip.label)])))
spOnly <- reroot(tree = spOnly, node.number = getMRCA(spOnly, outg), position = 
                   spOnly$edge.length[which(spOnly$edge[,2] ==  getMRCA(spOnly, outg))] / 2); plot(spOnly)
hptree <- reroot(tree = jenks2Trees_PCA_hptree, node.number = getMRCA(jenks2Trees_PCA_hptree, outg), position = 
                   jenks2Trees_PCA_hptree$edge.length[which(jenks2Trees_PCA_hptree$edge[,2] ==  getMRCA(jenks2Trees_PCA_hptree, outg))] / 2); plot(hptree)
coph_plot <- cophylo(mol_tree, spOnly, rotate = T)
assoc <- cbind(mol_tree$tip.label[mol_tree$tip.label != "Homo_sapiens"], mol_tree$tip.label[mol_tree$tip.label != "Homo_sapiens"])
assoc <- rbind(assoc, cbind(rep("Homo_sapiens", length(hptree$tip.label[grepl("Homo_sap", hptree$tip.label)])), hptree$tip.label[grepl("Homo_sap", hptree$tip.label)]))
coph_plot_hp <- cophylo(mol_tree, hptree, assoc = assoc, rotate = T)
plot_pop_tree <- T
for(i in 1:1){
  png(filename = "Documents/Harvati_Reanalysis_Manuscript/figures/figure1_jenks2_PCA.png", width = 1600, height = 800)
  if(plot_pop_tree){
    plot.cophylo(coph_plot_hp, mar=c(5,0,2,0), lwd = 4,  fsize = c(2,1.75), link.lwd=4)
  } else {
    plot.cophylo(coph_plot, mar=c(5,0,2,0), lwd = 4,  fsize = 2, link.lwd=4)
  }
  # df(); plot.cophylo(coph_plot)
  par(xpd=TRUE)
  
  if(!plot_pop_tree){
    col_left <- as.integer(internal_nodes_probs(coph_plot[[1]][[1]], jenks2Trees_PCA_nhp)*100)
  } else {
    col_left <- as.integer(internal_nodes_probs(coph_plot_hp[[1]][[1]], jenks2Trees_PCA_nhp)*100)
  }
  
  col_left[1] <- 100
  col_left[col_left == 0] <- 1
  col_left <- colfunc(100)[col_left]
  col_left[is.na(col_left)] <- colfunc(100)[1]
  nodelabels.cophylo(pch=21, frame="none", bg=col_left, cex=3.5, which = "left")
  
  if(!plot_pop_tree){
    col_right <- as.integer(internal_nodes_probs(coph_plot[[1]][[2]], jenks2Trees_PCA_nhp)*100)
  } else {
    col_right <- as.integer(internal_nodes_probs(coph_plot_hp[[1]][[2]], jenks2Trees_PCA_hp)*100)
  }
  col_right[1] <- 100
  col_right[col_right == 0] <- 1
  col_right <- colfunc(100)[col_right]
  col_right[is.na(col_right)] <- colfunc(100)[1]
  nodelabels.cophylo(pch=21, frame="none", bg=col_right, cex=3.5, which = "right")
  
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

#vs the JENKS-4 tree -- PCA transform
spOnly <- maxCladeCred(jenks4Trees_PCA_nhp)
jenks4Trees_PCA_hptree <- maxCladeCred(jenks4Trees_PCA_hp)
outg <- as.character(unlist(sapply(c("Macaca", "Mandrillus", "Papio"), function(genus) spOnly$tip.label[grepl(genus, spOnly$tip.label)])))
spOnly <- reroot(tree = spOnly, node.number = getMRCA(spOnly, outg), position = 
                   spOnly$edge.length[which(spOnly$edge[,2] ==  getMRCA(spOnly, outg))] / 2); plot(spOnly)
hptree <- reroot(tree = jenks4Trees_PCA_hptree, node.number = getMRCA(jenks4Trees_PCA_hptree, outg), position = 
                   jenks4Trees_PCA_hptree$edge.length[which(jenks4Trees_PCA_hptree$edge[,2] ==  getMRCA(jenks4Trees_PCA_hptree, outg))] / 2); plot(hptree)
coph_plot <- cophylo(mol_tree, spOnly, rotate = T)
assoc <- cbind(mol_tree$tip.label[mol_tree$tip.label != "Homo_sapiens"], mol_tree$tip.label[mol_tree$tip.label != "Homo_sapiens"])
assoc <- rbind(assoc, cbind(rep("Homo_sapiens", length(hptree$tip.label[grepl("Homo_sap", hptree$tip.label)])), hptree$tip.label[grepl("Homo_sap", hptree$tip.label)]))
coph_plot_hp <- cophylo(mol_tree, hptree, assoc = assoc, rotate = T)
plot_pop_tree <- T
for(i in 1:1){
  png(filename = "Documents/Harvati_Reanalysis_Manuscript/figures/figure1_jenks4_PCA.png", width = 1600, height = 800)
  if(plot_pop_tree){
    plot.cophylo(coph_plot_hp, mar=c(5,0,2,0), lwd = 4,  fsize = c(2,1.75), link.lwd=4)
  } else {
    plot.cophylo(coph_plot, mar=c(5,0,2,0), lwd = 4,  fsize = 2, link.lwd=4)
  }
  # df(); plot.cophylo(coph_plot)
  par(xpd=TRUE)
  
  if(!plot_pop_tree){
    col_left <- as.integer(internal_nodes_probs(coph_plot[[1]][[1]], jenks4Trees_PCA_nhp)*100)
  } else {
    col_left <- as.integer(internal_nodes_probs(coph_plot_hp[[1]][[1]], jenks4Trees_PCA_nhp)*100)
  }
  
  col_left[1] <- 100
  col_left[col_left == 0] <- 1
  col_left <- colfunc(100)[col_left]
  col_left[is.na(col_left)] <- colfunc(100)[1]
  nodelabels.cophylo(pch=21, frame="none", bg=col_left, cex=3.5, which = "left")
  
  if(!plot_pop_tree){
    col_right <- as.integer(internal_nodes_probs(coph_plot[[1]][[2]], jenks4Trees_PCA_nhp)*100)
  } else {
    col_right <- as.integer(internal_nodes_probs(coph_plot_hp[[1]][[2]], jenks4Trees_PCA_hp)*100)
  }
  col_right[1] <- 100
  col_right[col_right == 0] <- 1
  col_right <- colfunc(100)[col_right]
  col_right[is.na(col_right)] <- colfunc(100)[1]
  nodelabels.cophylo(pch=21, frame="none", bg=col_right, cex=3.5, which = "right")
  
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


#vs the distance trees
njTree <- read.tree(file = paste0("output/empirical_neighborJoining.txt"))
upgmaTree <- read.tree(file = paste0("output/empirical_upgma.txt"))
outg <- as.character(unlist(sapply(c("Macaca", "Mandrillus", "Papio"), function(genus) spOnly$tip.label[grepl(genus, spOnly$tip.label)])))
while(length(which(njTree$edge[,2] == getMRCA(njTree, outg))) == 0){
  njTree <- reroot(njTree, node.number = sample(size = 1, 1:njTree$Nnode))
}
hptree <- reroot(tree = njTree, node.number = getMRCA(njTree, outg), position = 
                   njTree$edge.length[which(njTree$edge[,2] ==  getMRCA(njTree, outg))] / 2); plot(hptree)
assoc <- cbind(mol_tree$tip.label[mol_tree$tip.label != "Homo_sapiens"], mol_tree$tip.label[mol_tree$tip.label != "Homo_sapiens"])
assoc <- rbind(assoc, cbind(rep("Homo_sapiens", length(hptree$tip.label[grepl("Homo_sap", hptree$tip.label)])), hptree$tip.label[grepl("Homo_sap", hptree$tip.label)]))
spOnly <- drop.tip(hptree, tip = assoc[startsWith(assoc[,2], "Homo_s"),2][-1]); 
spOnly$tip.label[spOnly$tip.label == assoc[startsWith(assoc[,2], "Homo_s"),2][1]] <- "Homo_sapiens"
coph_plot_hp <- cophylo(mol_tree, hptree, assoc = assoc, rotate = T)
plot_pop_tree <- T
for(i in 1:1){
  png(filename = "Documents/Harvati_Reanalysis_Manuscript/figures/figure1_nj_mahalanobis.png", width = 1600, height = 800)
  if(plot_pop_tree){
    plot.cophylo(coph_plot_hp, mar=c(0,0,2,0), lwd = 4,  fsize = c(2,1.75), link.lwd=4)
  } else {
    plot.cophylo(coph_plot, mar=c(0,0,2,0), lwd = 4,  fsize = 2, link.lwd=4)
  }
  # df(); plot.cophylo(coph_plot)
  par(xpd=TRUE)
  
  title("       Molecular                                          Morphological", cex.main = 3, line = -0.1)
  title(paste0("RF-Distance = ", RF.dist(mol_tree, spOnly)), cex.main = 1.75, line = -1.2)
  
  dev.off()
}
hptree <- upgmaTree
spOnly <- drop.tip(hptree, tip = assoc[startsWith(assoc[,2], "Homo_s"),2][-1]); 
spOnly$tip.label[spOnly$tip.label == assoc[startsWith(assoc[,2], "Homo_s"),2][1]] <- "Homo_sapiens"
coph_plot_hp <- cophylo(mol_tree, hptree, assoc = assoc, rotate = T)
for(i in 1:1){
  png(filename = "Documents/Harvati_Reanalysis_Manuscript/figures/figure1_upgma_mahalanobis.png", width = 1600, height = 800)
  if(plot_pop_tree){
    plot.cophylo(coph_plot_hp, mar=c(0,0,2,0), lwd = 4,  fsize = c(2,1.75), link.lwd=4)
  } else {
    plot.cophylo(coph_plot, mar=c(0,0,2,0), lwd = 4,  fsize = 2, link.lwd=4)
  }
  # df(); plot.cophylo(coph_plot)
  par(xpd=TRUE)
  
  title("       Molecular                                          Morphological", cex.main = 3, line = -0.1)
  title(paste0("RF-Distance = ", RF.dist(mol_tree, spOnly)), cex.main = 1.75, line = -1.2)
  
  dev.off()
}

#vs the parsimony trees
DC_PCA_tree <- maxCladeCred(read.tree(file = paste0("output/maximumParsimonyTree_empirical_PCA_divergenceCoding.txt")))
MVC_PCA_tree <- maxCladeCred(read.tree(file = paste0("output/maximumParsimonyTree_empirical_PCA_MVC.txt")))
DC_RAW_tree <- maxCladeCred(read.tree(file = paste0("output/maximumParsimonyTree_empirical_RAW_divergenceCoding.txt")))
MVC_RAW_tree <- maxCladeCred(read.tree(file = paste0("output/maximumParsimonyTree_empirical_RAW_MVC.txt")))
Jenks2_RAW_tree <- maxCladeCred(read.tree(file = paste0("output/maximumParsimonyTree_harvati_empirical_raw_jenks_2_expcats.txt")))
Jenks4_RAW_tree <- maxCladeCred(read.tree(file = paste0("output/maximumParsimonyTree_harvati_empirical_raw_jenks_4_expcats.txt")))
Jenks6_RAW_tree <- maxCladeCred(read.tree(file = paste0("output/maximumParsimonyTree_harvati_empirical_raw_jenks_6_expcats.txt")))
Jenks2_PCA_tree <- maxCladeCred(read.tree(file = paste0("output/maximumParsimonyTree_harvati_empirical_PCs_jenks_2_expcats.txt")))
Jenks4_PCA_tree <- maxCladeCred(read.tree(file = paste0("output/maximumParsimonyTree_harvati_empirical_PCs_jenks_4_expcats.txt")))
Jenks6_PCA_tree <- maxCladeCred(read.tree(file = paste0("output/maximumParsimonyTree_harvati_empirical_PCs_jenks_6_expcats.txt")))

parsTrees <- list(DC_PCA_tree, MVC_PCA_tree, DC_RAW_tree, MVC_RAW_tree, Jenks2_RAW_tree, Jenks4_RAW_tree, Jenks6_RAW_tree, Jenks2_PCA_tree, Jenks4_PCA_tree, Jenks6_PCA_tree)
parsTrees_names <- trimws(c("DC_PCA_tree", " MVC_PCA_tree", " DC_RAW_tree", " MVC_RAW_tree", " Jenks2_RAW_tree", " Jenks4_RAW_tree", " Jenks6_RAW_tree", " Jenks2_PCA_tree", " Jenks4_PCA_tree", " Jenks6_PCA_tree"))

#read in the stepwise parsimony distance matrices, i.e. euclidean distances on marginal traits
DC_PCA_dist <- readDist( "output/parsDistances_empirical_PCA_divergenceCoding.dist" )
MVC_PCA_dist <- readDist( "output/parsDistances_empirical_PCA_MVC.dist" )
DC_RAW_dist <- readDist( "output/parsDistances_empirical_RAW_divergenceCoding.dist" )
MVC_RAW_dist <- readDist( "output/parsDistances_empirical_RAW_MVC.dist" )
Jenks2_RAW_dist <- readDist( "output/jenksDistances_harvati_empirical_raw_jenks_2_expcats.dist")
Jenks4_RAW_dist <- readDist( "output/jenksDistances_harvati_empirical_raw_jenks_4_expcats.dist")
Jenks6_RAW_dist <- readDist( "output/jenksDistances_harvati_empirical_raw_jenks_6_expcats.dist")
Jenks2_PCA_dist <- readDist( "output/jenksDistances_harvati_empirical_PCs_jenks_2_expcats.dist")
Jenks4_PCA_dist <- readDist( "output/jenksDistances_harvati_empirical_PCs_jenks_4_expcats.dist")
Jenks6_PCA_dist <- readDist( "output/jenksDistances_harvati_empirical_PCs_jenks_6_expcats.dist")

parsTrees_dists <- list(DC_PCA_dist, MVC_PCA_dist, DC_RAW_dist, MVC_RAW_dist, Jenks2_RAW_dist, Jenks4_RAW_dist, Jenks6_RAW_dist, Jenks2_PCA_dist, Jenks4_PCA_dist, Jenks6_PCA_dist)
  
for(i in 1:length(parsTrees_dists)){
  parsTrees[[i]] <- nnls.phylo(parsTrees[[i]], parsTrees_dists[[i]])
  # parsTrees[[i]]$edge.length[parsTrees[[i]]$edge.length == 0] <- 0.1
}

for(tree_ind in 1:length(parsTrees_names)){
  tree_to_use <- parsTrees[[tree_ind]]
  tree_to_use_name <- parsTrees_names[[tree_ind]]
  tiplabs <- tree_to_use$tip.label
  outg <- as.character(unlist(sapply(c("Macaca", "Mandrillus", "Papio"), function(genus) tiplabs[grepl(genus, tiplabs)])))
  # hptree <- root(tree_to_use, outgroup = outg); plot(hptree)
  while(length(which(tree_to_use$edge[,2] == getMRCA(tree_to_use, outg))) == 0){
    tree_to_use <- reroot(tree_to_use, node.number = sample(size = 1, 1:tree_to_use$Nnode))
  }
  hptree <- reroot(tree = tree_to_use, node.number = getMRCA(tree_to_use, outg), position = 
                     tree_to_use$edge.length[which(tree_to_use$edge[,2] ==  getMRCA(tree_to_use, outg))] / 2); plot(hptree)
  assoc <- cbind(mol_tree$tip.label[mol_tree$tip.label != "Homo_sapiens"], mol_tree$tip.label[mol_tree$tip.label != "Homo_sapiens"])
  assoc <- rbind(assoc, cbind(rep("Homo_sapiens", length(hptree$tip.label[grepl("Homo_sap", hptree$tip.label)])), hptree$tip.label[grepl("Homo_sap", hptree$tip.label)]))
  spOnly <- drop.tip(hptree, tip = assoc[startsWith(assoc[,2], "Homo_s"),2][-1]); 
  spOnly$tip.label[spOnly$tip.label == assoc[startsWith(assoc[,2], "Homo_s"),2][1]] <- "Homo_sapiens"
  coph_plot_hp <- cophylo(mol_tree, hptree, assoc = assoc, rotate = T)
  plot_pop_tree <- T
  for(i in 1:1){
    png(filename = paste0("Documents/Harvati_Reanalysis_Manuscript/figures/figure1_parsimony_", tree_to_use_name, ".png"), width = 1600, height = 800)
    if(plot_pop_tree){
      plot.cophylo(coph_plot_hp, mar=c(0,0,2,0), lwd = 4,  fsize = c(2,1.75), link.lwd=4)
    } else {
      plot.cophylo(coph_plot, mar=c(0,0,2,0), lwd = 4,  fsize = 2, link.lwd=4)
    }
    # df(); plot.cophylo(coph_plot)
    par(xpd=TRUE)
    
    if(is.binary.phylo(spOnly)){
      distance_betw_trees <- RF.dist(mol_tree, spOnly)
    } else {
      distance_betw_trees <- mean(RF.dist(mol_tree, resolveAllNodes(spOnly)))
    }
    title("       Molecular                                          Morphological", cex.main = 3, line = -0.1)
    title(paste0("RF-Distance = ", distance_betw_trees), cex.main = 1.75, line = -1.2)
    
    dev.off()
  }
}

hptree <- upgmaTree
coph_plot_hp <- cophylo(mol_tree, hptree, assoc = assoc, rotate = T)
for(i in 1:1){
  png(filename = "Documents/Harvati_Reanalysis_Manuscript/figures/figure1_upgma_mahalanobis.png", width = 1600, height = 800)
  if(plot_pop_tree){
    plot.cophylo(coph_plot_hp, mar=c(0,0,2,0), lwd = 4,  fsize = c(2,1.75), link.lwd=4)
  } else {
    plot.cophylo(coph_plot, mar=c(0,0,2,0), lwd = 4,  fsize = 2, link.lwd=4)
  }
  # df(); plot.cophylo(coph_plot)
  par(xpd=TRUE)
  
  title("       Molecular                                          Morphological", cex.main = 3, line = -0.1)
  title(paste0("RF-Distance = ", RF.dist(mol_tree, spOnly)), cex.main = 1.75, line = -1.2)
  
  dev.off()
}


###########################################
## comparetrees plots for manuscript figure
###########################################

#compute mcc trees
tiplabs <- trees[[1]]$tip.label
dev.off()
mvBM_trees <- c(read.tree("output/harvati_noPCA_c1.trees"), read.tree("output/harvati_noPCA_c2.trees"))

mk_trees_pca <- c(read.tree("output/harvati2004_mkModel_nohomopops_PCA_c1.trees"), read.tree("output/harvati2004_mkModel_nohomopops_PCA_c2.trees"))

png(filename = "Documents/Harvati_Reanalysis_Manuscript/figures/figure2_compareTrees.png", width = 1080, height = 1080)

ptsiz <- 2
lnsiz <- 4
labsiz <- 2
axsiz <- 2
tsiz <- 2
lsiz = 2
textsiz = 4
par(mfrow =c(2,2))
par(mar = c(5,5,5,5))
# par(mfrow =c(2,2), mai=c(0.65,0.65,0.65,0.5), xpd = T)
tiplabs <- mvBM_trees[[1]]$tip.label
comparison_mvBM_Mol <- compareTrees(prop.part.df(tenk_trees_filtered), prop.part.df(mvBM_trees), tiplabs)
plot(comparison_mvBM_Mol, xlab = "molecular posterior probabilities", ylab = "morphological posterior probabilities", cex = ptsiz, cex.lab = labsiz, cex.axis = axsiz)
abline(0, 1, lwd = lnsiz, lty = 2, col = "darkgrey")
title(main = "mvBM vs. Molecular Compare-Trees Plot", cex.main = tsiz, line = 0.5)
box(which = "figure", lty = 3)
text(labels = "a)", cex = textsiz, x = -0.1, y = 1.1, xpd = T, font = 4)

hist(comparison_mvBM_Mol[comparison_mvBM_Mol[,1] > 0.5,2], breaks = 20, col = rgb(0,0,0,0.5), ylim = c(0,45), cex.lab = labsiz,
     xlab = "morphological probability", cex.axis = axsiz, main = "")
hist(comparison_mvBM_Mol[comparison_mvBM_Mol[,1] < 0.5,2], breaks = 20, add = T, col = rgb(red = 1,1,1,0.5))
title(main = "mvBM Molecular Probability Comparison", cex.main = tsiz, line = 0.5)
legend(x = "topright", legend = c("molecular probability < 0.5", "molecular probability > 0.5"), fill = c(rgb(red = 1,1,1,0.5), rgb(0,0,0,0.5)), cex = lsiz)
box(which = "plot")
box(which = "figure", lty = 3)
text(labels = "b)", cex = textsiz, x = -0.1, y = 49.5, xpd = T, font = 4)


comparison_mk_Mol <- compareTrees(prop.part.df(tenk_trees_filtered), prop.part.df(mk_trees_pca), tiplabs)
plot(comparison_mk_Mol, xlab = "molecular posterior probabilities", ylab = "morphological posterior probabilities", cex = ptsiz, cex.lab = labsiz, cex.axis = axsiz)
abline(0, 1, lwd = lnsiz, lty = 2, col = "darkgrey")
title(main = "Mk-model vs. Molecular Compare-Trees Plot", cex.main = tsiz, line = 0.5)
box(which = "figure", lty = 3)

hist(comparison_mk_Mol[comparison_mk_Mol[,1] > 0.5,2], breaks = 20, col = rgb(0,0,0,0.5), ylim = c(0,110), 
     xlab = "morphological probability", cex.lab = labsiz, cex.axis = axsiz, main = "")
hist(comparison_mk_Mol[comparison_mk_Mol[,1] < 0.5,2], breaks = 20, add = T, col = rgb(red = 1,1,1,0.5))
title(main = "Mk-model Molecular Probability Comparison", cex.main = tsiz, line = 0.5)
legend(x = "topright", legend = c("molecular probability < 0.5", "molecular probability > 0.5"), fill = c(rgb(red = 1,1,1,0.5), rgb(0,0,0,0.5)), cex = lsiz)
box(which = "plot")
box(which = "figure", lty = 3)

dev.off()


#######################################
#### TESTING OUT THE FORKING PATHS ####
#######################################

load("tenk_trees_filtered")
mol_tree_gct <- greedyCT(tenk_trees_filtered)
mol_tree_mcc <- maxCladeCred(tenk_trees_filtered)


trees1 <- read.tree("/Volumes/1TB/Harvati_Empirical/output/infPriorsMoreProps/harvati_PCA_bothSexes_noHomopops_mvBM_PCA99_c1.trees")
trees2 <- read.tree("/Volumes/1TB/Harvati_Empirical/output/infPriorsMoreProps/harvati_PCA_bothSexes_noHomopops_mvBM_PCA99_c2.trees")

mk_trees1 <- read.tree("/Volumes/1TB/Harvati_Empirical/output/empirical_42traits_jenks_nohomopops_4_PCA99_MK_trees_run_1.trees")
mk_trees2 <- read.tree("/Volumes/1TB/Harvati_Empirical/output/empirical_42traits_jenks_nohomopops_4_PCA99_MK_trees_run_2.trees")

freek_trees1 <- read.tree("/Volumes/1TB/Harvati_Empirical/output/empirical_42traits_jenks_nohomopops_4_PCA99_trees_run_1.trees")
freek_trees2 <- read.tree("/Volumes/1TB/Harvati_Empirical/output/empirical_42traits_jenks_nohomopops_4_PCA99_trees_run_2.trees")


# trees1 <- read.tree("/Volumes/1TB/Harvati_Empirical/output/harvati_PCA_noHomopops_fixRM_PCA99_c1.trees")
# trees2 <- read.tree("/Volumes/1TB/Harvati_Empirical/output/harvati_PCA_noHomopops_fixRM_PCA99_c2.trees")

DC_MP_tree <- read.tree("/Volumes/1TB/Harvati_Empirical/output/maximumParsimonyTree_empirical_noHomoPops_PCA99_divergenceCoding.txt")
nj_tree <- read.tree("/Volumes/1TB/Harvati_Empirical/output/empirical_PCA99_neighborJoining.txt")
upgma_tree <- read.tree("/Volumes/1TB/Harvati_Empirical/output/empirical_PCA99_upgma.txt")
# trees1 <- read.tree("/Volumes/1TB/Harvati_Empirical/output/harvati_PCA_noHomopops_uninfRM_PCA99_c1.trees")
# trees2 <- read.tree("/Volumes/1TB/Harvati_Empirical/output/harvati_PCA_noHomopops_uninfRM_PCA99_c2.trees")
# 
# trees1 <- read.tree("output/harvati_noPCA_c1.trees") 
# trees2 <- read.tree("output/harvati_noPCA_c2.trees")

trees <- mvBM_trees <- c(trees1, trees2)
mk_trees <- c(mk_trees1, mk_trees2)
freek_trees <- c(freek_trees1, freek_trees2)

tiplabs <- trees[[1]]$tip.label
comparison <- compareTrees(prop.part.df(mvBM_trees1), prop.part.df(mvBM_trees2), tiplabs)
plot(comparison, main = paste0("same analysis comparetrees, r^2 = ", round(cor(comparison)[1,2]^2, digits = 3)))

comparison <- compareTrees(prop.part.df(mk_trees1), prop.part.df(mk_trees2), tiplabs)
plot(comparison, main = paste0("same analysis comparetrees, r^2 = ", round(cor(comparison)[1,2]^2, digits = 3)))

comparison <- compareTrees(prop.part.df(freek_trees1), prop.part.df(freek_trees2), tiplabs)
plot(comparison, main = paste0("same analysis comparetrees, r^2 = ", round(cor(comparison)[1,2]^2, digits = 3)))

RF.dist(maxCladeCred(trees1), maxCladeCred(trees2))
RF.dist(greedyCT(mk_trees1), greedyCT(mk_trees2))
RF.dist(greedyCT(freek_trees1), greedyCT(freek_trees2))

mvBM_mcc <- maxCladeCred(trees) 
mvBM_gct <- greedyCT(trees)

mk_mcc <- maxCladeCred(mk_trees) 
mk_gct <- greedyCT(mk_trees)

freek_mcc <- maxCladeCred(freek_trees) 
freek_gct <- greedyCT(freek_trees)


outg <- as.character(unlist(sapply(c("Macaca", "Mandrillus", "Papio"), function(genus) test_mcc$tip.label[grepl(genus, test_mcc$tip.label)])))
test_mcc <- reroot(tree = test_mcc, node.number = getMRCA(test_mcc, outg), position = 
                     test_mcc$edge.length[which(test_mcc$edge[,2] ==  getMRCA(test_mcc, outg))] / 2)
test_gct <- reroot(tree = test_gct, node.number = getMRCA(test_gct, outg), position = 
                     test_gct$edge.length[which(test_gct$edge[,2] ==  getMRCA(test_gct, outg))] / 2)
sapply(1:21, function(x) nodeheight(mol_tree_mcc, x))

#write the 10k tree for mvBM rate matrix estimation
#but first add a bit of branch length to the kinda baboon and cynocephalus
# write.tree(mol_tree_mcc, file = "/Volumes/1TB/Harvati_Empirical/data2_neanF/mol_tree_mcc.nex")
# mol_tree_mcc_noNean <- drop.tip(mol_tree_mcc, "Homo_neanderthalensis")
# write.tree(mol_tree_mcc_noNean, file = "/Volumes/1TB/Harvati_Empirical/data2_neanF/mol_tree_mcc_noNean.nex")

plot(test_gct)
plot(test_mcc)

RF.dist(mvBM_gct, mol_tree_gct)
RF.dist(mvBM_mcc, mol_tree_gct)
RF.dist(mk_mcc, mol_tree_gct)
RF.dist(mk_gct, mol_tree_gct)
RF.dist(freek_mcc, mol_tree_gct)
RF.dist(freek_gct, mol_tree_gct)
RF.dist(upgma_tree, mol_tree_gct)
RF.dist(nj_tree, mol_tree_gct, normalize = F)
RF.dist(DC_MP_tree, mol_tree_gct, normalize = F)

all_candidate_trees <- c(molecular_tree = mol_tree_gct, mvBM = mvBM_mcc, mk = mk_mcc, freek = freek_mcc, nj = nj_tree, upgma = upgma_tree, max_parsimony = DC_MP_tree)

round(as.matrix(RF.dist(all_candidate_trees, normalize = T)), 2)

par(mfrow = c(3,1))
plot(table(sapply(1:length(trees), function(tree) 
  round(RF.dist(trees[tree], mol_tree_gct, normalize = T), 2))) / length(trees) * 100, xlim = c(0,1), xlab = "", ylab = "")
title("multivariate Brownian motion")
plot(table(sapply(1:length(mk_trees), function(tree) 
  round(RF.dist(mk_trees[tree], mol_tree_gct, normalize = T), 2))) / length(mk_trees) * 100, xlim = c(0,1), xlab = "", ylab = "percent", cex.lab = 1.5)
title("MK-model")
plot(table(sapply(1:length(mk_trees), function(tree) 
  round(RF.dist(freek_trees[tree], mol_tree_gct, normalize = T), 2))) / length(freek_trees) * 100, xlim = c(0,1), xlab = "normalized RF-distance", ylab = "", cex.lab = 1.5)
title("ordered free-K model")


table(sapply(1:length(trees), function(tree) RF.dist(trees[tree], mol_tree_gct)))

par(mfrow = c(1,1))
comparison <- compareTrees(prop.part.df(tenk_trees_filtered), prop.part.df(trees), tiplabs)
plot(comparison, main = "mvBM analysis comparetrees", xlab = "molecular", ylab = "morphological")

##############################################
### Getting species-only analyses in there ###
##############################################

# load("tenk_trees_filtered")
# tenk_trees_filtered_spOnly <- tenk_trees_filtered
# for(i in 1:length(tenk_trees_filtered_spOnly)){
#   tr <- tenk_trees_filtered_spOnly[[i]]
#   tips_to_drop <- c("Gorilla_gorilla_beringei", "Pan_troglodytes_schweinfurthii", "Pan_troglodytes_verus", "Papio_hamadryas_anubis", "Papio_hamadryas_cynocephalus",
#                     "Papio_hamadryas_papio", "Papio_hamadryas_ursinus", "Papio_hamadryas_kindae")      
#   tr <- drop.tip(tr, tips_to_drop)
#   tip_names_subsp <- tr$tip.label
#   tip_names_subsp <- as.vector(sapply(tip_names_subsp, function(tipname) paste0(strsplit(tipname, "_")[[1]][1:2], collapse = "_")))
#   tr$tip.label <- tip_names_subsp
#   tenk_trees_filtered_spOnly[[i]] <- tr
# }
# mol_tree_mcc_spOnly <- maxCladeCred(tenk_trees_filtered_spOnly)
# mol_tree_mcc_spOnly_noNean <- drop.tip(mol_tree_mcc_spOnly, "Homo_neanderthalensis")
# sapply(1:21, function(x) nodeheight(mol_tree_mcc_spOnly, x))[1:length(mol_tree_mcc_spOnly$tip.label)]
# # write.tree(mol_tree_mcc_spOnly, file = "/Volumes/1TB/Harvati_Empirical/data2_neanF/mol_tree_mcc_spOnly.nex")
# # write.tree(mol_tree_mcc_spOnly_noNean, file = "/Volumes/1TB/Harvati_Empirical/data2_neanF/mol_tree_mcc_spOnly_noNean.nex")
# save(tenk_trees_filtered_spOnly, file = "tenk_trees_filtered_spOnly")

load("tenk_trees_filtered_spOnly")
mol_tree_mcc_spOnly <- maxCladeCred(tenk_trees_filtered_spOnly)

uvBM_trees1 <- read.tree("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_bothSexes_noHomopops_justSpecies_uvBM_PCA99_c1.trees")
uvBM_trees2 <- read.tree("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_bothSexes_noHomopops_justSpecies_uvBM_PCA99_c2.trees")

mvBM_trees1 <- read.tree("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_bothSexes_noHomopops_justSpecies_inf_mvBM_PCA99_c1.trees")
mvBM_trees2 <- read.tree("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_bothSexes_noHomopops_justSpecies_inf_mvBM_PCA99_c2.trees")

mvBM_trees <- c(mvBM_trees1, mvBM_trees2)
uvBM_trees <- c(uvBM_trees1, uvBM_trees2)
mvBM_tree_mcc <- (maxCladeCred(mvBM_trees))
uvBM_tree_mcc <- (maxCladeCred(uvBM_trees))
# mvBM_tree_gct <- (greedyCT(mvBM_trees))
# uvBM_tree_gct <- (greedyCT(uvBM_trees))
# plot(mcc_tree)
# plot(mol_tree_mcc)
RF.dist(mol_tree_mcc_spOnly, mvBM_mcc_tree, normalize = T)
comparison <- compareTrees(prop.part.df(tenk_trees_filtered_spOnly), prop.part.df(trees), tiplabs)
plot(comparison, main = "mvBM analysis comparetrees", xlab = "molecular", ylab = "morphological")
table(sapply(1:length(trees), function(tree) RF.dist(trees[tree], mol_tree_mcc_spOnly, normalize = T)))

mk_trees1 <- read.tree("/Volumes/1TB/Harvati_Empirical/output/empirical_39traits_jenks_noHomoPops_justSpecies_4_PCA99_MK_trees_run_1.trees")
mk_trees2 <- read.tree("/Volumes/1TB/Harvati_Empirical/output/empirical_39traits_jenks_noHomoPops_justSpecies_4_PCA99_MK_trees_run_2.trees")
mk_trees <- c(mk_trees1, mk_trees2)
mk_tree_mcc <- maxCladeCred(mk_trees)

freek_trees1 <- read.tree("/Volumes/1TB/Harvati_Empirical/output/empirical_39traits_jenks_noHomoPops_justSpecies_4_PCA99_trees_run_1.trees")
freek_trees2 <- read.tree("/Volumes/1TB/Harvati_Empirical/output/empirical_39traits_jenks_noHomoPops_justSpecies_4_PCA99_trees_run_2.trees")
freek_trees <- c(freek_trees1, freek_trees2)
freek_tree_mcc <- maxCladeCred(freek_trees)

DC_MP_tree <- read.tree("/Volumes/1TB/Harvati_Empirical/output/justSpecies/maximumParsimonyTree_empirical_noHomoPops_justSpecies_PCA99_divergenceCoding.txt")
jenks_MP_tree <- read.tree("/Volumes/1TB/Harvati_Empirical/output/justSpecies/maximumParsimonyTree_empirical_noHomoPops_justSpecies_PCA99_jenksCoding.txt")
nj_tree <- read.tree("/Volumes/1TB/Harvati_Empirical/output/justSpecies/empirical_PCA99_neighborJoining.txt")
upgma_tree <- read.tree("/Volumes/1TB/Harvati_Empirical/output/justSpecies/empirical_PCA99_upgma.txt")

# mvBM_tree_gct <- (greedyCT(trees))
# mk_tree_gct <- greedyCT(mk_trees)
# freek_tree_gct <- greedyCT(freek_trees)

#reroot all unrooted trees
outg <- as.character(unlist(sapply(c("Macaca", "Mandrillus", "Papio"), function(genus) mvBM_tree_mcc$tip.label[grepl(genus, mvBM_tree_mcc$tip.label)])))
mvBM_tree_mcc <- reroot(tree = mvBM_tree_mcc, node.number = getMRCA(mvBM_tree_mcc, outg), position = 
                          mvBM_tree_mcc$edge.length[which(mvBM_tree_mcc$edge[,2] ==  getMRCA(mvBM_tree_mcc, outg))] / 2); plot(mvBM_tree_mcc)
uvBM_tree_mcc <- reroot(tree = uvBM_tree_mcc, node.number = getMRCA(uvBM_tree_mcc, outg), position = 
                          uvBM_tree_mcc$edge.length[which(uvBM_tree_mcc$edge[,2] ==  getMRCA(uvBM_tree_mcc, outg))] / 2); plot(uvBM_tree_mcc)
mk_tree_mcc <- reroot(tree = mk_tree_mcc, node.number = getMRCA(mk_tree_mcc, outg), position = 
                          mk_tree_mcc$edge.length[which(mk_tree_mcc$edge[,2] ==  getMRCA(mk_tree_mcc, outg))] / 2); plot(mk_tree_mcc)
freek_tree_mcc <- reroot(tree = freek_tree_mcc, node.number = getMRCA(freek_tree_mcc, outg), position = 
                          freek_tree_mcc$edge.length[which(freek_tree_mcc$edge[,2] ==  getMRCA(freek_tree_mcc, outg))] / 2); plot(freek_tree_mcc)


nj_tree$edge.length[nj_tree$edge.length < 0] <- 0
while(length(which(nj_tree$edge[,2] == getMRCA(nj_tree, outg))) == 0){
  nj_tree <- reroot(nj_tree, node.number = sample(size = 1, 1:nj_tree$Nnode))
}
nj_tree <- reroot(tree = nj_tree, node.number = getMRCA(nj_tree, outg), position = 
                          nj_tree$edge.length[which(nj_tree$edge[,2] ==  getMRCA(nj_tree, outg))] / 2); plot(nj_tree)
 

DC_MP_tree <- reroot(tree = DC_MP_tree, node.number = getMRCA(DC_MP_tree, outg), position = 
                          DC_MP_tree$edge.length[which(DC_MP_tree$edge[,2] ==  getMRCA(DC_MP_tree, outg))] / 2); plot(DC_MP_tree)
jenks_MP_tree <- reroot(tree = jenks_MP_tree, node.number = getMRCA(jenks_MP_tree, outg), position = 
                          jenks_MP_tree$edge.length[which(jenks_MP_tree$edge[,2] ==  getMRCA(jenks_MP_tree, outg))] / 2); plot(jenks_MP_tree)



mvBM_RFdists <- sapply(1:length(mvBM_trees), function(tree) round(RF.dist(mvBM_trees[tree], mol_tree_mcc_spOnly, normalize = T), 2))
mvBM_table <- table(mvBM_RFdists) / length(mvBM_trees) * 100
uvBM_RFdists <- sapply(1:length(uvBM_trees), function(tree) round(RF.dist(uvBM_trees[tree], mol_tree_mcc_spOnly, normalize = T), 2))
uvBM_table <- table(uvBM_RFdists) / length(uvBM_trees) * 100
mk_RFdists <- sapply(1:length(mk_trees), function(tree) round(RF.dist(mk_trees[tree], mol_tree_mcc_spOnly, normalize = T), 2))
mk_table <- table(mk_RFdists) / length(mk_trees) * 100
freek_RFdists <- sapply(1:length(mk_trees), function(tree) round(RF.dist(freek_trees[tree], mol_tree_mcc_spOnly, normalize = T), 2))
freek_table <- table(freek_RFdists) / length(freek_trees) * 100

par(mfrow = c(4,1))
plot(mvBM_table, xlim = c(0,1), xlab = "normalized RF-Distance", ylab = "percent", ylim = c(0,50), cex.lab = 1.5)
legend(x = "topright", legend = paste0("µ = ", round(mean(mvBM_RFdists), 3), ", σ = ", round(sd(mvBM_RFdists), 3)))
title("multivariate Brownian motion")
plot(uvBM_table, xlim = c(0,1), xlab = "", ylab = "percent", ylim = c(0,50), cex.lab = 1.5)
title("univariate Brownian motion")
plot(mk_table, xlim = c(0,1), xlab = "", ylab = "percent", cex.lab = 1.5, ylim = c(0,50), cex.lab = 1.5)
title("MK-model")
plot(freek_table, xlim = c(0,1), xlab = "normalized RF-distance", ylab = "percent", cex.lab = 1.5, ylim = c(0,50))
title("ordered free-K model")

sum(as.numeric(mvBM_table) / 100 * as.numeric(attr(mvBM_table, "dimnames")[[1]]))
sum(as.numeric(mk_table) / 100 * as.numeric(attr(mk_table, "dimnames")[[1]]))
sum(as.numeric(freek_table) / 100 * as.numeric(attr(freek_table, "dimnames")[[1]]))

all_candidate_trees_mcc <- c(molecular_tree = mol_tree_mcc_spOnly, mvBM = mvBM_tree_mcc, uvBM = uvBM_tree_mcc, mk = mk_tree_mcc, freek = freek_tree_mcc, nj = nj_tree, upgma = upgma_tree, max_parsimony_DC = DC_MP_tree, max_parsimony_jenks = jenks_MP_tree)
all_candidate_trees_gct <- c(molecular_tree = mol_tree_mcc_spOnly, mvBM = mvBM_tree_mcc, uvBM = uvBM_tree_gct, mk = mk_tree_gct, freek = freek_tree_gct, nj = nj_tree, upgma = upgma_tree, max_parsimony_DC = DC_MP_tree, max_parsimony_jenks = jenks_MP_tree)

round(as.matrix(RF.dist(all_candidate_trees_mcc, normalize = T)), 2)
round(as.matrix(RF.dist(all_candidate_trees_gct, normalize = T)), 2)

#######################################
### MAKING THE FINAL SET OF FIGURES ###
#######################################

colfunc <- colorRampPalette(c("white", "black"))
discretizeContData <- function(data, range = c(0,1), bins = 10){
  bintervals <- seq(from = range[1], to = range[2], length.out = bins+1)
  binMeans <- bintervals[-1] + diff(bintervals) / 2 - bintervals[2]
  discrData <- sapply(data, function(num) sum(num - bintervals > 0))
  discrFreqs <- sapply(1:bins, function(bin) sum(discrData == bin))
  discrFreqs <- discrFreqs / sum(discrFreqs)
  return(rbind(freqs = discrFreqs, binLocations = binMeans))
}

tiplabs <- mol_tree_mcc_spOnly$tip.label
coph_plot_mvBM <- cophylo(mol_tree_mcc_spOnly, mvBM_tree_mcc, rotate = T)
mol_tree_mcc_spOnly <- (coph_plot_mvBM$trees[[1]])
coph_plot_uvBM <- cophylo(mol_tree_mcc_spOnly, uvBM_tree_mcc, rotate = T)
coph_plot_mk <- cophylo(mol_tree_mcc_spOnly, mk_tree_mcc, rotate = T)
coph_plot_freek <- cophylo(mol_tree_mcc_spOnly, freek_tree_mcc, rotate = T)

comparison_mvBM <- compareTrees(prop.part.df(mvBM_trees), prop.part.df(tenk_trees_filtered_spOnly), tiplabs)
comparison_uvBM <- compareTrees(prop.part.df(uvBM_trees), prop.part.df(tenk_trees_filtered_spOnly), tiplabs)
comparison_mk <- compareTrees(prop.part.df(mk_trees), prop.part.df(tenk_trees_filtered_spOnly), tiplabs)
comparison_freek <- compareTrees(prop.part.df(freek_trees), prop.part.df(tenk_trees_filtered_spOnly), tiplabs)

ptsiz <- 2
lnsiz <- 4
labsiz <- 2
axsiz <- 2
tsiz <- 2
lsiz = 2
textsiz = 4

########################################
########################################
############### FIGURE 1 ###############
########################################
########################################

png(filename = "Documents/Harvati_Reanalysis_Manuscript/figures/figure1_final.png", width = 2800, height = 3000)
# par(mfrow = c(4, 3))


layout(matrix(c(1,1,1,1,1,2,2,2,3,3,3,4,4,4,4,4,5,5,5,6,6,6,7,7,7,7,7,8,8,8,9,9,9,10,10,10,10,10,11,11,11,12,12,12), nrow = 4, ncol = 11, byrow = TRUE))

##########################
#### vs the mvBM tree ####
##########################

#cophylo
plot.cophylo(coph_plot_mvBM, mar=c(1,1,4,1), lwd = 4.5, fsize = 3.5, link.lwd=4)

par(xpd=TRUE)
col_left <- as.integer(internal_nodes_probs(coph_plot_mvBM[[1]][[1]], mvBM_trees)*100)
col_left[1] <- 100
col_left[col_left == 0] <- 1
col_left <- colfunc(100)[col_left]
col_left[is.na(col_left)] <- colfunc(100)[1]
nodelabels.cophylo(pch=21, frame="none", bg=col_left, cex=3.5, which = "left")

col_right <- as.integer(internal_nodes_probs(coph_plot_mvBM[[1]][[2]], mvBM_trees)*100)
col_right[1] <- 100
col_right[col_right == 0] <- 1
col_right <- colfunc(100)[col_right]
col_right[is.na(col_right)] <- colfunc(100)[1]
nodelabels.cophylo(pch=21, frame="none", bg=col_right, cex=3.5, which = "right")

title("       Molecular                         Morphological", cex.main = 4.5, line = -0.1)
title(paste0("RF-Distance = ", RF.dist(mol_tree_mcc_spOnly, mvBM_tree_mcc, normalize = T)), cex.main = 2.75, line = -1.2)

box(which = "figure", lty = 5)
text(labels = "a)", cex = textsiz, x = 0.53, y = 1.065, xpd = T, font = 4, col = "darkred")
text(labels = "mvBM", cex = textsiz, x = -0.5, y = 1, xpd = T, font = 4, srt = 90, col = "darkred")
segments(x0 = -0.485, x1 = -0.485, y0 = 0.925, y1 = 1.075, lwd = 3)

#comparetrees
par(mar=c(8,9,8,2))
plot(comparison_mvBM[,1], comparison_mvBM[,2] + 1/3, main = "Compare-Trees Plot and Histogram", xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
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

highprobs <- discretizeContData(comparison_mvBM[comparison_mvBM[,2] > 0.5,1])
lowprobs <- discretizeContData(comparison_mvBM[comparison_mvBM[,2] < 0.5,1])
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

text(labels = "b)", cex = textsiz, x = 1.0365, y = 1.51, xpd = T, font = 4, col = "darkred")
box(which = "figure", lty = 5)

#RF-dist of entire posterior
plot(mvBM_table, xlim = c(0,1), xlab = "", ylab = "", ylim = c(0,50), 
     cex.lab = 4, cex.axis = 3, lwd = 4.5, xaxt = "n", yaxt = "n", main = "RF-Distances from Molecular Tree", cex.main = 4.5)
title(ylab = "percent", xlab = "normalized RF-Distance", cex.lab = 4, line = 5.75)
axis(1, at = 1:7/10, labels = rep("", 7), lwd = 4, cex.axis = 3, tck = -0.015)
mtext(text = 1:7/10, side = 1, at = 1:7/10, cex = 2, line = 2.25)
axis(2, at = 0:5*10, labels =  rep("", 6), lwd = 4, cex.axis = 3, tck = -0.015)
mtext(text = 0:5*10, side = 2, at = 0:5*10, cex = 2, line = 2)
box(lwd=3)
legend(x = "topright", legend = paste0("µ = ", "0.460", ", σ = ", round(sd(mvBM_RFdists), 3)), cex = 3, bty = "n") #round(mean(mvBM_RFdists), 3) chops off 0
text(labels = "c)", cex = textsiz, x = 1.0365, y = 56.5, xpd = T, font = 4, col = "darkred")
box(which = "figure", lty = 5)



##########################
#### vs the uvBM tree ####
##########################

#cophylo
par(xpd=TRUE)
plot.cophylo(coph_plot_uvBM, mar=c(1,1,4,1), lwd = 4.5, fsize = 3.5, link.lwd=4)

col_left <- as.integer(internal_nodes_probs(coph_plot_uvBM[[1]][[1]], uvBM_trees)*100)
col_left[1] <- 100
col_left[col_left == 0] <- 1
col_left <- colfunc(100)[col_left]
col_left[is.na(col_left)] <- colfunc(100)[1]
nodelabels.cophylo(pch=21, frame="none", bg=col_left, cex=3.5, which = "left")

col_right <- as.integer(internal_nodes_probs(coph_plot_uvBM[[1]][[2]], uvBM_trees)*100)
col_right[1] <- 100
col_right[col_right == 0] <- 1
col_right <- colfunc(100)[col_right]
col_right[is.na(col_right)] <- colfunc(100)[1]
nodelabels.cophylo(pch=21, frame="none", bg=col_right, cex=3.5, which = "right")

title("       Molecular                         Morphological", cex.main = 4.5, line = -0.1)
title(paste0("RF-Distance = ", RF.dist(mol_tree_mcc_spOnly, uvBM_tree_mcc, normalize = T)), cex.main = 2.75, line = -1.2)

box(which = "figure", lty = 5)
text(labels = "uvBM", cex = textsiz, x = -0.5, y = 1, xpd = T, font = 4, srt = 90, col = "darkred")
segments(x0 = -0.485, x1 = -0.485, y0 = 0.925, y1 = 1.075, lwd = 3)

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

#RF-dist of entire posterior
plot(uvBM_table, xlim = c(0,1), xlab = "", ylab = "", ylim = c(0,50), 
     cex.lab = 4, cex.axis = 3, lwd = 4.5, xaxt = "n", yaxt = "n", main = "RF-Distances from Molecular Tree", cex.main = 4.5)
title(ylab = "percent", xlab = "normalized RF-Distance", cex.lab = 4, line = 5.75)
axis(1, at = 1:7/10, labels = rep("", 7), lwd = 4, cex.axis = 3, tck = -0.015)
mtext(text = 1:7/10, side = 1, at = 1:7/10, cex = 2, line = 2.25)
axis(2, at = 0:5*10, labels =  rep("", 6), lwd = 4, cex.axis = 3, tck = -0.015)
mtext(text = 0:5*10, side = 2, at = 0:5*10, cex = 2, line = 2)
box(lwd=3)
legend(x = "topright", legend = paste0("µ = ", round(mean(uvBM_RFdists), 3), ", σ = ", round(sd(uvBM_RFdists), 3)), cex = 3, bty = "n")
box(which = "figure", lty = 5)



##########################
#### vs the mk tree ####
##########################

#cophylo
par(xpd=TRUE)
plot.cophylo(coph_plot_mk, mar=c(1,1,4,1), lwd = 4.5, fsize = 3.5, link.lwd=4)

col_left <- as.integer(internal_nodes_probs(coph_plot_mk[[1]][[1]], mk_trees)*100)
col_left[1] <- 100
col_left[col_left == 0] <- 1
col_left <- colfunc(100)[col_left]
col_left[is.na(col_left)] <- colfunc(100)[1]
nodelabels.cophylo(pch=21, frame="none", bg=col_left, cex=3.5, which = "left")

col_right <- as.integer(internal_nodes_probs(coph_plot_mk[[1]][[2]], mk_trees)*100)
col_right[1] <- 100
col_right[col_right == 0] <- 1
col_right <- colfunc(100)[col_right]
col_right[is.na(col_right)] <- colfunc(100)[1]
nodelabels.cophylo(pch=21, frame="none", bg=col_right, cex=3.5, which = "right")

title("       Molecular                         Morphological", cex.main = 4.5, line = -0.1)
title(paste0("RF-Distance = ", RF.dist(mol_tree_mcc_spOnly, mk_tree_mcc, normalize = T)), cex.main = 2.75, line = -1.2)

box(which = "figure", lty = 5)
text(labels = "Mk", cex = textsiz, x = -0.5, y = 1, xpd = T, font = 4, srt = 90, col = "darkred")
segments(x0 = -0.485, x1 = -0.485, y0 = 0.95, y1 = 1.05, lwd = 3)

#comparetrees
par(mar=c(8,9,8,2))
plot(comparison_mk[,1], comparison_mk[,2] + 1/3, main = "Compare-Trees Plot and Histogram", xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
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

highprobs <- discretizeContData(comparison_mk[comparison_mk[,2] > 0.5,1])
lowprobs <- discretizeContData(comparison_mk[comparison_mk[,2] < 0.5,1])
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

#RF-dist of entire posterior
plot(mk_table, xlim = c(0,1), xlab = "", ylab = "", ylim = c(0,50), 
     cex.lab = 4, cex.axis = 3, lwd = 4.5, xaxt = "n", yaxt = "n", main = "RF-Distances from Molecular Tree", cex.main = 4.5)
title(ylab = "percent", xlab = "normalized RF-Distance", cex.lab = 4, line = 5.75)
axis(1, at = 2:10/10, labels = rep("", 9), lwd = 4, cex.axis = 3, tck = -0.015)
mtext(text = 2:10/10, side = 1, at = 2:10/10, cex = 2, line = 2.25)
axis(2, at = 0:5*10, labels =  rep("", 6), lwd = 4, cex.axis = 3, tck = -0.015)
mtext(text = 0:5*10, side = 2, at = 0:5*10, cex = 2, line = 2)
box(lwd=3)
legend(x = "topright", legend = paste0("µ = ", round(mean(mk_RFdists), 3), ", σ = ", round(sd(mk_RFdists), 3)), cex = 3, bty = "n")
box(which = "figure", lty = 5)



##########################
#### vs the freek tree ####
##########################

#cophylo
par(xpd=TRUE)
plot.cophylo(coph_plot_freek, mar=c(8,1,4,1), lwd = 4.5, fsize = 3.5, link.lwd=4)

col_left <- as.integer(internal_nodes_probs(coph_plot_freek[[1]][[1]], freek_trees)*100)
col_left[1] <- 100
col_left[col_left == 0] <- 1
col_left <- colfunc(100)[col_left]
col_left[is.na(col_left)] <- colfunc(100)[1]
nodelabels.cophylo(pch=21, frame="none", bg=col_left, cex=3.5, which = "left")

col_right <- as.integer(internal_nodes_probs(coph_plot_freek[[1]][[2]], freek_trees)*100)
col_right[1] <- 100
col_right[col_right == 0] <- 1
col_right <- colfunc(100)[col_right]
col_right[is.na(col_right)] <- colfunc(100)[1]
nodelabels.cophylo(pch=21, frame="none", bg=col_right, cex=3.5, which = "right")

title("       Molecular                         Morphological", cex.main = 4.5, line = -0.1)
title(paste0("RF-Distance = ", RF.dist(mol_tree_mcc_spOnly, freek_tree_mcc, normalize = T)), cex.main = 2.75, line = -1.2)

box(which = "figure", lty = 5)
text(labels = "ordered CTMC", cex = textsiz, x = -0.5, y = 0.9, xpd = T, font = 4, srt = 90, col = "darkred")
segments(x0 = -0.485, x1 = -0.485, y0 = 0.725, y1 = 1.075, lwd = 3)


xl <- -0.3; yb <- -0.115; xr <- 0.3; yt <- -0.065;
rect(
  head(seq(xl,xr,(xr-xl)/100),-1),
  yb,
  tail(seq(xl,xr,(xr-xl)/100),-1),
  yt,
  col=colfunc(100)
)
text(labels = seq(from = 0, to = 1, length.out = 11), y = yb-(yt-yb)/2.5, x = seq(xl,xr,length.out = 11), las=2, cex=2)
text(labels = "posterior probability", x = (xl+xr)/2, y = yt+(yt-yb)/2.5, font = 2, cex = 2.5)

# clade_prob(c("Macaca_fa", "Macaca_mu"), trees)
# clade_prob(c("Homo", "Pan"), trees)
# clade_prob(c("Pan"), trees)

#comparetrees
par(mar=c(8,9,8,2))
plot(comparison_freek[,1], comparison_freek[,2] + 1/3, main = "Compare-Trees Plot and Histogram", xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
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

highprobs <- discretizeContData(comparison_freek[comparison_freek[,2] > 0.5,1])
lowprobs <- discretizeContData(comparison_freek[comparison_freek[,2] < 0.5,1])
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

#RF-dist of entire posterior
plot(freek_table, xlim = c(0,1), xlab = "", ylab = "", ylim = c(0,50), 
     cex.lab = 4, cex.axis = 3, lwd = 4.5, xaxt = "n", yaxt = "n", main = "RF-Distances from Molecular Tree", cex.main = 4.5)
title(ylab = "percent", xlab = "normalized RF-Distance", cex.lab = 4, line = 5.75)
axis(1, at = 1:8/10, labels = rep("", 8), lwd = 4, cex.axis = 3, tck = -0.015)
mtext(text = 1:8/10, side = 1, at = 1:8/10, cex = 2, line = 2.25)
axis(2, at = 0:5*10, labels =  rep("", 6), lwd = 4, cex.axis = 3, tck = -0.015)
mtext(text = 0:5*10, side = 2, at = 0:5*10, cex = 2, line = 2)
box(lwd=3)
legend(x = "topright", legend = paste0("µ = ", round(mean(freek_RFdists), 3), ", σ = ", round(sd(freek_RFdists), 3)), cex = 3, bty = "n")
box(which = "figure", lty = 5)

dev.off()


########################################
########################################
############### FIGURE 2 ###############
########################################
########################################


jenks_MP_tree
DC_MP_tree

#read in the stepwise parsimony distance matrices, i.e. euclidean distances on marginal traits
DC_dist <- readDist( "/Volumes/1TB/Harvati_Empirical/output/justSpecies/parsDistances_empirical_noHomoPops_justSpecies_PCA99_divergenceCoding.dist" )
Jenks_dist <- readDist( "/Volumes/1TB/Harvati_Empirical/output/justSpecies/parsDistances_empirical_noHomoPops_justSpecies_PCA99_jenksCoding.dist" )

#get BLs on trees
DC_MP_tree <- nnls.phylo(DC_MP_tree, DC_dist)
jenks_MP_tree <- nnls.phylo(jenks_MP_tree, Jenks_dist)

#reroot with cercopithecoids
DC_MP_tree <- reroot(tree = DC_MP_tree, node.number = getMRCA(DC_MP_tree, outg), position = 
                       DC_MP_tree$edge.length[which(DC_MP_tree$edge[,2] ==  getMRCA(DC_MP_tree, outg))] / 2); plot(DC_MP_tree)
jenks_MP_tree <- reroot(tree = jenks_MP_tree, node.number = getMRCA(jenks_MP_tree, outg), position = 
                          jenks_MP_tree$edge.length[which(jenks_MP_tree$edge[,2] ==  getMRCA(jenks_MP_tree, outg))] / 2); plot(jenks_MP_tree)


tiplabs <- mol_tree_mcc_spOnly$tip.label
coph_plot_DCMP <- cophylo(mol_tree_mcc_spOnly, DC_MP_tree, rotate = T)
coph_plot_JenksMP <- cophylo(mol_tree_mcc_spOnly, jenks_MP_tree, rotate = T)
coph_plot_nj <- cophylo(mol_tree_mcc_spOnly, nj_tree, rotate = T)
coph_plot_upgma <- cophylo(mol_tree_mcc_spOnly, upgma_tree, rotate = T)

ptsiz <- 2
lnsiz <- 4
labsiz <- 2
axsiz <- 2
tsiz <- 2
lsiz = 2
textsiz = 4




















#making the actual figre
png(filename = "Documents/Harvati_Reanalysis_Manuscript/figures/figure2_final.png", width = 2800, height = 1500)
par(mfrow = c(2, 2))

#vs. nj tree
plot.cophylo(coph_plot_nj, mar=c(1,1,4,1), lwd = 4.5, fsize = 3.5, link.lwd=4)

par(xpd=TRUE)
col_left <- as.integer(internal_nodes_probs(coph_plot_nj[[1]][[1]], nj_tree)*100)
col_left[1] <- 100
col_left[col_left == 0] <- 1
col_left <- colfunc(100)[col_left]
col_left[is.na(col_left)] <- colfunc(100)[1]
nodelabels.cophylo(pch=21, frame="none", bg=col_left, cex=3.5, which = "left")

col_right <- as.integer(internal_nodes_probs(coph_plot_nj[[1]][[2]], mol_tree_mcc_spOnly)*100)
col_right[1] <- 100
col_right[col_right == 0] <- 1
col_right <- colfunc(100)[col_right]
col_right[is.na(col_right)] <- colfunc(100)[1]
nodelabels.cophylo(pch=21, frame="none", bg=col_right, cex=3.5, which = "right")

title("       Molecular                         Morphological", cex.main = 4.5, line = 0)
title(paste0("RF-Distance = ", RF.dist(mol_tree_mcc_spOnly, nj_tree, normalize = T)), cex.main = 2.75, line = -1.2)

box(which = "figure", lty = 5)
text(labels = "a)", cex = textsiz, x = 0.53, y = 1.065, xpd = T, font = 4, col = "darkred")
text(labels = "NJ", cex = textsiz, x = -0.5, y = 1, xpd = T, font = 4, srt = 90, col = "darkred")
segments(x0 = -0.485, x1 = -0.485, y0 = 0.95, y1 = 1.05, lwd = 3)

#vs. upgma tree
plot.cophylo(coph_plot_upgma, mar=c(1,1,4,1), lwd = 4.5, fsize = 3.5, link.lwd=4)

par(xpd=TRUE)
col_left <- as.integer(internal_nodes_probs(coph_plot_upgma[[1]][[1]], upgma_tree)*100)
col_left[1] <- 100
col_left[col_left == 0] <- 1
col_left <- colfunc(100)[col_left]
col_left[is.na(col_left)] <- colfunc(100)[1]
nodelabels.cophylo(pch=21, frame="none", bg=col_left, cex=3.5, which = "left")

col_right <- as.integer(internal_nodes_probs(coph_plot_upgma[[1]][[2]], mol_tree_mcc_spOnly)*100)
col_right[1] <- 100
col_right[col_right == 0] <- 1
col_right <- colfunc(100)[col_right]
col_right[is.na(col_right)] <- colfunc(100)[1]
nodelabels.cophylo(pch=21, frame="none", bg=col_right, cex=3.5, which = "right")

title("       Molecular                         Morphological", cex.main = 4.5, line = 0)
title(paste0("RF-Distance = ", RF.dist(mol_tree_mcc_spOnly, upgma_tree, normalize = T)), cex.main = 2.75, line = -1.2)

box(which = "figure", lty = 5)
text(labels = "b)", cex = textsiz, x = 0.53, y = 1.065, xpd = T, font = 4, col = "darkred")
text(labels = "UPGMA", cex = textsiz, x = -0.5, y = 0.95, xpd = T, font = 4, srt = 90, col = "darkred")
segments(x0 = -0.485, x1 = -0.485, y0 = 0.825, y1 = 1.075, lwd = 3)

legend(x = -0.59, y = 0.15, legend = c("present", "absent"), col = c("black","black"), pch=c(19, 1), cex = 3.5, xpd = NA)


#vs. divergence coding max pars tree
plot.cophylo(coph_plot_DCMP, mar=c(1,1,4,1), lwd = 4.5, fsize = 3.5, link.lwd=4)

par(xpd=TRUE)
col_left <- as.integer(internal_nodes_probs(coph_plot_DCMP[[1]][[1]], DC_MP_tree)*100)
col_left[1] <- 100
col_left[col_left == 0] <- 1
col_left <- colfunc(100)[col_left]
col_left[is.na(col_left)] <- colfunc(100)[1]
nodelabels.cophylo(pch=21, frame="none", bg=col_left, cex=3.5, which = "left")

col_right <- as.integer(internal_nodes_probs(coph_plot_DCMP[[1]][[2]], mol_tree_mcc_spOnly)*100)
col_right[1] <- 100
col_right[col_right == 0] <- 1
col_right <- colfunc(100)[col_right]
col_right[is.na(col_right)] <- colfunc(100)[1]
nodelabels.cophylo(pch=21, frame="none", bg=col_right, cex=3.5, which = "right")

title("       Molecular                         Morphological", cex.main = 4.5, line = 0)
title(paste0("RF-Distance = ", RF.dist(mol_tree_mcc_spOnly, DC_MP_tree, normalize = T)), cex.main = 2.75, line = -1.2)

box(which = "figure", lty = 5)
text(labels = "c)", cex = textsiz, x = 0.53, y = 1.065, xpd = T, font = 4, col = "darkred")
text(labels = "MP-DC", cex = textsiz, x = -0.5, y = 0.95, xpd = T, font = 4, srt = 90, col = "darkred")
segments(x0 = -0.485, x1 = -0.485, y0 = 0.835, y1 = 1.065, lwd = 3)

#vs. jenks-4 coding max pars tree
plot.cophylo(coph_plot_DCMP, mar=c(1,1,4,1), lwd = 4.5, fsize = 3.5, link.lwd=4)

par(xpd=TRUE)
col_left <- as.integer(internal_nodes_probs(coph_plot_JenksMP[[1]][[1]], jenks_MP_tree)*100)
col_left[1] <- 100
col_left[col_left == 0] <- 1
col_left <- colfunc(100)[col_left]
col_left[is.na(col_left)] <- colfunc(100)[1]
nodelabels.cophylo(pch=21, frame="none", bg=col_left, cex=3.5, which = "left")

col_right <- as.integer(internal_nodes_probs(coph_plot_JenksMP[[1]][[2]], mol_tree_mcc_spOnly)*100)
col_right[1] <- 100
col_right[col_right == 0] <- 1
col_right <- colfunc(100)[col_right]
col_right[is.na(col_right)] <- colfunc(100)[1]
nodelabels.cophylo(pch=21, frame="none", bg=col_right, cex=3.5, which = "right")

title("       Molecular                         Morphological", cex.main = 4.5, line = 0)
title(paste0("RF-Distance = ", RF.dist(mol_tree_mcc_spOnly, jenks_MP_tree, normalize = T)), cex.main = 2.75, line = -1.2)

box(which = "figure", lty = 5)
text(labels = "d)", cex = textsiz, x = 0.53, y = 1.065, xpd = T, font = 4, col = "darkred")
text(labels = "MP-Jenks", cex = textsiz, x = -0.5, y = 0.95, xpd = T, font = 4, srt = 90, col = "darkred")
segments(x0 = -0.485, x1 = -0.485, y0 = 0.81, y1 = 1.08, lwd = 3)

dev.off()



