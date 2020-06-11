library(geomorph)
library(Matrix)
library(phytools)
library(phangorn)

normalize <- T
noPCA <- T
PCA_covBG <- F
collapseHomo <- F
collapseSubsp <- F

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
if(collapseSubsp){
  codes_linnaen[,4] <- sapply(1:length(codes_linnaen[,4]), function(x) paste(strsplit(codes_linnaen[,4], " ")[[x]][1:2], collapse = " "))
}
species <- unique(codes_linnaen[,4])

#what is the 5th letter? subspecies? probably not, since TT is troglodytes troglodytes
# coded_genus <- substr(codes, 1, 2)
# genus_key <- cbind(c("HM", "MC", "MN", "PP", "PN", "GO"), c("Homo", "Macaca", "Mandrillus", "Papio", "Pan", "Gorilla"))
# sex <- substr(codes, 6, 6)
# linnaen <- as.character(unique(taxa[,1]))
# genus <- sapply(1:length(linnaen), function(x) strsplit(linnaen[x], split = " ")[[1]][1])
# species <- sapply(1:length(linnaen), function(x) strsplit(linnaen[x], split = " ")[[1]][2])
# subspecies <- sapply(1:length(linnaen), function(x) strsplit(linnaen[x], split = " ")[[1]][3])

#load in the molecular tree

#mean center all landmarks
for(i in 1:15){for(j in 1:3){data[i,j,] <- data[i,j,] - mean(data[i,j,])}}
t.data <- gpagen(data)$data
t.data$logCsize <- log(t.data$Csize)
t.data <- as.matrix(subset(x = t.data, select = colnames(t.data)[colnames(t.data) != "Csize"]))
covTotal <- cov(t.data) #total covariance matrix
if(PCA_covBG == T){
  d_PCA <- as.data.frame(t.data)
  ntraits_PCA <- length(d_PCA[1,])
  Population <- rep("species", nrow(d_PCA))
  for(i in 1:nrow(codes_linnaen)){Population[startsWith(rownames(d_PCA), codes_linnaen[,1][i])] <- codes_linnaen[,4][i]}
  d_PCA$Population <- Population
  pops_PCA <- sapply(1:length(unique(d_PCA$Population)), function (x) as.matrix((t(d_PCA[d_PCA$Population == unique(d_PCA$Population)[x],-(47)]))))
  
  #get traits
  traits_PCA <- list()
  for(i in 1:ntraits_PCA){ 
    print(i)
    trait <- list()
    trait[[1]] <- as.numeric(pops_PCA[[1]][i,])
    for (j in 2:27){
      trait[[j]] <- as.numeric(pops_PCA[[j]][i,])
    } 
    traits_PCA[[i]] <- trait
  }
  
  #make cov matrix
  covBG <- matrix(nrow = ntraits_PCA, ncol = ntraits_PCA)
  for (i in 1:ntraits_PCA){
    print(i)
    for(j in 1:ntraits_PCA){
      trait1 <- traits_PCA[[i]]
      trait2 <- traits_PCA[[j]]
      trait1means <- sapply(1:length(trait1), function(x) mean(trait1[[x]]))
      trait2means <- sapply(1:length(trait2), function(x) mean(trait2[[x]]))
      trait1deviations <- sapply(1:length(trait1), function(x) trait1[[x]] - trait1means[x])
      trait2deviations <- sapply(1:length(trait2), function(x) trait2[[x]] - trait2means[x])
      deviationsProducts <- sapply(1:length(trait1), function(x) trait1deviations[[x]] * trait2deviations[[x]])
      sumALL <- 0
      for (k in (1:length(deviationsProducts))){
        sumALL <- sumALL + sum(deviationsProducts[[k]])
      }
      indiv <- sum(sapply(1:length(trait1), function(x) length(trait1[[x]])))
      covariance <- sumALL/(indiv - length(deviationsProducts))
      covBG[i,j] <- covariance
    }
  }
}


#do PCA on procrustes-transformed data using total covariance matrix
if(PCA_covBG != T){
  V <- eigen(covTotal)$vectors
  L <- eigen(covTotal)$values
} else{
  V <- eigen(covBG)$vectors
  L <- eigen(covBG)$values
}
VarExpl <- round(L / sum(L), 4) * 100
PC_Scores <- t(V) %*% t(t.data) 
if(normalize){
  PC_Scores <- PC_Scores / L^0.5 #set sd to 1
  PC_Scores <- PC_Scores - apply(PC_Scores, 1, mean) #mean center
}
PC_Scores <- PC_Scores[1:(nrow(PC_Scores)-7),]

if(!noPCA){
  d <- as.data.frame(t(PC_Scores))
} else {
  d <- as.data.frame(t.data)
}

ntraits <- length(d[1,])


Population <- rep("species", nrow(d))
for(i in 1:nrow(codes_linnaen)){Population[startsWith(rownames(d), codes_linnaen[,1][i])] <- codes_linnaen[,4][i]}
d$Population <- Population
if(!noPCA){
  pops <- sapply(1:length(unique(d$Population)), function (x) as.matrix((t(d[d$Population == unique(d$Population)[x],-(41)]))))
} else {
  pops <- sapply(1:length(unique(d$Population)), function (x) as.matrix((t(d[d$Population == unique(d$Population)[x],-(47)]))))
}
#get traits
traits <- list()
npop <- length(pops)
for(i in 1:ntraits){ 
  print(i)
  trait <- list()
  trait[[1]] <- as.numeric(pops[[1]][i,])
  for (j in 2:npop){
    trait[[j]] <- as.numeric(pops[[j]][i,])
  } 
  traits[[i]] <- trait
}

#make cov matrix
cov <- matrix(nrow = ntraits, ncol = ntraits)
for (i in 1:ntraits){
  print(i)
  for(j in 1:ntraits){
    trait1 <- traits[[i]]
    trait2 <- traits[[j]]
    trait1means <- sapply(1:length(trait1), function(x) mean(trait1[[x]]))
    trait2means <- sapply(1:length(trait2), function(x) mean(trait2[[x]]))
    trait1deviations <- sapply(1:length(trait1), function(x) trait1[[x]] - trait1means[x])
    trait2deviations <- sapply(1:length(trait2), function(x) trait2[[x]] - trait2means[x])
    deviationsProducts <- sapply(1:length(trait1), function(x) trait1deviations[[x]] * trait2deviations[[x]])
    sumALL <- 0
    for (k in (1:length(deviationsProducts))){
      sumALL <- sumALL + sum(deviationsProducts[[k]])
    }
    indiv <- sum(sapply(1:length(trait1), function(x) length(trait1[[x]])))
    covariance <- sumALL/(indiv - length(deviationsProducts))
    cov[i,j] <- covariance
  }
}
if(noPCA){cov <- nearPD(cov)$mat}

#compute matrix of population means
traitsHARVATI <- matrix(nrow=npop, ncol=ntraits)
for(i in 1:npop){
  for(j in 1:ntraits){
    traitsHARVATI[i,j] <- mean(as.numeric(pops[[i]][j,]))
  }
}
rownames(traitsHARVATI) <- gsub(x = gsub(x = gsub(x = unique(d$Population), 
                                         pattern =  " ", replacement =  "_"), 
                                         pattern = ")", replacement = ""),
                                         pattern = "\\(", replacement = "")

colnames(traitsHARVATI) <- colnames(d)[-(ntraits+1)]

rownames(cov) <- colnames(cov) <- colnames(d)[-(ntraits+1)]

##########################################
############ USEFUL FUNCTIONS ############
##########################################

#write means to file function
meansNexus <- function (characterMatrix, outputFilePath, traitsOnColumns = T){
  if (traitsOnColumns) {ntax = dim(characterMatrix)[1]; nchar = dim(characterMatrix)[2]}
  else {ntax = dim(characterMatrix)[2]; nchar = dim(characterMatrix)[1]}
  sink(outputFilePath, append=F)
  cat("#NEXUS\n\nBegin data;\n\tDimensions ntax=", ntax, "nchar=", nchar, ";\n\tFormat datatype=Continuous missing=?;\n\tMatrix\n\n")
  for(i in 1:length(characterMatrix[,1])) {cat(rownames(characterMatrix)[i], characterMatrix[i,], "\n")}
  cat(";\nEnd;")
  sink()
}

ratesTSV <- function (rateMatrix, outputFilePath) {
  if(!isSymmetric(rateMatrix)) {stop("Error: This matrix is not symmetric.")}
  sink(outputFilePath, append=F)
  for(i in 1:length(rateMatrix[,1])) {
    cat(paste0(as.vector(rateMatrix[i, (1:(length(rateMatrix[i,])-1))])), sep = "\t")
    cat("\t")
    cat(rateMatrix[i, (length(rateMatrix[i,]))], "\n\n", sep = "")
  }
  sink()
}

if(!noPCA){
  if(PCA_covBG){
    meansNexus(traitsHARVATI, paste0("/Users/nikolai/data/Harvati_", ifelse(normalize, "normalized_", ""), ifelse(collapseHomo, "collapseHomo_", ""), ifelse(collapseSubsp, "collapseSubsp_", ""), "bgCovPCA.nex"))
    ratesTSV(cov, paste0("/Users/nikolai/data/Harvati_", ifelse(normalize, "normalized_", ""), ifelse(collapseHomo, "collapseHomo_", ""), ifelse(collapseSubsp, "collapseSubsp_", ""), "_bgCovPCAbgRM.tsv"))
    ratesTSV(cov2cor(cov), paste0("/Users/nikolai/data/Harvati_", ifelse(normalize, "normalized_", ""), ifelse(collapseHomo, "collapseHomo_", ""), ifelse(collapseSubsp, "collapseSubsp_", ""), "_bgCovPCAbgRM_corr.tsv"))
  } else {
    meansNexus(traitsHARVATI, paste0("/Users/nikolai/data/Harvati_", ifelse(normalize, "normalized_", ""), ifelse(collapseHomo, "collapseHomo_", ""), ifelse(collapseSubsp, "collapseSubsp_", ""), "totalCovPCA.nex"))
    ratesTSV(cov, paste0("/Users/nikolai/data/Harvati_", ifelse(normalize, "normalized_", ""), ifelse(collapseHomo, "collapseHomo_", ""), ifelse(collapseSubsp, "collapseSubsp_", ""), "totalCovPCAbgRM.tsv"))
    ratesTSV(cov2cor(cov), paste0("/Users/nikolai/data/Harvati_", ifelse(normalize, "normalized_", ""), ifelse(collapseHomo, "collapseHomo_", ""), ifelse(collapseSubsp, "collapseSubsp_", ""), "totalCovPCAbgRM_corr.tsv"))
  }
} else {
  meansNexus(traitsHARVATI, paste0("/Users/nikolai/data/Harvati_noPCA", ifelse(collapseHomo, "_homocollapse", "_homopops"), ".nex"))
  ratesTSV(cov, paste0("/Users/nikolai/data/Harvati_noPCA", ifelse(collapseHomo, "_homocollapse", "_homopops"), ".tsv"))
  ratesTSV(cov2cor(cov), paste0("/Users/nikolai/data/Harvati_noPCA_corr", ifelse(collapseHomo, "_homocollapse", "_homopops"), ".tsv"))
}

library(ape)
library(phangorn)
library(phytools)

if(PCA_covBG){
  trees1 <- read.tree("output/bgCovPCA_fixRM_c1.trees")
  trees2 <- read.tree("output/bgCovPCA_fixRM_c2.trees")
} else {
  trees1 <- read.tree("output/totalCovPCA_fixRM_c1.trees")
  trees2 <- read.tree("output/totalCovPCA_fixRM_c2.trees")
}
if(normalize){
  trees1 <- read.tree("output/totalCovPCA_fixRM_normalized_c1.trees")
  trees2 <- read.tree("output/totalCovPCA_fixRM_normalized_c2.trees")
}
if(originalTest){
  trees1 <- read.tree("output/testHarvati.trees")
  trees2 <- read.tree("output/testHarvati_2.trees")
}
if(collapseHomo){
  trees1 <- read.tree("/Users/nikolai/output/totalCovPCA_fixRM_normalized_collapseHomo_c2.trees")
  trees <- trees1
}
if(collapseSubsp){
  trees1 <- read.tree("output/totalCovPCA_fixRM_normalized_collapseSubsp_c2.trees")
  trees <- trees1
}
if(noPCA){
  trees1 <- read.tree("output/harvati_noPCA_c1.trees")
  trees2 <- read.tree("output/harvati_noPCA_c2.trees")
  njTree <- read.tree(file = paste0("output/empirical_neighborJoining.txt"))
  upgmaTree <- read.tree(file = paste0("output/empirical_upgma.txt"))
  mp_tree <- read.tree(file = "maximumParsimonyTree_empirical.txt")
  mp_tree_PCA <- read.tree(file = "maximumParsimonyTree_empirical_PCA.txt")
}
trees <- c(trees1, trees2)

tree <- maxCladeCred(trees)
tree <- consensus(trees, p = 0.5)
startsWith2 <- function(x, pre){apply(sapply(1:length(pre), function(n) startsWith(x, pre[n])), 1, any)}
tree <- root(tree, outgroup = tree$tip.label[startsWith2(tree$tip.label, c("Mac", "Pap", "Man"))], resolve.root = T, edgelabel = T)
njTree <- root(njTree, outgroup = njTree$tip.label[startsWith2(njTree$tip.label, c("Mac", "Pap", "Man"))], resolve.root = T, edgelabel = T)
mp_tree <- root(mp_tree, outgroup = mp_tree$tip.label[startsWith2(mp_tree$tip.label, c("Hom", "Pan", "Gor"))], resolve.root = T, edgelabel = T)
mp_tree_PCA <- root(mp_tree_PCA, outgroup = mp_tree_PCA$tip.label[startsWith2(mp_tree_PCA$tip.label, c("Mac", "Pap", "Man"))], resolve.root = T, edgelabel = T)

pp <- prop.clades(tree, trees, rooted = F)/length(trees)

par(mfrow = c(1,1))

# png(filename = "howellsMCC.png", width = 1500, height = 800)
plot(tree)
edgelabels(round(pp, 2), sapply(1:tree$Nnode + length(tree$tip.label), 
                                match, tree$edge[, 2]), bg = 1, frame = "none", adj = c(0,1.5))
plot(upgmaTree); title("upgma tree")
plot(njTree); title("neighbor joining tree")
dev.off()
par(mfrow = c(1,2))
plot(mp_tree); title("maximum parsimony tree")
plot(mp_tree_PCA); title("maximum parsimony tree with PCA transform")

plot(cophylo(mp_tree, mp_tree_PCA)); title("\nmax-pars tree                                                                                                                                                                                                 max-pars tree w/ pca")

# tree_norm <- tree
# plot(cophylo(tree, tree_norm))


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

clade_monophyly <- prop.part.df(trees, cutoff = 0.0001)
cladeName <- c("Macaca_fa", "Macaca_syl")
clade <- list(trees[[1]]$tip.label[startsWith2(trees[[1]]$tip.label, cladeName)], trees[[1]]$tip.label[!startsWith2(trees[[1]]$tip.label, cladeName)])
clade_ind <- which(sapply(1:length(clade_monophyly$cladeNames), function(x) isTRUE(all.equal(sort(clade_monophyly$cladeNames[[x]]), sort(clade[[1]]))) |
                                                       isTRUE(all.equal(sort(clade_monophyly$cladeNames[[x]]), sort(clade[[2]])))))
clade_monophyly$postProbs[clade_ind]
