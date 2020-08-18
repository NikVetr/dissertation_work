setwd("/Volumes/1TB/Harvati_Empirical/")
library(geomorph)
library(Matrix)

meanValueCodingDo <- T
numStartingTrees <- 50 #number of independent random starting trees
unchangedRun <- 500 #how long to keep the ratchet running without any improvement
traceOutput <- 1 #should the parsimony ratchet output stuff to screen?
l <- 1
HomoPops_justSpecies <- F
collapseHomo <- !HomoPops_justSpecies

library(abind)
library(ape)
library(coda)
library(rwty)
library(phangorn)
library(lattice)
library(latticeExtra)
library(MCMCpack)
library(readtext)
library(StatMatch)
library(mvtnorm)

maha <- function(x1, x2, vcvm, squared = FALSE){
  squaredDist <- t(x1-x2) %*% solve(vcvm) %*% (x1-x2)
  if(!squared){
    squaredDist
  } else {
    squaredDist^0.5
  }
} 
mahaMatrix <- function(data, vcvm, squared = FALSE){
  mat <- matrix (nrow = nrow(data), ncol = nrow(data), 0)
  for(i in 1:nrow(mat)){
    for(j in i:nrow(mat)){
      mat[i,j] <- maha(data[i,], data[j,], vcvm, squared = squared)
    }
  }
  mat <- mat + t(mat)
  rownames(mat) <- colnames(mat) <-rownames(data)
  mat
}
convNum2Str <- function(nums, key){
  sapply(1:length(nums), function(x) key[nums[x]])
}
wagnerCostMatrix <- function(dim){
  mat <- matrix(0, nrow = dim, ncol = dim)
  for(i in 1:(dim-1)){
    for(j in (i+1):dim){
      mat[i,j] <- abs(i-j)
    }
  }
  mat + t(mat)
}
undoMeansNexus <- function(fileString){
  nexdata <- readLines(con = fileString, warn = F)
  nexdata <- nexdata[(which(nexdata == "\tMatrix")+2):(which(nexdata == "End;")-2)]
  names <- unlist(lapply(nexdata, function(l) strsplit(l, split = " ")[[1]][1]))
  data <- do.call(rbind, lapply(nexdata, function(l) as.numeric(strsplit(l, split = " ")[[1]][-1])))
  rownames(data) <- names
  return(data)
}

#get true tree
tips <- rownames(traitsHARVATI)

#get traits
traits <- traitsHARVATI

#get discrete traits via divergence coding
data_subset <- c(95,99,100)[2]
sex <- "males"
collapseHomo = T
load(paste0("/Users/nikolai/data/Harvati_", ifelse(collapseHomo, "no", ""),"HomoPops_justSpecies_", sex,"_PCA", data_subset,"_PopData"))
numPops <- length(pops)
ntraits <- nrow(pops$Homo_sapiens)

divergenceCodedTraits <- matrix(data = 0, nrow = numPops, ncol = ntraits)
colnames(divergenceCodedTraits) <- paste0("PC", 1:ntraits)
rownames(divergenceCodedTraits) <- names(pops)


pval <- 0.05

#correct small diffs... lol
noPCA = F
csd <- ifelse(noPCA, F, T)
for(trait in 1:ntraits){
  divergenceMatrix <- matrix(data = 0, nrow = numPops, ncol = numPops)
  for(poprow in 2:numPops){
    for(popcol in 1:(poprow - 1)){
      poprowindiv <- as.numeric(pops[[poprow]][trait,])
      popcolindiv <- as.numeric(pops[[popcol]][trait,])
      if(csd){
        if(length(popcolindiv) > 1){
          popcolindiv <- popcolindiv - poprowindiv[1]
        }
        if(length(poprowindiv) > 1){
          poprowindiv <- poprowindiv - poprowindiv[1]
        }
      }
      #test for significance
      if(length(poprowindiv) > 1 & length(popcolindiv) > 1){
        significance <- t.test(poprowindiv, popcolindiv)$p.value < pval
      } else {
        if(length(poprowindiv) == 1){
          significance <- abs(poprowindiv - mean(popcolindiv)) / sd(popcolindiv) > 2
        }
        if(length(popcolindiv) == 1){ #this is really just for the female papio hamadryas papio
          significance <- abs(popcolindiv - mean(poprowindiv)) / sd(poprowindiv) > 2
        }
      }
      # print(t.test(poprowindiv, popcolindiv)$p.value)
      if(significance){
        if(mean(popcolindiv) > mean(poprowindiv)){
          divergenceMatrix[poprow, popcol] <- 1
        } else {
          divergenceMatrix[poprow, popcol] <- -1
        }
      } else {
        divergenceMatrix[poprow, popcol] <- 0
      }
    }
  }
  divergenceMatrix <- divergenceMatrix - t(divergenceMatrix)
  divergenceCodedTraits[,trait] <- sapply(1:numPops, function(x) sum(divergenceMatrix[,x]))
}

for(foo in 1:length(divergenceCodedTraits[1,])){
  divergenceCodedTraits[,foo] <- divergenceCodedTraits[,foo] + abs(min(divergenceCodedTraits[,foo])) + 1
}

save(divergenceCodedTraits, file = paste0("/Users/nikolai/data/Harvati_", ifelse(collapseHomo, "no", ""),"HomoPops_justSpecies_", sex,"_PCA", data_subset,"_DivCodedTraits"))

load(paste0("/Users/nikolai/data/Harvati_", ifelse(collapseHomo, "no", ""),"HomoPops_justSpecies_", "males","_PCA", data_subset,"_DivCodedTraits"))
male_divergence_coded_traits <- divergenceCodedTraits
load(paste0("/Users/nikolai/data/Harvati_", ifelse(collapseHomo, "no", ""),"HomoPops_justSpecies_", "females","_PCA", data_subset,"_DivCodedTraits"))
female_divergence_coded_traits <- divergenceCodedTraits

addFemaleNeandertal = T
if(addFemaleNeandertal){
  fictitious_female_neandertal <- female_divergence_coded_traits["Homo_sapiens",] 
  female_divergence_coded_traits <- rbind(Homo_sapiens = female_divergence_coded_traits[1,], Homo_neanderthalensis = fictitious_female_neandertal, female_divergence_coded_traits[2:nrow(female_divergence_coded_traits),])
}

divergenceCodedTraits <- cbind(female_divergence_coded_traits, male_divergence_coded_traits)

#write distance matrix to file
eucl <- F; univSlide <- T
d <- divergenceCodedTraits
if(eucl){distmat <- stats::dist(d)}
if(univSlide){
  distmat_univ <- lapply(1:length(d[1,]), function(tr) stats::dist(d[,tr]))
  distmat <- distmat_univ[[1]]
  for(i in 2:length(distmat_univ)){
    distmat <- distmat + distmat_univ[[i]]
  }
}

writeDist(distmat, file = paste0("output/justSpecies/parsDistances_empirical_",
                                 ifelse(collapseHomo, "no", ""), "HomoPops_justSpecies_", "PCA", data_subset, "_divergenceCoding.dist"))

#conver to phyDat format
DCT.pD <- phyDat(divergenceCodedTraits, type = "USER", levels = 1:(numPops*2-1))

#use the parsimony ratchet to search for most parsimonious trees a number of times
mpTrees <- list()
for(sepStart in 1:numStartingTrees){
  mpTree <- pratchet(DCT.pD, maxit = 1e5, method = "sankoff", k = unchangedRun, start = rtree(n = length(tips), tip.label = sample(tips)), trace = traceOutput, all = T, cost = wagnerCostMatrix(numPops*2-1))
  mpTreeScore <- parsimony(tree = mpTree, data = DCT.pD, method = "sankoff", cost = wagnerCostMatrix(numPops*2-1))
  mpTrees[[sepStart]] <- c(mpTree, mpTreeScore) 
}

#find unique most parsimonious trees
mpScores <- sapply(1:length(mpTrees), function(x) mpTrees[[x]][[2]][1])
mostPars <- which(mpScores == min(mpScores))
bestTreesRough <- lapply(1:length(mostPars), function(x) mpTrees[[mostPars[x]]][[1]])
bestTrees <- bestTreesRough[[1]]
if(length(bestTreesRough) > 1){for(treeGroup in 2:length(bestTreesRough)){bestTrees <- c(bestTrees, bestTreesRough[[treeGroup]])}}
if(class(bestTrees) != "phylo"){bestTrees <- unique(bestTrees)}
write.tree(bestTrees, paste0("output/justSpecies/maximumParsimonyTree_empirical_",
                             ifelse(collapseHomo, "no", ""), "HomoPops_justSpecies_", "PCA", data_subset, "_divergenceCoding.txt"))

#calculate average RF distance from the true tree
# if(class(bestTrees) == "multiPhylo"){avgRFdistDC <- mean(sapply(1:length(bestTrees), function(x) dist.topo(unroot(trueTree), bestTrees[[x]])[1]))
# } else {avgRFdistDC <- dist.topo(unroot(trueTree), bestTrees)[1]}

#get mahalonibis distance matrix
#read in the data
percPCA <- "99"
females <- undoMeansNexus(fileString = paste0("data2_neanF/Harvati_noHomoPops_justSpecies_females_PCA", percPCA, ".nex"))
males <- undoMeansNexus(fileString = paste0("data2_neanF/Harvati_noHomoPops_justSpecies_males_PCA", percPCA, ".nex"))


if(all.equal(rownames(females), rownames(males))){
  traitsHARVATI <- cbind(females, males)
}

mahDistMat <- mahaMatrix(traitsHARVATI, diag(ncol(traitsHARVATI)), squared = F)
head(mahDistMat)

#compute trees
njTree <- nj(mahDistMat)
upgmaTree <- upgma(mahDistMat)
# rfDistNJ <- dist.topo(x = trueTree, y = njTree, method = "PH85")
# rfDistUPGMA <- dist.topo(x = trueTree, y = upgmaTree, method = "PH85")

write.tree(njTree, file = paste0("output/justSpecies/empirical_PCA", percPCA ,"_neighborJoining.txt"))
write.tree(upgmaTree, file = paste0("output/justSpecies/empirical_PCA", percPCA ,"_upgma.txt"))

