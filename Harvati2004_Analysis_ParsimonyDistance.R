library(geomorph)
library(Matrix)

normalize <- T
noPCA <- F
PCA_covWG <- T
collapseHomo <- F
collapseSubsp <- F
addCentroidSize <- T
removeLastPCs <- F

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
if(collapseSubsp){
  codes_linnaen[,4] <- sapply(1:length(codes_linnaen[,4]), function(x) paste(strsplit(codes_linnaen[,4], " ")[[x]][1:2], collapse = " "))
}
#what is the 5th letter? subspecies? probably not, since TT is troglodytes troglodytes
# coded_genus <- substr(codes, 1, 2)
# genus_key <- cbind(c("HM", "MC", "MN", "PP", "PN", "GO"), c("Homo", "Macaca", "Mandrillus", "Papio", "Pan", "Gorilla"))
# sex <- substr(codes, 6, 6)
# linnaen <- as.character(unique(taxa[,1]))
# genus <- sapply(1:length(linnaen), function(x) strsplit(linnaen[x], split = " ")[[1]][1])
# species <- sapply(1:length(linnaen), function(x) strsplit(linnaen[x], split = " ")[[1]][2])
# subspecies <- sapply(1:length(linnaen), function(x) strsplit(linnaen[x], split = " ")[[1]][3])


#mean center all landmarks
for(i in 1:15){for(j in 1:3){data[i,j,] <- data[i,j,] - mean(data[i,j,])}}
t.data <- gpagen(data)$data
if(addCentroidSize){t.data$logCsize <- log(t.data$Csize)}
t.data <- as.matrix(subset(x = t.data, select = colnames(t.data)[colnames(t.data) != "Csize"]))
covTotal <- cov(t.data) #total covariance matrix
if(PCA_covWG == T){
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
  covWG <- matrix(nrow = ntraits_PCA, ncol = ntraits_PCA)
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
      covWG[i,j] <- covariance
    }
  }
}


#do PCA on procrustes-transformed data using total covariance matrix
if(PCA_covWG != T){
  V <- eigen(covTotal)$vectors
  L <- eigen(covTotal)$values
} else{
  V <- eigen(covWG)$vectors
  L <- eigen(covWG)$values
}
VarExpl <- round(L / sum(L), 4) * 100
PC_Scores <- t(V) %*% t(t.data) 
if(normalize){
  PC_Scores <- PC_Scores / L^0.5 #set sd to 1
  PC_Scores <- PC_Scores - apply(PC_Scores, 1, mean) #mean center
}
if(removeLastPCs){PC_Scores <- PC_Scores[1:(nrow(PC_Scores)-7),]}

if(!noPCA){
  d <- as.data.frame(t(PC_Scores))
} else {
  d <- as.data.frame(t.data)
}

ntraits <- length(d[1,])


Population <- rep("species", nrow(d))
for(i in 1:nrow(codes_linnaen)){Population[startsWith(rownames(d), codes_linnaen[,1][i])] <- codes_linnaen[,4][i]}
d$Population <- Population
if(!noPCA & removeLastPCs){
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




meanValueCodingDo <- T
numStartingTrees <- 50 #number of independent random starting trees
unchangedRun <- 500 #how long to keep the ratchet running without any improvement
traceOutput <- 1 #should the parsimony ratchet output stuff to screen?
l <- 1

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


#get true tree
tips <- rownames(traitsHARVATI)

#get traits
traits <- traitsHARVATI

#get pooled within-group covariance matrix
wgCM <- cov <- as.matrix(cov)
#wgCM <- as.matrix(read.table(file = paste0("rates/", fileName, "_rates.tsv")))

#get simulated individuals      
numPops <- length(pops)

#get discrete traits via divergence coding
divergenceCodedTraits <- matrix(data = 0, nrow = numPops, ncol = ntraits)
rownames(divergenceCodedTraits) <- rownames(traits)
pval <- 0.05

for(trait in 1:ntraits){
  divergenceMatrix <- matrix(data = 0, nrow = numPops, ncol = numPops)
  for(poprow in 2:numPops){
    for(popcol in 1:(poprow - 1)){
      poprowindiv <- as.numeric(pops[[poprow]][trait,])
      popcolindiv <- as.numeric(pops[[popcol]][trait,])
      significance <- t.test(poprowindiv, popcolindiv)$p.value < pval
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
write.tree(bestTrees, paste0("output/maximumParsimonyTree_empirical", ifelse(noPCA, "_RAW", "_PCA"), "_divergenceCoding.txt"))

#calculate average RF distance from the true tree
# if(class(bestTrees) == "multiPhylo"){avgRFdistDC <- mean(sapply(1:length(bestTrees), function(x) dist.topo(unroot(trueTree), bestTrees[[x]])[1]))
# } else {avgRFdistDC <- dist.topo(unroot(trueTree), bestTrees)[1]}

if(meanValueCodingDo){
  #do mean value coding
  traitMeans <- sapply(1:length(traits[1,]), function(x) mean(traits[,x]))
  meanBinarizedTraits <- t(sapply(1:length(traits[,1]), function(x) traits[x,] > traitMeans))
  rownames(meanBinarizedTraits) <- rownames(traits)
  meanBinarizedTraits[meanBinarizedTraits] <- "L"
  meanBinarizedTraits[meanBinarizedTraits == "FALSE"] <- "S"
  mBT.pD <- phyDat(meanBinarizedTraits, type = "USER", levels = c("L", "S"))
  mpTrees <- list()
  for(sepStart in 1:numStartingTrees){
    mpTree <- pratchet(mBT.pD, maxit = 1e5, method = "sankoff", k = unchangedRun, start = rtree(n = length(tips), tip.label = sample(tips)), trace = traceOutput, all = T)
    mpTreeScore <- parsimony(tree = mpTree, data = mBT.pD, method = "sankoff")
    mpTrees[[sepStart]] <- c(mpTree, mpTreeScore) 
  }
  mpScores <- sapply(1:length(mpTrees), function(x) mpTrees[[x]][[2]][1])
  mostPars <- which(mpScores == min(mpScores))
  bestTreesRough_MVC <- lapply(1:length(mostPars), function(x) mpTrees[[mostPars[x]]][[1]])
  bestTrees_MVC <- bestTreesRough_MVC[[1]]
  if(length(bestTreesRough_MVC) > 1){for(treeGroup in 2:length(bestTreesRough_MVC)){bestTrees_MVC <- c(bestTrees_MVC, bestTreesRough_MVC[[treeGroup]])}}
  if(class(bestTrees_MVC) != "phylo"){bestTrees_MVC <- unique(bestTrees_MVC)}
  # if(class(bestTrees) == "multiPhylo"){avgRFdistMVC <- mean(sapply(1:length(bestTrees), function(x) dist.topo(unroot(trueTree), bestTrees[[x]])[1]))
  # } else {avgRFdistMVC <- dist.topo(unroot(trueTree), bestTrees)[1]}
}
write.tree(bestTrees, paste0("output/maximumParsimonyTree_empirical", ifelse(noPCA, "_RAW", "_PCA"), "_MVC.txt"))

#get mahalonibis distance matrix
mahDistMat <- mahaMatrix(traits, wgCM, squared = F)
head(mahDistMat)

#compute trees
njTree <- nj(mahDistMat)
upgmaTree <- upgma(mahDistMat)
# rfDistNJ <- dist.topo(x = trueTree, y = njTree, method = "PH85")
# rfDistUPGMA <- dist.topo(x = trueTree, y = upgmaTree, method = "PH85")

write.tree(bestTrees, file = paste0("output/empirical_divergenceCodingParsimony.txt"))
write.tree(njTree, file = paste0("output/empirical_neighborJoining.txt"))
write.tree(upgmaTree, file = paste0("output/empirical_upgma.txt"))

