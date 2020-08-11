library(geomorph)
library(Matrix)
library(phytools)
library(phangorn)
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
library(classInt)

ssqm <- function(nums){sum((nums-mean(nums))^2)}
closest <- function(num, nums, eps = 1E-6){
  diffs <- abs(nums-num)
  which(diffs < (min(diffs) + eps) )
}

convNum2Str <- function(nums, key){
  sapply(1:length(nums), function(x) key[nums[x]])
}
meansNexus <- function (characterMatrix, outputFilePath, traitsOnColumns = T, dataType = "Continuous"){
  if (traitsOnColumns) {ntax = dim(characterMatrix)[1]; nchar = dim(characterMatrix)[2]}
  else {ntax = dim(characterMatrix)[2]; nchar = dim(characterMatrix)[1]}
  sink(outputFilePath, append=F)
  cat(paste0("#NEXUS\n\nBegin data;\n\tDimensions ntax=", ntax, " nchar=", nchar, ";\n\tFormat datatype=", 
             dataType, " missing=?", ifelse(dataType == "Standard", " SYMBOLS=\"01\"", "") ,";\n\tMatrix\n\n"))
  for(i in 1:length(characterMatrix[,1])) {cat(rownames(characterMatrix)[i], characterMatrix[i,], "\n")}
  cat(";\nEnd;")
  sink()
}
writeTSV <- function (object, outputFilePath) {
  rwn <- rownames(object)
  sink(outputFilePath, append=F)
  for(i in 1:length(object[,1])) {
    cat(rwn[i], "\t")
    if(length(object[i,]) > 1){
      cat(paste0(as.vector(object[i, (1:(length(object[i,])-1))])), sep = "\t")
      cat("\t")
    }
    cat(object[i, (length(object[i,]))], "\n\n", sep = "")
  }
  sink()
}
make_symmetric <- function(mat, method = c("upper.tri", "lower.tri", "mean")[1]){
  mat_diag <- diag(diag(mat))
  if(method == "upper.tri"){
    mat_upper <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
    mat_upper[upper.tri(mat_upper)] <- mat[upper.tri(mat)]
    mat_symm <- mat_diag + mat_upper + t(mat_upper)
  } else if (method == "lower.tri"){
    mat_lower <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
    mat_lower[lower.tri(mat_lower)] <- mat[lower.tri(mat)]
    mat_symm <- mat_diag + mat_lower + t(mat_lower)
  } else if (method == "mean"){
    mat_symm <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
    mat_symm[upper.tri(mat_symm)] <- mat[upper.tri(mat)]
    mat_symm[lower.tri(mat_symm)] <- mat[lower.tri(mat)]
    mat_symm <- mat_diag + (mat_symm + t(mat_symm))/2
  } else {stop("you did not specify a valid method")}
  rownames(mat_symm) <- rownames(mat)
  colnames(mat_symm) <- colnames(mat)
  if(all(mat_symm - t(mat_symm) == 0)){
    return(mat_symm)
  } else {stop("sorry something went wrong")}
}
setwd("~")


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

convNum2Str <- function(nums, key){
  sapply(1:length(nums), function(x) key[nums[x]])
}
meansNexus <- function (characterMatrix, outputFilePath, traitsOnColumns = T, dataType = "Continuous"){
  if (traitsOnColumns) {ntax = dim(characterMatrix)[1]; nchar = dim(characterMatrix)[2]}
  else {ntax = dim(characterMatrix)[2]; nchar = dim(characterMatrix)[1]}
  sink(outputFilePath, append=F)
  cat(paste0("#NEXUS\n\nBegin data;\n\tDimensions ntax=", ntax, " nchar=", nchar, ";\n\tFormat datatype=", 
             dataType, " missing=?", ifelse(dataType == "Standard", " SYMBOLS=\"01\"", "") ,";\n\tMatrix\n\n"))
  for(i in 1:length(characterMatrix[,1])) {cat(rownames(characterMatrix)[i], characterMatrix[i,], "\n")}
  cat(";\nEnd;")
  sink()
}
writeTSV <- function (object, outputFilePath) {
  rwn <- rownames(object)
  sink(outputFilePath, append=F)
  for(i in 1:length(object[,1])) {
    cat(rwn[i], "\t")
    cat(paste0(as.vector(object[i, (1:(length(object[i,])-1))])), sep = "\t")
    cat("\t")
    cat(object[i, (length(object[i,]))], "\n\n", sep = "")
  }
  sink()
}
make_symmetric <- function(mat, method = c("upper.tri", "lower.tri", "mean")[1]){
  mat_diag <- diag(diag(mat))
  if(method == "upper.tri"){
    mat_upper <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
    mat_upper[upper.tri(mat_upper)] <- mat[upper.tri(mat)]
    mat_symm <- mat_diag + mat_upper + t(mat_upper)
  } else if (method == "lower.tri"){
    mat_lower <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
    mat_lower[lower.tri(mat_lower)] <- mat[lower.tri(mat)]
    mat_symm <- mat_diag + mat_lower + t(mat_lower)
  } else if (method == "mean"){
    mat_symm <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
    mat_symm[upper.tri(mat_symm)] <- mat[upper.tri(mat)]
    mat_symm[lower.tri(mat_symm)] <- mat[lower.tri(mat)]
    mat_symm <- mat_diag + (mat_symm + t(mat_symm))/2
  } else {stop("you did not specify a valid method")}
  rownames(mat_symm) <- rownames(mat)
  colnames(mat_symm) <- colnames(mat)
  if(all(mat_symm - t(mat_symm) == 0)){
    return(mat_symm)
  } else {stop("sorry something went wrong")}
}



traits <- traitsHARVATI

V <- eigen(cov)$vectors
L <- eigen(cov)$values
VarExpl <- round(L / sum(L), 4) * 100
PC_Scores <- t(V) %*% t(traits) 
traits_PCA <- t(PC_Scores)

fileName <- paste0("harvati_empirical", ifelse(collapseHomo, "_nohomopops", ""))
allTraits <- abind::abind(traits, traits_PCA, along = 3)
ntraits <- 46

for(type in 1:2){
  print(paste0("traits of the ", type, " kind"))
  traits <- allTraits[,,type]
  njenks_range <- 1:ntips #specify range of break points from none to everyone in a unique category
  all_jenks_traits <- array(0, dim = c(dim(traits), length(njenks_range)))
  all_GVFs <- matrix(0, nrow = length(njenks_range), ncol = ntraits)
  
  #iterate over all numbers of categories
  for(njenks in njenks_range){ 
    
    cat(paste0(njenks, " "))
    
    if(njenks > 1){
      
      traitMeans <- sapply(1:ntraits, function(trait) mean(traits[,trait]))
      traitMags <- sapply(1:ntraits, function(trait) mean(abs(traits[,trait])))
      underflowHack <- 10^round(log(abs(1/traitMags)) / log(10))
      SDAM <- sapply(1:ntraits, function(trait) ssqm(traits[,trait] * underflowHack[trait]))
      
      #having weird behavior with PC 41... underflow error? implement ad hoc solution... also 42 & 43
      jenksCodedTraits <- sapply(1:ntraits, function(trait) findCols(classInt::classIntervals(var = traits[,trait] - traits[1,trait], style = "jenks", n = njenks)))
      
      rownames(jenksCodedTraits) <- rownames(traits)
      
      SDCM <- sapply(1:ntraits, function(trait) 
        sum(sapply(1:njenks, function(jcat) ssqm(traits[jenksCodedTraits[,trait] == jcat,trait] * underflowHack[trait]))))
      
      #compute goodness of variance fit
      GVF <- (SDAM - SDCM) / SDAM
      
    } else if(njenks == 1){
      
      jenksCodedTraits <- matrix(1, nrow(traits), ncol(traits), dimnames = list(rownames(traits)))
      GVF <- rep(0, ntraits)
      
    }
    
    all_jenks_traits[,,njenks] <- jenksCodedTraits - 1
    all_GVFs[njenks,] <- GVF
    
  }
  
  #find the expected number of categories at each GVF threshold
  GVF_range <- 1:999/1000
  exp_ncats <- sapply(GVF_range, function(GFVt) mean(sapply(1:ntraits, function(trait) min(which(GFVt < all_GVFs[,trait])))))
  GVF_thresholds <- sapply(jenks_exp, function(ncats) GVF_range[median(closest(ncats, exp_ncats))])
  
  ncats_per_exp <- sapply(GVF_thresholds, function(GVFt) sapply(1:ntraits, function(trait) min(which(all_GVFs[,trait] > GVFt))))
  
  print("\nnow to find write the data...\n\n")
  for(jenks_exp_cats in 1:length(jenks_exp)){
    print(jenks_exp_cats)
    which_ncats <- ncats_per_exp[,jenks_exp_cats]
    unique_ncats <- unique(which_ncats)
    sink(paste0("data/", fileName, ifelse(type == 1, "_raw", "_PCs"), "_jenks_", jenks_exp[jenks_exp_cats], "_expcats.txt"))
    cat(unique_ncats)
    sink()
    
    for(uncats in unique_ncats){
      cat(paste0(uncats, " "))
      jenks_traits <- as.matrix(all_jenks_traits[,which_ncats == uncats,uncats])
      rownames(jenks_traits) <- rownames(traits)
      
      traitsPath_jenksCoding <- paste0("data/", fileName, ifelse(type == 1, "_raw", "_PCs"), 
                                       "_jenks_", jenks_exp[jenks_exp_cats], "_expcats_", uncats, "_cats_traits.tsv")
      writeTSV(jenks_traits, traitsPath_jenksCoding)
    }
  }
}

nrun <- 2
for(transData in c("raw", "PCs")){
  
  for(jenks_exp_cats in 1:length(jenks_exp)){
    
    ncats <- as.numeric(strsplit(readLines(paste0("data/", fileName, "_", transData, "_jenks_", jenks_exp[jenks_exp_cats], "_expcats.txt")), split = " ")[[1]])
    ncats <- sort(ncats)
    
    for (l in 1:nrun) {
      
      tempScript <- readLines("~/scripts/metascript_discrete_jenks.rev")
      tempScript[1] <- paste0("setwd(\"", getwd(), "\") ")
      tempScript[2] <- paste0(tempScript[2], l)
      tempScript[3] <- paste0("replicateNum <- ", "\"noReplicate\"")
      tempScript[5] <- paste0("coding <- \"", "jenks", ifelse(collapseHomo, "_nohomopops", ""), "\"")
      tempScript[6] <- paste0("dataFileName = ", paste0("v(", paste0("\"", fileName, "_", transData, "_jenks_", jenks_exp[jenks_exp_cats], "_expcats_", ncats, "_cats_traits.tsv\"", collapse = ", "), ")"))
      tempScript[7] <- paste0("transform = \"", transData, "\"")
      tempScript[8] <- paste0("ncats = v(", paste0(ncats, collapse = ", "), ")")
      tempScript[9] <- paste0("jenks_exp_cats = ", jenks_exp[jenks_exp_cats])
      tempScript[10] <- paste0("total_ntraits = ", ntraits)
      
      
      writeLines(tempScript, paste0("scripts/", fileName, "_", transData, "_jenks_", jenks_exp[jenks_exp_cats], "_expcats_script_run_", l, ".Rev"))
    }
  }
}

 
