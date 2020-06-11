setwd(dir = "/Volumes/1TB/Harvati/")

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
library(BAMMtools)

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

nreps <- 100
traitNumIter <- c((46-1)*2^(0:4)+1)
jenks_exp <- c(2,4,8)

ssqm <- function(nums){sum((nums-mean(nums))^2)}
closest <- function(num, nums, eps = 1E-6){
  diffs <- abs(nums-num)
  which(diffs < (min(diffs) + eps) )
}

for(i in 1:nreps){
  
  cat(paste0("\n", i, ":"))
  
  for(j in 1:length(traitNumIter)){
    
    cat(paste0(" ", j, " "))
    
    ntraits <- traitNumIter[j]
    replicateNum <- i
    
    fileName <- paste0("simulations_", ntraits, "traits_", "replicate_", replicateNum)
    
    #get true tree
    trueTree <- read.nexus(file = paste0("trees/", fileName, "_tree.nex"))
    tips <- trueTree$tip.label
    ntips <- length(tips)
    
    #get traits
    traits <- strsplit(readLines(paste0("data/", fileName, "_traits.nex")), split = "\t")
    traits <- traits[8:(length(traits)-2)]
    rwn <- sapply(1:length(traits), function(x) strsplit(traits[[x]], split = " ")[[1]][1])
    traits <- t(as.matrix(sapply(1:length(traits), function(x) as.numeric(strsplit(traits[[x]], split = " ")[[1]][-1]))))
    rownames(traits) <- rwn
    
    #get pooled within-group covariance matrix
    cov <- as.matrix(read.table(paste0("rates/", fileName ,"_rates.tsv")))
    cov <- make_symmetric(cov)
    
    #apply PC transform to simulated pop means
    
    V <- eigen(cov)$vectors
    L <- eigen(cov)$values
    VarExpl <- round(L / sum(L), 4) * 100
    PC_Scores <- t(V) %*% t(traits) 
    traits_PCA <- t(PC_Scores)
    
    allTraits <- abind(traits, traits_PCA, along = 3)
    
    for(type in 1:2){
      
      traits <- allTraits[,,type]
      njenks_range <- 1:ntips #specify range of break points from none to everyone in a unique category
      all_jenks_traits <- array(0, dim = c(dim(traits), length(njenks_range)))
      all_GVFs <- matrix(0, nrow = length(njenks_range), ncol = ntraits)
      
      #iterate over all numbers of categories
      for(njenks in njenks_range){ 
        
        if(njenks > 1){
          #observation: first break is always the set min, so ignore it
        jenksBreaks <- sapply(1:ntraits, function(trait) getJenksBreaks(traits[,trait], k = njenks)[-1])
        traitMeans <- sapply(1:ntraits, function(trait) mean(traits[,trait]))
        traitMags <- sapply(1:ntraits, function(trait) mean(abs(traits[,trait])))
        underflowHack <- 10^round(log(abs(1/traitMags)) / log(10))
        SDAM <- sapply(1:ntraits, function(trait) ssqm(traits[,trait] * underflowHack[trait]))
        
        if(njenks > 2){
          jenksCodedTraits <- sapply(1:ntraits, function(trait) sapply(1:ntips, function(tip) sum(traits[tip,trait] >= jenksBreaks[,trait]))) + 1
        } else {
          jenksCodedTraits <- sapply(1:ntraits, function(trait) sapply(1:ntips, function(tip) sum(traits[tip,trait] >= jenksBreaks[trait]))) + 1
        }
        rownames(jenksCodedTraits) <- rownames(traits)
        
        SDCM <- sapply(1:ntraits, function(trait) 
          sum(sapply(1:njenks, function(jcat) ssqm(traits[jenksCodedTraits[,trait] == jcat,trait] * underflowHack[trait]))))
        
        #compute goodness of variance fit
        GVF <- (SDAM - SDCM) / SDAM
        
        } else if(njenks == 1){
        
            jenksCodedTraits <- matrix(1, nrow(traits), ncol(traits), dimnames = list(rownames(traits)))
            GVF <- rep(0, ntraits)
        
        }
        
        all_jenks_traits[,,njenks] <- jenksCodedTraits
        all_GVFs[njenks,] <- GVF
        
      }
      
      #find the expected number of categories at each GVF threshold
      GVF_range <- 1:999/1000
      exp_ncats <- sapply(GVF_range, function(GFVt) mean(sapply(1:ntraits, function(trait) min(which(GFVt < all_GVFs[,trait])))))
      GVF_thresholds <- sapply(jenks_exp, function(ncats) GVF_range[median(closest(ncats, exp_ncats))])
      
      ncats_per_exp <- sapply(GVF_thresholds, function(GVFt) sapply(1:ntraits, function(trait) min(which(all_GVFs[,trait] > GVFt))))
      
      for(jenks_exp_cats in 1:length(jenks_exp)){
        
        which_ncats <- ncats_per_exp[,jenks_exp_cats]
        unique_ncats <- unique(which_ncats)
        sink(paste0("data_discrete/", fileName, ifelse(type == 1, "_raw", "_PCs"), "_jenks_", jenks_exp[jenks_exp_cats], "_expcats.txt"))
        cat(unique_ncats)
        sink()
        
        for(uncats in unique_ncats){
          jenks_traits <- as.matrix(all_jenks_traits[,which_ncats == uncats,uncats])
          rownames(jenks_traits) <- rownames(traits)

          traitsPath_jenksCoding <- paste0("data_discrete/", fileName, ifelse(type == 1, "_raw", "_PCs"), 
                                           "_jenks_", jenks_exp[jenks_exp_cats], "_expcats_", uncats, "_cats_traits.tsv")
          writeTSV(jenks_traits, traitsPath_jenksCoding)
        }
        
      }
      

    }
  }
}
