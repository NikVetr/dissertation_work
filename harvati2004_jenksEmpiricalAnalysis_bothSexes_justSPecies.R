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
undoMeansNexus <- function(fileString){
  nexdata <- readLines(con = fileString, warn = F)
  nexdata <- nexdata[(which(nexdata == "\tMatrix")+2):(which(nexdata == "End;")-2)]
  names <- unlist(lapply(nexdata, function(l) strsplit(l, split = " ")[[1]][1]))
  data <- do.call(rbind, lapply(nexdata, function(l) as.numeric(strsplit(l, split = " ")[[1]][-1])))
  rownames(data) <- names
  return(data)
}
ssqm <- function(nums){sum((nums-mean(nums))^2)}
closest <- function(num, nums, eps = 1E-6){
  diffs <- abs(nums-num)
  which(diffs < (min(diffs) + eps) )
}

setwd("/Volumes/1TB/Harvati_Empirical/")

homopops = F
collapseHomo = T
percPCA <- "99"
if(homopops){
  females <- undoMeansNexus(fileString = paste0("data2_neanF/Harvati_HomoPops_females_PCA", percPCA, ".nex"))
  males <- undoMeansNexus(fileString = paste0("data2_neanF/Harvati_HomoPops_males_PCA", percPCA, ".nex"))
} else {
  females <- undoMeansNexus(fileString = paste0("data2_neanF/Harvati_noHomoPops_justSpecies_females_PCA", percPCA, ".nex"))
  males <- undoMeansNexus(fileString = paste0("data2_neanF/Harvati_noHomoPops_justSpecies_males_PCA", percPCA, ".nex"))
}

if(all.equal(rownames(females), rownames(males))){
  traitsHARVATI <- cbind(females, males)
}

traits <- traitsHARVATI
ntips = nrow(traits)
ntraits = ncol(traits)

fileName <- paste0("harvati_empirical", ifelse(collapseHomo, "_noHomoPops_justSpecies", ""))
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
jenks_exp = 4
GVF_range <- 1:999/1000
exp_ncats <- sapply(GVF_range, function(GFVt) mean(sapply(1:ntraits, function(trait) min(which(GFVt < all_GVFs[,trait])))))
GVF_thresholds <- sapply(jenks_exp, function(ncats) GVF_range[median(closest(ncats, exp_ncats))])

ncats_per_exp <- sapply(GVF_thresholds, function(GVFt) sapply(1:ntraits, function(trait) min(which(all_GVFs[,trait] > GVFt))))

print("\nnow to find write the data...\n\n")
jenks_exp_cats = 1
print(jenks_exp_cats)
which_ncats <- ncats_per_exp[,jenks_exp_cats]
unique_ncats <- sort(unique(which_ncats))

sink(paste0("data/", fileName, "_PCA", percPCA, "_jenks_", jenks_exp[jenks_exp_cats], "_expcats.txt"))
cat(unique_ncats)
sink()

for(uncats in unique_ncats){
  cat(paste0(uncats, " "))
  jenks_traits <- as.matrix(all_jenks_traits[,which_ncats == uncats,uncats])
  rownames(jenks_traits) <- rownames(traits)
  
  traitsPath_jenksCoding <- paste0("data/", fileName, "_PCA", percPCA, 
                                   "_jenks_", jenks_exp[jenks_exp_cats], "_expcats_", uncats, "_cats_traits.tsv")
  writeTSV(jenks_traits, traitsPath_jenksCoding)
}



nrun <- 2
for(transData in c("PCA99")){
  
  for(jenks_exp_cats in 1:length(jenks_exp)){
    
    ncats <- as.numeric(strsplit(readLines(paste0("data/", fileName, "_PCA", percPCA, "_jenks_", jenks_exp[jenks_exp_cats], "_expcats.txt"), warn = F), split = " ")[[1]])
    ncats <- sort(ncats)
    
    for (l in 1:nrun) {
      
      tempScript <- readLines("~/scripts/metascript_discrete_jenks.rev")
      tempScript[1] <- paste0("setwd(\"", getwd(), "\") ")
      tempScript[2] <- paste0(tempScript[2], l)
      tempScript[3] <- paste0("replicateNum <- ", "\"noReplicate\"")
      tempScript[5] <- paste0("coding <- \"", "jenks", ifelse(collapseHomo, "_noHomoPops_justSpecies", ""), "\"")
      tempScript[6] <- paste0("dataFileName = ", paste0("v(", paste0("\"", fileName, "_", transData, "_jenks_", jenks_exp[jenks_exp_cats], "_expcats_", ncats, "_cats_traits.tsv\"", collapse = ", "), ")"))
      tempScript[7] <- paste0("transform = \"", transData, "\"")
      tempScript[8] <- paste0("ncats = v(", paste0(ncats, collapse = ", "), ")")
      tempScript[9] <- paste0("jenks_exp_cats = ", jenks_exp[jenks_exp_cats])
      tempScript[10] <- paste0("total_ntraits = ", ntraits)
      
      
      writeLines(tempScript, paste0("scripts/", fileName, "_", transData, "_jenks_", jenks_exp[jenks_exp_cats], "_expcats_script_run_", l, ".Rev"))
    }
  }
}

 
