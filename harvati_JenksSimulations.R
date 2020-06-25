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

nreps <- 100
traitNumIter <- c((46-1)*2^(0:4)+1)
jenks_exp <- c(2,4,6,8)

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
          #observation: first break is always the set min, so ignore it... but it only happens some of the time lol, so check when that is
        # jenksBreaks <- sapply(1:ntraits, function(trait) getJenksBreaks(traits[,trait], k = njenks-1))
          
        # jenksBreaks <- sapply(1:ntraits, function(trait) getJenksBreaks(traits[,trait], k = njenks-1))
        # jenksBreaks <- sapply(1:ntraits, function(trait) as.numeric(classInt::classIntervals(var = traits[,trait], style = "jenks", n = njenks)$brks[2:njenks]))
        # jenksBreaks2 <- sapply(1:ntraits, function(trait) getJenksBreaks(traits[,trait], k = njenks)[-1])
        # jenksBreaks3 <- sapply(1:ntraits, function(trait) getJenksBreaks(traits[,trait], k = njenks+1)[-c(1:2)])
        # jenksBreaks4 <- sapply(1:ntraits, function(trait) getJenksBreaks(traits[,trait], k = njenks+2)[-c(1:3)]) #never need this one, but just in case...
        # 
        ### actually we don't need to do this adhoc, can just compare GFVs ###
        # minTraits <- sapply(1:ntraits, function(trait) min(traits[,trait]))
        # minBreaks <- sapply(1:ntraits, function(trait) min(ifelse(njenks > 2, jenksBreaks[,trait],jenksBreaks[trait])))
        # if(njenks == 2){
        #   jenksBreaks[minTraits == minBreaks] <- jenksBreaksP1[minTraits == minBreaks]
        #   minBreaks <- sapply(1:ntraits, function(trait) min(ifelse(njenks > 2, jenksBreaks[,trait],jenksBreaks[trait])))
        #   jenksBreaks[minTraits == minBreaks] <- jenksBreaksP2[minTraits == minBreaks]
        # }
        # if(njenks > 2){
        #   #sometimes the first two breaks are also equal, so check for that, too
        #   jenksBreaks[,jenksBreaks[1,] == jenksBreaks[2,]] <- jenksBreaksP1[,jenksBreaks[1,] == jenksBreaks[2,]]
        #   minBreaks <- sapply(1:ntraits, function(trait) min(ifelse(njenks > 2, jenksBreaks[,trait],jenksBreaks[trait])))
        #   #I also think it's unlikely that we're parseling off the set min if njenks < e.g. 10, so let's adjust for that too
        #   if(njenks < 10){
        #     jenksBreaks[,minTraits == minBreaks] <- jenksBreaksP1[,minTraits == minBreaks]
        #     minBreaks <- sapply(1:ntraits, function(trait) min(ifelse(njenks > 2, jenksBreaks[,trait],jenksBreaks[trait])))
        #     jenksBreaks[,minTraits == minBreaks] <- jenksBreaksP2[,minTraits == minBreaks]
        #   }
        # }
        
        #to avoid breaks being on boundaries, shift all closer to 0
        #this ensures that the lowest trait gets score = 0, and the highest trait gets score = max, when the thresholds are on those boundaries
        # jenksBreaks <- jenksBreaks * (1-10^-5)
        # jenksBreaks2 <- jenksBreaks2 * (1-10^-5)
        # jenksBreaks3 <- jenksBreaks3 * (1-10^-5)
        # jenksBreaks4 <- jenksBreaks4 * (1-10^-5)
        traitMeans <- sapply(1:ntraits, function(trait) mean(traits[,trait]))
        traitMags <- sapply(1:ntraits, function(trait) mean(abs(traits[,trait])))
        underflowHack <- 10^round(log(abs(1/traitMags)) / log(10))
        SDAM <- sapply(1:ntraits, function(trait) ssqm(traits[,trait] * underflowHack[trait]))
        
        if(njenks > 2){
          # jenksCodedTraits <- sapply(1:ntraits, function(trait) sapply(1:ntips, function(tip) sum(traits[tip,trait] >= jenksBreaks[,trait]))) + 1
          jenksCodedTraits <- sapply(1:ntraits, function(trait) findCols(classInt::classIntervals(var = traits[,trait] - traits[1,trait], style = "jenks", n = njenks)))
          # jenksCodedTraits2 <- sapply(1:ntraits, function(trait) sapply(1:ntips, function(tip) sum(traits[tip,trait] >= jenksBreaks2[,trait]))) + 1
          # jenksCodedTraits3 <- sapply(1:ntraits, function(trait) sapply(1:ntips, function(tip) sum(traits[tip,trait] >= jenksBreaks3[,trait]))) + 1
          # jenksCodedTraits4 <- sapply(1:ntraits, function(trait) sapply(1:ntips, function(tip) sum(traits[tip,trait] >= jenksBreaks4[,trait]))) + 1
        } else {
          jenksCodedTraits <- sapply(1:ntraits, function(trait) findCols(classInt::classIntervals(var = traits[,trait] - traits[1,trait], style = "jenks", n = njenks)))
          # jenksCodedTraits2 <- sapply(1:ntraits, function(trait) sapply(1:ntips, function(tip) sum(traits[tip,trait] >= jenksBreaks2[trait]))) + 1
          # jenksCodedTraits3 <- sapply(1:ntraits, function(trait) sapply(1:ntips, function(tip) sum(traits[tip,trait] >= jenksBreaks3[trait]))) + 1
          # jenksCodedTraits4 <- sapply(1:ntraits, function(trait) sapply(1:ntips, function(tip) sum(traits[tip,trait] >= jenksBreaks4[trait]))) + 1
        }
        rownames(jenksCodedTraits) <- rownames(traits)
        # rownames(jenksCodedTraits2) <- rownames(traits)
        # rownames(jenksCodedTraits3) <- rownames(traits)
        # rownames(jenksCodedTraits4) <- rownames(traits)
        
        SDCM <- sapply(1:ntraits, function(trait) 
          sum(sapply(1:njenks, function(jcat) ssqm(traits[jenksCodedTraits[,trait] == jcat,trait] * underflowHack[trait]))))
        # SDCM2 <- sapply(1:ntraits, function(trait) 
        #   sum(sapply(1:njenks, function(jcat) ssqm(traits[jenksCodedTraits2[,trait] == jcat,trait] * underflowHack[trait]))))
        # SDCM3 <- sapply(1:ntraits, function(trait) 
        #   sum(sapply(1:njenks, function(jcat) ssqm(traits[jenksCodedTraits3[,trait] == jcat,trait] * underflowHack[trait]))))
        # SDCM4 <- sapply(1:ntraits, function(trait) 
        #   sum(sapply(1:njenks, function(jcat) ssqm(traits[jenksCodedTraits4[,trait] == jcat,trait] * underflowHack[trait]))))
        
        #compute goodness of variance fit
        GVF <- (SDAM - SDCM) / SDAM
        # GVF2 <- (SDAM - SDCM2) / SDAM
        # GVF3 <- (SDAM - SDCM3) / SDAM
        # GVF4 <- (SDAM - SDCM4) / SDAM
        
        # GVFs <- cbind(GVF, GVF2, GVF3, GVF4)
        # bestGVF <- sapply(1:ntraits, function(trait) which.max(GVFs[trait,]))
        
        # jenksCodedTraits[,bestGVF == 2] <- jenksCodedTraits2[,bestGVF == 2]
        # jenksCodedTraits[,bestGVF == 3] <- jenksCodedTraits3[,bestGVF == 3]
        # jenksCodedTraits[,bestGVF == 4] <- jenksCodedTraits4[,bestGVF == 4]
        # GVF[bestGVF == 2] <- GVF2[bestGVF == 2]
        # GVF[bestGVF == 3] <- GVF3[bestGVF == 3]
        # GVF[bestGVF == 4] <- GVF4[bestGVF == 4]
        
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

#convert from 0 as starting trait index to 1 and back again in jenks' coding
#because revbayes is very inflexible lol

nreps <- 100
traitNumIter <- c((46-1)*2^(0:4)+1)
incr <- T
decr <- !incr
for(i in 1:nreps){
  
  cat(paste0("\n", i, ":"))
  
  for(j in 1:length(traitNumIter)){
    
    cat(paste0(" ", j, " "))
    
    ntraits <- traitNumIter[j]
    replicateNum <- i
    
    fileName <- paste0("simulations_", ntraits, "traits_", "replicate_", replicateNum)
    
    for(transData in c("raw", "PCs")){
      
      unique_ncats <- as.numeric(strsplit(readLines(paste0("data_discrete/", fileName, "_", transData, "_jenks_", jenks_exp[jenks_exp_cats], "_expcats.txt")), split = " ")[[1]])
      
      for(uncats in unique_ncats){
        traitsPath_jenksCoding <- paste0("data_discrete/", fileName, "_", transData, 
                                         "_jenks_", jenks_exp[jenks_exp_cats], "_expcats_", uncats, "_cats_traits.tsv")
        jenksCodedTraits <- read.table(traitsPath_jenksCoding, header = F, sep = "\t", colClasses = "character")
        rwn <- trimws(jenksCodedTraits[,1]); jenksCodedTraits[] <- as.numeric(as.matrix(jenksCodedTraits)); rownames(jenksCodedTraits) <- rwn; 
        jenksCodedTraits <- jenksCodedTraits[,-1]
        
        for(foo in 1:length(jenksCodedTraits[1,])){
          jenksCodedTraits[,foo] <- jenksCodedTraits[,foo] + ifelse(incr, 1, ifelse(decr, -1, 0))
        }
        
        writeTSV(jenksCodedTraits, traitsPath_jenksCoding)
        
      }
    }
  }
}


#now construct the revbayes scripts to analyze these data
nreps <- 100
traitNumIter <- c((46-1)*2^(0:4)+1)
jenks_exp <- c(2,4,6,8)
jenks_exp <- c(2)
nrun <- 2

counter <- 0
for(i in 1:nreps){
  
  cat(paste0("\n", i, ":"))
  
  for(j in 1:length(traitNumIter)){
    
    cat(paste0(" ", j, " "))
    
    ntraits <- traitNumIter[j]
    replicateNum <- i
    
    fileName <- paste0("simulations_", ntraits, "traits_", "replicate_", replicateNum)
    
    #get true tree just in case lol
    trueTree <- read.nexus(file = paste0("trees/", fileName, "_tree.nex"))
    tips <- trueTree$tip.label
    ntips <- length(tips)
    
    counter <- counter + 1
    if(counter == 1){
      file.remove("jobs_rb_freek_jenks.txt")
    }
    
    for(transData in c("raw", "PCs")){
      
      for(jenks_exp_cats in 1:length(jenks_exp)){
        
        ncats <- as.numeric(strsplit(readLines(paste0("data_discrete/", fileName, "_", transData, "_jenks_", jenks_exp[jenks_exp_cats], "_expcats.txt")), split = " ")[[1]])
        ncats <- sort(ncats)
        
        for (l in 1:nrun) {
          
          tempScript <- readLines("/Volumes/1TB/Harvati/metascript_discrete_jenks.rev")
          tempScript[1] <- paste0("setwd(\"", getwd(), "/\") ")
          tempScript[2] <- paste0(tempScript[2], l)
          tempScript[3] <- paste0("replicateNum <- ", replicateNum)
          tempScript[5] <- paste0("coding <- \"", "jenks", "\"")
          tempScript[6] <- paste0("dataFileName = ", paste0("v(", paste0("\"", fileName, "_", transData, "_jenks_", jenks_exp[jenks_exp_cats], "_expcats_", ncats, "_cats_traits.tsv\"", collapse = ", "), ")"))
          tempScript[7] <- paste0("transform = \"", transData, "\"")
          tempScript[8] <- paste0("ncats = v(", paste0(ncats, collapse = ", "), ")")
          tempScript[9] <- paste0("jenks_exp_cats = ", jenks_exp[jenks_exp_cats])
          tempScript[10] <- paste0("total_ntraits = ", ntraits)
          
          
          writeLines(tempScript, paste0("scripts_discrete/", fileName, "_", transData, "_jenks_", jenks_exp[jenks_exp_cats], "_expcats_script_run_", l, ".Rev"))
          
          scriptName <- "jobs_rb_freek_jenks.txt"
          sink(file = scriptName, append = T)
          cat(paste0("cd ", getwd(), "; ",  "/Users/nikolai/repos/revbayes-development/projects/cmake/rb2 ", 
                     paste0("scripts_discrete/", fileName, "_", transData, "_jenks_", jenks_exp[jenks_exp_cats], "_expcats_script_run_", l, ".Rev\n")))
          sink()
        }
      }   
    }
  }
}

#freek constructor
path <- "orderedFreekScripts/orderedFreekNUM.Rev"
temp <- readLines(path)
for(i in 1:500){
  cat(paste0(i, " "))
  newScript <- gsub(x = temp, pattern = "NUM", replacement = i)
  writeLines(newScript, con = gsub(x = path, pattern = "NUM", replacement = i))
}

#terminal command; can use system() but screen output is messy
cat(paste0("cd ", getwd(), "; parallel --jobs 8 --sshloginfile /Volumes/macOS/Users/nikolai/instances_noDogCat --eta --results terminal_output_jenks --files < jobs_rb_freek_jenks.txt"))
