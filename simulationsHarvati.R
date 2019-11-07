setwd("/Volumes/1TB/Harvati")

library(ape)
library(coda)
library(rwty)
library(phangorn)
library(lattice)
library(latticeExtra)
library(MCMCpack)
library(mvMORPH)
library(phytools)

# cov <- readLines("/Users/nikolai/data/Harvati_noPCA.tsv")
# cov <- strsplit(cov, split = "\t")
# ntraits <- length(cov[[1]])
# cov <- as.matrix(sapply(1:ntraits, function(x) as.numeric(cov[[(1:ntraits*2-1)[x]]])))
load("rateMatrix")

origTraits <- strsplit(readLines("/Users/nikolai/data/Harvati_noPCA.nex"), split = "\t")
origTraits <- origTraits[8:(length(origTraits)-2)]
rwn <- sapply(1:length(origTraits), function(x) strsplit(origTraits[[x]], split = " ")[[1]][1])
origTraits <- t(as.matrix(sapply(1:length(origTraits), function(x) as.numeric(strsplit(origTraits[[x]], split = " ")[[1]][-1]))))
rownames(origTraits) <- rwn


trees1 <- read.tree("/Users/nikolai/output/harvati_noPCA_c1.trees")
trees2 <- read.tree("/Users/nikolai/output/harvati_noPCA_c2.trees")
trees <- c(trees1, trees2)
startsWith2 <- function(x, pre){apply(sapply(1:length(pre), function(n) startsWith(x, pre[n])), 1, any)}
trees <- lapply(1:length(trees), function(x)
  root(trees[[x]], outgroup = tree$tip.label[startsWith2(tree$tip.label, c("Mac", "Pap", "Man"))], resolve.root = T, edgelabel = T))

ntips <- length(trees[[1]]$tip.label)

meansNexus <- function (characterMatrix, outputFilePath, traitsOnColumns = T){
  if (traitsOnColumns) {ntax = dim(characterMatrix)[1]; nchar = dim(characterMatrix)[2]}
  else {ntax = dim(characterMatrix)[2]; nchar = dim(characterMatrix)[1]}
  sink(outputFilePath, append=F)
  cat("#NEXUS\n\nBegin data;\n\tDimensions ntax=", ntax, "nchar=", nchar, ";\n\tFormat datatype=Continuous missing=?;\n\tMatrix\n\n")
  for(i in 1:length(characterMatrix[,1])) {cat(rownames(characterMatrix)[i], characterMatrix[i,], "\n")}
  cat(";\nEnd;")
  sink()
}


ratesPhy <- function (rateMatrix, outputFilePath) {
  if(!isSymmetric(rateMatrix)) {stop("Error: This matrix is not symmetric.")}
  sink(outputFilePath, append=F)
  cat(dim(rateMatrix)[1],"\n")
  for(i in 1:length(rateMatrix[,1])) {
    cat(rownames(rateMatrix)[i], paste0(as.vector(rateMatrix[i, (1:(length(rateMatrix[i,])-1))])))
    cat(" ")
    cat(rateMatrix[i, (length(rateMatrix[i,]))], "\n", sep = "")
  }
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

traitNumIter <- c(ntraits)

ntips

#remove old files if necessary
file.remove(paste0("trees/", list.files("trees/")))
dir.create("trees")
file.remove(paste0("data/", list.files("data/")))
dir.create("data")
file.remove(paste0("rates/", list.files("rates/")))
dir.create("rates")
file.remove(paste0("scripts/", list.files("scripts/")))
dir.create("scripts")

nrun <- 2
nrep <- 500

# i iterates over replicates
for(i in 1:nrep){
  
    # j iterates over trait numbers
    for (j in 1:length(traitNumIter)) {
      
      
      print(c(i, j))
      
      ntraits <- traitNumIter[j]
      replicateNum <- i
      
      fileName <- paste0("simulations_", ntraits, "traits_", "replicate_", replicateNum)
      
      # write the covariance matrices
      ratesPathPhy <- paste0("rates/", fileName, "_rates.Phy")
      ratesPhy(cov, ratesPathPhy)
      
      # simulate the traits and write to file
      traitsPath <- paste0("data/", fileName, "_traits.nex")
      
      #use Harvati trees
      tree <- trees[[sample(size = 1, x = 1:length(trees))]]
      
      traits <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = dim(cov)[1], sigma=cov, mu= rep(0, times = dim(cov)[1])))
      meansNexus(traits, traitsPath)
      
      treePath <- paste0("trees/", fileName, "_tree.nex")
      write.nexus(tree, file = treePath, translate = T)
      
      if(i == 1 & j == 1){
        file.remove("jobs.txt")
      }
      
      # write the master scripts
      for (l in 1:nrun) {
        
        tempScript <- readLines("metascript.txt")
        tempScript[1] <- paste0("setwd(\"", getwd(), "/\") ")
        tempScript[2] <- paste0(tempScript[2], l)
        tempScript[5] <- paste0("replicateNum <- ", replicateNum)
        tempScript[8] <- paste0("means <- readContinuousCharacterData(\"data/", fileName, "_traits.nex\")")
        tempScript[11] <- paste0("rate = readDistanceMatrix(\"rates/", fileName, "_rates.Phy\")")
        writeLines(tempScript, paste0("scripts/", fileName, "_script_run_", l, ".Rev"))
        
        sink(file = "jobs.txt", append = T)
        cat(paste0("cd ", getwd(), "; ",  "/Users/nikolai/repos/revbayes-development/projects/cmake/rb scripts/", fileName, "_script_run_", l, ".Rev\n"))
        sink()
      }
        
  }
      
}
    
#terminal command; can use system() but screen output is messy
cat(paste0("cd ", getwd(), "; parallel --jobs 7 --sshloginfile /Volumes/macOS/Users/nikolai/instances_all --eta --results terminal_output --files < jobs.txt"))
