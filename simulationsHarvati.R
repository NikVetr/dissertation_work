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


#MCMC diagnostics
metaScript <- readLines("diagTemp.R", warn = F)
dir.create("diagnosticScripts")
nrep <- 500
for(i in 1:nrep){
  tempScript <- metaScript
  tempScript[1] <- paste0("setwd(dir = \"", getwd(), "/\") ")
  tempScript[2] <- paste0(tempScript[2], ntraits)
  tempScript[3] <- paste0(tempScript[3], i)
  writeLines(tempScript, paste0("diagnosticScripts/", i, ".R"))
}

sink(file = "diagnosticJobs.txt", append = F)

for(i in 1:nrep){
  cat(paste0("cd ", getwd(), "; ",  "/usr/local/bin/RScript ", "diagnosticScripts/", i, ".R\n"))
}   

sink()

#terminal command; can use system() but screen output is messy
cat(paste0("cd ", getwd(), "; parallel --jobs 7 --sshloginfile /Volumes/macOS/Users/nikolai/instances_all --eta --results terminal_output_diags --files < diagnosticJobs.txt"))


# instantiate matrices to read diagnostic output back in
traitNumIter <- c(46)

diagnostics <- matrix(ncol = 5 + 376, nrow = nrep*14)
colnames(diagnostics) <- c("diagnostic", "chain", "ntraits", "replicateNum","Posterior","Likelihood","TL",
                           "collessIndex","pair1","pair2","pair3","pair4","pair5","pair6","pair7","pair8","pair9","pair10","pair11",
                           "pair12","pair13","pair14","pair15","pair16","pair17","pair18","pair19","pair20","pair21","pair22","pair23",
                           "pair24","pair25","pair26","pair27","pair28","pair29","pair30","pair31","pair32","pair33","pair34","pair35",
                           "pair36","pair37","pair38","pair39","pair40","pair41","pair42","pair43","pair44","pair45","pair46","pair47",
                           "pair48","pair49","pair50","pair51","pair52","pair53","pair54","pair55","pair56","pair57","pair58","pair59",
                           "pair60","pair61","pair62","pair63","pair64","pair65","pair66","pair67","pair68","pair69","pair70","pair71",
                           "pair72","pair73","pair74","pair75","pair76","pair77","pair78","pair79","pair80","pair81","pair82","pair83",
                           "pair84","pair85","pair86","pair87","pair88","pair89","pair90","pair91","pair92","pair93","pair94","pair95",
                           "pair96","pair97","pair98","pair99","pair100","pair101","pair102","pair103","pair104","pair105","pair106",
                           "pair107","pair108","pair109","pair110","pair111","pair112","pair113","pair114","pair115","pair116","pair117",
                           "pair118","pair119","pair120","pair121","pair122","pair123","pair124","pair125","pair126","pair127","pair128",
                           "pair129","pair130","pair131","pair132","pair133","pair134","pair135","pair136","pair137","pair138","pair139",
                           "pair140","pair141","pair142","pair143","pair144","pair145","pair146","pair147","pair148","pair149","pair150",
                           "pair151","pair152","pair153","pair154","pair155","pair156","pair157","pair158","pair159","pair160","pair161",
                           "pair162","pair163","pair164","pair165","pair166","pair167","pair168","pair169","pair170","pair171","pair172",
                           "pair173","pair174","pair175","pair176","pair177","pair178","pair179","pair180","pair181","pair182","pair183",
                           "pair184","pair185","pair186","pair187","pair188","pair189","pair190","pair191","pair192","pair193","pair194",
                           "pair195","pair196","pair197","pair198","pair199","pair200","pair201","pair202","pair203","pair204","pair205",
                           "pair206","pair207","pair208","pair209","pair210","pair211","pair212","pair213","pair214","pair215","pair216",
                           "pair217","pair218","pair219","pair220","pair221","pair222","pair223","pair224","pair225","pair226","pair227",
                           "pair228","pair229","pair230","pair231","pair232","pair233","pair234","pair235","pair236","pair237","pair238",
                           "pair239","pair240","pair241","pair242","pair243","pair244","pair245","pair246","pair247","pair248","pair249",
                           "pair250","pair251","pair252","pair253","pair254","pair255","pair256","pair257","pair258","pair259","pair260",
                           "pair261","pair262","pair263","pair264","pair265","pair266","pair267","pair268","pair269","pair270","pair271",
                           "pair272","pair273","pair274","pair275","pair276","pair277","pair278","pair279","pair280","pair281","pair282",
                           "pair283","pair284","pair285","pair286","pair287","pair288","pair289","pair290","pair291","pair292","pair293",
                           "pair294","pair295","pair296","pair297","pair298","pair299","pair300","pair301","pair302","pair303","pair304",
                           "pair305","pair306","pair307","pair308","pair309","pair310","pair311","pair312","pair313","pair314","pair315",
                           "pair316","pair317","pair318","pair319","pair320","pair321","pair322","pair323","pair324","pair325","pair326",
                           "pair327","pair328","pair329","pair330","pair331","pair332","pair333","pair334","pair335","pair336","pair337",
                           "pair338","pair339","pair340","pair341","pair342","pair343","pair344","pair345","pair346","pair347","pair348",
                           "pair349","pair350","pair351","rfDist1",
                           "kfDist1","clade1","clade2","clade3","clade4","clade5","clade6","clade7","clade8","clade9","clade10","clade11",
                           "clade12","clade13","clade14","clade15","clade16","clade17","clade18","clade19","clade20")

comptrees <- matrix(ncol = 4, nrow = nrep)
colnames(comptrees) <- c("ntraits", "replicateNum", "correlation", "p.value")

counter <- 0
for(i in 1:nrep){
  print(i)    
  
  counter <- counter + 1
  
  ntraits <- 46
  replicateNum <- i
  
  fileName <- paste0("simulations_", ntraits, "traits_", "replicate_", replicateNum)
  
  load(file = paste0("diagnostics/geweke/", fileName, "_c1.txt"))
  load(file = paste0("diagnostics/geweke/", fileName, "_c2.txt"))
  load(file = paste0("diagnostics/geweke/", fileName, "_cb.txt"))
  
  load(file = paste0("diagnostics/ess/", fileName, "_c1.txt"))
  load(file = paste0("diagnostics/ess/", fileName, "_c2.txt"))
  load(file = paste0("diagnostics/ess/", fileName, "_cb.txt"))
  
  load(file = paste0("diagnostics/heiwel/", fileName, "_c1.txt"))
  load(file = paste0("diagnostics/heiwel/", fileName, "_c2.txt"))
  load(file = paste0("diagnostics/heiwel/", fileName, "_cb.txt"))
  
  if(is.vector(heiwel)){heiwel <- rbind(heiwel, heiwel)}; if(is.vector(heiwel1)){heiwel1 <- rbind(heiwel1, heiwel1)}; if(is.vector(heiwel2)){heiwel2 <- rbind(heiwel2, heiwel2)}
  
  load(file = paste0("diagnostics/psrf/", fileName, "_cb.txt"))    
  if(is.vector(psrf)){psrf <- rbind(psrf, psrf); print("changing...")}
  
  diagnostics[(((counter-1)*14+1):(counter*14)),] <- rbind(
    c(deparse(substitute(geweke1)),1,ntraits,replicateNum, geweke1[1:377]),
    c(deparse(substitute(ess1)),1,ntraits,replicateNum, ess1[1:377]),
    c(paste0(deparse(substitute(heiwel1)), "_CVM.p.val"),1,ntraits,replicateNum, heiwel1[1,1:377]),
    c(paste0(deparse(substitute(heiwel1)), "_halfwidth.passed"),1,ntraits,replicateNum, heiwel1[2,1:377]),
    
    c(deparse(substitute(geweke2)),2,ntraits,replicateNum, geweke2[1:377]),
    c(deparse(substitute(ess2)),2,ntraits,replicateNum, ess2[1:377]),
    c(paste0(deparse(substitute(heiwel2)), "_CVM.p.val"),2,ntraits,replicateNum, heiwel2[1,1:377]),
    c(paste0(deparse(substitute(heiwel2)), "_halfwidth.passed"),2,ntraits,replicateNum, heiwel2[2,1:377]),
    
    c(deparse(substitute(psrf_point)),"both",ntraits,replicateNum, psrf[1,1:377]),
    c(deparse(substitute(psrf_upper)),"both",ntraits,replicateNum, psrf[2,1:377]),
    c(deparse(substitute(geweke)),"both",ntraits,replicateNum, geweke[1:377]),
    c(deparse(substitute(ess)),"both",ntraits,replicateNum, ess[1:377]),
    c(paste0(deparse(substitute(heiwel)), "_CVM.p.val"),"both",ntraits,replicateNum, heiwel[1,1:377]),
    c(paste0(deparse(substitute(heiwel)), "_halfwidth.passed"),"both",ntraits,replicateNum, heiwel[2,1:377])
  )
  
  comptrees[counter,] <- unlist(read.table(file = paste0("diagnostics/compareTree/comptrees_for_", fileName, ".txt")))
}

diagnostics <- as.data.frame(diagnostics)
diagnostics[,6:381] <- sapply(6:381, function(x) as.numeric(as.character(diagnostics[,x])))
save(diagnostics, file = "diagnostics.txt")
save(comptrees, file = "CompareTreesCorrsPs")


#loading diagnostics back in

print(paste0("# of original failures: ", length(read.delim("job_MCMC_Fail.txt", header = F)[[1]])/2))
load(file = "diagnostics.txt")
load(file = "CompareTreesCorrsPs")

#look at ESS and R2
essBothChains <- diagnostics[diagnostics$diagnostic == "ess",]
thresholdESS <- 1000
thresholdR2 <- 0.95
failures <- c(which(as.numeric(as.character(essBothChains$Posterior)) < thresholdESS), 
              which(as.numeric(as.character(essBothChains$Likelihood)) < thresholdESS),
              which(as.numeric(as.character(essBothChains$TL)) < thresholdESS),
              which(as.numeric(as.character(essBothChains$rfDist1)) < thresholdESS),
              which(as.numeric(as.character(essBothChains$kfDist1)) < thresholdESS),
              which(as.numeric(as.character(essBothChains$kfDist1)) < thresholdESS),
              unique(unlist(sapply(1:length(which(startsWith(colnames(essBothChains), "pair"))), 
                                   function(x) which(as.numeric(as.character(essBothChains[,which(startsWith(colnames(essBothChains), "pair"))[x]])) < thresholdESS)))),
              which(comptrees[,3]^2 < thresholdR2))
failures <- sort(unique(failures))
failures <- essBothChains[failures,3:5]
failures <- failures[2*(1:(length(failures[,1])/2)),]
summary(failures$rateMatrixCorrelation)
