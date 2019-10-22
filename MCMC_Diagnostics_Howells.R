setwd(dir = "/Volumes/2TB/uvBM_sims_Howells//")
dir.create("diagnostics")
dir.create("diagnostics/rwty")
dir.create("diagnostics/compareTree")

library(abind)
library(parallel)
library(foreach)
library(doParallel)
library(ape)
library(coda)
library(rwty)
library(phangorn)
library(lattice)
library(latticeExtra)
library(MCMCpack)

#find unfinished files
unfinishedFiles <- setdiff(1:100, sort(as.numeric(gsub(gsub(list.files("diagnostics/")[startsWith(list.files("diagnostics/"), "diag")],
                                         pattern = "diagnostics", replacement = ""), pattern = ".txt", replacement = ""))))

convNum2Str <- function(nums, key){
  sapply(1:length(nums), function(x) key[nums[x]])
}

prop.part.df <- function(trees, cutoff = 0.01){
  numTrees <- length(trees)
  out <- prop.part(trees)
  outList <- as.list.data.frame(out)
  pps <- attributes(outList)$number/numTrees
  props <- data.frame(pps, as.matrix(outList)); colnames(props) <- c("postProbs", "clade")
  props <- props[order(-pps),]
  props <- props[props[,1] > cutoff,]
  rownames(props) <- 1:nrow(props)
  props$cladeNames <- sapply(1:length(props[,1]), function(x) sort(convNum2Str(props$clade[[x]], attributes(outList)$labels)))
  props <- props[,c(1,3)]
  props[-1,]
}

compareTrees <- function(biparts1, biparts2, tipLabs){
  matchProbs <- matrix(0, nrow = sum(length(biparts1[,1]), length(biparts2[,1])), ncol = 2)
  counter <- 1
  biparts2$notSeen <- 1
  for(clade in 1:length(biparts1[,2])){
    cladeName <- biparts1[,2][[clade]]
    isThere <- sapply(1:length(biparts2[,1]), function(x) identical(cladeName, biparts2[x,2][[1]]))
    altCladeName <- sort(setdiff(tipLabs, cladeName))
    orIsThere <- sapply(1:length(biparts2[,1]), function(x) identical(altCladeName, biparts2[x,2][[1]]))
    if(any(isThere)){
      biparts2$notSeen[isThere] <- 0
      matchProbs[counter,] <- c(biparts1[,1][clade], biparts2[,1][isThere])
    } else if (any(orIsThere)) {
      biparts2$notSeen[orIsThere] <- 0
      matchProbs[counter,] <- c(biparts1[,1][clade], biparts2[,1][orIsThere])
    } else {
      matchProbs[counter,] <- c(biparts1[,1][[clade]], 0)
    }
    counter <- counter + 1
  }
  if(sum(biparts2$notSeen) > 0){
    matchProbs[counter:(counter+sum(biparts2$notSeen)-1),] <- cbind(0, biparts2[biparts2$notSeen == 1,1])
  }
  matchProbs <- matchProbs[rowSums(matchProbs == c(0,0)) != 2,]
  matchProbs
}

clade.PA <- function(trees, numClades = NA){
  biparts <- prop.part.df(trees)
  if(is.na(numClades)){
    numClades <- length(biparts[,1])
  } else if(numClades > length(biparts[,1])){
    numClades <- length(biparts[,1])
  } 
  tipLabs <- trees[[1]]$tip.label
  cladePA <- matrix(0, nrow = length(trees), ncol = numClades)
  colnames(cladePA) <- sapply(1:numClades, function(x) paste0(biparts[,2][[x]], collapse = "-"))
  for(clade in 1:numClades){
    #print(clade)
    for(treeNum in 1:length(trees)){
      treeBipart <- prop.part.df(trees[treeNum])[,2]
      treeBipart <- sapply(1:length(treeBipart), function(x) paste0(treeBipart[[x]], collapse = "-"))
      cladeName <- biparts[,2][[clade]]
      altCladeName <- paste0(sort(setdiff(tipLabs, cladeName)), collapse = "-")
      cladeName <- paste0(cladeName, collapse = "-")
      cladePA[treeNum, clade] <- as.numeric(any(treeBipart == cladeName | treeBipart == altCladeName))
    }
  }
  cladePA
}

cl <- makeCluster(4, outfile="")
registerDoParallel(cl)
getDoParWorkers()
# stopCluster(cl)

# foreach(m=1:100, .packages = c("ape", "coda", "rwty", "phangorn", "lattice", "latticeExtra", "MCMCpack")) %dopar% {
foreach(unfinishedTODO=1:length(unfinishedFiles), .packages = c("ape", "coda", "rwty", "phangorn", "lattice", "latticeExtra", "MCMCpack")) %dopar% { m <- unfinishedFiles[unfinishedTODO]
    
  diagnostics <- matrix(ncol = 5 + 461)[-1,]
  colnames(diagnostics) <- c("diagnostic", "chain", "ntraits", "rateMatrixMisspecification", "replicateNum","Posterior","Likelihood","TL",
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
                             "pair349","pair350","pair351","pair352","pair353","pair354","pair355","pair356","pair357","pair358","pair359",
                             "pair360","pair361","pair362","pair363","pair364","pair365","pair366","pair367","pair368","pair369","pair370",
                             "pair371","pair372","pair373","pair374","pair375","pair376","pair377","pair378","pair379","pair380","pair381",
                             "pair382","pair383","pair384","pair385","pair386","pair387","pair388","pair389","pair390","pair391","pair392",
                             "pair393","pair394","pair395","pair396","pair397","pair398","pair399","pair400","pair401","pair402","pair403",
                             "pair404","pair405","pair406","pair407","pair408","pair409","pair410","pair411","pair412","pair413","pair414",
                             "pair415","pair416","pair417","pair418","pair419","pair420","pair421","pair422","pair423","pair424","pair425",
                             "pair426","pair427","pair428","pair429","pair430","pair431","pair432","pair433","pair434","pair435","rfDist1",
                             "kfDist1","clade1","clade2","clade3","clade4","clade5","clade6","clade7","clade8","clade9","clade10","clade11",
                             "clade12","clade13","clade14","clade15","clade16","clade17","clade18","clade19","clade20")
  
  comptrees <- matrix(ncol = 5)[-1,]
  colnames(comptrees) <- c("ntraits", "rateMatrixMisspecification", "replicateNum", "correlation", "p.value")
  
  #Also output marginal histograms and trace plots for continuous parameters, split	frequencies	&	presence/absence for trees
  
  degreeMisspecification <- c("perfect", "interregional", "chimpHuman")
  degreeMisspecificationNames <- c("perfect", "interregional", "interspecific")
  traitNumIter <- c(2,4,8,16,32,57)
  rwty <- F
  
  # k iterates over degree of misspecification
  for (k in 1:length(degreeMisspecification)) {
    
    for (j in 1:6) { #number of traits
      
      for(i in m:m){
        
        cat(c(i, j, k, "... \n"))
        
        ntraits <- traitNumIter[j]
        rateMatrixMisspecification <- degreeMisspecification[k]
        replicateNum <- i
        
        fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixMisspecification, "RateMatrix_replicate_", replicateNum)
        trueTree <- read.nexus(file = paste0("trees/", fileName, "_tree.nex"))
        
        #read in mcmc output
        log1 <- read.table(file = paste0("output/", fileName, "_params_run_1.log"), header = T)
        log2 <- read.table(file = paste0("output/", fileName, "_params_run_2.log"), header = T)
        trees1 <- read.tree(file = paste0("output/", fileName, "_trees_run_1.trees"))
        trees2 <- read.tree(file = paste0("output/", fileName, "_trees_run_2.trees"))
        
        #compute tree length
        log1$TL <- sapply(1:length(log1[,1]), function(x) sum(log1[x,5:62]))
        log1$collessIndex <- 100
        log2$TL <- sapply(1:length(log2[,1]), function(x) sum(log2[x,5:62]))
        log2$collessIndex <- 100
        
        #compute patristic distances
        npairs <- ( length(trees1[[1]]$tip.label) * length(trees1[[1]]$tip.label) - length(trees1[[1]]$tip.label) ) / 2
        patristicDists1 <- matrix(0, nrow = length(trees1), ncol = npairs)
        colnames(patristicDists1) <- paste0("pair", 1:npairs)
        for(treeNum in 1:length(trees1)){
          patristicDists <- cophenetic.phylo(trees1[[treeNum]])
          patristicDists <- patristicDists[sort(rownames(patristicDists)),sort(rownames(patristicDists))]
          patristicDists <- patristicDists[upper.tri(patristicDists)]
          patristicDists1[treeNum,] <- patristicDists
        }
        
        npairs <- ( length(trees2[[1]]$tip.label) * length(trees2[[1]]$tip.label) - length(trees2[[1]]$tip.label) ) / 2
        patristicDists2 <- matrix(0, nrow = length(trees2), ncol = npairs)
        colnames(patristicDists2) <- paste0("pair", 1:npairs)
        for(treeNum in 1:length(trees2)){
          patristicDists <- cophenetic.phylo(trees2[[treeNum]])
          patristicDists <- patristicDists[sort(rownames(patristicDists)),sort(rownames(patristicDists))]
          patristicDists <- patristicDists[upper.tri(patristicDists)]
          patristicDists2[treeNum,] <- patristicDists
        }
        
        #compute rf dists from true tree 
        rfDist1 <- sapply(1:length(trees1), function(x) dist.topo(x = unroot(trueTree), y = unroot(trees1[[x]]), method = "PH85")[1])
        rfDist2 <- sapply(1:length(trees2), function(x) dist.topo(x = unroot(trueTree), y = unroot(trees2[[x]]), method = "PH85")[1])
        
        #kf dists from true tree
        kfDist1 <- sapply(1:length(trees1), function(x) dist.topo(x = unroot(trueTree), y = unroot(trees1[[x]]), method = "score")[1])
        kfDist2 <- sapply(1:length(trees2), function(x) dist.topo(x = unroot(trueTree), y = unroot(trees2[[x]]), method = "score")[1])
        
        #compare trees plots
        biparts1 <- prop.part.df(trees1, cutoff = 0.005)
        biparts2 <- prop.part.df(trees2, cutoff = 0.005)
        compareTreeProbs <- compareTrees(biparts1, biparts2, trueTree$tip.label)
        compareTreesCorr <- cor(compareTreeProbs)[1,2]
        png(paste0("diagnostics/compareTree/", fileName, "_compareTreePlot.png"), width = 500, height = 500)
        plot(compareTreeProbs, xlim = c(0,1), ylim = c(0,1)); abline(0,1, col = 2, lwd = 2); text(0.95,0, paste0("r = ", round(compareTreesCorr, 4)))
        dev.off()
        
        #clade P/A, 20 most probable splits
        cladePAtrees1 <- clade.PA(trees1, 20)
        cladePAtrees2 <- clade.PA(trees2, 20)
        # plot(cladePAtrees1[,3], type = "l")
        
        #put together master log file 
        log1 <- cbind(log1[,c("Posterior", "Likelihood", "TL", "collessIndex")], patristicDists1, rfDist1, kfDist1, cladePAtrees1)
        log2 <- cbind(log2[,c("Posterior", "Likelihood", "TL", "collessIndex")], patristicDists2, rfDist2, kfDist2, cladePAtrees2)
        
        #single-chain diagnostics
        geweke1 <- sapply(1:length(log1[1,]), function(x) geweke.diag(log1[,x])[1]$z[[1]])
        ess1 <- sapply(1:length(log1[1,]), function(x) effectiveSize(log1[,x])[[1]])
        heiwel1 <- sapply(1:length(log1[1,]), function(x) heidel.diag(log1[,x])[c(3,4)])
        
        geweke2 <- sapply(1:length(log2[1,]), function(x) geweke.diag(log2[,x])[1]$z[[1]])
        ess2 <- sapply(1:length(log2[1,]), function(x) effectiveSize(log2[,x])[[1]])
        heiwel2 <- sapply(1:length(log2[1,]), function(x) heidel.diag(log2[,x])[c(3,4)])
        
        #multi-chain diagnostics
        psrf <- sapply(1:length(log1[1,]), function(x) gelman.diag(mcmc.list(mcmc(log1[,x]), mcmc(log2[,x])))$psrf)
        comptrecor <- cor.test(compareTreeProbs[,1], compareTreeProbs[,2])[c("p.value", "estimate")]
        geweke <- sapply(1:length(log2[1,]), function(x) geweke.diag(c(log1[,x],log2[,x]))[1]$z[[1]])
        ess <- sapply(1:length(log2[1,]), function(x) effectiveSize(c(log1[,x],log2[,x]))[[1]])
        heiwel <- sapply(1:length(log2[1,]), function(x) heidel.diag(c(log1[,x],log2[,x]))[c(3,4)])
        
        if(rwty){
          ### rwty stuff ###
          #write log files
          write.table(x = rbind(colnames(log1[,1:4]),log1[,1:4]), file = paste0("diagnostics/rwty/log1rwty", m), row.names = F, col.names = T)
          write.table(x = rbind(colnames(log2[,1:4]),log2[,1:4]), file = paste0("diagnostics/rwty/log2rwty", m), row.names = F, col.names = T)
          
          #read in trees with log files
          trees1rwty <- load.trees(file = paste0("output/", fileName, "_trees_run_1.trees"), logfile = paste0("diagnostics/rwty/log1rwty", m), type = "revbayes", gens.per.tree = 1)
          trees2rwty <- load.trees(file = paste0("output/", fileName, "_trees_run_2.trees"), logfile = paste0("diagnostics/rwty/log2rwty", m), type = "revbayes", gens.per.tree = 1)
          
          #perform rwty analysis    
          trees.rwty <- analyze.rwty(list(trees1rwty,trees2rwty))
        }
        
        
        #add diagnostics to diagnostics matrix        
        #first five cols are c("diagnostic", "chain", "ntraits", "rateMatrixMisspecification", "replicateNum")
        
        diagnostics <- rbind(
          diagnostics,
          c(deparse(substitute(geweke1)),1,ntraits,rateMatrixMisspecification,replicateNum, geweke1), 
          c(deparse(substitute(ess1)),1,ntraits,rateMatrixMisspecification,replicateNum, ess1),
          c(paste0(deparse(substitute(heiwel1)), "_CVM.p.val"),1,ntraits,rateMatrixMisspecification,replicateNum, heiwel1[1,]),
          c(paste0(deparse(substitute(heiwel1)), "_halfwidth.passed"),1,ntraits,rateMatrixMisspecification,replicateNum, heiwel1[2,]), 
          
          c(deparse(substitute(geweke2)),2,ntraits,rateMatrixMisspecification,replicateNum, geweke2), 
          c(deparse(substitute(ess2)),2,ntraits,rateMatrixMisspecification,replicateNum, ess2), 
          c(paste0(deparse(substitute(heiwel2)), "_CVM.p.val"),2,ntraits,rateMatrixMisspecification,replicateNum, heiwel2[1,]), 
          c(paste0(deparse(substitute(heiwel1)), "_halfwidth.passed"),2,ntraits,rateMatrixMisspecification,replicateNum, heiwel2[2,]), 
          
          c(deparse(substitute(psrf_point)),"both",ntraits,rateMatrixMisspecification,replicateNum, psrf[1,]), 
          c(deparse(substitute(psrf_upper)),"both",ntraits,rateMatrixMisspecification,replicateNum, psrf[2,]), 
          c(deparse(substitute(geweke)),"both",ntraits,rateMatrixMisspecification,replicateNum, geweke), 
          c(deparse(substitute(ess)),"both",ntraits,rateMatrixMisspecification,replicateNum, ess), 
          c(paste0(deparse(substitute(heiwel)), "_CVM.p.val"),"both",ntraits,rateMatrixMisspecification,replicateNum, heiwel[1,]), 
          c(paste0(deparse(substitute(heiwel)), "_halfwidth.passed"),"both",ntraits,rateMatrixMisspecification,replicateNum, heiwel[2,])
        )
        
        comptrees <- rbind(
          comptrees,
          c(ntraits,rateMatrixMisspecification,replicateNum, comptrecor$estimate, comptrecor$p.value)
        )
        
      }
    }
  }
  
  #for single run
  cat("writing data to disk...\n")
  
  write.table(x = diagnostics, file = paste0("diagnostics/diagnostics", m, ".txt"))
  write.table(x = comptrees, file = paste0("diagnostics/comptrees", m, ".txt"))
  
}

diagnostics <- read.table(file = paste0("diagnostics/diagnostics", 1, ".txt"))
for(i in 2:100){
  print(i)
  if(file.exists(paste0("diagnostics/diagnostics", i, ".txt"))){
    diagnostics <- rbind(diagnostics, read.table(file = paste0("diagnostics/diagnostics", i, ".txt")))
  } else (print(paste0("file ", i, " was not found???")))
}

# max(gelrub, na.rm = T)
# sum(geweke < qnorm(0.025) | geweke > qnorm(0.975), na.rm = T) / (length(geweke) - sum(is.na(geweke)))

# 
# #read back in the files
# MCMC_Diagnostics_PSRF <- read.table(file = "analyses/MCMC_Diagnostics_PSRF1")[-1,]
# MCMC_Diagnostics_Geweke <- read.table(file = "analyses/MCMC_Diagnostics_Geweke1")[-1,]
# MCMC_Diagnostics_ESS <- read.table(file = "analyses/MCMC_Diagnostics_ESS1")[-1,]
# MCMC_Diagnostics_NodalFreqCorr <- read.table(file = "analyses/MCMC_Diagnostics_NodalFreqCorr1")[-1,]
# Tree_Dist <- read.table(file = "analyses/Tree_Dist1")[-1,]
# MJCT_Compatibility <- read.table(file = "analyses/MJCT_Compatibility1")[-1,]
# 
# for (i in 2:100){
#   print(i)
#   MCMC_Diagnostics_PSRF <- rbind(MCMC_Diagnostics_PSRF, read.table(file = paste0("analyses/MCMC_Diagnostics_PSRF", i))[-1,])
#   MCMC_Diagnostics_Geweke <- rbind(MCMC_Diagnostics_Geweke, read.table(file = paste0("analyses/MCMC_Diagnostics_Geweke", i))[-1,])
#   MCMC_Diagnostics_ESS <- rbind(MCMC_Diagnostics_ESS, read.table(file = paste0("analyses/MCMC_Diagnostics_ESS", i))[-1,])
#   MCMC_Diagnostics_NodalFreqCorr <- rbind(MCMC_Diagnostics_NodalFreqCorr, read.table(file = paste0("analyses/MCMC_Diagnostics_NodalFreqCorr", i))[-1,])
#   Tree_Dist <- rbind(Tree_Dist, read.table(file = paste0("analyses/Tree_Dist", i))[-1,])
#   MJCT_Compatibility <- rbind(MJCT_Compatibility, read.table(file = paste0("analyses/MJCT_Compatibility", i))[-1,])
# }
