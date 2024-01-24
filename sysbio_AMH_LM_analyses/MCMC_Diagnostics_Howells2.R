setwd(dir = "/Volumes/2TB/500/full/uvBM_sims_Howells_500//")
noisyness <- "noiseless"
dir.create("diagnostics")
dir.create("bipartitionPosteriors")
dir.create("diagnostics/rwty")
dir.create("diagnostics/compareTree")
dir.create("diagnostics/geweke")
dir.create("diagnostics/ess")
dir.create("diagnostics/heiwel")
dir.create("diagnostics/psrf")


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
library(apTreeshape)

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

clade.PA <- function(trees, numClades = NA, allTrees = NA, cladesOfInterest = NA){ 
  
  #unroot rooted strictly bifurcating trees, or else face weird bug in prop.part
  if(trees[[1]]$Nnode == (length(trees[[1]]$tip.label) - 1)){
    trees <- unroot(trees) 
  }
  
  tipLabs <- trees[[1]]$tip.label
  
  if(any(is.na(cladesOfInterest))){
    if(any(is.na(allTrees))){
      if(is.na(numClades)){
        biparts <- prop.part.df(trees)
        numClades <- length(biparts[,1])
      }
    } else {
      biparts <- prop.part.df(allTrees)
    }
    cladeNames <- sapply(1:numClades, function(x) biparts[,2][[x]])
  } else {
    cladeNames <- cladesOfInterest
    numClades <- length(cladeNames)
  }
  
  
  #faster vectorized version
  
  iterBiparts <- lapply(1:length(trees), function(x) prop.part.df(trees[x])[,2])
  iterBiparts <- lapply(1:length(iterBiparts), function(x) sapply(1:length(iterBiparts[[x]]), function(y) paste0(iterBiparts[[x]][[y]], collapse = "-")))
  altCladeNames <- sapply(1:numClades, function(x) paste0(sort(setdiff(tipLabs, cladeNames[[x]])), collapse = "-"))
  cladeNames <- sapply(1:numClades, function(x) paste0(cladeNames[[x]], collapse = "-"))
  cladePA <- t(sapply(1:length(trees), function(x) sapply(1:numClades, function(y) as.numeric(any(iterBiparts[[x]] == cladeNames[[y]] | iterBiparts[[x]] == altCladeNames[[y]])))))
  colnames(cladePA) <- cladeNames
  
  #for loop version
  # cladePA <- matrix(0, nrow = length(trees), ncol = numClades)
  # colnames(cladePA) <- sapply(1:numClades, function(x) paste0(biparts[,2][[x]], collapse = "-"))
  # for(clade in 1:numClades){
  #   print(clade)
  #   for(treeNum in 1:length(trees)){
  #     treeBipart <- prop.part.df(trees[treeNum])[,2]
  #     treeBipart <- sapply(1:length(treeBipart), function(x) paste0(treeBipart[[x]], collapse = "-"))
  #     cladeName <- biparts[,2][[clade]]
  #     altCladeName <- paste0(sort(setdiff(tipLabs, cladeName)), collapse = "-")
  #     cladeName <- paste0(cladeName, collapse = "-")
  #     cladePA[treeNum, clade] <- as.numeric(any(treeBipart == cladeName | treeBipart == altCladeName))
  #   }
  # }
  
  
  cladePA
}


cl <- makeCluster(1, outfile="")
registerDoParallel(cl)
getDoParWorkers()
# stopCluster(cl)

foreach(m=1:500, .packages = c("ape", "coda", "rwty", "phangorn", "lattice", "latticeExtra", "MCMCpack")) %dopar% {

  #Also output marginal histograms and trace plots for continuous parameters, split	frequencies	&	presence/absence for trees
  
  traitNumIter <- c(2,4,8,16,32,57)
  degreeMisspecification <- c("perfect", "interregional", "interspecific")
  degreeMisspecificationNames <- c("perfect", "interregional", "interspecific")
  rwty <- F
  savePostBips <- F
  
  # k iterates over degree of misspecification
  for (k in 1:length(degreeMisspecification)) {
  
    for (j in 1:length(traitNumIter)) { #number of traits
      
      for(i in m:m){
        
        #instantiate matrices
        # diagnostics <- matrix(ncol = 5 + 461)[-1,]
        # colnames(diagnostics) <- c("diagnostic", "chain", "ntraits", "rateMatrixMisspecification", "replicateNum","Posterior","Likelihood","TL",
        #                            "collessIndex","pair1","pair2","pair3","pair4","pair5","pair6","pair7","pair8","pair9","pair10","pair11",
        #                            "pair12","pair13","pair14","pair15","pair16","pair17","pair18","pair19","pair20","pair21","pair22","pair23",
        #                            "pair24","pair25","pair26","pair27","pair28","pair29","pair30","pair31","pair32","pair33","pair34","pair35",
        #                            "pair36","pair37","pair38","pair39","pair40","pair41","pair42","pair43","pair44","pair45","pair46","pair47",
        #                            "pair48","pair49","pair50","pair51","pair52","pair53","pair54","pair55","pair56","pair57","pair58","pair59",
        #                            "pair60","pair61","pair62","pair63","pair64","pair65","pair66","pair67","pair68","pair69","pair70","pair71",
        #                            "pair72","pair73","pair74","pair75","pair76","pair77","pair78","pair79","pair80","pair81","pair82","pair83",
        #                            "pair84","pair85","pair86","pair87","pair88","pair89","pair90","pair91","pair92","pair93","pair94","pair95",
        #                            "pair96","pair97","pair98","pair99","pair100","pair101","pair102","pair103","pair104","pair105","pair106",
        #                            "pair107","pair108","pair109","pair110","pair111","pair112","pair113","pair114","pair115","pair116","pair117",
        #                            "pair118","pair119","pair120","pair121","pair122","pair123","pair124","pair125","pair126","pair127","pair128",
        #                            "pair129","pair130","pair131","pair132","pair133","pair134","pair135","pair136","pair137","pair138","pair139",
        #                            "pair140","pair141","pair142","pair143","pair144","pair145","pair146","pair147","pair148","pair149","pair150",
        #                            "pair151","pair152","pair153","pair154","pair155","pair156","pair157","pair158","pair159","pair160","pair161",
        #                            "pair162","pair163","pair164","pair165","pair166","pair167","pair168","pair169","pair170","pair171","pair172",
        #                            "pair173","pair174","pair175","pair176","pair177","pair178","pair179","pair180","pair181","pair182","pair183",
        #                            "pair184","pair185","pair186","pair187","pair188","pair189","pair190","pair191","pair192","pair193","pair194",
        #                            "pair195","pair196","pair197","pair198","pair199","pair200","pair201","pair202","pair203","pair204","pair205",
        #                            "pair206","pair207","pair208","pair209","pair210","pair211","pair212","pair213","pair214","pair215","pair216",
        #                            "pair217","pair218","pair219","pair220","pair221","pair222","pair223","pair224","pair225","pair226","pair227",
        #                            "pair228","pair229","pair230","pair231","pair232","pair233","pair234","pair235","pair236","pair237","pair238",
        #                            "pair239","pair240","pair241","pair242","pair243","pair244","pair245","pair246","pair247","pair248","pair249",
        #                            "pair250","pair251","pair252","pair253","pair254","pair255","pair256","pair257","pair258","pair259","pair260",
        #                            "pair261","pair262","pair263","pair264","pair265","pair266","pair267","pair268","pair269","pair270","pair271",
        #                            "pair272","pair273","pair274","pair275","pair276","pair277","pair278","pair279","pair280","pair281","pair282",
        #                            "pair283","pair284","pair285","pair286","pair287","pair288","pair289","pair290","pair291","pair292","pair293",
        #                            "pair294","pair295","pair296","pair297","pair298","pair299","pair300","pair301","pair302","pair303","pair304",
        #                            "pair305","pair306","pair307","pair308","pair309","pair310","pair311","pair312","pair313","pair314","pair315",
        #                            "pair316","pair317","pair318","pair319","pair320","pair321","pair322","pair323","pair324","pair325","pair326",
        #                            "pair327","pair328","pair329","pair330","pair331","pair332","pair333","pair334","pair335","pair336","pair337",
        #                            "pair338","pair339","pair340","pair341","pair342","pair343","pair344","pair345","pair346","pair347","pair348",
        #                            "pair349","pair350","pair351","pair352","pair353","pair354","pair355","pair356","pair357","pair358","pair359",
        #                            "pair360","pair361","pair362","pair363","pair364","pair365","pair366","pair367","pair368","pair369","pair370",
        #                            "pair371","pair372","pair373","pair374","pair375","pair376","pair377","pair378","pair379","pair380","pair381",
        #                            "pair382","pair383","pair384","pair385","pair386","pair387","pair388","pair389","pair390","pair391","pair392",
        #                            "pair393","pair394","pair395","pair396","pair397","pair398","pair399","pair400","pair401","pair402","pair403",
        #                            "pair404","pair405","pair406","pair407","pair408","pair409","pair410","pair411","pair412","pair413","pair414",
        #                            "pair415","pair416","pair417","pair418","pair419","pair420","pair421","pair422","pair423","pair424","pair425",
        #                            "pair426","pair427","pair428","pair429","pair430","pair431","pair432","pair433","pair434","pair435","rfDist1",
        #                            "kfDist1","clade1","clade2","clade3","clade4","clade5","clade6","clade7","clade8","clade9","clade10","clade11",
        #                            "clade12","clade13","clade14","clade15","clade16","clade17","clade18","clade19","clade20")
        
        comptrees <- matrix(ncol = 5)[-1,]
        colnames(comptrees) <- c("ntraits", "rateMatrixMisspecification", "replicateNum", "Misspecification", "p.value")
        
        
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
        # log1$TL <- sapply(1:length(log1[,1]), function(x) sum(log1[x,startsWith(colnames(log1), "BL")]))
        # log2$TL <- sapply(1:length(log2[,1]), function(x) sum(log2[x,startsWith(colnames(log2), "BL")]))
        log1$TL <- sapply(1:length(trees1), function(x) sum(trees1[[x]]$edge.length))
        log2$TL <- sapply(1:length(trees2), function(x) sum(trees2[[x]]$edge.length))
        
        #colless index filler -- unrooted trees have no colless index, so let's root with the Bushman outgroup
        log1$collessIndex <- sapply(1:length(trees1), function(x) colless(as.treeshape(root(trees1[[x]], outgroup = "BUSHMAN", resolve.root = T))))
        log2$collessIndex <- sapply(1:length(trees2), function(x) colless(as.treeshape(root(trees2[[x]], outgroup = "BUSHMAN", resolve.root = T))))
        
    
        #compute patristic distances
        npairs <- ( length(trees1[[1]]$tip.label) * length(trees1[[1]]$tip.label) - length(trees1[[1]]$tip.label) ) / 2
        patristicDists1 <- matrix(0, nrow = length(trees1), ncol = npairs)
        colnames(patristicDists1) <- paste0("pair", 1:npairs)
        pairNames <- sort(trees1[[1]]$tip.label)
        pairNames <- unlist(sapply(2:length(pairNames), function(x) sapply(1:(x-1), function(y) paste0(c(pairNames[x], pairNames[y]), collapse = "-"))))
        colnames(patristicDists1) <- pairNames
        for(treeNum in 1:length(trees1)){
          patristicDists <- cophenetic.phylo(trees1[[treeNum]])
          patristicDists <- patristicDists[sort(rownames(patristicDists)),sort(rownames(patristicDists))]
          patristicDists <- patristicDists[upper.tri(patristicDists)]
          patristicDists1[treeNum,] <- patristicDists
        }
    
        npairs <- ( length(trees2[[1]]$tip.label) * length(trees2[[1]]$tip.label) - length(trees2[[1]]$tip.label) ) / 2
        patristicDists2 <- matrix(0, nrow = length(trees2), ncol = npairs)
        colnames(patristicDists2) <- paste0("pair", 1:npairs)
        pairNames <- sort(trees2[[1]]$tip.label)
        pairNames <- unlist(sapply(2:length(pairNames), function(x) sapply(1:(x-1), function(y) paste0(c(pairNames[x], pairNames[y]), collapse = "-"))))
        colnames(patristicDists1) <- pairNames
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
        biparts <- prop.part.df(c(trees1, trees2), cutoff = 0.00005)
        cladeNames <- sapply(1:20, function(x) biparts[,2][[x]])
        cladePAtrees1 <- clade.PA(trees1, cladesOfInterest = cladeNames)
        cladePAtrees2 <- clade.PA(trees2, cladesOfInterest = cladeNames)
        # plot(cladePAtrees1[,3], type = "l")
        
        if(savePostBips){
          cat("processing trees...\n")
          biparts$trueState <- 0
          trueClades <- prop.part.df(unroot(trueTree))[,2]
          
          for(clade in 1:length(trueClades)){
            cladeName <- trueClades[[clade]]
            isThere <- sapply(1:length(biparts[,1]), function(x) identical(cladeName, biparts[x,2][[1]]))
            altCladeName <- sort(setdiff(trueTree$tip.label, cladeName))
            orIsThere <- sapply(1:length(biparts[,1]), function(x) identical(altCladeName, biparts[x,2][[1]]))
            if(any(isThere)){
              biparts$trueState[isThere] <- 1
            } else if (any(orIsThere)) {
              biparts$trueState[orIsThere] <- 1
            } else {
              biparts <- rbind(biparts, c(0, paste(cladeName, collapse = " "), 1))
            }
          }
          save(x = biparts, file = paste0("bipartitionPosteriors/", fileName, "_bipartitionProbs"))
        }
   
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
        
        # diagnostics <- rbind(
        #   diagnostics,
        #   c(deparse(substitute(geweke1)),1,ntraits,rateMatrixMisspecification,replicateNum, geweke1), 
        #   c(deparse(substitute(ess1)),1,ntraits,rateMatrixMisspecification,replicateNum, ess1),
        #   c(paste0(deparse(substitute(heiwel1)), "_CVM.p.val"),1,ntraits,rateMatrixMisspecification,replicateNum, heiwel1[1,]),
        #   c(paste0(deparse(substitute(heiwel1)), "_halfwidth.passed"),1,ntraits,rateMatrixMisspecification,replicateNum, heiwel1[2,]), 
        #   
        #   c(deparse(substitute(geweke2)),2,ntraits,rateMatrixMisspecification,replicateNum, geweke2), 
        #   c(deparse(substitute(ess2)),2,ntraits,rateMatrixMisspecification,replicateNum, ess2), 
        #   c(paste0(deparse(substitute(heiwel2)), "_CVM.p.val"),2,ntraits,rateMatrixMisspecification,replicateNum, heiwel2[1,]), 
        #   c(paste0(deparse(substitute(heiwel1)), "_halfwidth.passed"),2,ntraits,rateMatrixMisspecification,replicateNum, heiwel2[2,]), 
        # 
        #   c(deparse(substitute(psrf_point)),"both",ntraits,rateMatrixMisspecification,replicateNum, psrf[1,]), 
        #   c(deparse(substitute(psrf_upper)),"both",ntraits,rateMatrixMisspecification,replicateNum, psrf[2,]), 
        #   c(deparse(substitute(geweke)),"both",ntraits,rateMatrixMisspecification,replicateNum, geweke), 
        #   c(deparse(substitute(ess)),"both",ntraits,rateMatrixMisspecification,replicateNum, ess), 
        #   c(paste0(deparse(substitute(heiwel)), "_CVM.p.val"),"both",ntraits,rateMatrixMisspecification,replicateNum, heiwel[1,]), 
        #   c(paste0(deparse(substitute(heiwel)), "_halfwidth.passed"),"both",ntraits,rateMatrixMisspecification,replicateNum, heiwel[2,])
        #   )
        
        comptrees <- rbind(
          comptrees,
          c(ntraits,rateMatrixMisspecification,replicateNum, comptrecor$estimate, comptrecor$p.value)
        )
        
        #for single pair of analyses
        cat("writing data to disk...\n")
        
        # write.table(x = diagnostics, file = paste0("diagnostics/diagnostics_for_", fileName, ".txt"))
        write.table(x = comptrees, file = paste0("diagnostics/compareTree/comptrees_for_", fileName, ".txt"))
        
        save(geweke1, file = paste0("diagnostics/geweke/", fileName, "_c1.txt"))
        save(geweke2, file = paste0("diagnostics/geweke/", fileName, "_c2.txt"))
        save(geweke, file = paste0("diagnostics/geweke/", fileName, "_cb.txt"))
        
        save(ess1, file = paste0("diagnostics/ess/", fileName, "_c1.txt"))
        save(ess2, file = paste0("diagnostics/ess/", fileName, "_c2.txt"))
        save(ess, file = paste0("diagnostics/ess/", fileName, "_cb.txt"))
        
        save(heiwel1, file = paste0("diagnostics/heiwel/", fileName, "_c1.txt"))
        save(heiwel2, file = paste0("diagnostics/heiwel/", fileName, "_c2.txt"))
        save(heiwel, file = paste0("diagnostics/heiwel/", fileName, "_cb.txt"))
        
        save(psrf, file = paste0("diagnostics/psrf/", fileName, "_cb.txt"))
        
      }
    }
  }
}

# instantiate matrices
traitNumIter <- c(2,4,8,16,32,57)
degreeMisspecification <- c("perfect", "interregional", "chimpHuman")

diagnostics <- matrix(ncol = 5 + 461, nrow = length(traitNumIter)*length(degreeMisspecification)*500*14)
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

comptrees <- matrix(ncol = 5, nrow = length(traitNumIter)*length(degreeMisspecification)*500)
colnames(comptrees) <- c("ntraits", "rateMatrixMisspecification", "replicateNum", "correlation", "p.value")

counter <- 0
for(i in 1:500){
  
  # cat(c(i, "... \n"))
  
  for (k in 1:length(degreeMisspecification)) {
    
    for (j in 1:length(traitNumIter)) { #number of traits
      
      counter <- counter + 1
      cat(c(i, j, k, "... \n"))
      
      ntraits <- traitNumIter[j]
      rateMatrixMisspecification <- degreeMisspecification[k]
      replicateNum <- i
      
      fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixMisspecification, "RateMatrix_replicate_", replicateNum)
      
      load(file = paste0("diagnostics/geweke/", fileName, "_c1.txt"))
      load(file = paste0("diagnostics/geweke/", fileName, "_c2.txt"))
      load(file = paste0("diagnostics/geweke/", fileName, "_cb.txt"))
      
      load(file = paste0("diagnostics/ess/", fileName, "_c1.txt"))
      load(file = paste0("diagnostics/ess/", fileName, "_c2.txt"))
      load(file = paste0("diagnostics/ess/", fileName, "_cb.txt"))
      
      load(file = paste0("diagnostics/heiwel/", fileName, "_c1.txt"))
      load(file = paste0("diagnostics/heiwel/", fileName, "_c2.txt"))
      load(file = paste0("diagnostics/heiwel/", fileName, "_cb.txt"))
      
      load(file = paste0("diagnostics/psrf/", fileName, "_cb.txt"))    
      
      diagnostics[(((counter-1)*14+1):(counter*14)),] <- rbind(
        c(deparse(substitute(geweke1)),1,ntraits,rateMatrixMisspecification,replicateNum, geweke1[1:461]),
        c(deparse(substitute(ess1)),1,ntraits,rateMatrixMisspecification,replicateNum, ess1[1:461]),
        c(paste0(deparse(substitute(heiwel1)), "_CVM.p.val"),1,ntraits,rateMatrixMisspecification,replicateNum, heiwel1[1,1:461]),
        c(paste0(deparse(substitute(heiwel1)), "_halfwidth.passed"),1,ntraits,rateMatrixMisspecification,replicateNum, heiwel1[2,1:461]),
        
        c(deparse(substitute(geweke2)),2,ntraits,rateMatrixMisspecification,replicateNum, geweke2[1:461]),
        c(deparse(substitute(ess2)),2,ntraits,rateMatrixMisspecification,replicateNum, ess2[1:461]),
        c(paste0(deparse(substitute(heiwel2)), "_CVM.p.val"),2,ntraits,rateMatrixMisspecification,replicateNum, heiwel2[1,1:461]),
        c(paste0(deparse(substitute(heiwel2)), "_halfwidth.passed"),2,ntraits,rateMatrixMisspecification,replicateNum, heiwel2[2,1:461]),
        
        c(deparse(substitute(psrf_point)),"both",ntraits,rateMatrixMisspecification,replicateNum, psrf[1,1:461]),
        c(deparse(substitute(psrf_upper)),"both",ntraits,rateMatrixMisspecification,replicateNum, psrf[2,1:461]),
        c(deparse(substitute(geweke)),"both",ntraits,rateMatrixMisspecification,replicateNum, geweke[1:461]),
        c(deparse(substitute(ess)),"both",ntraits,rateMatrixMisspecification,replicateNum, ess[1:461]),
        c(paste0(deparse(substitute(heiwel)), "_CVM.p.val"),"both",ntraits,rateMatrixMisspecification,replicateNum, heiwel[1,1:461]),
        c(paste0(deparse(substitute(heiwel)), "_halfwidth.passed"),"both",ntraits,rateMatrixMisspecification,replicateNum, heiwel[2,1:461])
      )
      
      comptrees_instance <- read.table(file = paste0("diagnostics/compareTree/comptrees_for_", fileName, ".txt"))
      comptrees[counter,] <- as.matrix(comptrees_instance)
      
    }
  }
}

#redoing comptrees whoopsie
comptrees <- matrix(ncol = 5, nrow = length(traitNumIter)*length(degreeMisspecification)*500); counter <- 0
colnames(comptrees) <- c("ntraits", "rateMatrixMisspecification", "replicateNum", "correlation", "p.value")
for(i in 1:500){ for (k in 1:length(degreeMisspecification)) {for (j in 1:length(traitNumIter)) { 
      counter <- counter + 1; cat(c(i, j, k, "... \n"))
      ntraits <- traitNumIter[j]; rateMatrixMisspecification <- degreeMisspecification[k]; replicateNum <- i
      fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixMisspecification, "RateMatrix_replicate_", replicateNum)
      comptrees[counter,] <- as.matrix(read.table(file = paste0("diagnostics/compareTree/comptrees_for_", fileName, ".txt")))}}}

diagnostics <- as.data.frame(diagnostics)
diagnostics[,6:466] <- sapply(6:466, function(x) as.numeric(as.character(diagnostics[,x])))
save(diagnostics, file = "allDiagnostics")
save(comptrees, file = "compareTreesCorrsPs")

# examine diagnostic output
setwd(dir = "/Volumes/2TB/500/full/mvBM_sims_Howells_500//")
load(file = "allDiagnostics")
load(file = "compareTreesCorrsPs")

#look at ESS and R2
essBothChains <- diagnostics[diagnostics$diagnostic == "ess",]
thresholdESS <- 1000
thresholdR2 <- 0.95
failures <- c(which(essBothChains$Posterior < thresholdESS), 
              which(essBothChains$Likelihood < thresholdESS),
              which(essBothChains$TL < thresholdESS),
              which(essBothChains$collessIndex < thresholdESS),
              which(essBothChains$rfDist1 < thresholdESS),
              which(essBothChains$kfDist1 < thresholdESS),
              which(essBothChains$kfDist1 < thresholdESS),
              unique(unlist(sapply(1:length(which(startsWith(colnames(essBothChains), "pair"))), 
                                   function(x) which(essBothChains[,which(startsWith(colnames(essBothChains), "pair"))[x]] < thresholdESS)))),
              which(as.numeric(comptrees[,4])^2 < thresholdR2))
failures <- sort(unique(failures))
failures <- essBothChains[failures,3:5]
summary(failures$rateMatrixMisspecification)

#write the jobs that need re-doing
counter <- 0
nrun <- 2
for (q in 1:length(failures[,1])){
  for (l in 1:nrun) {
    counter <- counter + 1
    
    ntraits <- failures[q,1]
    rateMatrixMisspecification <- failures[q,2]
    replicateNum <- failures[q,3]
    
    fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixMisspecification, "RateMatrix_", "replicate_", replicateNum)
    
    jobName <- "job_MCMC_Fail.txt"
    if(counter == 1 & l == 1){
      scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
      sink(jobName, append = F)
      cat(scriptPath, "\n", sep = "")
      sink()
    } else {
      scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
      sink(jobName, append = T)
      cat(scriptPath, "\n", sep = "")
      sink()
    }
  }
}


#############################################
#############################################
# re-analyze the jobs that need re-analysis #
#############################################
#############################################

setwd(dir = "/Volumes/2TB/500/full/mvBM_sims_Howells_500/")

traitNumIter <- c(2,4,8,16,32,57)
degreeMisspecification <- c("perfect", "interregional", "chimpHuman")
redosScript <- read.delim("job_MCMC_Fail.txt", header = F)
redosScript <- sapply(1:length(redosScript$V1), function(x) trimws(redosScript$V1[x]))
if(getwd() == "/Volumes/2TB/500/full/mvBM_sims_PB_noiseless_500"){redosScript <- sapply(1:length(redosScript), function(x) strsplit(redosScript[x], " ")[[1]][4])}
nt <- sapply(1:length(redosScript), function(x) strsplit(redosScript[x], split = "_")[[1]][2]); nt <- as.numeric(substr(nt, 1, nchar(nt) - 6))
dm <- substr(sapply(1:length(redosScript), function(x) strsplit(redosScript[x], split = "_")[[1]][3]), 1, 1); 
for(i in 1:length(degreeMisspecification)){dm[startsWith(degreeMisspecification[i], dm)] <- degreeMisspecification[i]}
rn <- as.numeric(sapply(1:length(redosScript), function(x) strsplit(redosScript[x], split = "_")[[1]][5]))
redos <- cbind(nt, dm, rn)
redos <- redos[2*(1:(length(redos[,1])/2)),]

# #move output files over if need be
# redone <- list.files("output2/output")
# file.exists(paste0("scripts/", redone[1]))


#####################################################################
# alternatively, let's try running these diagnostics on the cluster #
#####################################################################

metaScript <- readLines("/Volumes/2TB/500/full/mvBM_sims_Howells_500/diagTemp.R", warn = F)
dir.create("diagnosticScripts")
for(i in 1:length(redos[,1])){
  tempScript <- metaScript
  tempScript[1] <- paste0(tempScript[1], "\"", getwd(), "\"")
  tempScript[2] <- paste0(tempScript[2], redos[i,1])
  tempScript[3] <- paste0(tempScript[3], "\"", redos[i,2], "\"")
  tempScript[4] <- paste0(tempScript[4], redos[i,3])
  writeLines(tempScript, paste0("diagnosticScripts/", i, ".R"))
}

sink(file = "diagnosticJobs.txt")

for(i in 1:length(redos[,1])){
  cat(paste0("cd ", getwd(), "; ",  "/usr/local/bin/RScript ", "diagnosticScripts/", i, ".R\n"))
}   

sink()

#terminal command; can use system() but screen output is messy
cat(paste0("cd ", getwd(), "; parallel --jobs 8 --sshloginfile /Volumes/macOS/Users/nikolai/instances_all --eta --results terminal_output_diags --files < diagnosticJobs.txt"))

#####################################################################
# alternatively... we have now run these diagnostics on the cluster #
#####################################################################

# instantiate matrices to read diagnostic output from redone runs back in


diagnostics <- matrix(ncol = 5 + 461, nrow = length(redos[,1])*14)
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

comptrees <- matrix(ncol = 5, nrow = length(redos[,1]))
colnames(comptrees) <- c("ntraits", "rateMatrixMisspecification", "replicateNum", "correlation", "p.value")

counter <- 0
for(i in 1:length(redos[,1])){
  print(i)    
  
  counter <- counter + 1
  
  ntraits <- as.numeric(redos[i,1])
  rateMatrixMisspecification <- redos[i,2]
  replicateNum <- as.numeric(redos[i,3])
  
  fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixMisspecification, "RateMatrix_replicate_", replicateNum)
  
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
    c(deparse(substitute(geweke1)),1,ntraits,rateMatrixMisspecification,replicateNum, geweke1[1:461]),
    c(deparse(substitute(ess1)),1,ntraits,rateMatrixMisspecification,replicateNum, ess1[1:461]),
    c(paste0(deparse(substitute(heiwel1)), "_CVM.p.val"),1,ntraits,rateMatrixMisspecification,replicateNum, heiwel1[1,1:461]),
    c(paste0(deparse(substitute(heiwel1)), "_halfwidth.passed"),1,ntraits,rateMatrixMisspecification,replicateNum, heiwel1[2,1:461]),
    
    c(deparse(substitute(geweke2)),2,ntraits,rateMatrixMisspecification,replicateNum, geweke2[1:461]),
    c(deparse(substitute(ess2)),2,ntraits,rateMatrixMisspecification,replicateNum, ess2[1:461]),
    c(paste0(deparse(substitute(heiwel2)), "_CVM.p.val"),2,ntraits,rateMatrixMisspecification,replicateNum, heiwel2[1,1:461]),
    c(paste0(deparse(substitute(heiwel2)), "_halfwidth.passed"),2,ntraits,rateMatrixMisspecification,replicateNum, heiwel2[2,1:461]),
    
    c(deparse(substitute(psrf_point)),"both",ntraits,rateMatrixMisspecification,replicateNum, psrf[1,1:461]),
    c(deparse(substitute(psrf_upper)),"both",ntraits,rateMatrixMisspecification,replicateNum, psrf[2,1:461]),
    c(deparse(substitute(geweke)),"both",ntraits,rateMatrixMisspecification,replicateNum, geweke[1:461]),
    c(deparse(substitute(ess)),"both",ntraits,rateMatrixMisspecification,replicateNum, ess[1:461]),
    c(paste0(deparse(substitute(heiwel)), "_CVM.p.val"),"both",ntraits,rateMatrixMisspecification,replicateNum, heiwel[1,1:461]),
    c(paste0(deparse(substitute(heiwel)), "_halfwidth.passed"),"both",ntraits,rateMatrixMisspecification,replicateNum, heiwel[2,1:461])
  )
  
  comptrees[counter,] <- as.matrix(read.table(file = paste0("diagnostics/compareTree/comptrees_for_", fileName, ".txt")))
}

diagnostics <- as.data.frame(diagnostics)
diagnostics[,6:466] <- sapply(6:466, function(x) as.numeric(as.character(diagnostics[,x])))
save(diagnostics, file = "redoDiagnostics")
save(comptrees, file = "redoCompareTreesCorrsPs")


#loading diagnostics back in

setwd(dir = "/Volumes/2TB/500/full/mvBM_sims_Howells_500/") 
setwd(dir = "/Volumes/2TB/500/full/uvBM_sims_Howells_500/") 
print(paste0("# of original failures: ", length(read.delim("job_MCMC_Fail.txt", header = F)[[1]])/2))
load(file = "redoDiagnostics")
load(file = "redoCompareTreesCorrsPs")

#look at ESS and R2
essBothChains <- diagnostics[diagnostics$diagnostic == "ess",]
thresholdESS <- 1000
thresholdR2 <- 0.95
failures <- c(which(essBothChains$Posterior < thresholdESS), 
              which(essBothChains$Likelihood < thresholdESS),
              which(essBothChains$TL < thresholdESS),
              which(essBothChains$collessIndex < thresholdESS),
              which(essBothChains$rfDist1 < thresholdESS),
              which(essBothChains$kfDist1 < thresholdESS),
              which(essBothChains$kfDist1 < thresholdESS),
              unique(unlist(sapply(1:length(which(startsWith(colnames(essBothChains), "pair"))), 
                                   function(x) which(essBothChains[,which(startsWith(colnames(essBothChains), "pair"))[x]] < thresholdESS)))),
              which(as.numeric(comptrees[,4])^2 < thresholdR2))
failures <- sort(unique(failures))
failures <- essBothChains[failures,3:5]
summary(failures$rateMatrixMisspecification)
print(paste0("# of original failures: ", length(read.delim("job_MCMC_Fail.txt", header = F)[[1]])/2))






### OKKKK LET'S JUST RUN COMPLETELY NEW ANALYSES FOR THE ONES THAT ARE STILL FAILING ###
### OKKKK LET'S JUST RUN COMPLETELY NEW ANALYSES FOR THE ONES THAT ARE STILL FAILING ###
### OKKKK LET'S JUST RUN COMPLETELY NEW ANALYSES FOR THE ONES THAT ARE STILL FAILING ###
### OKKKK LET'S JUST RUN COMPLETELY NEW ANALYSES FOR THE ONES THAT ARE STILL FAILING ###
### OKKKK LET'S JUST RUN COMPLETELY NEW ANALYSES FOR THE ONES THAT ARE STILL FAILING ###



d.orig <- read.csv("Howell.csv")
d <- d.orig
# str(d)
# d[d == 0] <- NA
# d <- d[complete.cases(d),]
men <- d[d$Sex == "M",]
#men <- d
#head(as.matrix(men))
#unique(men$Population)[1]
pops <- sapply(1:length(unique(men$Population)), function (x) as.matrix((t(men[men$Population == unique(men$Population)[x],-(1:4)]))))

#get traits
traits <- list()
for(i in 1:82){ 
  print(i)
  trait <- list()
  trait[[1]] <- pops[[1]][i,]
  for (j in 2:30){
    trait[[j]] <- pops[[j]][i,]
  } 
  traits[[i]] <- trait
}

#make cov matrix
cov <- matrix(nrow = 82, ncol = 82)
for (i in 1:82){
  print(i)
  for(j in 1:82){
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

npop <- length(pops)
#compute matrix of population means
traitsHOWELL <- matrix(nrow=npop, ncol=82)
for(i in 1:npop){
  for(j in 1:82){
    traitsHOWELL[i,j] <- mean(pops[[i]][j,])
  }
}
rownames(traitsHOWELL) <- unique(men$Population)
colnames(traitsHOWELL) <- colnames(men)[-(1:4)]

rownames(cov) <- colnames(cov) <- colnames(men)[-(1:4)]
covHuman <- cov

#taken from roseman and weaver 2004, only BPC -> BPL (due to a typo? it's the basion-prosthion length )
linMeasNames <- c("GOL", "NOL", "BNL", "BBH", "XCB", "XFB", "STB", "ZYB", "AUB", "WCB", "ASB", "BPL", "NPH", "NLH", "OBH", "OBB", "JUB", "NLB", "MAB", "MDH", "MDB", "ZMB", "SSS", "FMB", "NAS", "EKB", "DKB", "NDS", "WNB", "SIS", "IML", "XML", "MLS", "WMH", "FOL", "FRC", "FRS", "FRF", "PAC", "PAS", "PAF", "OCC", "OCS", "OCF", "VRR", "NAR", "SSR", "PRR", "DKR", "ZOR", "FMR", "EKR", "ZMR", "AVR", "DKS", "SOS", "GLS")
linMeasTraitsHuman <- traitsHOWELL[,linMeasNames]
covFULL <- linMeasCovHuman <- covHuman[linMeasNames, linMeasNames]
targetRegionalMBD <- c(0.3016502, 0.5492471, 1.0165738, 1.9375714, 3.8145961, 7.8122877) #(2,4,8,16,32,64)
targetRegionalMBD <- c(0.3016502, 0.5492471, 1.0165738, 1.9375714, 3.8145961, 6.939882, 7.8122877) #incl 57
targetRegionalMBD <- c(0.3016502, 0.5492471, 1.0165738, 1.9375714, 3.8145961, 6.939882, 7.8122877) #(2,4,8,16,32,57)

targetChimpHumanMBD <- c(0.7828651, 1.2172830, 1.9157832, 3.1369000, 5.3419662, 7.6690310, 8.0803367) #(2,4,8,16,32,57,64)
targetChimpHumanMBD <- c(0.7828651, 1.2172830, 1.9157832, 3.1369000, 5.3419662, 7.6690310) #(2,4,8,16,32,57)


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


remote <- TRUE

# CREATING REVBAYES REPLICATE SCRIPTS

library(ape)
library(coda)
library(rwty)
library(phangorn)
library(lattice)
library(latticeExtra)
library(MCMCpack)
library(mvMORPH)
library(phytools)
library(TESS)

failedRedos <- failures
failedRedos[,3] <- as.numeric(as.character(failedRedos[,3]))
failedRedos_temp <- failedRedos; failedRedos_temp[,3] <- as.numeric(as.character(failedRedos_temp[,3])) + 500; 
failedRedos <- rbind(failedRedos, failedRedos_temp); rm(failedRedos_temp)

#incl 57 but not 64
dfMatrix <- matrix(c(0,0,0,0,0,0,58,67,76,86,98,113,10,16,26,41,62,101), nrow = 6, ncol = 3)


#getting the trees
howellsTrees1 <- read.tree(file = "~/output/howells_nonlog.trees")

howellsTrees2 <- read.tree(file = "~/output/howells_nonlog2.trees")

trees <- c(howellsTrees1, howellsTrees2)


traitNumIter <- c(2,4,8,16,32,57)
degreeMisspecification <- c("perfect", "interregional", "chimpHuman")
ntips <- 30
nrun <- 2
noisyness <- "noiseless"

counter <- 0
for(i in 1:length(failedRedos[,1])){
  

  print(c(i, i / length(failedRedos[,1]) * 100)); counter <- counter + 1
    
  ntraits <- as.numeric(as.character(failedRedos[i,1]))
  rateMatrixMisspecification <- as.character(failedRedos[i,2])
  replicateNum <- 500 + as.numeric(as.character(failedRedos[i,3]))
  
  
  if (ntraits < 58){
    traitsToUse <- sample(linMeasNames, size = ntraits, replace = F)
  } else {
    traitsToUse <- linMeasNames
  }
  cov <- covFULL[traitsToUse, traitsToUse]
  
  # generate more traits if needed

  covOrig <- cov
  
  # wiggle traits if needed
  
  if (rateMatrixMisspecification != "perfect"){
    dfWish <- dfMatrix[which(traitNumIter == ntraits),which(degreeMisspecification == rateMatrixMisspecification)]
    cov <- rwish(dfWish, cov)/dfWish
  }
  
  if(getwd() == "/Volumes/2TB/500/full/uvBM_sims_Howells_500"){
    cov <- diag(diag(cov))
  }
  
  fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixMisspecification, "RateMatrix_", "replicate_", replicateNum)
  
  # write the covariance matrices
  ratesPath <- paste0("rates/", fileName, "_rates.Phy")
  ratesPhy(round(covOrig, 4), ratesPath)
  
  ratesPathTSV <- paste0("rates/", fileName, "_rates.tsv")
  ratesTSV(round(covOrig, 4), ratesPathTSV)
  
  # simulate the traits and write to file
  traitsPath <- paste0("data/", fileName, "_traits.nex")
  
  #use howells trees
  tree <- trees[[sample(size = 1, x = 1:length(trees))]]
  if(ntips < 30){
    tipsToKeep <- sample(tree$tip.label, ntips, replace = F)
    tipsToDrop <- setdiff(tree$tip.label, tipsToKeep)
    tree <- drop.tip(phy = tree, tip = tipsToDrop)
    tree <- root(tree, outgroup = tree$tip.label[1], resolve.root = T)
  } else {tree <- root(tree, outgroup = "BUSHMAN", resolve.root = T)}
  
  traits <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = dim(cov)[1], sigma=round(cov, 4), mu= rep(0, times = dim(cov)[1])))
  meansNexus(traits, traitsPath)
  treePath <- paste0("trees/", fileName, "_tree.nex")
  write.nexus(tree, file = treePath, translate = T)
  
  # write the master scripts
  for (l in 1:nrun) {
    scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
    sink(scriptPath, append = F)
    
    cat("## Analysis of Simulated Data:", ntraits, "traits,", rateMatrixMisspecification, "rate matrix specification, replication number", replicateNum, "##\n\n")
    cat("replicateNum <-", replicateNum, "\n")
    cat("runNum <-", l, "\n")
    cat(paste0("rateMatrixSpecification <- \"", rateMatrixMisspecification, "\"\n\n"))
    cat(paste0("noisyness = \"", noisyness, "\"\n\n"))
    cat("mi = 0\nseed(runNum)\n\n")
    
    cat("# read in the data\n")
    cat(paste0("means <- readContinuousCharacterData(\"data/", fileName, "_traits.nex\")\n\n"))
    
    cat("# read in the rate matrix\n")
    cat(paste0("rate = readDistanceMatrix(\"rates/", fileName, "_rates.Phy\")\n"))
    cat("rateMatrix = rate.symmetricMatrix()\n\n")
    
    cat("# get some useful variables from the data\nnumTips = means.ntaxa()\nnumBranches = 2 * numTips - 3\nntraits <- means.nchar()\ntaxa = means.taxa()\nnumNodes = numTips * 2 - 1\n\n")
    
    cat("# specify tree topology prior and moves\nsource(\"simulationsTreeTopology.Rev\")\n\n")
    
    cat("# specify branch length priors and moves\nsource(\"simulationsBranchLengths.Rev\")\n\n")
    
    cat("# assemble the tree\nphylogeny := treeAssembly(topology, br_lens)\n\n")
    
    cat("# put it all together\nmvBNdata ~ dnPhyloMultivariateBrownianREML(phylogeny, rep(1.0, numBranches), rateMatrix, nSites=ntraits, method=\"transform\")\nmvBNdata.clamp(means)\n\n")
    
    cat("# specify monitors\nsource(\"simulationsMonitors_orig.Rev\")\n\n")
    
    cat("# define model\nmymodel = model(phylogeny)\n\n")
    
    cat("# run MCMC\nsource(\"simulationsMCMC_orig.Rev\")\n\n")
    
    cat("# analyze mcmc output\nsource(\"simulationsAnalysis.Rev\")\n\n")
    
    cat("q()")
    
    sink()
    
    # add master script name to list
    
    jobName <- "job_pb_redoFailures.txt"
    if(counter == 1 & l == 1){
      scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
      if(remote){scriptPath <- paste0("cd ", getwd(), "; /Volumes/macOS/Users/nikolai/repos/revbayes-development/projects/cmake/rb ", scriptPath)}
      sink(jobName, append = F)
      cat(scriptPath, "\n", sep = "")
      sink()
    } else {
      scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
      if(remote){scriptPath <- paste0("cd ", getwd(), "; /Volumes/macOS/Users/nikolai/repos/revbayes-development/projects/cmake/rb ", scriptPath)}
      sink(jobName, append = T)
      cat(scriptPath, "\n", sep = "")
      sink()
    }
  }
}    

cat(paste0("cd ", getwd(), "; parallel --jobs 7 --sshloginfile /Volumes/macOS/Users/nikolai/instances_all --eta --results terminal_output_redoFailures --files < job_pb_redoFailures.txt"))


