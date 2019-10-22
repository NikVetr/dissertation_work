library(ape)
library(RColorBrewer)

setwd(dir = "/Volumes/2TB/500/full/mvBM_sims_Howells_500/")
setwd(dir = "/Volumes/2TB/500/full/uvBM_sims_Howells_500/")
traitNumIter <- c(2,4,8,16,32,57)
degreeMisspecification <- c("perfect", "interregional", "chimpHuman")
noisyness <- "noiseless"


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



probs <- seq(0.5, 0.999, by = .002)
posteriorBipartsCAR <- matrix(data = 0, ncol = length(probs)*2+3, nrow = 500*length(traitNumIter)*length(degreeMisspecification))
colnames(posteriorBipartsCAR) <- c("ntraits", "correlation", "replicate", paste0(probs, " ID"), paste0(probs, " Correct"))

counter <- 0
for (k in 1:length(degreeMisspecification)) { #rateMatrix
  
  for (j in 1:length(traitNumIter)) { #trait number
    
    for(i in 1:500){ #replicate
      counter <- counter + 1
      cat(c(k,j,i, "... \n"))
      
      ntraits <- traitNumIter[j]
      rateMatrixMisspecification <- degreeMisspecification[k]
      replicateNum <- i 
      
      fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixMisspecification, "RateMatrix_replicate_", replicateNum)
      load(file = paste0("bipartitionPosteriors/", fileName, "_bipartitionProbs"))
      posteriorBiparts <- posteriorBiparts[-1,]
      
      #get internality measure
      internality <- sapply(1:length(posteriorBiparts[,1]), function(x) length(posteriorBiparts[x,2][[1]]))
      internality[internality > 15] <- 30 - internality[internality > 15]
      posteriorBiparts$internality <- internality
      posteriorBiparts$ntraits <- ntraits
      posteriorBiparts$CMS <- rateMatrixMisspecification
      posteriorBiparts$rep <- replicateNum
      
      nNodes <- sapply(1:length(probs), function(x){
        above <- posteriorBiparts[posteriorBiparts[,1] > probs[x],]
        c(length(above[,1]), sum(as.numeric(above$trueState)))
      })
      posteriorBipartsCAR[counter,] <-c(ntraits, k, replicateNum, nNodes[1,], nNodes[2,])
      
    }
  }
}
write.table(x = posteriorBipartsCAR, file = "analyses/PosteriorBipartsCAR_Howells")
# write.table(x = posteriorBipartsCAR, file = "analyses/PosteriorBipartsCAR_PB_NoError")
# write.table(x = posteriorBipartsCAR, file = "analyses/PosteriorBipartsCAR_PB_Error")
setwd(dir = "/Volumes/2TB/500/full/mvBM_sims_Howells_500/")
setwd(dir = "/Volumes/2TB/500/full/uvBM_sims_Howells_500/")
posteriorBipartsCAR <- read.table(file = "analyses/PosteriorBipartsCAR_Howells")


probs <- seq(0.5, 0.999, by = .002)
degreeMisspecification <- c("perfect", "interregional", "chimpHuman")
degreeMisspecificationNames <- c("Perfect", "Interregional", "Interspecific")
traitNumIter <- c(2,4,8,16,32,57)
maxUnrootedResolution <- 27
par(mfrow = c(1,1))
for(q in 1:3){
  cols <- brewer.pal(6, "Dark2")
  CMS <-q
  posteriorBipartsCAR_filtered <- posteriorBipartsCAR[posteriorBipartsCAR[,1] ==  traitNumIter[1] & posteriorBipartsCAR[,2] == CMS,]
  avgNumNodes <- sapply(1:(length(probs)*2), function(x) mean(posteriorBipartsCAR_filtered[,x+3]))
  plot(probs, avgNumNodes[1:(length(avgNumNodes)/2)] / maxUnrootedResolution, type = "l", lwd = 2, xlim = c(0.5,1), ylim = c(0,1), yaxt = "n", 
       main = paste0("Cumulative Average Resolution Curve (Howells' Trees)\nMisspecification Condition: ", degreeMisspecificationNames[CMS], ", Model: ", strsplit(strsplit(getwd(), split = "/")[[1]][4], split = "_")[[1]][1]), 
       xlab = "Minimum Probability", ylab = paste0("Proportion of Nodes Resolved (of ", maxUnrootedResolution, ")"), col = cols[1])
  axis(side = 2 ,at = seq(0,1,0.1))
  lines(probs, avgNumNodes[(1+(length(avgNumNodes)/2)):length(avgNumNodes)] / maxUnrootedResolution, col = cols[1], lty = 3)
  for (k in CMS:CMS) { #misspecification
    
    for (j in 2:6) { #trait number
      
      ntraits <- traitNumIter[j]
      rateMatrixMisspecification <- degreeMisspecification[k]
      posteriorBipartsCAR_filtered <- posteriorBipartsCAR[posteriorBipartsCAR[,1] ==  ntraits & posteriorBipartsCAR[,2] == k,]
      avgNumNodes <- sapply(1:(length(probs)*2), function(x) mean(posteriorBipartsCAR_filtered[,x+3]))
      
      lines(probs, avgNumNodes[1:(length(avgNumNodes)/2)] / maxUnrootedResolution, type = "l", xlim = c(0.5,1), ylim = c(0,1), lwd = 2, col = cols[j])
      lines(probs, avgNumNodes[(1+(length(avgNumNodes)/2)):length(avgNumNodes)] / maxUnrootedResolution, col = cols[j], lty = 3)
      
    }
  }
  legend(x = "topright",legend = traitNumIter,col=cols,lwd = 2, ncol = 2, title = "nTraits")
}


#getting 3x2 cumulative average resolution curves down
directories <- c("/Volumes/2TB/500/full/mvBM_sims_Howells_500/", "/Volumes/2TB/500/full/uvBM_sims_Howells_500/")
staircase <- T; if(staircase){modernHumans_M_nonlog <- read.tree("/Volumes/macOS/Users/nikolai/output/howells_nonlog.trees")}

#colored version
png(filename = "~/mvbm_manuscript/figures/CARcurvesHowells.png", width = 800, height = 400)
par(mfrow = c(2,3), mar=c(2,2,3,2)) # mar is c(bottom, left, top, right)
for(condition in 1:2){
  setwd(dir = directories[condition])
  posteriorBipartsCAR <- read.table(file = "analyses/PosteriorBipartsCAR_Howells")
  degreeMisspecification <- c("perfect", "interregional", "chimpHuman")
  degreeMisspecificationNames <- c("Perfect", "Interregional", "Interspecific")
  traitNumIter <- c(2,4,8,16,32,57)
  maxUnrootedResolution <- 27
  cols <- brewer.pal(6, "Dark2")
  for(q in 1:3){
    CMS <-q
    posteriorBipartsCAR_filtered <- posteriorBipartsCAR[posteriorBipartsCAR[,1] ==  traitNumIter[1] & posteriorBipartsCAR[,2] == CMS,]
    avgNumNodes <- sapply(1:(length(probs)*2), function(x) mean(posteriorBipartsCAR_filtered[,x+3]))
    plot(probs, avgNumNodes[1:(length(avgNumNodes)/2)] / maxUnrootedResolution, type = "l", lwd = 2, xlim = c(0.5,1), ylim = c(0,1), yaxt = "n", 
         main = paste0("Cumulative Average Resolution Curve (Howells' Trees)\nMisspecification Condition: ", degreeMisspecificationNames[CMS], ", Model: ", strsplit(strsplit(getwd(), split = "/")[[1]][6], split = "_")[[1]][1]), 
         xlab = "Minimum Probability", ylab = paste0("Proportion of Nodes Resolved (of ", maxUnrootedResolution, ")"), col = cols[1])
    axis(side = 2 ,at = seq(0,1,0.1))
    lines(probs, avgNumNodes[(1+(length(avgNumNodes)/2)):length(avgNumNodes)] / maxUnrootedResolution, col = cols[1], lty = 3)
      for (k in CMS:CMS) { #misspecification
      
      for (j in 2:6) { #trait number
        
        ntraits <- traitNumIter[j]
        rateMatrixMisspecification <- degreeMisspecification[k]
        posteriorBipartsCAR_filtered <- posteriorBipartsCAR[posteriorBipartsCAR[,1] ==  ntraits & posteriorBipartsCAR[,2] == k,]
        avgNumNodes <- sapply(1:(length(probs)*2), function(x) mean(posteriorBipartsCAR_filtered[,x+3]))
        
        lines(probs, avgNumNodes[1:(length(avgNumNodes)/2)] / maxUnrootedResolution, type = "l", xlim = c(0.5,1), ylim = c(0,1), lwd = 2, col = cols[j])
        lines(probs, avgNumNodes[(1+(length(avgNumNodes)/2)):length(avgNumNodes)] / maxUnrootedResolution, col = cols[j], lty = 3)
        
        if(staircase){
          #get stairstep CAR curve
          biparts <- prop.part.df(modernHumans_M_nonlog)$postProbs
          stairs <- sapply(500:1000, function(x) sum(x/1000 < biparts)/27)
          lines((500:1000)/1000, stairs)
          
        }
        
      }
    }
    if(condition == 1 & q == 3){
      if(!staircase){
      legend(x = "topright",legend = c(traitNumIter, "true"), col=c(cols,1), lwd = 2, ncol = 2, title = "# Traits", lty = c(rep(1,6),3))
      } else {
      legend(x = "topright",legend = c(traitNumIter, "true", "A-M"), col=c(cols,1,1), lwd = 2, ncol = 2, title = "# Traits", lty = c(rep(1,6),3,1))
      }
    }
  }

}
dev.off()

#b/w version
directories <- c("/Volumes/2TB/500/full/mvBM_sims_Howells_500/", "/Volumes/2TB/500/full/uvBM_sims_Howells_500/")
staircase <- T; if(staircase){modernHumans_M_nonlog <- read.tree("/Volumes/macOS/Users/nikolai/output/howells_nonlog.trees")}

png(filename = "~/mvbm_manuscript/figures/CARcurvesHowells_BW.png", width = 800, height = 400)
par(mfrow = c(2,3), mar=c(2,2,2,1), oma=c(1.75,1.75,0,0)) # mar is c(bottom, left, top, right)
q <- 2
n_adj <- -0.005
cex_adj <- 1.2
probs <- seq(0.5, 0.999, by = .002)
for(condition in 1:2){
  setwd(dir = directories[condition])
  posteriorBipartsCAR <- read.table(file = "analyses/PosteriorBipartsCAR_Howells")
  degreeMisspecification <- c("perfect", "interregional", "chimpHuman")
  degreeMisspecificationNames <- c("Perfect", "Interregional", "Interspecific")
  traitNumIter <- c(2,4,8,16,32,57)
  maxUnrootedResolution <- 27
  cols <- rep(1,6)
  for(q in 1:3){
    CMS <-q
    posteriorBipartsCAR_filtered <- posteriorBipartsCAR[posteriorBipartsCAR[,1] ==  traitNumIter[1] & posteriorBipartsCAR[,2] == CMS,]
    avgNumNodes <- sapply(1:(length(probs)*2), function(x) mean(posteriorBipartsCAR_filtered[,x+3]))
    plot(probs, avgNumNodes[1:(length(avgNumNodes)/2)] / maxUnrootedResolution, type = "l", lwd = 2, xlim = c(0.485,1), ylim = c(0,1), yaxt = "n", 
         main = paste0(degreeMisspecificationNames[CMS], ", ", strsplit(strsplit(getwd(), split = "/")[[1]][6], split = "_")[[1]][1]), 
         xlab = "Minimum Probability", ylab = paste0("Proportion of Nodes Resolved (of ", maxUnrootedResolution, ")"), col = cols[1])
    axis(side = 2 ,at = seq(0,1,0.1))
    lines(probs, avgNumNodes[(1+(length(avgNumNodes)/2)):length(avgNumNodes)] / maxUnrootedResolution, col = cols[1], lty = 3)
    text(2, x = 0.49 + n_adj, y = avgNumNodes[1] / maxUnrootedResolution, font = 2, cex = cex_adj)
    
    for (k in CMS:CMS) { #misspecification
      
      for (j in 2:6) { #trait number
        
        ntraits <- traitNumIter[j]
        rateMatrixMisspecification <- degreeMisspecification[k]
        posteriorBipartsCAR_filtered <- posteriorBipartsCAR[posteriorBipartsCAR[,1] ==  ntraits & posteriorBipartsCAR[,2] == k,]
        avgNumNodes <- sapply(1:(length(probs)*2), function(x) mean(posteriorBipartsCAR_filtered[,x+3]))
        
        lines(probs, avgNumNodes[1:(length(avgNumNodes)/2)] / maxUnrootedResolution, type = "l", xlim = c(0.5,1), ylim = c(0,1), lwd = 2, col = cols[j])
        lines(probs, avgNumNodes[(1+(length(avgNumNodes)/2)):length(avgNumNodes)] / maxUnrootedResolution, col = cols[j], lty = 3)
        text(ntraits, x = (if (ntraits == 128) {0.4875} else if (ntraits > 9) {0.49} else {0.49}) + n_adj, 
             y = avgNumNodes[1] / maxUnrootedResolution, font = 2, cex = cex_adj)
        
        if(staircase){
          #get stairstep CAR curve
          biparts <- prop.part.df(modernHumans_M_nonlog)$postProbs
          stairs <- sapply(500:1000, function(x) sum(x/1000 < biparts)/27)
          lines((500:1000)/1000, stairs, lty = 5)
          
        }
        
      }
    }
    if(condition == 1 & q == 3){
      if(!staircase){
        legend(x = "topright",legend = c("proportion of total nodes inferred", "proportion of true nodes inferred"), 
               col=c(1,1), lwd = 2, ncol = 1, lty = c(2,1))
      } else {
        legend(x = "topright",legend = c("proportion of total nodes inferred", "proportion of true nodes inferred", "empirical"), 
               col=c(1,1,1), lwd = 2, ncol = 1, lty = c(1,3,2))
      }
    }
  }
}
mtext(text="Probability Threshold",side=1,line=.5,outer=TRUE)
mtext(text="Proportion of Nodes Inferred",side=2,line=.5,outer=TRUE)
dev.off()
