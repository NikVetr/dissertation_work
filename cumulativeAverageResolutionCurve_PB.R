setwd(dir = "/Volumes/2TB/500/full/mvBM_sims_PB_noiseless_500/")
setwd(dir = "/Volumes/2TB/500/full/uvBM_sims_PB_noiseless_500/")
setwd(dir = "/Volumes/2TB/500/full/mvBM_sims_PB_noisy_500/")
setwd(dir = "/Volumes/2TB/500/full/uvBM_sims_PB_noisy_500/")


cl <- makeCluster(5, outfile="")
registerDoParallel(cl)
getDoParWorkers()

# traitNumIter <- c(2,4,8,16,32,64,128)
# degreeCorrelation <- c("2", "4", "6", "8", "9")
# noisyness <- strsplit(getwd(), "_")[[1]][startsWith(strsplit(getwd(), "_")[[1]], "nois")]
# 
# probs <- seq(0.5, 0.999, by = .002)
# bipartsCAR <- matrix(data = 0, ncol = length(probs)*2+3, nrow = 500*length(traitNumIter)*length(degreeCorrelation))
# colnames(bipartsCAR) <- c("ntraits", "correlation", "replicate", paste0(probs, " ID"), paste0(probs, " Correct"))
#   
# for (k in 1:length(degreeCorrelation)) { #rateMatrix
degreeCorrelation <- c("2", "4", "6", "8", "9")
if(!dir.exists("analyses")){dir.create("analyses")}
foreach(k=1:length(degreeCorrelation)) %dopar% {
  
  counter <- 0
  traitNumIter <- c(2,4,8,16,32,64,128)
  degreeCorrelation <- c("2", "4", "6", "8", "9")
  noisyness <- strsplit(getwd(), "_")[[1]][startsWith(strsplit(getwd(), "_")[[1]], "nois")]
  
  probs <- seq(0.5, 0.999, by = .002)
  bipartsCAR <- matrix(data = 0, ncol = length(probs)*2+3, nrow = 500*length(traitNumIter))
  colnames(bipartsCAR) <- c("ntraits", "correlation", "replicate", paste0(probs, " ID"), paste0(probs, " Correct"))
  
    
  for (j in 1:length(traitNumIter)) { #trait number
    
    for(i in 1:500){ #replicate
      counter <- counter + 1
      cat(c(k,j,i, "... \n"))
      
      ntraits <- traitNumIter[j]
      rateMatrixCorrelation <- degreeCorrelation[k]
      replicateNum <- i 
      
      fileName <- paste0("simulations_", ntraits, "traits_perfectRateMatrix_", rateMatrixCorrelation, "correlationStrength_", noisyness, "_replicate_", replicateNum)
      load(file = paste0("bipartitionPosteriors/", fileName, "_bipartitionProbs"))
      # biparts <- biparts[-1,] #do only if the "everything:nothing" bipartition was left in there
      
      #get internality measure
      # internality <- sapply(1:length(biparts[,1]), function(x) length(biparts[x,2][[1]]))
      # internality[internality > 15] <- 30 - internality[internality > 15]
      # biparts$internality <- internality
      biparts$internality <- 0
      
      biparts$ntraits <- ntraits
      biparts$CMS <- rateMatrixCorrelation
      biparts$rep <- replicateNum
      
      nNodes <- sapply(1:length(probs), function(x){
        above <- biparts[biparts[,1] > probs[x],]
        c(length(above[,1]), sum(as.numeric(above$trueState)))
      })
      bipartsCAR[counter,] <-c(ntraits, k, replicateNum, nNodes[1,], nNodes[2,])
      
    }
  }
  if(noisyness == "noiseless"){
    write.table(x = bipartsCAR, file = paste0("analyses/bipartsCAR_PB_NoError_", rateMatrixCorrelation))
  } else if(noisyness == "noisy"){
    write.table(x = bipartsCAR, file = paste0("analyses/bipartsCAR_PB_Error_", rateMatrixCorrelation))
  }
}

bipartsCAR <- matrix(data = 0, ncol = length(probs)*2+3, nrow = 500*length(traitNumIter)*length(degreeCorrelation))
colnames(bipartsCAR) <- c("ntraits", "correlation", "replicate", paste0(probs, " ID"), paste0(probs, " Correct"))

noisyness <- strsplit(getwd(), "_")[[1]][startsWith(strsplit(getwd(), "_")[[1]], "nois")]
for(i in 1:length(degreeCorrelation)){
  print(i)
  rateMatrixCorrelation <- degreeCorrelation[i]
  if(noisyness == "noiseless"){
    sub_bipartsCAR <- as.matrix(read.table(file = paste0("analyses/bipartsCAR_PB_NoError_", rateMatrixCorrelation)))
  } else if(noisyness == "noisy"){
    sub_bipartsCAR <- as.matrix(read.table(file = paste0("analyses/bipartsCAR_PB_Error_", rateMatrixCorrelation)))
  }
  bipartsCAR[(500*length(traitNumIter)*(i-1)+1):(i*500*length(traitNumIter)),] <- sub_bipartsCAR
}
if(noisyness == "noiseless"){
  write.table(x = bipartsCAR, file = paste0("analyses/bipartsCAR_PB_NoError"))
} else if(noisyness == "noisy"){
  write.table(x = bipartsCAR, file = paste0("analyses/bipartsCAR_PB_Error"))
}


setwd(dir = "/Volumes/2TB/500/full/mvBM_sims_PB_noiseless_500/")
setwd(dir = "/Volumes/2TB/500/full/uvBM_sims_PB_noiseless_500/")
bipartsCAR <- read.table(file = "analyses/bipartsCAR_PB_NoError")

setwd(dir = "/Volumes/2TB/500/full/mvBM_sims_PB_noisy_500/")
setwd(dir = "/Volumes/2TB/500/full/uvBM_sims_PB_noisy_500/")
bipartsCAR <- read.table(file = "analyses/bipartsCAR_PB_Error")



degreeCorrelation <- c("2", "4", "6", "8", "9")
traitNumIter <- c(2,4,8,16,32,64,128)
maxUnrootedResolution <- 27
par(mfrow = c(1,1))
for(q in 1:5){
cols = c(1,2,3,4,5,6,7)
cols <- brewer.pal(7, "Dark2")
CMS <- q
bipartsCAR_filtered <- bipartsCAR[bipartsCAR[,1] ==  traitNumIter[1] & bipartsCAR[,2] == CMS,]
avgNumNodes <- sapply(1:(length(probs)*2), function(x) mean(bipartsCAR_filtered[,x+3]))
plot(probs, avgNumNodes[1:(length(avgNumNodes)/2)] / maxUnrootedResolution, type = "l", lwd = 2, xlim = c(0.5,1), ylim = c(0,1), yaxt = "n", 
     main = paste0("Cumulative Average Resolution Curve (Pure-Birth Tree)\nCorrelation Strength: 0.", degreeCorrelation[CMS]), xlab = "Minimum Probability", ylab = paste0("Proportion of Nodes Resolved (of ", maxUnrootedResolution, ")"), col = cols[1])
     # main = paste0("Cumulative Average Resolution Curve (Pure-Birth Tree)\nNo Measurement Error"), xlab = "Minimum Probability", ylab = paste0("Proportion of Nodes Resolved (of ", maxUnrootedResolution, ")"), col = cols[1])
axis(side = 2 ,at = seq(0,1,0.1))
lines(probs, avgNumNodes[(1+(length(avgNumNodes)/2)):length(avgNumNodes)] / maxUnrootedResolution, col = cols[1], lty = 3)
for (k in CMS:CMS) { #misspecification
  
  for (j in 2:7) { #trait number
    
    ntraits <- traitNumIter[j]
    rateMatrixCorrelation <- degreeCorrelation[k]
    bipartsCAR_filtered <- bipartsCAR[bipartsCAR[,1] ==  ntraits & bipartsCAR[,2] == k,]
    avgNumNodes <- sapply(1:(length(probs)*2), function(x) mean(bipartsCAR_filtered[,x+3]))
    
    lines(probs, avgNumNodes[1:(length(avgNumNodes)/2)] / maxUnrootedResolution, type = "l", xlim = c(0.5,1), ylim = c(0,1), lwd = 2, col = cols[j])
    lines(probs, avgNumNodes[(1+(length(avgNumNodes)/2)):length(avgNumNodes)] / maxUnrootedResolution, col = cols[j], lty = 3)
    
    }
}
legend(x = "topright",legend = traitNumIter,col=cols,lwd = 2, ncol = 2, title = "nTraits")
}


#getting 2x2 cumulative average resolution curves down
directories <- c("/Volumes/2TB/500/full/mvBM_sims_PB_noiseless_500/",
                 "/Volumes/2TB/500/full/mvBM_sims_PB_noisy_500/",
                 "/Volumes/2TB/500/full/uvBM_sims_PB_noiseless_500/",
                 "/Volumes/2TB/500/full/uvBM_sims_PB_noisy_500/")
par(mfrow = c(2,2), mar=c(3,2,2,1))
for(q in 1:5){
  for(condition in 1:4){
    setwd(dir = directories[condition])
    noisyness <- strsplit(directories[condition], split = "_")[[1]][4]
    if(noisyness == "noiseless"){
      bipartsCAR <- read.table(file = "analyses/bipartsCAR_PB_NoError")
    } else if (noisyness == "noisy") {
      bipartsCAR <- read.table(file = "analyses/bipartsCAR_PB_Error")
    }
    degreeCorrelation <- c("2", "4", "6", "8", "99")
    traitNumIter <- c(2,4,8,16,32,64,128)
    maxUnrootedResolution <- 27
    cols <- brewer.pal(7, "Dark2")
    CMS <-q
    bipartsCAR_filtered <- bipartsCAR[bipartsCAR[,1] ==  traitNumIter[1] & bipartsCAR[,2] == CMS,]
    avgNumNodes <- sapply(1:(length(probs)*2), function(x) mean(bipartsCAR_filtered[,x+3]))
    
    
    plot(probs, avgNumNodes[1:(length(avgNumNodes)/2)] / maxUnrootedResolution, type = "l", lwd = 2, xlim = c(0.5,1), ylim = c(0,1), yaxt = "n", 
    main = paste0(strsplit(strsplit(directories[condition], split = "_")[[1]][1], "/")[[1]][4], ", ", noisyness), 
    xlab = "Minimum Probability", ylab = paste0("Proportion of Nodes Resolved (of ", maxUnrootedResolution, ")"), col = cols[1])

    
    # main = paste0("Cumulative Average Resolution Curve (Pure-Birth Tree)\nCorrelation Strength: 0.", degreeCorrelation[CMS]), xlab = "Minimum Probability", ylab = paste0("Proportion of Nodes Resolved (of ", maxUnrootedResolution, ")"), col = cols[1])
    # main = paste0("Cumulative Average Resolution Curve (Pure-Birth Tree)\nNo Measurement Error"), xlab = "Minimum Probability", ylab = paste0("Proportion of Nodes Resolved (of ", maxUnrootedResolution, ")"), col = cols[1])
    axis(side = 2 ,at = seq(0,1,0.1))
    lines(probs, avgNumNodes[(1+(length(avgNumNodes)/2)):length(avgNumNodes)] / maxUnrootedResolution, col = cols[1], lty = 3)
    for (k in CMS:CMS) { #misspecification
      
      for (j in 2:7) { #trait number
        
        ntraits <- traitNumIter[j]
        rateMatrixCorrelation <- degreeCorrelation[k]
        bipartsCAR_filtered <- bipartsCAR[bipartsCAR[,1] ==  ntraits & bipartsCAR[,2] == k,]
        avgNumNodes <- sapply(1:(length(probs)*2), function(x) mean(bipartsCAR_filtered[,x+3]))
        
        lines(probs, avgNumNodes[1:(length(avgNumNodes)/2)] / maxUnrootedResolution, type = "l", xlim = c(0.5,1), ylim = c(0,1), lwd = 2, col = cols[j])
        lines(probs, avgNumNodes[(1+(length(avgNumNodes)/2)):length(avgNumNodes)] / maxUnrootedResolution, col = cols[j], lty = 3)
        
      }
    }
    if(condition == 2){
      legend(x = "topright",legend = c(traitNumIter, "true"), col=c(cols,1), lwd = 2, ncol = 2, title = "# Traits", lty = c(rep(1,7),3))
    }
  }
  mtext(paste0("Correlation Strength: 0.", degreeCorrelation[CMS]), side = 3, line = -27, outer = TRUE, cex = 1.75, font = 2)
  
}


#getting 0.4 black and white CAR curve for main text
#getting 2x2 cumulative average resolution curves down
directories <- c("/Volumes/2TB/500/full/mvBM_sims_PB_noiseless_500/",
                 "/Volumes/2TB/500/full/mvBM_sims_PB_noisy_500/",
                 "/Volumes/2TB/500/full/uvBM_sims_PB_noiseless_500/",
                 "/Volumes/2TB/500/full/uvBM_sims_PB_noisy_500/")
png(filename = "~/mvbm_manuscript/figures/CARcurvesPBT.png", width = 1000, height = 600)
par(mfrow = c(2,2), mar=c(2,2,2,1), oma=c(1.75,1.75,0,0))
q <- 2
n_adj <- -0.005
cex_adj <- 1.2
probs <- seq(0.5, 0.999, by = .002)
for(condition in 1:4){
  setwd(dir = directories[condition])
  noisyness <- strsplit(directories[condition], split = "_")[[1]][4]
  if(noisyness == "noiseless"){
    bipartsCAR <- read.table(file = "analyses/bipartsCAR_PB_NoError")
  } else if (noisyness == "noisy") {
    bipartsCAR <- read.table(file = "analyses/bipartsCAR_PB_Error")
  }
  degreeCorrelation <- c("2", "4", "6", "8", "99")
  traitNumIter <- c(2,4,8,16,32,64,128)
  maxUnrootedResolution <- 27
  cols <- rep("black", 7)
  CMS <-q
  bipartsCAR_filtered <- bipartsCAR[bipartsCAR[,1] ==  traitNumIter[1] & bipartsCAR[,2] == CMS,]
  avgNumNodes <- sapply(1:(length(probs)*2), function(x) mean(bipartsCAR_filtered[,x+3]))
  
  
  plot(probs, avgNumNodes[1:(length(avgNumNodes)/2)] / maxUnrootedResolution, type = "l", lwd = 2, xlim = c(0.485,1), ylim = c(0,1), yaxt = "n", 
       main = paste0(strsplit(strsplit(directories[condition], split = "_")[[1]][1], "/")[[1]][6], ", ", noisyness), 
       xlab = "Minimum Probability", ylab = paste0("Proportion of Nodes Resolved (of ", maxUnrootedResolution, ")"), col = cols[1])
  
  axis(side = 2 ,at = seq(0,1,0.1))
  text(2, x = 0.4925 + n_adj, y = avgNumNodes[1] / maxUnrootedResolution, font = 2, cex = cex_adj)
  lines(probs, avgNumNodes[(1+(length(avgNumNodes)/2)):length(avgNumNodes)] / maxUnrootedResolution, col = cols[1], lty = 2)
  for (k in CMS:CMS) { #misspecification
    
    for (j in 2:7) { #trait number
      
      ntraits <- traitNumIter[j]
      rateMatrixCorrelation <- degreeCorrelation[k]
      bipartsCAR_filtered <- bipartsCAR[bipartsCAR[,1] ==  ntraits & bipartsCAR[,2] == k,]
      avgNumNodes <- sapply(1:(length(probs)*2), function(x) mean(bipartsCAR_filtered[,x+3]))
      
      lines(probs, avgNumNodes[1:(length(avgNumNodes)/2)] / maxUnrootedResolution, type = "l", xlim = c(0.5,1), ylim = c(0,1), lwd = 2, col = cols[j])
      text(ntraits, x = (if (ntraits == 128) {0.4875} else if (ntraits > 9) {0.49} else {0.4925}) + n_adj, 
           y = avgNumNodes[1] / maxUnrootedResolution, font = 2, cex = cex_adj)
      lines(probs, avgNumNodes[(1+(length(avgNumNodes)/2)):length(avgNumNodes)] / maxUnrootedResolution, col = cols[j], lty = 2)
      
    }
  }
  if(condition == 2){
    legend(x = "topright",legend = c("proportion of total nodes inferred", "proportion of true nodes inferred"), 
           col=c(1,1), lwd = 2, ncol = 1, lty = c(1,2))
  }
}
mtext(text="Probability Threshold",side=1,line=.5,outer=TRUE)
mtext(text="Proportion of Nodes Inferred",side=2,line=.5,outer=TRUE)
mtext(paste0("r = 0.", degreeCorrelation[CMS]), side = 3, line = -26, outer = TRUE, cex = 1.75, font = 2)
dev.off()



#getting 0.4 black and white CAR curve for main text
#swap solid and dashed lines
directories <- c("/Volumes/2TB/500/full/mvBM_sims_PB_noiseless_500/",
                 "/Volumes/2TB/500/full/mvBM_sims_PB_noisy_500/",
                 "/Volumes/2TB/500/full/uvBM_sims_PB_noiseless_500/",
                 "/Volumes/2TB/500/full/uvBM_sims_PB_noisy_500/")
png(filename = "~/mvbm_manuscript/figures/CARcurvesPBT_swapLineStyle.png", width = 1000, height = 600)
par(mfrow = c(2,2), mar=c(2,2,2,1), oma=c(1.75,1.75,0,0))
q <- 2
n_adj <- -0.005
cex_adj <- 1.2
probs <- seq(0.5, 0.999, by = .002)
for(condition in 1:4){
  setwd(dir = directories[condition])
  noisyness <- strsplit(directories[condition], split = "_")[[1]][4]
  if(noisyness == "noiseless"){
    bipartsCAR <- read.table(file = "analyses/bipartsCAR_PB_NoError")
  } else if (noisyness == "noisy") {
    bipartsCAR <- read.table(file = "analyses/bipartsCAR_PB_Error")
  }
  degreeCorrelation <- c("2", "4", "6", "8", "99")
  traitNumIter <- c(2,4,8,16,32,64,128)
  maxUnrootedResolution <- 27
  cols <- rep("black", 7)
  CMS <-q
  bipartsCAR_filtered <- bipartsCAR[bipartsCAR[,1] ==  traitNumIter[1] & bipartsCAR[,2] == CMS,]
  avgNumNodes <- sapply(1:(length(probs)*2), function(x) mean(bipartsCAR_filtered[,x+3]))
  
  
  plot(probs, avgNumNodes[1:(length(avgNumNodes)/2)] / maxUnrootedResolution, type = "l", lwd = 1, lty = 2, xlim = c(0.485,1), ylim = c(0,1), yaxt = "n", 
       main = paste0(strsplit(strsplit(directories[condition], split = "_")[[1]][1], "/")[[1]][6], ", ", noisyness), 
       xlab = "Minimum Probability", ylab = paste0("Proportion of Nodes Resolved (of ", maxUnrootedResolution, ")"), col = cols[1])
  
  axis(side = 2 ,at = seq(0,1,0.1))
  text(2, x = 0.4925 + n_adj, y = avgNumNodes[1] / maxUnrootedResolution, font = 2, cex = cex_adj)
  lines(probs, avgNumNodes[(1+(length(avgNumNodes)/2)):length(avgNumNodes)] / maxUnrootedResolution, col = cols[1], lty = 1, lwd = 2)
  for (k in CMS:CMS) { #misspecification
    
    for (j in 2:7) { #trait number
      
      ntraits <- traitNumIter[j]
      rateMatrixCorrelation <- degreeCorrelation[k]
      bipartsCAR_filtered <- bipartsCAR[bipartsCAR[,1] ==  ntraits & bipartsCAR[,2] == k,]
      avgNumNodes <- sapply(1:(length(probs)*2), function(x) mean(bipartsCAR_filtered[,x+3]))
      
      lines(probs, avgNumNodes[1:(length(avgNumNodes)/2)] / maxUnrootedResolution, type = "l", xlim = c(0.5,1), ylim = c(0,1), lwd = 1, col = cols[j], lty = 2)
      text(ntraits, x = (if (ntraits == 128) {0.4875} else if (ntraits > 9) {0.49} else {0.4925}) + n_adj, 
           y = avgNumNodes[1] / maxUnrootedResolution, font = 2, cex = cex_adj)
      lines(probs, avgNumNodes[(1+(length(avgNumNodes)/2)):length(avgNumNodes)] / maxUnrootedResolution, col = cols[j], lty = 1, lwd = 2)
      
    }
  }
  if(condition == 2){
    legend(x = "topright",legend = c("proportion of total nodes inferred", "proportion of true nodes inferred"), 
           col=c(1,1), lwd = 2, ncol = 1, lty = c(2,1))
  }
}
mtext(text="Probability Threshold",side=1,line=.5,outer=TRUE)
mtext(text="Proportion of Nodes Inferred",side=2,line=.5,outer=TRUE)
mtext(paste0("r = 0.", degreeCorrelation[CMS]), side = 3, line = -26, outer = TRUE, cex = 1.75, font = 2)
dev.off()


