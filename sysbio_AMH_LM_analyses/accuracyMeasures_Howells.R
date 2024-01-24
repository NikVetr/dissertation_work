setwd(dir = "/Volumes/2TB/mvBM_sims_Howells")
# setwd(dir = "/Volumes/2TB/uvBM_sims_Howells")

degreeMisspecification <- c("perfect", "interregional", "chimpHuman")
degreeMisspecificationNames <- c("perfect", "interregional", "interspecific")
traitNumIter <- c(2,4,8,16,32,57)

#direct distance matrix methods
distMatrixUPGMADists <- read.table(file = "analyses/distMatrixUPGMADists.txt"); 
distMatrixUPGMADists$rateMatrixMisspecification <- as.factor(sapply(1:length(distMatrixUPGMADists$rateMatrixMisspecification), function(x) which(distMatrixUPGMADists$rateMatrixMisspecification[x] == degreeMisspecification)))
distMatrixNJDists <- read.table(file = "analyses/distMatrixNJDists.txt")
distMatrixNJDists$rateMatrixMisspecification <- as.factor(sapply(1:length(distMatrixNJDists$rateMatrixMisspecification), function(x) which(distMatrixNJDists$rateMatrixMisspecification[x] == degreeMisspecification)))

#MCC trees
mccDists.mvBM <- read.table(file = "analyses/mccDists.txt"); mccDists.mvBM$type <- "MCC-mvBM"; 
mccDists.mvBM$rateMatrixMisspecification <- mccDists.mvBM$rateMatrixMisspecification + 2
mccDists.mvBM$rateMatrixMisspecification[mccDists.mvBM$rateMatrixMisspecification == 4] <- 2; mccDists.mvBM$rateMatrixMisspecification[mccDists.mvBM$rateMatrixMisspecification == 5] <- 1
mccDists.uvBM <- read.table(file = paste0(gsub(x = getwd(), pattern = "mvBM", replacement = "uvBM"), "/analyses/mccDists.txt")); mccDists.uvBM$type <- "MCC-uvBM"
mccDists.uvBM$rateMatrixMisspecification <- mccDists.uvBM$rateMatrixMisspecification + 2
mccDists.uvBM$rateMatrixMisspecification[mccDists.uvBM$rateMatrixMisspecification == 4] <- 2; mccDists.uvBM$rateMatrixMisspecification[mccDists.uvBM$rateMatrixMisspecification == 5] <- 1

#Maximum Parsimony Trees
maxPars <- read.table("parsRatchetRF/mpRFdists1.txt"); for(i in 2:100){maxPars <- rbind(maxPars, read.table(paste0("parsRatchetRF/mpRFdists", i, ".txt")))}
maxPars <- cbind(maxPars[,c(1,2)], rep("max-Pars", length(maxPars[,1])), maxPars[,c(3,4)], rep(NA, length(maxPars[,1]))); colnames(maxPars) <- colnames(mccDists.mvBM)
maxPars$rateMatrixMisspecification <- match(maxPars$rateMatrixMisspecification, degreeMisspecification)

#all trees
allDists <- rbind(rbind(rbind(rbind(distMatrixUPGMADists, distMatrixNJDists), mccDists.mvBM), mccDists.uvBM), maxPars)
allDistsMeans <- mccDists.mvBM[1:(3*length(degreeMisspecification)*length(traitNumIter)),]

types <- c("NJ", "UPGMA", "max-Pars", "MCC-uvBM", "MCC-mvBM")
counter <- 0
for (i in 1:length(traitNumIter)){
  for(j in 1:length(degreeMisspecification)){
    for(k in 1:length(types)){
      counter <- counter + 1
      cat(paste0(counter, "."))
      
      ntraits <- traitNumIter[i]
      # rateMatrixMisspecification <- degreeMisspecification[j]
      rateMatrixMisspecification <- j
      type <- types[k]
      
      dists <- allDists[allDists$ntraits == ntraits & allDists$rateMatrixMisspecification == rateMatrixMisspecification & allDists$type == type,]
      rfDist <- mean(dists$rfDist)
      kfDist <- mean(dists$kfDist)
      
      allDistsMeans[counter,] <- c(ntraits, rateMatrixMisspecification, type, "all", rfDist, kfDist)
    }    
  }
}
maxRFdist <- 54
allDistsMeans$rfDist <- as.numeric(allDistsMeans$rfDist) / maxRFdist
allDistsMeans$kfDist <- as.numeric(allDistsMeans$kfDist)
allDistsMeans$ntraits <- as.numeric(allDistsMeans$ntraits)

par(mfrow = c(1,3))
cols <- brewer.pal(length(types), "Dark2")
for(j in 1:length(degreeMisspecification)){
  rateMatrixMisspecification <- j
  data <- allDistsMeans[allDistsMeans$type == types[1] & allDistsMeans$rateMatrixMisspecification == rateMatrixMisspecification,c(1,5)]
  data <- data[order(data[,1]),]; data$ntraits <- 1:length(data$ntraits)
  plot(data, type = "l", lwd = 1,
       xlim = c(1,7), ylim = c(0,1), col = cols[1], ylab = "Relative Robinson-Fould's Distance (/54)", xlab = "Number of Characters", xaxt = "n")
  axis(1, at = 1:length(data$ntraits), labels = traitNumIter)
  title(paste0("Howells' trees, ", degreeMisspecificationNames[j] ," Rate Matrix Specification"))
  for(i in 2:length(types)){
    type = types[i]
    data <- allDistsMeans[allDistsMeans$type == type & allDistsMeans$rateMatrixMisspecification == rateMatrixMisspecification,c(1,5)]
    data <- data[order(data[,1]),]; data$ntraits <- 1:length(data$ntraits)
    if(type == "MCC-mvBM") {lines(data, col = cols[i], lwd = 2)} else {lines(data, col = cols[i], lwd = 1)}
  }
  legend(x = "bottomleft",legend = types, col=cols,lwd = 2, ncol = 1, title = "Tree Algorithm")
}


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################


# #modal distance matrix methods
# treeDistsNBLs <- read.table(file = "analyses/univariateModalNBLstrees.txt")
# treeDists <- read.table(file = "analyses/univariateModalvsMCCtrees.txt")
# 
# #average mcmc tree
# load("analyses/MCMC_Diagnostics_Distances_All_inclPrior_thinned")
# par(mfrow = c(1,1))
# means <- MCMC_Diagnostics_Distances[,,1]
# means[,5] <- sapply(1:2200, function (x) mean(as.numeric(MCMC_Diagnostics_Distances[x,5,])))
# means[,6] <- sapply(1:2200, function (x) mean(as.numeric(MCMC_Diagnostics_Distances[x,6,])))
# means <- as.data.frame(means)
# colnames(means) <- c("ntraits", "rateMatrixSpecification", "replicateNum", "idk", "rfDist", "kfDist")
# means$rfDist <- as.numeric(as.character(means$rfDist))
# means$kfDist <- as.numeric(as.character(means$kfDist))
# means