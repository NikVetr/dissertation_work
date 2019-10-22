setwd(dir = "/Volumes/2TB/mvBM_sims_all")
traitNumIter <- c(2,4,8,16,32,57,64)
degreeMisspecification <- c("perfect", "interregional", "chimpHuman")


probs <- seq(0.5, 0.999, by = .002)
posteriorBipartsCAR <- matrix(data = 0, ncol = length(probs)*2+3, nrow = 2100)
colnames(posteriorBipartsCAR) <- c("ntraits", "misspecification", "replicate", paste0(probs, " ID"), paste0(probs, " Correct"))
  
counter <- 0
for (k in 1:3) { #misspecification
  
  for (j in 1:6) { #trait number
    
    for(i in 1:100){ #replicate
      counter <- counter + 1
      cat(c(k,j,i, "... \n"))
      
      ntraits <- traitNumIter[j]
      rateMatrixSpecification <- degreeMisspecification[k]
      replicateNum <- i 
      
      fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_", "replicate_", replicateNum)
      load(file = paste0("bipartitionPosteriors/", fileName, "_bipartitionProbs"))
      posteriorBiparts <- posteriorBiparts[-1,]
      
      #get internality measure
      internality <- sapply(1:length(posteriorBiparts[,1]), function(x) length(posteriorBiparts[x,2][[1]]))
      internality[internality > 15] <- 30 - internality[internality > 15]
      posteriorBiparts$internality <- internality
      posteriorBiparts$ntraits <- ntraits
      posteriorBiparts$RMS <- rateMatrixSpecification
      posteriorBiparts$rep <- replicateNum
      
      nNodes <- sapply(1:length(probs), function(x){
        above <- posteriorBiparts[posteriorBiparts[,1] > probs[x],]
        c(length(above[,1]), sum(as.numeric(above$trueState)))
      })
      posteriorBipartsCAR[counter,] <-c(ntraits, k, replicateNum, nNodes[1,], nNodes[2,])
      
    }
  }
}
# write.table(x = posteriorBipartsCAR, file = "analyses/PosteriorBipartsCAR")
posteriorBipartsCAR <- read.table(file = "analyses/PosteriorBipartsCAR")
posteriorBipartsCAR <- posteriorBipartsCAR[posteriorBipartsCAR[,1] != 0,]

posteriorBipartsCAR_filtered <- posteriorBipartsCAR[posteriorBipartsCAR[,1] ==  64 & posteriorBipartsCAR[,2] == 1,]
avgNumNodes <- sapply(1:(length(probs)*2), function(x) mean(posteriorBipartsCAR_filtered[,x+3]))
plot(probs, avgNumNodes[1:(length(avgNumNodes)/2)], type = "l", xlim = c(0.5,1), ylim = c(0,27), 
     main = "Cumulative Average Resolution Curve", xlab = "Minimum Probability", ylab = "Number of Nodes Resolved")
lines(probs, avgNumNodes[(1+(length(avgNumNodes)/2)):length(avgNumNodes)], col = 3)
lines(probs, avgNumNodes[1:(length(avgNumNodes)/2)] - avgNumNodes[(1+(length(avgNumNodes)/2)):length(avgNumNodes)], col = 2)
posteriorBipartsCAR


maxUnrootedResolution <- 27
cols <- brewer.pal(6, "Dark2")
degreeMisspecification <- c("perfect", "interregional", "chimp-human")
for(q in 1:3){
RMS <-q
posteriorBipartsCAR_filtered <- posteriorBipartsCAR[posteriorBipartsCAR[,1] ==  traitNumIter[1] & posteriorBipartsCAR[,2] == RMS,]
avgNumNodes <- sapply(1:(length(probs)*2), function(x) mean(posteriorBipartsCAR_filtered[,x+3]))
plot(probs, avgNumNodes[1:(length(avgNumNodes)/2)] / maxUnrootedResolution, type = "l", lwd = 2, xlim = c(0.5,1), ylim = c(0,1), yaxt = "n", 
     main = paste0("Cumulative Average Resolution Curve (Howells' Posterior Trees)\n", toupper(degreeMisspecification[RMS]), " Rate Matrix"), xlab = "Minimum Probability", ylab = paste0("Proportion of Nodes Resolved (of ", maxUnrootedResolution, ")"), col = cols[1])
# main = paste0("Cumulative Average Resolution Curve (Pure-Birth Tree)\nNo Measurement Error"), xlab = "Minimum Probability", ylab = paste0("Proportion of Nodes Resolved (of ", maxUnrootedResolution, ")"), col = cols[1])
axis(side = 2 ,at = seq(0,1,0.1))
lines(probs, avgNumNodes[(1+(length(avgNumNodes)/2)):length(avgNumNodes)] / maxUnrootedResolution, col = cols[1], lty = 3)
for (k in RMS:RMS) { #misspecification
  
  for (j in 2:6) { #trait number
    
    ntraits <- traitNumIter[j]
    rateMatrixSpecification <- degreeMisspecification[k]
    posteriorBipartsCAR_filtered <- posteriorBipartsCAR[posteriorBipartsCAR[,1] ==  ntraits & posteriorBipartsCAR[,2] == k,]
    avgNumNodes <- sapply(1:(length(probs)*2), function(x) mean(posteriorBipartsCAR_filtered[,x+3]))
    
    lines(probs, avgNumNodes[1:(length(avgNumNodes)/2)] / maxUnrootedResolution, type = "l", xlim = c(0.5,1), ylim = c(0,27), lwd = 2, 
         main = "Cumulative Average Resolution Curve", xlab = "Minimum Probability", ylab = "Number of Nodes Resolved", col = cols[j])
    lines(probs, avgNumNodes[(1+(length(avgNumNodes)/2)):length(avgNumNodes)] / maxUnrootedResolution, col = cols[j], lty = 3)
    
    }
}
legend(0.956, 0.99,traitNumIter[1:6],col=cols,lwd = 2, title = "nTraits")
}


par(mfrow = c(1,1))
intervalSize <- 0.05
breaks <- seq(0,1,intervalSize)

# for all the rate matrix probs
expected <- rep(0, length(breaks)-1)
observed <- rep(0, length(breaks)-1)
for(i in 1:length(breaks)-1){
  inInterval <- posteriorBipartsReduced[posteriorBipartsReduced[,1] > breaks[i] & posteriorBipartsReduced[,1] < breaks[i+1],]
  expected[i] <- mean(as.numeric(inInterval[,1]))
  observed[i] <- mean(as.numeric(inInterval[,3]))
}
plot(expected, observed, ylab = "Observed Frequency", xlab = "Inferred Probability", 
     main = paste0("Calibration Curve for Binned (", intervalSize, ") Nodal Posterior Probabilities\nPerfect Rate Matrix Specifixation, All Traits Analyses"))
abline(a = 0, b = 1, col = 2, lwd = 2)
