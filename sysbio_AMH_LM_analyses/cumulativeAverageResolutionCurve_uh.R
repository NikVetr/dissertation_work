setwd(dir = "/Users/nikolai/mvBM_sims_underhand")
traitNumIter <- c(2,4,8,16,32,64)
degreeCorrelation <- c("2", "4", "6", "8", "9")
degreeCorrelation <- c("2", "4", "6", "8")


probs <- seq(0.5, 0.999, by = .002)
posteriorBipartsCAR <- matrix(data = 0, ncol = length(probs)*2+3, nrow = 2400)
colnames(posteriorBipartsCAR) <- c("ntraits", "correlation", "replicate", paste0(probs, " ID"), paste0(probs, " Correct"))
  
counter <- 0
for (k in 5:5) { #rateMatrix
  
  for (j in 1:6) { #trait number
    
    for(i in 201:207){ #replicate
      counter <- counter + 1
      cat(c(k,j,i, "... \n"))
      
      ntraits <- traitNumIter[j]
      rateMatrixCorrelation <- degreeCorrelation[k]
      replicateNum <- i 
      
      fileName <- paste0("simulations_", ntraits, "traits_perfectRateMatrix_", rateMatrixCorrelation, "correlationStrength_", "replicate_", replicateNum)
      load(file = paste0("bipartitionPosteriors/", fileName, "_bipartitionProbs"))
      posteriorBiparts <- posteriorBiparts[-1,]
      
      #get internality measure
      internality <- sapply(1:length(posteriorBiparts[,1]), function(x) length(posteriorBiparts[x,2][[1]]))
      internality[internality > 15] <- 30 - internality[internality > 15]
      posteriorBiparts$internality <- internality
      posteriorBiparts$ntraits <- ntraits
      posteriorBiparts$CMS <- rateMatrixCorrelation
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
posteriorBipartsCAR <- posteriorBipartsCAR[posteriorBipartsCAR[,1] != 0,]

posteriorBipartsCAR_filtered <- posteriorBipartsCAR[posteriorBipartsCAR[,1] ==  64 & posteriorBipartsCAR[,2] == 1,]
avgNumNodes <- sapply(1:(length(probs)*2), function(x) mean(posteriorBipartsCAR_filtered[,x+3]))
plot(probs, avgNumNodes[1:(length(avgNumNodes)/2)], type = "l", xlim = c(0.5,1), ylim = c(0,27), 
     main = "Cumulative Average Resolution Curve", xlab = "Minimum Probability", ylab = "Number of Nodes Resolved")
lines(probs, avgNumNodes[(1+(length(avgNumNodes)/2)):length(avgNumNodes)], col = 3)
lines(probs, avgNumNodes[1:(length(avgNumNodes)/2)] - avgNumNodes[(1+(length(avgNumNodes)/2)):length(avgNumNodes)], col = 2)
posteriorBipartsCAR

par(mfrow = c(1,1))
for(q in 1:5){
cols = c(1,2,3,4,5,6)
CMS <-q
posteriorBipartsCAR_filtered <- posteriorBipartsCAR[posteriorBipartsCAR[,1] ==  traitNumIter[1] & posteriorBipartsCAR[,2] == CMS,]
avgNumNodes <- sapply(1:(length(probs)*2), function(x) mean(posteriorBipartsCAR_filtered[,x+3]))
plot(probs, avgNumNodes[1:(length(avgNumNodes)/2)], type = "l", lwd = 2, xlim = c(0.5,1), ylim = c(0,27), 
     main = paste0("Cumulative Average Resolution Curve (Symmetric Tree)\nCorrelation Strength: 0.", degreeCorrelation[CMS]), xlab = "Minimum Probability", ylab = "Number of Nodes Resolved", col = cols[1])
lines(probs, avgNumNodes[(1+(length(avgNumNodes)/2)):length(avgNumNodes)], col = cols[1], lty = 3)
for (k in CMS:CMS) { #misspecification
  
  for (j in 2:6) { #trait number
    
    ntraits <- traitNumIter[j]
    rateMatrixCorrelation <- degreeCorrelation[k]
    posteriorBipartsCAR_filtered <- posteriorBipartsCAR[posteriorBipartsCAR[,1] ==  ntraits & posteriorBipartsCAR[,2] == k,]
    avgNumNodes <- sapply(1:(length(probs)*2), function(x) mean(posteriorBipartsCAR_filtered[,x+3]))
    
    lines(probs, avgNumNodes[1:(length(avgNumNodes)/2)], type = "l", xlim = c(0.5,1), ylim = c(0,27), lwd = 2, col = cols[j])
    lines(probs, avgNumNodes[(1+(length(avgNumNodes)/2)):length(avgNumNodes)], col = cols[j], lty = 3)
    
    }
}
legend(0.942, 27,traitNumIter,col=cols,lwd = 2, title = "nTraits")
}
