library(geiger)
library(phytools)

ntips <- 450
ntraits <- 1
delta2 <- -1
lambda <- 4

t <- 0:(100*treeHeight)/100
r <- sapply(1:length(t), function(x) if(t[x] > burstLocation){1 + delta2*exp(-lambda*(t[x]-burstLocation))} else{1})
plot(t,r,type = "l", ylim = c(0, 1+qexp(p = 0.99, rate = 1)), lwd = 3, col = "red")

opacity <- 0.3
nsamp <- 500

disparities <- matrix(0, nrow = nsamp, ncol = ntips)
times <- matrix(0, nrow = nsamp, ncol = ntips)

for(samp in 1:nsamp){
  print(samp)
  
  #simulate tree 
  tree <- pbtree(n = ntips)
  edgeStartEndLength <- nodeHeights(tree)
  edgeStartEndLength <- cbind(edgeStartEndLength, edgeStartEndLength[,2] - edgeStartEndLength[,1])
  treeHeight <- nodeheight(tree, 1)
  burstLocation <- treeHeight / 2
  rootState <- rep(0, ntraits)
  baseRates <- diag(ntraits)
  
  #compute branch rates
  branchRates <- rep(1, length(tree$edge.length))
  for(i in 1:length(branchRates)){
    start <- edgeStartEndLength[i,1] 
    end <- edgeStartEndLength[i,2]
    if(start > burstLocation){
      branchRates[i] <- branchRates[i] + 
        (delta2 / lambda * (exp(-lambda*(start - burstLocation)) - exp(-lambda*(end - burstLocation))))/edgeStartEndLength[i,3]
    } else if(end > burstLocation) {
      branchRates[i] <- branchRates[i] + 
        (delta2 / lambda * (1 - exp(-lambda*(end - burstLocation))))/edgeStartEndLength[i,3]
    }
  }
  
  #simulate traits
  rescTree <- tree
  rescTree$edge.length <- rescTree$edge.length*branchRates
  treeCov <- vcv(rescTree)
  traits <- as.vector(rmvnorm(n = 1, mean = rep(rootState, ntips), sigma = kronecker(treeCov, baseRates)))
  traits <- as.matrix(sapply(1:ntips, function(x) traits[(((x-1)*ntraits)+1):(x*ntraits)]))
  rownames(traits) <- rownames(treeCov)
  
  #compute disparity through time graph
  dtt.nums <- dtt(phy = tree, data = traits, index = c("avg.sq", "avg.manhattan")[1], plot = F)
  if(samp == 1){
    plot(dtt.nums$times, dtt.nums$dtt, lwd = 2, col = rgb(0,0,0,opacity), xlim = c(0,1), ylim = c(0,1.25), 
         type = "l", xlab = "time", ylab = "average subclade disparity", 
         main = paste0("ntips = ", ntips, ", ntraits = ", ntraits, ", delta = ", delta2, ", base rate = ", baseRates[1,1], ", lambda = ", lambda))
    segments(x0 = burstLocation / treeHeight, x1 = burstLocation / treeHeight, y0 = 0.03, y1 = 1.25, lwd = 4, col = "red", lty = 2)
    text("BURST", x = burstLocation / treeHeight, y = 0, font = 2)
    } else {
    lines(dtt.nums$times, dtt.nums$dtt, lwd = 2, col = rgb(0,0,0,opacity))
  }
  if(samp == nsamp){
    segments(x0 = burstLocation / treeHeight, x1 = burstLocation / treeHeight, y0 = 0.03, y1 = 1.25, lwd = 4, col = "red", lty = 2)
    text("BURST", x = burstLocation / treeHeight, y = 0, font = 2)
  }
  disparities[samp,] <- dtt.nums$dtt
  times[samp,] <- dtt.nums$times
  
}

mean.disparity <- sapply(1:ntips, function(x) mean(disparities[,x]))
mean.time <- sapply(1:ntips, function(x) mean(times[,x]))

lines(mean.time, y = mean.disparity, lwd = 4, col = 3)
legend(x = 0.85, y = 1.285, col = 3, lwd = 3, legend = "mean disparity")
