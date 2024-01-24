library(geiger)
library(phytools)

ntips <- 100
ntraits <- 1
delta2 <- 2
lambda <- 3

tree <- pbtree(n = ntips)
edgeStartEndLength <- nodeHeights(tree)
edgeStartEndLength <- cbind(edgeStartEndLength, edgeStartEndLength[,2] - edgeStartEndLength[,1])
treeHeight <- nodeheight(tree, 1)
burstLocation <- treeHeight / 2
rootState <- rep(0, ntraits)
baseRates <- diag(ntraits)

treeHeight <- 20
t <- 0:(100*treeHeight)/100
r <- sapply(1:length(t), function(x) if(t[x] > burstLocation){1 + delta2*exp(-lambda*(t[x]-burstLocation))} else{1})
plot(t,r,type = "l", ylim = c(0, 1+qexp(p = 0.99, rate = 1)), lwd = 3, col = "red", ylab = "evolutionary rate", xlab = "time", cex.lab=1.25)
# title(paste0("ntips = ", ntips, ", ntraits = ", ntraits, ", delta = ", delta2, ", base rate = ", baseRates[1,1], ", lambda = ", lambda))

opacity <- 0.3
nsamp <- 200

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
         type = "l", xlab = "time", ylab = "average subclade disparity", cex.lab = 1.25,
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
legend(x = 0.79, y = 1.285, col = 3, lwd = 3, legend = "mean disparity")

################################################################
# comparing mean disparities across different burst magnitudes #
################################################################

nsamp <- 1000
ntips <- 300
ntraits <- 1
delta2s <- c(-1,-1/2,0,1,2,5,20,100,1000)
mean.disparities <- matrix(0, nrow = length(delta2s), ncol = ntips)
mean.times <- matrix(0, nrow = length(delta2s), ncol = ntips)
lambda <- 4

for(burstsize in 1:length(delta2s)){
  print(burstsize)
  delta2 <- delta2s[burstsize]
  disparities <- matrix(0, nrow = nsamp, ncol = ntips)
  times <- matrix(0, nrow = nsamp, ncol = ntips)
  
  for(samp in 1:nsamp){
    if(samp %% 100 == 0){cat(paste0(" ", samp))}
    
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
    disparities[samp,] <- dtt.nums$dtt
    times[samp,] <- dtt.nums$times
    
  }
  
  mean.disparity <- sapply(1:ntips, function(x) mean(disparities[,x]))
  mean.time <- sapply(1:ntips, function(x) mean(times[,x]))
  
  mean.disparities[burstsize,] <- mean.disparity
  mean.times[burstsize,] <- mean.time

}


pdf(file = "dissertation/figures/expected_disparity_through_time.pdf", width = 14, height = 7)
par(mar = c(4,4.5,3,2), mfrow = c(1,2), xpd = TRUE)

r <- sapply(0:500, function(x) if(x / 100 > 2.5){1 + delta2*exp(-lambda*(x/100-2.5))} else{1})
plot(0:500/500,r,type = "l", ylim = c(0, 4), lwd = 3, col = "red", ylab = "Evolutionary Rate", xlab = "Standardized Time", 
     cex.lab=1.25, cex.lab = 1.5, cex.axis = 1.25, main = "Evolutionary Rate Through Time", cex.main = 2, yaxt = "n")
text("a)", x = -0.15, y = 4.4, cex = 2, font = 3)
axis(side = 2, at = c(1,3), labels = c(TeX("$\\sigma^2_0$"), TeX("$\\sigma^2_0 + \\Delta$")), las = 1, cex.axis = 1.25)
cols <- viridis::plasma(length(delta2s), begin = 0, end = 0.8)
plot(mean.times[1,], y = mean.disparities[1,], lwd = 2, col = cols[1], type = "l", 
     xlab = "Standardized Time", ylab = "Standardized Expected Disparity", cex.lab = 1.5, cex.axis = 1.25, main = "Expected Disparity Through Time", cex.main = 2)
text("b)", x = -0.15, y = 1.1, cex = 2, font = 3)
segments(x0 = 0.5, x1 = 0.5, y0 = 0.025, y1 = 0.99, lty = 2, lwd = 2.25, col = rgb(0,0,0,0.5))
text(x = 0.5, y = 0, labels = TeX("$\\Delta Burst$"), cex = 1.5, col = rgb(0,0,0,0.5))
for(i in 1:length(delta2s)){
  lines(mean.times[i,], y = mean.disparities[i,], lwd = 2, col = cols[i])
}
text(x = 0.615, y = 1.015, labels = TeX(paste0("n$_{tips}$ = ", ntips, ", $\\sigma^2_0$ = ", baseRates[1,1], ", $\\lambda = ", lambda)))
legend(x = "topright", col = cols, lwd = 2, legend = sapply(1:length(delta2s), function(d2) TeX(paste0("$\\Delta$ = ", delta2s[d2]))), box.lty = 2)
dev.off()
