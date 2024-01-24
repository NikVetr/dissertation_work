library(phytools)

t <- 100  # total time
n <- 4  # total taxa
sig2 <- 0.01

b <- (log(n) - log(2))/t
tree <- pbtree(b = b, n = n, t = t, type = "discrete")
plot(tree)
## simulate evolution along each edge
X <- lapply(tree$edge.length, function(x) c(0, cumsum(rnorm(n = x, sd = sqrt(sig2)))))
## reorder the edges of the tree for pre-order traversal
cw <- reorder(tree)
## now simulate on the tree
ll <- tree$edge.length + 1
for (i in 1:nrow(cw$edge)) {
  pp <- which(cw$edge[, 2] == cw$edge[i, 1])
  if (length(pp) > 0) 
    X[[i]] <- X[[i]] + X[[pp]][ll[pp]] else X[[i]] <- X[[i]] + X[[1]][1]
}
## get the starting and ending points of each edge for plotting
H <- nodeHeights(tree)
## plot the simulation
plot(H[1, 1], X[[1]][1], ylim = range(X), xlim = range(H), xlab = "time", ylab = "phenotype")
for (i in 1:(length(X))) lines(H[i, 1]:H[i, 2], X[[i]])
## add tip labels if desired
yy <- sapply(1:length(tree$tip.label), function(x, y) which(x == y), y = tree$edge[, 
                                                                                   2])
yy <- sapply(yy, function(x, y) y[[x]][length(y[[x]])], y = X)
text(x = max(H), y = yy, tree$tip.label)

## simulate a tree
tree<-pbtree(n=5,tip.label=LETTERS)
plotTree(tree)
x<-fastBM(tree)
x
phenogram(tree,x,spread.labels=TRUE)
phenogram(tree,x,spread.labels=TRUE,spread.cost=c(1,0))
fancyTree(tree,type="phenogram95",x=x,spread.cost=c(1,0))

X<-fastBM(tree,nsim=3, theta = c(1,2,3)) ## simulate 3 characters
phylomorphospace(tree,X[,c(1,2)],xlab="trait 1",ylab="trait 2")
phylomorphospace3d(tree,X)

# ----------



#Brownian Motion Gif

startValue <- 5
simL <- 100

branch1 <- startValue
branch1[2:simL] <- branch1[1] + cumsum(rnorm(n=simL-1, sd=1, mean=0))

branch2 <- startValue
branch2[2:round(simL/2)] <- branch2[1] + cumsum(rnorm(n=round(simL/2)-1, sd=1, mean=0))

branch3 <- branch2[length(branch2)] 
branch3[2:(simL - round(simL/2) + 1)] <- branch3[1] + cumsum(rnorm(n=(simL - round(simL/2)), sd=1, mean=0))

branch4 <- branch2[length(branch2)] 
branch4[2:(simL - round(simL/2) + 1)] <- branch4[1] + cumsum(rnorm(n=(simL - round(simL/2)), sd=1, mean=0))

for(i in 1:simL) {
  print(i)
  ymax <- ceiling(max(branch1, branch2, branch3, branch4))
  ymin <- floor(min(branch1, branch2, branch3, branch4))
  
  png(paste0("testplot", i, ".png"))
  plot(branch1[1:i], type = "l", xlim = c(1,simL), ylim = c(ymin,ymax), xlab = "Time", 
       ylab = "Trait Value", main = "Brownian Motion", lwd = 2)
  lines(branch2[1:min(i, round(simL/2))], col = "red", lwd = 2)
  if(i > round(simL/2)) {
    lines(round(simL/2):i, branch3[1:(i-(round(simL/2) - 1))], col = "orange", lwd = 2)
    lines(round(simL/2):i, branch4[1:(i-(round(simL/2) - 1))], col = "purple", lwd = 2)
  }
  dev.off()
}

# ----------

#BM with trend Gif

# setwd(choose.dir())
startValue <- 5
simL <- 100
trend <- .4

branch1 <- startValue
branch1[2:simL] <- branch1[1] + cumsum(rnorm(n=simL-1, sd=1, mean=trend))

branch2 <- startValue
branch2[2:round(simL/2)] <- branch2[1] + cumsum(rnorm(n=round(simL/2)-1, sd=1, mean=trend))

branch3 <- branch2[length(branch2)] 
branch3[2:(simL - round(simL/2) + 1)] <- branch3[1] + cumsum(rnorm(n=(simL - round(simL/2)), sd=1, mean=trend))

branch4 <- branch2[length(branch2)] 
branch4[2:(simL - round(simL/2) + 1)] <- branch4[1] + cumsum(rnorm(n=(simL - round(simL/2)), sd=1, mean=trend))

for(i in 1:simL) {
  print(i)
  ymax <- ceiling(max(branch1, branch2, branch3, branch4))
  ymin <- floor(min(branch1, branch2, branch3, branch4))
  
  png(paste0("testplot", i, ".png"))
  plot(branch1[1:i], type = "l", xlim = c(1,simL), ylim = c(ymin,ymax), xlab = "Time", 
       ylab = "Trait Value", main = "Brownian Motion w/ Trend", lwd = 2)
  lines(branch2[1:min(i, round(simL/2))], col = "red", lwd = 2)
  if(i > round(simL/2)) {
    lines(round(simL/2):i, branch3[1:(i-(round(simL/2) - 1))], col = "orange", lwd = 2)
    lines(round(simL/2):i, branch4[1:(i-(round(simL/2) - 1))], col = "purple", lwd = 2)
  }
  dev.off()
}


# -----------
  

#Ornstein-Uhlenbeck Gif

OUsim <- function(time,theta,alpha,sigma2,x0){
  dw  <- rnorm(time, 0)
  x <- c(x0) # the first trait value is the ancestral value given by x0
  for (i in 2:(time+1)) {
    x[i]  <-  x[i-1] + alpha*(theta-x[i-1]) + sigma2*dw[i-1] # this then calculates the trait value for each subsequent generation 
  }
  return(x);
}

# setwd(choose.dir())
startValue <- 5
simL <- 100
alpha = .1
theta = 25

branch1 <- OUsim(time = simL, theta = theta, alpha = alpha, sigma2 = 1, x0 = startValue)
branch2 <- OUsim(time = round(simL/2) - 1, theta = theta, alpha = alpha, sigma2 = 1, x0 = startValue)
branch3 <- OUsim(time = simL - round(simL/2), theta = theta, alpha = alpha, sigma2 = 1, x0 = branch2[length(branch2)])  
branch4 <- OUsim(time = simL - round(simL/2), theta = theta, alpha = alpha, sigma2 = 1, x0 = branch2[length(branch2)])

for(i in 1:simL) {
  print(i)
  ymax <- ceiling(max(branch1, branch2, branch3, branch4))
  ymin <- floor(min(branch1, branch2, branch3, branch4))
  
  png(paste0("testplot", i, ".png"))
  plot(branch1[1:i], type = "l", xlim = c(1,simL), ylim = c(ymin,ymax), xlab = "Time", 
       ylab = "Trait Value", main = "Ornstein-Uhlenbeck, One Optimum", lwd = 2)
  lines(branch2[1:min(i, round(simL/2))], col = "red", lwd = 2)
  if(i > round(simL/2)) {
    lines(round(simL/2):i, branch3[1:(i-(round(simL/2) - 1))], col = "orange", lwd = 2)
    lines(round(simL/2):i, branch4[1:(i-(round(simL/2) - 1))], col = "purple", lwd = 2)
  }
  if (i < (.9*simL)) text (labels = "Optimum", x = simL - simL/20, y = theta)
  dev.off()
}

# ----------

#Ornstein-Uhlenbeck 2 optima Gif

OUsim <- function(time,theta,alpha,sigma2,x0){
  dw  <- rnorm(time, 0)
  x <- c(x0) # the first trait value is the ancestral value given by x0
  for (i in 2:(time+1)) {
    x[i]  <-  x[i-1] + alpha*(theta-x[i-1]) + sigma2*dw[i-1] # this then calculates the trait value for each subsequent generation 
  }
  return(x);
}

# setwd(choose.dir())
startValue <- 5
simL <- 100
alpha = .1
theta1 = 2
theta2 = 12
sigma2 = .5

branch1 <- OUsim(time = simL, theta = theta1, alpha = alpha, sigma2 = sigma2, x0 = startValue)
branch2 <- OUsim(time = round(simL/2) - 1, theta = theta2, alpha = alpha, sigma2 = sigma2, x0 = startValue)
branch3 <- OUsim(time = simL - round(simL/2), theta = theta2, alpha = alpha, sigma2 = sigma2, x0 = branch2[length(branch2)])  
branch4 <- OUsim(time = simL - round(simL/2), theta = theta2, alpha = alpha, sigma2 = sigma2, x0 = branch2[length(branch2)])

for(i in 1:simL) {
  print(i)
  ymax <- ceiling(max(branch1, branch2, branch3, branch4))
  ymin <- floor(min(branch1, branch2, branch3, branch4))
  
  png(paste0("testplot", i, ".png"))
  plot(branch1[1:i], type = "l", xlim = c(1,simL), ylim = c(ymin,ymax), xlab = "Time", 
       ylab = "Trait Value", main = "Ornstein-Uhlenbeck, Two Optima", lwd = 2)
  lines(branch2[1:min(i, round(simL/2))], col = "red", lwd = 2)
  if(i > round(simL/2)) {
    lines(round(simL/2):i, branch3[1:(i-(round(simL/2) - 1))], col = "orange", lwd = 2)
    lines(round(simL/2):i, branch4[1:(i-(round(simL/2) - 1))], col = "purple", lwd = 2)
  }
  if (i < (.9*simL)) text (labels = "Optimum 1", x = simL - simL/20, y = theta1)
  if (i < (.9*simL)) text (labels = "Optimum 2", col = "red", x = simL - simL/20, y = theta2)
  dev.off()
}


# -----------------

#Early or Late Burst

### a function that requires r which determines the exponential change in the rate over time, time the length of the simulation and the initial rate of 
#Brownian motion, it assumes the starting ancestral trait value is 0

eblb<-function(time, r, sigma0){ 
  displace<-vector(,length=time) # set up an empty vector of the same length as time
  for(i in 1:time){
    e<-exp(1)
    displace[i]<-rnorm(1, mean=0, sd=sigma0*e^(-r*i)) # sample from normal distribution with mean of zero and standard deviation determined by sigma0*e^(-r*i)			
  }
  eblbtraj<-cumsum(displace) # calculates the trajectory over time
  return(eblbtraj) 
}


# setwd(choose.dir())
startValue <- 5
simL <- 100
r = -.04
sigma0 = 1

branch1 <- startValue + eblb(time = simL, r=r, sigma0 = sigma0)
branch1 <- append(startValue, branch1)

branch2 <- startValue + eblb(time = round(simL/2) - 1, r = r, sigma0 = sigma0)
branch2 <- append(startValue, branch2)

branch3 <- branch2[length(branch2)] + 
  eblb(time = simL - round(simL/2), r = r, sigma0 = sigma0*exp(1)^(-r*round(simL/2)))  
branch3 <- append(branch2[length(branch2)], branch3)

branch4 <- branch2[length(branch2)] + 
  eblb(time = simL - round(simL/2), r = r, sigma0 = sigma0*exp(1)^(-r*round(simL/2)))
branch4 <- append(branch2[length(branch2)], branch4)
for(i in 1:simL) {
  print(i)
  ymax <- ceiling(max(branch1, branch2, branch3, branch4))
  ymin <- floor(min(branch1, branch2, branch3, branch4))
  
  png(paste0("testplot", i, ".png"))
  plot(branch1[1:i], type = "l", xlim = c(1,simL), ylim = c(ymin,ymax), xlab = "Time", 
       ylab = "Trait Value", main = "Late Burst", lwd = 2)
  lines(branch2[1:min(i, round(simL/2))], col = "red", lwd = 2)
  if(i > round(simL/2)) {
    lines(round(simL/2):i, branch3[1:(i-(round(simL/2) - 1))], col = "orange", lwd = 2)
    lines(round(simL/2):i, branch4[1:(i-(round(simL/2) - 1))], col = "purple", lwd = 2)
  }
  dev.off()
}


# -------------
  
  # tree

tree <- pbtree(n=3)
str(tree)
tree$edge.length <- c(100, 50, 50, 50)
tree$tip.label <- c("Black sp.", "Orange sp.", "Purple sp.")
plot.phylo(tree, edge.width = 2, edge.color = c("black", "red", "orange", "purple"), 
           type = "cla", font = 2, label.offset = 2, tip.color = c("black", "orange", "purple"))
axisPhylo()
title("Phylogeny Used in Simulations")






#### animation of hill climbing algorithm


#let's try the hill-climbing algorithm out in a bunch of simulations but change the proposal mechanism

#5 tips, 30 traits, correct 64% of the time across 250 simulations
#5 tips, 15 traits, correct 51% of the time across 250 simulations
#80 traits, 128/171 correct; 61/79
#make starting tree a not pbtree? rtree from ape

iterOptim <- 5000
numTips <- 6
numTraits <- 60
change <- rep(0,19)
trueTree <- (pbtree(n=numTips))
trueTree$tip.label <- c("A", "B", "C", "D", "E", "F")
plot(trueTree)
bestTree <- unroot(rtree(n=numTips))
bestTree$tip.label <- c("A", "B", "F", "D", "C", "E")
plot(bestTree)

#bestTree$edge.length <- rep(1, (2*numTips-3))
sigmaBM <- CovMatrix(numTraits)
sigmaBM <- as.matrix(nearPD(sigmaBM)$mat)

  
  
traits <- mvSIM(tree = trueTree, nsim = 1, model = "BM1", param = list(ntraits = numTraits, sigma=sigmaBM, mu= rep(0, numTraits)))
#traits <- traits*0
trueTree <- unroot(trueTree)
bestLL <- PlantQuasiLL(traits, sigmaBM, bestTree, tip=1)
trueLL <- PlantQuasiLL(traits, sigmaBM, trueTree, tip=1)

change <- rep(0,19)


phylos <- list()
logLiks <- vector()
RFDist <- vector()
pathDiff <- vector()
branchDist <- vector()
quadDiff <- vector()
changes <- list()

phylos[[1]] <- bestTree
logLiks[1] <- bestLL
RFDist[1] <- RF.dist(bestTree, trueTree)
pathDiff[1] <- treedist(bestTree, trueTree)[3]
branchDist[1] <- treedist(bestTree, trueTree)[2]
quadDiff[1] <- treedist(bestTree, trueTree)[4]
changes[[1]] <- paste("Tree Rearrangements = ", sum(change[1:12]), 
                      "; Branch Length Changes = ", sum(change[13:19]))

for (i in 1:iterOptim){
  prop <- sample(x = 1:19, size = 1)
  if(prop <= 6){testTree <- rNNI(bestTree)} else
    if(prop > 6 & prop < 13){testTree <- rSPR(bestTree)} else
      if(prop > 12 & prop < 17){testTree <- bestTree
      edge <- sample(1:length(testTree$edge.length), size=1)
      testTree$edge.length[edge] <- testTree$edge.length[edge]*runif(n = 1, min = .8, max = 1.2)} else
        if(prop==18){testTree <- bestTree
        edge <- sample(1:length(testTree$edge.length), size=1)
        testTree$edge.length[edge] <- testTree$edge.length[edge]*runif(n = 1, min = .8, max = 1.2)
        edge <- sample(1:length(testTree$edge.length), size=1)
        testTree$edge.length[edge] <- testTree$edge.length[edge]*runif(n = 1, min = .8, max = 1.2)} else
          if(prop==19){testTree <- bestTree
          testTree$edge.length <- testTree$edge.length*runif(n = 1, min = .8, max = 1.2)} 
  
  testLL <- PlantQuasiLL(traits, sigmaBM, testTree, tip=1)
  if(testLL > bestLL){
    bestLL <- testLL
    bestTree <- testTree
    change[prop] <- change[prop] + 1
    
    phylos[[length(phylos) + 1]] <- bestTree
    logLiks[length(logLiks) + 1] <- bestLL
    RFDist[length(RFDist) + 1] <- RF.dist(bestTree, trueTree)
    pathDiff[length(pathDiff) + 1] <- treedist(bestTree, trueTree)[3]
    branchDist[length(branchDist) + 1] <- treedist(bestTree, trueTree)[2]
    quadDiff[length(quadDiff) + 1] <- treedist(bestTree, trueTree)[4]
    changes[[length(changes) + 1]] <- paste("Tree Rearrangements = ", sum(change[1:12]), 
                          "; Branch Length Changes = ", sum(change[13:19]))
    
    
  }
  if(i %% 100 == 0){print(c(i, paste0("trueLL=",trueLL), paste0("bestLL=",bestLL)))}
  if(i %% 1000 == 0){print(paste("Tree Rearrangements = ", sum(change[1:12]), 
                                 "; Branch Length Changes = ", sum(change[13:19])))
    print(paste0("Robinson-Foulds Distance = ", RF.dist(trueTree, bestTree)))}
}




phylos <- phylos[1:148]
logLiks <- logLiks[1:148]
RFDist <- RFDist[1:148]
quadDiff <- quadDiff[1:148]
branchDist <- branchDist[1:148]
pathDiff <- pathDiff[1:148]
changes <- changes[1:148]

phylo.diff.alter(trueTree, phylos[[148]], numTraits, trueLL, logLiks[148])
simL <- 148

#Likelihood plot with types of changes
par(mfrow = c(1,1))
ymax <- ceiling(max(logLiks))
ymin <- floor(min(logLiks))
plot(logLiks, type = "l", ylab = "Log-Likelihood", lwd = 2, xlab = "Successful Perturbation Number", xlim = c(1,simL), ylim = c(ymin,ymax))
abline(h = trueLL, col = "blue", lwd = 2, lty = 2)
text(labels = "Data-Generating Tree's Log-Likelihood", x = 41, y = -67, col = "navy")
text(labels = changes[[1]], x = 75, y = -144)
title("Climbing Mt. Log-Likelihood")


for(i in 1:148) {
  print(i)
  png(filename = paste0("logLiks", i, ".png"), width = 720, height = 720)
  par(mfrow = c(1,1))
  ymax <- ceiling(max(logLiks))
  ymin <- floor(min(logLiks))
  plot(logLiks[1:i], type = "l", ylab = "Log-Likelihood", lwd = 2, xlab = "Successful Perturbation Number", xlim = c(1,simL), ylim = c(ymin,ymax))
  abline(h = trueLL, col = "blue", lwd = 2, lty = 2)
  text(labels = "Data-Generating Tree's Log-Likelihood", x = 41, y = -67, col = "navy")
  text(labels = changes[[i]], x = 75, y = -144)
  title("Climbing Mt. Log-Likelihood")
  dev.off()
}


#Branch Distance Plot
ymax <- ceiling(max(branchDist)*10)/10
ymin <- floor(min(branchDist))
plot(branchDist, type = "l", ylab = "Branch Score Distance", lwd = 2, xlab = "Successful Perturbation Number", xlim = c(1,simL), ylim = c(ymin,ymax), 
     pin = c(1,1))
title("Robinson-Foulds Distance = ", RFDist[i])


for(i in 1:148) {
  print(i)
  png(filename = paste0("Distances", i, ".png"), width = 360, height = 360)
  par(mfrow = c(1,1))
  ymax <- ceiling(max(branchDist)*10)/10
  ymin <- floor(min(branchDist))
  plot(branchDist[1:i], type = "l", ylab = "Branch Score Distance", lwd = 2, xlab = "", xlim = c(1,simL), ylim = c(ymin,ymax))
  title(paste0("Robinson-Foulds Distance = ", RFDist[i]))
  dev.off()
}


#Plot Both Trees
phylo.diff.alter(trueTree, phylos[[148]], numTraits, trueLL, logLiks[148])

for(i in 1:148) {
  print(i)
  png(paste0("phylocompare", i, ".png"), width = 720, height = 500)
  par(mfrow = c(2,1))
  phylo.diff.alter(trueTree, phylos[[i]], numTraits, trueLL, logLiks[i])
  dev.off()
}





phylo.diff.alter <- function(x, y, numTraits, trueLL, bestLL, ...){
  uniqT1 <- distinct.edges(x, y)
  uniqT2 <- distinct.edges(y, x)
  treeA.cs <- rep("black", dim(x$edge)[1])
  treeA.cs[uniqT1] = "red"
  treeB.cs <- rep("black", dim(y$edge)[1])
  treeB.cs[uniqT2] = "red"
  par(mfrow = c(1, 2))
  plot(x, edge.color = treeA.cs, ...)
  title(c("True, Data-Generating Tree", paste(numTraits, "traits")), paste("Log Likelihood =", round(trueLL, digits = 4)))
  axisPhylo()
  plot(y, edge.color = treeB.cs, ...)
  title(c("Current Best Tree", paste(numTraits, "traits")), paste("Log Likelihood =", round(bestLL, digits = 4)))
  axisPhylo()
  invisible()
}

