library(phytools)
library(mvtnorm)

if(!dir.exists("~/middleburst")){dir.create("middleburst")}
setwd("~/middleBurst")

#################################################################
#################### Specify Single Analysis ####################
#################################################################

ntips <- 100
tree <- pbtree(n = ntips)
plot(tree); edgelabels()
edgeStartEndLength <- nodeHeights(tree)
edgeStartEndLength <- cbind(edgeStartEndLength, edgeStartEndLength[,2] - edgeStartEndLength[,1])
treeHeight <- nodeheight(tree, 1)
burstLocation <- treeHeight / 2
ntraits <- 2
rootState <- rep(0, ntraits)
baseRates <- diag(ntraits)

delta2 <- rexp(1, rate = 1)
lambda <- 10^runif(n = 1, 0,1)

delta2 <- 5
lambda <- 2

t <- 0:(100*treeHeight)/100
r <- sapply(1:length(t), function(x) if(t[x] > burstLocation){1 + delta2*exp(-lambda*(t[x]-burstLocation))} else{1})
plot(t,r,type = "l", ylim = c(0, 1+qexp(p = 0.99, rate = 1)), lwd = 3, col = "red")
for(i in 1:2e3){
  prior_delta2 <- rexp(1, rate = 1)
  prior_lambda <- 10^runif(n = 1, 0,1)
  prior_r <- sapply(1:length(t), function(x) if(t[x] > burstLocation){1 + prior_delta2*exp(-prior_lambda*(t[x]-burstLocation))} else{1})
  lines(t,prior_r,type = "l", ylim = c(0, 1+delta2), col = rgb(0,0,0,0.05))
}
lines(t,r, ylim = c(0, 1+delta2), lwd = 3, col = "red")

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

plotBranchbyTrait(tree, x = branchRates, mode = "edges", palette = "rainbow")
segments(x0 = burstLocation, x1 = burstLocation, y0 = 0, y1 = 200, lwd = 2, col = "red", lty = 2)
text("BURST", x = burstLocation, y = -5, font = 2)
plotBranchbyTrait(tree, x = branchRates, mode = "edges", palette = "gray")
segments(x0 = burstLocation, x1 = burstLocation, y0 = 0, y1 = 200, lwd = 2, col = "red", lty = 2)
text("BURST", x = burstLocation, y = -5, font = 2)

rescTree <- tree
rescTree$edge.length <- rescTree$edge.length*branchRates
treeCov <- vcv(rescTree)
traits <- as.vector(rmvnorm(n = 1, mean = rep(rootState, ntips), sigma = kronecker(treeCov, baseRates)))
traits <- t(sapply(1:ntips, function(x) traits[(((x-1)*ntraits)+1):(x*ntraits)]))
rownames(traits) <- rownames(treeCov)

#compute disparity through time graph
library(geiger)
dtt(phy = tree, data = traits[,2])

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

meansNexus(traits, "traits.nex")
ratesTSV(baseRates, "rates.tsv")
write.nexus(tree, file = "tree.nex", translate = T)

params <- data.frame(delta2 = delta2, lambda = lambda, burstLocation = burstLocation)
write.table(params, file = "modelParameters.txt")

sink("script.rev", append = F)
  
  cat("# read in the data\n")
  cat(paste0("means = readContinuousCharacterData(\"traits.nex\")\n\n"))
  
  cat("# read in the rate matrix\n")
  cat(paste0("rateMatrix = readMatrix(\"rates.tsv\")\n\n"))
  
  cat("# read in the tree\ntree <- readTrees(\"tree.nex\")[1]\n\n")
  cat("# specify age of burst\nburstAge =", treeHeight, "-", burstLocation, "\n\n")
  
  cat("# get some useful variables from the data\nnumTips = means.ntaxa()\nnumBranches = 2 * numTips - 2\nntraits <- means.nchar()\ntaxa = means.taxa()\nnumNodes = numTips * 2 - 1\nrootNode = numBranches + 1\ntreeHeight = tree.nodeAge(rootNode)\n\n")
  
  cat("#specify priors on the relaxed clock parameters\ndelta2 ~ dnExp(1)\n")
  cat("lambda_log ~ dnUniform(0,1)\nlambda := 10^lambda_log\n\n")
  
  cat("#specify moves on the relaxed clock parameters\n")
  cat("mi = 0\n")
  cat("moves[++mi] = mvSlide(delta2, delta = 0.25, tune = true, weight = 5)\n")
  cat("moves[++mi] = mvScale(delta2, lambda = 0.25, tune = true, weight = 5)\n")
  cat("moves[++mi] = mvSlide(lambda_log, delta = 0.25, tune = true, weight = 5)\n")
  cat("moves[++mi] = mvScale(lambda_log, lambda = 0.25, tune = true, weight = 5)\n\n")
  
  cat("# iterate over branches and, using deterministic operators, compute branch rates\n")
  cat("baseBranchRates <- rep(1, numBranches)\n") 
  cat("for(i in 1:numBranches){
      branchRates[i] := sqrt(abs( baseBranchRates[i] + 
      ifelse((tree.nodeAge(i)  + tree.branchLength(i)) < burstAge, 
      (delta2 / lambda * (exp(-lambda*(burstAge - (tree.nodeAge(i) + tree.branchLength(i)))) - exp(-lambda*(burstAge - tree.nodeAge(i))))) / tree.branchLength(i), 0) + 
      ifelse(tree.nodeAge(i) < burstAge & (tree.nodeAge(i) + tree.branchLength(i)) > burstAge, 
      (delta2 / lambda * (1 - exp(-lambda*(burstAge - tree.nodeAge(i))))) / tree.branchLength(i), 0) )) 
      \n}\n\n")
  
  cat("# put it all together\nmvBNdata ~ dnPhyloMultivariateBrownianREML(tree, abs(branchRates), rateMatrix, nSites=ntraits, method=\"transform\")\nmvBNdata.clamp(means)\n\n")
  
  cat("# specify monitors for the simulations model\n")
  cat("monitors[1] = mnModel(filename=\"params.log\", printgen=5, separator = TAB)\n")
  cat("monitors[2] = mnFile(tree, filename=\"trees.trees\", printgen=5, separator = TAB)\n")
  cat("monitors[3] = mnScreen(printgen=1000)\n\n")
  
  cat("# define model\nmymodel = model(tree)\n\n")
  
  cat("# define and run the mcmc, with burnin\nmymcmc = mcmc(mymodel, monitors, moves, nruns = 1)\n")
  cat("mymcmc.burnin(generations=400,tuningInterval=40)\n")
  cat("mymcmc.operatorSummary()\n")
  cat("mymcmc.run(generations=1200)\n")
  cat("mymcmc.operatorSummary()\n\n")
  
  cat("q()")

sink()

###################################################################
#################### Specify Multiple Analyses ####################
###################################################################

n_analyses <- 10
simulateData <- T

meansNexus <- function (characterMatrix, outputFilePath, traitsOnColumns = T){
  if (traitsOnColumns) {ntax = dim(characterMatrix)[1]; nchar = dim(characterMatrix)[2]}
  else {ntax = dim(characterMatrix)[2]; nchar = dim(characterMatrix)[1]}
  sink(outputFilePath, append=F)
  cat("#NEXUS\n\nBegin data;\n\tDimensions ntax=", ntax, "nchar=", nchar, ";\n\tFormat datatype=Continuous missing=?;\n\tMatrix\n\n")
  for(i in 1:length(characterMatrix[,1])) {cat(rownames(characterMatrix)[i], characterMatrix[i,], "\n")}
  cat(";\nEnd;")
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

for(num in 1:n_analyses){
  print(num)
  if(simulateData){
    ntips <- 200
    tree <- pbtree(n = ntips)
    edgeStartEndLength <- nodeHeights(tree)
    edgeStartEndLength <- cbind(edgeStartEndLength, edgeStartEndLength[,2] - edgeStartEndLength[,1])
    treeHeight <- nodeheight(tree, 1)
    burstLocation <- treeHeight / 2
    ntraits <- 8
    rootState <- rep(0, ntraits)
    baseRates <- diag(ntraits)
    
    delta2 <- rexp(1, rate = 1)
    lambda <- 10^runif(n = 1, 0,1)
    
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
    
    rescTree <- tree
    rescTree$edge.length <- rescTree$edge.length*branchRates
    treeCov <- vcv(rescTree)
    traits <- as.vector(rmvnorm(n = 1, mean = rep(rootState, ntips), sigma = kronecker(treeCov, baseRates)))
    traits <- t(sapply(1:ntips, function(x) traits[(((x-1)*ntraits)+1):(x*ntraits)]))
    rownames(traits) <- rownames(treeCov)
    
    meansNexus(traits, paste0("traits_", num, ".nex"))
    ratesTSV(baseRates, paste0("rates_", num, ".tsv"))
    write.nexus(tree, file = paste0("tree_", num, ".nex"), translate = T)
    
    params <- data.frame(delta2 = delta2, lambda = lambda, burstLocation = burstLocation)
    write.table(params, file = paste0("modelParameters_", num, ".txt"))
  }
  
  #WRITE THE .REV SCRIPT
  sink(paste0("script_", num, ".Rev"), append = F)
  
  cat("# read in the data\n")
  cat(paste0("means = readContinuousCharacterData(\"", paste0("traits_", num, ".nex"), "\")\n\n"))
  
  cat("# read in the rate matrix\n")
  cat(paste0("rateMatrix = readMatrix(\"", paste0("rates_", num, ".tsv"), "\")\n\n"))
  
  cat(paste0("# read in the tree\ntree <- readTrees(\"", paste0("tree_", num, ".nex"), "\")[1]\n\n"))
  cat("# specify age of burst\nburstAge =", treeHeight, "-", burstLocation, "\n\n")
  
  cat("# get some useful variables from the data\nnumTips = means.ntaxa()\nnumBranches = 2 * numTips - 2\nntraits <- means.nchar()\ntaxa = means.taxa()\nnumNodes = numTips * 2 - 1\nrootNode = numBranches + 1\ntreeHeight = tree.nodeAge(rootNode)\n\n")
  
  cat("#specify priors on the relaxed clock parameters\ndelta2 ~ dnExp(1)\n")
  cat("lambda_log ~ dnUniform(0,1)\nlambda := 10^lambda_log\n\n")
  
  cat("#specify moves on the relaxed clock parameters\n")
  cat("mi = 0\n")
  cat("moves[++mi] = mvSlide(delta2, delta = 0.25, tune = true, weight = 5)\n")
  cat("moves[++mi] = mvScale(delta2, lambda = 0.25, tune = true, weight = 5)\n")
  cat("moves[++mi] = mvSlide(lambda_log, delta = 0.25, tune = true, weight = 5)\n")
  cat("moves[++mi] = mvScale(lambda_log, lambda = 0.25, tune = true, weight = 5)\n\n")
  
  cat("# iterate over branches and, using deterministic operators, compute branch rates\n")
  cat("baseBranchRates <- rep(1, numBranches)\n") 
  cat("for(i in 1:numBranches){
      branchRates[i] := sqrt(abs( baseBranchRates[i] + 
      ifelse((tree.nodeAge(i)  + tree.branchLength(i)) < burstAge, 
      (delta2 / lambda * (exp(-lambda*(burstAge - (tree.nodeAge(i) + tree.branchLength(i)))) - exp(-lambda*(burstAge - tree.nodeAge(i))))) / tree.branchLength(i), 0) + 
      ifelse(tree.nodeAge(i) < burstAge & (tree.nodeAge(i) + tree.branchLength(i)) > burstAge, 
      (delta2 / lambda * (1 - exp(-lambda*(burstAge - tree.nodeAge(i))))) / tree.branchLength(i), 0) )) 
      \n}\n\n")
  
  cat("# put it all together\nmvBNdata ~ dnPhyloMultivariateBrownianREML(tree, abs(branchRates), rateMatrix, nSites=ntraits, method=\"transform\")\nmvBNdata.clamp(means)\n\n")
  
  cat("# specify monitors for the simulations model\n")
  cat(paste0("monitors[1] = mnModel(filename=\"", paste0("params_", num, ".log"), "\", printgen=5, separator = TAB)\n"))
  cat(paste0("monitors[2] = mnFile(tree, filename=\"", paste0("trees_", num, ".trees"), "\", printgen=5, separator = TAB)\n"))
  cat("monitors[3] = mnScreen(printgen=1000)\n\n")
  
  cat("# define model\nmymodel = model(tree)\n\n")
  
  cat("# define and run the mcmc, with burnin\nmymcmc = mcmc(mymodel, monitors, moves, nruns = 1)\n")
  cat("mymcmc.burnin(generations=1000,tuningInterval=100)\n")
  cat("mymcmc.operatorSummary()\n")
  cat("mymcmc.run(generations=4000)\n")
  cat("mymcmc.operatorSummary()\n\n")
  
  cat("q()")
  
  sink()
  
  #specify job file for ppss
  jobName <- "job.txt"
  if(num == 1){
    scriptPath <- paste0("script_", num, ".Rev")
    sink(jobName, append = F)
    cat(scriptPath, "\n", sep = "")
    sink()
  } else {
    scriptPath <- paste0("script_", num, ".Rev")
    sink(jobName, append = T)
    cat(scriptPath, "\n", sep = "")
    sink()
  }
  
}


#############################################################
#################### Analyze MCMC Output ####################
#############################################################

#read in single MCMC analysis
# log <- read.table(file = paste0("params.log"), header = T)[c("lambda", "delta2")]
# log.branchRates <- read.table(file = paste0("params.log"), header = T); log.branchRates <- log.branchRates[startsWith(colnames(log.branchRates), "branchRates")]
# lambda <- read.table(file = paste0("modelParameters.txt"), header = T)["lambda"][1,1]
# delta2 <- read.table(file = paste0("modelParameters.txt"), header = T)["delta2"][1,1]
# burstLocation <- read.table(file = paste0("modelParameters.txt"), header = T)["burstLocation"][1,1]
# tree <- read.nexus(file = paste0("tree.nex"))
# treeHeight <- nodeheight(tree, 1)

#read in indexed multiple analysis

printToFile <- T
if(printToFile){
  if(!dir.exists("~/middleburst/plots")){dir.create("~/middleburst/plots")}
}
setwd("~/middleBurst")

for(num in 1:n_analyses){
  # num <- 2
  print(num)
  log <- read.table(file = paste0("params_", num, ".log"), header = T)[c("lambda", "delta2")]
  log.branchRates <- read.table(file = paste0("params_", num, ".log"), header = T); log.branchRates <- log.branchRates[startsWith(colnames(log.branchRates), "branchRates")]
  lambda <- read.table(file = paste0("modelParameters_", num, ".txt"), header = T)["lambda"][1,1]
  delta2 <- read.table(file = paste0("modelParameters_", num, ".txt"), header = T)["delta2"][1,1]
  burstLocation <- read.table(file = paste0("modelParameters_", num, ".txt"), header = T)["burstLocation"][1,1]
  tree <- read.nexus(file = paste0("tree_", num, ".nex"))
  treeHeight <- nodeheight(tree, 1)
  
  t <- 0:(100*treeHeight)/100
  r <- sapply(1:length(t), function(x) if(t[x] > burstLocation){1 + delta2*exp(-lambda*(t[x]-burstLocation))} else{1})
  
  if(printToFile){png(filename = paste0("plots/plot_", num, ".png"), width = 1600, height = 400)}
  
  par(mfrow = c(1,3))
  plot(t,r,type = "l", ylim = c(0, 1+qexp(p = 0.99, rate = 1)), lwd = 3, col = "black")
  
  #plot prior
  for(i in 1:2e3){
    prior_delta2 <- rexp(1, rate = 1)
    prior_lambda <- 10^runif(n = 1, 0,1)
    prior_r <- sapply(1:length(t), function(x) if(t[x] > burstLocation){1 + prior_delta2*exp(-prior_lambda*(t[x]-burstLocation))} else{1})
    lines(t,prior_r,type = "l", ylim = c(0, 1+delta2), col = rgb(1,0,0,0.05))
  }

  #plot posterior
  for(i in 1:length(log$lambda)){
    post_delta2 <- log$delta2[i]
    post_lambda <- log$lambda[i]
    post_r <- sapply(1:length(t), function(x) if(t[x] > burstLocation){1 + post_delta2*exp(-post_lambda*(t[x]-burstLocation))} else{1})
    lines(t,post_r,type = "l", ylim = c(0, 1+delta2), col = rgb(0,0,1,0.075))
  }
  lines(t,r, ylim = c(0, 1+delta2), lwd = 3, col = "black")
  
  legend(x = 0, y = 1+qexp(p = 0.99, rate = 1), legend = c("prior", "posterior", "true value"), lwd = c(2,2,4), col = c("red", "blue", "black"))
  
  #visualize rate heterogeneity
  edgeStartEndLength <- nodeHeights(tree)
  edgeStartEndLength <- cbind(edgeStartEndLength, edgeStartEndLength[,2] - edgeStartEndLength[,1])
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
  
  plotBranchbyTrait(tree, x = branchRates, mode = "edges", palette = "rainbow")
  segments(x0 = burstLocation, x1 = burstLocation, y0 = 0, y1 = 200, lwd = 2, col = "red", lty = 2)
  text("BURST", x = burstLocation, y = -5, font = 2)
  plotBranchbyTrait(tree, x = branchRates, mode = "edges", palette = "gray")
  segments(x0 = burstLocation, x1 = burstLocation, y0 = 0, y1 = 200, lwd = 2, col = "red", lty = 2)
  text("BURST", x = burstLocation, y = -5, font = 2)
  
  if(printToFile){dev.off()}
  
}

############################################################################################
#################### Specify Multiple Analyses with Bernoulli Indicator ####################
############################################################################################

if(!dir.exists("~/middleburst_bernoulli")){dir.create("~/middleburst_bernoulli")}
setwd("~/middleburst_bernoulli/")

n_analyses <- 12
simulateData <- T

meansNexus <- function (characterMatrix, outputFilePath, traitsOnColumns = T){
  if (traitsOnColumns) {ntax = dim(characterMatrix)[1]; nchar = dim(characterMatrix)[2]}
  else {ntax = dim(characterMatrix)[2]; nchar = dim(characterMatrix)[1]}
  sink(outputFilePath, append=F)
  cat("#NEXUS\n\nBegin data;\n\tDimensions ntax=", ntax, "nchar=", nchar, ";\n\tFormat datatype=Continuous missing=?;\n\tMatrix\n\n")
  for(i in 1:length(characterMatrix[,1])) {cat(rownames(characterMatrix)[i], characterMatrix[i,], "\n")}
  cat(";\nEnd;")
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

for(num in 1:n_analyses){
  print(num)
  if(simulateData){
    ntips <- 200
    tree <- pbtree(n = ntips)
    edgeStartEndLength <- nodeHeights(tree)
    edgeStartEndLength <- cbind(edgeStartEndLength, edgeStartEndLength[,2] - edgeStartEndLength[,1])
    treeHeight <- nodeheight(tree, 1)
    burstLocation <- treeHeight / 2
    ntraits <- 8
    rootState <- rep(0, ntraits)
    baseRates <- diag(ntraits)
    
    burstPA <- rbinom(1,1,0.5)
    delta2 <- rexp(1, rate = 1)*burstPA
    bernBurst <- rbinom(1,1,0.5)
    lambda <- 10^runif(n = 1, 0,1)
    
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
    
    rescTree <- tree
    rescTree$edge.length <- rescTree$edge.length*branchRates
    treeCov <- vcv(rescTree)
    traits <- as.vector(rmvnorm(n = 1, mean = rep(rootState, ntips), sigma = kronecker(treeCov, baseRates)))
    traits <- t(sapply(1:ntips, function(x) traits[(((x-1)*ntraits)+1):(x*ntraits)]))
    rownames(traits) <- rownames(treeCov)
    
    meansNexus(traits, paste0("traits_", num, ".nex"))
    ratesTSV(baseRates, paste0("rates_", num, ".tsv"))
    write.nexus(tree, file = paste0("tree_", num, ".nex"), translate = T)
    
    params <- data.frame(delta2 = delta2, lambda = lambda, burstLocation = burstLocation)
    write.table(params, file = paste0("modelParameters_", num, ".txt"))
  }
  
  #WRITE THE .REV SCRIPT
  sink(paste0("script_", num, ".Rev"), append = F)
  
  cat("# read in the data\n")
  cat(paste0("means = readContinuousCharacterData(\"", paste0("traits_", num, ".nex"), "\")\n\n"))
  
  cat("# read in the rate matrix\n")
  cat(paste0("rateMatrix = readMatrix(\"", paste0("rates_", num, ".tsv"), "\")\n\n"))
  
  cat(paste0("# read in the tree\ntree <- readTrees(\"", paste0("tree_", num, ".nex"), "\")[1]\n\n"))
  cat("# specify age of burst\nburstAge =", treeHeight, "-", burstLocation, "\n\n")
  
  cat("# get some useful variables from the data\nnumTips = means.ntaxa()\nnumBranches = 2 * numTips - 2\nntraits <- means.nchar()\ntaxa = means.taxa()\nnumNodes = numTips * 2 - 1\nrootNode = numBranches + 1\ntreeHeight = tree.nodeAge(rootNode)\n\n")
  
  cat("#specify priors on the relaxed clock parameters\ndelta2 ~ dnExp(1)\n")
  cat("burstPA ~ dnBernoulli(0.5)\n")
  cat("lambda_log ~ dnUniform(0,1)\nlambda := 10^lambda_log\n\n")
  
  cat("#specify moves on the relaxed clock parameters\n")
  cat("mi = 0\n")
  cat("moves[++mi] = mvBinarySwitch(burstPA, weight = 10)\n")
  cat("moves[++mi] = mvSlide(delta2, delta = 0.25, tune = true, weight = 5)\n")
  cat("moves[++mi] = mvScale(delta2, lambda = 0.25, tune = true, weight = 5)\n")
  cat("moves[++mi] = mvSlide(lambda_log, delta = 0.25, tune = true, weight = 5)\n")
  cat("moves[++mi] = mvScale(lambda_log, lambda = 0.25, tune = true, weight = 5)\n\n")
  
  cat("# iterate over branches and, using deterministic operators, compute branch rates\n")
  cat("baseBranchRates <- rep(1, numBranches)\n") 
  cat("for(i in 1:numBranches){
      branchRates[i] := sqrt(abs( baseBranchRates[i] + burstPA * (
      ifelse((tree.nodeAge(i)  + tree.branchLength(i)) < burstAge, 
      (delta2 / lambda * (exp(-lambda*(burstAge - (tree.nodeAge(i) + tree.branchLength(i)))) - exp(-lambda*(burstAge - tree.nodeAge(i))))) / tree.branchLength(i), 0) + 
      ifelse(tree.nodeAge(i) < burstAge & (tree.nodeAge(i) + tree.branchLength(i)) > burstAge, 
      (delta2 / lambda * (1 - exp(-lambda*(burstAge - tree.nodeAge(i))))) / tree.branchLength(i), 0)) )) 
      \n}\n\n")
  
  cat("# put it all together\nmvBNdata ~ dnPhyloMultivariateBrownianREML(tree, abs(branchRates), rateMatrix, nSites=ntraits, method=\"transform\")\nmvBNdata.clamp(means)\n\n")
  
  cat("# specify monitors for the simulations model\n")
  cat(paste0("monitors[1] = mnModel(filename=\"", paste0("params_", num, ".log"), "\", printgen=5, separator = TAB)\n"))
  cat(paste0("monitors[2] = mnFile(tree, filename=\"", paste0("trees_", num, ".trees"), "\", printgen=5, separator = TAB)\n"))
  cat("monitors[3] = mnScreen(printgen=1000)\n\n")
  
  cat("# define model\nmymodel = model(tree)\n\n")
  
  cat("# define and run the mcmc, with burnin\nmymcmc = mcmc(mymodel, monitors, moves, nruns = 1)\n")
  cat("mymcmc.burnin(generations=1000,tuningInterval=100)\n")
  cat("mymcmc.operatorSummary()\n")
  cat("mymcmc.run(generations=4000)\n")
  cat("mymcmc.operatorSummary()\n\n")
  
  cat("q()")
  
  sink()
  
  #specify job file for ppss
  jobName <- "job.txt"
  if(num == 1){
    scriptPath <- paste0("script_", num, ".Rev")
    sink(jobName, append = F)
    cat(scriptPath, "\n", sep = "")
    sink()
  } else {
    scriptPath <- paste0("script_", num, ".Rev")
    sink(jobName, append = T)
    cat(scriptPath, "\n", sep = "")
    sink()
  }
  
}

