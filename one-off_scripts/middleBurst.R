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

for(num in 10:12){
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

############################################
# make single pretty figure for dissertation
############################################

setwd("~/middleburst_450/")
num = 8
log <- read.table(file = paste0("params_", num, ".log"), header = T)[c("lambda", "delta2")]
log.branchRates <- read.table(file = paste0("params_", num, ".log"), header = T); log.branchRates <- log.branchRates[startsWith(colnames(log.branchRates), "branchRates")]
lambda <- read.table(file = paste0("modelParameters_", num, ".txt"), header = T)["lambda"][1,1]
delta2 <- read.table(file = paste0("modelParameters_", num, ".txt"), header = T)["delta2"][1,1]
burstLocation <- read.table(file = paste0("modelParameters_", num, ".txt"), header = T)["burstLocation"][1,1]
tree <- read.nexus(file = paste0("tree_", num, ".nex"))
treeHeight <- nodeheight(tree, 1)

t <- 0:(100*treeHeight)/100
r <- sapply(1:length(t), function(x) if(t[x] > burstLocation){1 + delta2*exp(-lambda*(t[x]-burstLocation))} else{1})


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

fig_label <- function(text, region="figure", pos="topleft", cex=NULL, shrinkX = 0.95, shrinkY = 0.95, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1*shrinkX, y1*shrinkY, text, cex=cex, ...)
  return(invisible(c(x,y)))
}

#actually make the figure

grDevices::cairo_pdf(filename = "/Volumes/macOS/Users/nikolai/dissertation/figures/middle_burst.pdf", width = 2250 / 72, height = 750 / 72)

par(mfrow = c(1,3), mar = c(8.5,8.5,6,2))
plot(t,r,type = "l", ylim = c(0, 1+qexp(p = 0.999, rate = 1)), lwd = 3, col = "navy", main = "", xaxt = "n", yaxt = "n", frame.plot = F, xlab = "", ylab = "")
axis(side = 1, tck = -0.02, lwd = 3, lwd.ticks = 3, cex.axis = 3, labels = rep("", 6), at = 0:5)
mtext(text = 0:5, side = 1, line = 3, at = 0:5, cex = 2)
axis(side = 2, tck = -0.02, lwd = 3, lwd.ticks = 3, cex.axis = 3, labels = rep("", 5), at = 0:4 * 2)
mtext(text = 0:4 * 2, side = 2, line = 1.75, at = 0:4 * 2, cex = 2)

box("plot", lwd = 3, lty = 1)
title(xlab = latex2exp::TeX(paste0("Time")), cex.lab = 3, line = 6.5)
title(ylab = latex2exp::TeX(paste0("Rate ($\\sigma^2$)")), cex.lab = 3, line = 4.5)
box(which = "figure", lwd = 2, lty = 2)
fig_label(text = paste0("a", ")"), cex = 5, shrinkX = 0.85, shrinkY = 0.975)
title(main = latex2exp::TeX(paste0("Rate Variation ($\\sigma^2$) Through Time")), cex.main = 3.75)
legend(x = "topleft", y = 1+qexp(p = 0.99, rate = 1), legend = c("prior", "posterior", "true value"), box.lty = 2, box.lwd = 4, cex = 3,
       lwd = c(3,3,6), col = c(wesanderson::wes_palette("FantasticFox1", 5)[4], wesanderson::wes_palette("FantasticFox1", 5)[3], "navy"))

#plot prior
for(i in 1:5e3){
  prior_delta2 <- rexp(1, rate = 1)
  prior_lambda <- 10^runif(n = 1, 0,1)
  prior_r <- sapply(1:length(t), function(x) if(t[x] > burstLocation){1 + prior_delta2*exp(-prior_lambda*(t[x]-burstLocation))} else{1})
  lines(t,prior_r,type = "l", ylim = c(0, 1+delta2), col = make.transparent(wesanderson::wes_palette("FantasticFox1", 5)[4], 0.05))
}

#plot posterior
for(i in 1:length(log$lambda)){
  post_delta2 <- log$delta2[i]
  post_lambda <- log$lambda[i]
  post_r <- sapply(1:length(t), function(x) if(t[x] > burstLocation){1 + post_delta2*exp(-post_lambda*(t[x]-burstLocation))} else{1})
  lines(t,post_r,type = "l", ylim = c(0, 1+delta2), col = make.transparent(wesanderson::wes_palette("FantasticFox1", 5)[3], 0.075))
}
lines(t,r, ylim = c(0, 1+delta2), lwd = 4, col = "black")

par(mar = c(8.5,8.5,8.5,2))
colfunc <- colorRampPalette(viridis::magma(5))
plotBranchbyTrait(tree, x = branchRates, mode = "edges", palette = colfunc, show.tip.label = F, legend = F)

segments(x0 = burstLocation, x1 = burstLocation, y0 = 5, y1 = 455, lwd = 5, col = rgb(1,0,0,0.6), lty = 2)
text("BURST", x = burstLocation, y = -5, font = 2, cex = 3)
box(which = "figure", lwd = 2, lty = 2)
fig_label(text = paste0("b", ")"), cex = 5, shrinkX = -0.65, shrinkY = 0.975)

xl <- 0.1; yb <- 70; xr <- 0.3; yt <- 300; 
rect(
  xl,
  head(seq(yb,yt,(yt-yb)/100),-1),
  xr,
  yt = tail(seq(yb,yt,(yt-yb)/100),-1),
  col=colfunc(100), lwd = 0.15
)
text(labels = c("1.0", round(seq(from = 1, to = 4.4, length.out = 12), 1)[-1]), x = xr+(xr-xl)/2, y = seq(yb,yt,length.out = 12), las=2, cex=1.75)
text(labels = latex2exp::TeX(paste0("$\\sigma^2$")), x = xl + (xr-xl)/1.25, y = yt+(yt-yb)/20, font = 2, cex = 3.5)

colfunc <- colorRampPalette(c("white", "black"))
plotBranchbyTrait(tree, x = branchRates, mode = "edges", palette = colfunc, show.tip.label = F, legend = F)
segments(x0 = burstLocation, x1 = burstLocation, y0 = 5, y1 = 455, lwd = 5, col = rgb(1,0,0,0.6), lty = 2)
text("BURST", x = burstLocation, y = -5, font = 2, cex = 3)
box(which = "figure", lwd = 2, lty = 2)
fig_label(text = paste0("c", ")"), cex = 5, shrinkX = -0.5, shrinkY = 0.975)

xl <- 0.1; yb <- 70; xr <- 0.3; yt <- 300; 
rect(
  xl,
  head(seq(yb,yt,(yt-yb)/100),-1),
  xr,
  yt = tail(seq(yb,yt,(yt-yb)/100),-1),
  col=colfunc(100), lwd = 0.15
)
text(labels = c("1.0", round(seq(from = 1, to = 4.4, length.out = 12), 1)[-1]), x = xr+(xr-xl)/2, y = seq(yb,yt,length.out = 12), las=2, cex=1.75)
text(labels = latex2exp::TeX(paste0("$\\sigma^2$")), x = xl + (xr-xl)/1.25, y = yt+(yt-yb)/20, font = 2, cex = 3.5)


dev.off()

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

