library(Matrix)
library(matrixcalc)
library(adephylo)
library(phytools)
library(phytools)
library(mvMORPH)
library(mvtnorm)
library(mnormt)
library(ggplot2)
library(phangorn)
library(distory)
library(psych)
library(abind)
library(sirt)
library(irtoys)
library(mvProbit)
library(miscTools)
library(bayesm)

trees <- rmtree(N = 1000, n = 10)

convNum2Str <- function(nums, key){
  sapply(1:length(nums), function(x) key[nums[x]])
}

prop.part.df <- function(trees, cutoff = 0.01){
  numTrees <- length(trees)
  out <- prop.part(trees)
  outList <- as.list.data.frame(out)
  pps <- attributes(outList)$number/numTrees
  props <- data.frame(pps, as.matrix(outList)); colnames(props) <- c("postProbs", "clade")
  props <- props[order(-pps),]
  props <- props[props[,1] > cutoff,]
  rownames(props) <- 1:nrow(props)
  props$cladeNames <- sapply(1:length(props[,1]), function(x) paste(sort(convNum2Str(props$clade[[x]], attributes(outList)$labels)), collapse = " "))
  props <- props[,c(1,3)]
  props
}
