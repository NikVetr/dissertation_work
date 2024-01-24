library(gtools)
library(ape)

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

ntips <- 6
tr <- unroot(rtree(n = ntips))
tips <- sort(tr$tip.label)

2^(ntips-1) - ntips - 1

clades <- combinations(ntips,r=2,v=tips)
clades <- sapply(1:length(clades[,1]), function(x) paste(clades[x,], collapse = " "))
for(i in 3:(ntips-2)){
  newClades <- combinations(ntips,r=i,v=tips)
  newClades <- sapply(1:length(newClades[,1]), function(x) paste(newClades[x,], collapse = " "))
  clades <- c(clades, newClades)
}

cladeFreq <- data.frame(cladeNames = clades, cladeFreqs = rep(0, length(clades)))

clades <- paste(sort(tips), collapse = " ")
nIter <- 1e5
for(i in 1:nIter){
  if(i %% 1000 == 0){
    print(paste0(i/nIter*100,"%"))
  }
  tr <- unroot(rtree(n = ntips))
  tr$tip.label <- sample(length(tips), x = tips, replace = F)
  randClades <- prop.part.df(tr)$cladeNames[-1]
  altClades <- sapply(1:length(randClades), function(x) paste(sort(setdiff(tips, strsplit(randClades[[x]], split = " ")[[1]])), collapse = " "))
  redundantClades <- c(randClades, altClades)
  for(j in 1:length(redundantClades)){
    corrMatch <- which(redundantClades[j] == cladeFreq[,1])
    cladeFreq[corrMatch,2] <- cladeFreq[corrMatch,2] + 1
  }
}

cladeFreq[,2] <- cladeFreq[,2]/nIter

