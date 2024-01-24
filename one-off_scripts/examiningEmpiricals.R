setwd("/Volumes/macOS/Users/nikolai/output/")
# 
# library(Matrix)
# library(matrixcalc)
library(adephylo)
library(phytools)
# library(mvMORPH)
# library(mvtnorm)
# library(mnormt)
# library(ggplot2)
library(phangorn)
# library(distory)
# library(psych)
# library(abind)
# library(sirt)
# library(irtoys)
# library(mvProbit)
# library(miscTools)
# library(bayesm)
# library(phyloTop)
# library(parallel)
# library(doParallel)
# library(foreach)
# library(rwty)
# library(latticeExtra)
# library(bonsai)
library(smacof)
library(MASS)

convNum2Str <- function(nums, key){
  sapply(1:length(nums), function(x) key[nums[x]])
}


prop.part.df <- function(trees, cutoff = 0.01, bs = T){
  if(class(trees) == "multiPhylo"){
    if(trees[[1]]$Nnode == (length(trees[[1]]$tip.label) - 1)){
      trees <- unroot(trees) #unroot rooted trees
    }
    tipLabs <- trees[[1]]$tip.label
    numTrees <- length(trees)
    if(bs) {
      out <- as.prop.part(bitsplits(trees))
    } else {
      out <- prop.part(trees)
    }
    outList <- as.list.data.frame(out)
    pps <- attributes(outList)$number/numTrees
    props <- data.frame(pps, as.matrix(outList)); colnames(props) <- c("postProbs", "clade")
    props <- props[order(-pps),]
    props <- props[props[,1] > cutoff,]
    rownames(props) <- 1:nrow(props)
    props$cladeNames <- sapply(1:length(props[,1]), function(x) sort(convNum2Str(props$clade[[x]], attributes(outList)$labels)))
    props <- props[,c(1,3)]
    if(!bs) {
      props <- props[-1,]
    }
    allClades <- c(props$cladeNames, lapply(1:length(props$postProbs), function(x) sort(setdiff(tipLabs, props$cladeNames[[x]]))))
    if(any(duplicated(allClades))){
      # print("duplicate clades found")
      dupClades <- allClades[duplicated(allClades)]
      dupCladesComp <- lapply(1:length(dupClades), function(x) sort(setdiff(tipLabs, dupClades[[x]])))
      matchedDupes <- cbind(1:length(dupClades), sapply(1:length(dupClades), function(x) which(sapply(1:length(dupCladesComp), function(y) setequal(dupClades[[x]], dupCladesComp[[y]])))))
      dupes <- matrix(data = 0, nrow = length(matchedDupes[,1])/2, ncol = 2)
      for(i in 1:length(matchedDupes[,1])){
        pair <- sort(matchedDupes[i,])
        if(!any(sapply(1:length(dupes[,1]), function(x) pair == dupes[x,]))){
          dupes[min(which(apply(dupes, 1, sum) == 0)),] <- pair
        }
      }
      for(i in 1:length(dupes[,1])){
        clade1 <- dupClades[dupes[i,1]][[1]]
        clade2 <- dupClades[dupes[i,2]][[1]]
        propsInd1 <- which(sapply(1:length(props[,1]), function(x) setequal(clade1, props$cladeNames[[x]])))
        propsInd2 <- which(sapply(1:length(props[,1]), function(x) setequal(clade2, props$cladeNames[[x]])))
        props[propsInd1,1] <- props[propsInd1,1] + props[propsInd2,1]
        props <- props[-propsInd2,]
      }
      props <- props[order(props[,1], decreasing = T),]
    }
    props
  } else if (class(trees) == "phylo"){
    if(trees$Nnode == (length(trees$tip.label) - 1)){
      trees <- unroot(trees) #unroot rooted trees
    }
    tipLabs <- trees$tip.label
    numTrees <- 1
    if(bs) {
      out <- as.prop.part(bitsplits(trees))
    } else {
      out <- prop.part(trees)
    }
    outList <- as.list.data.frame(out)
    pps <- attributes(outList)$number/numTrees
    props <- data.frame(pps, as.matrix(outList)); colnames(props) <- c("postProbs", "clade")
    props <- props[order(-pps),]
    props <- props[props[,1] > cutoff,]
    rownames(props) <- 1:nrow(props)
    props$cladeNames <- sapply(1:length(props[,1]), function(x) sort(convNum2Str(props$clade[[x]], attributes(outList)$labels)))
    props <- props[,c(1,3)]
    if(!bs) {
      props <- props[-1,]
    }
    props
  }
}

compareTrees <- function(biparts1, biparts2, tipLabs){
  matchProbs <- matrix(0, nrow = sum(length(biparts1[,1]), length(biparts2[,1])), ncol = 2)
  counter <- 1
  biparts2$notSeen <- 1
  for(clade in 1:length(biparts1[,2])){
    cladeName <- biparts1[,2][[clade]]
    isThere <- sapply(1:length(biparts2[,1]), function(x) identical(cladeName, biparts2[x,2][[1]]))
    altCladeName <- sort(setdiff(tipLabs, cladeName))
    orIsThere <- sapply(1:length(biparts2[,1]), function(x) identical(altCladeName, biparts2[x,2][[1]]))
    if(any(isThere)){
      biparts2$notSeen[isThere] <- 0
      matchProbs[counter,] <- c(biparts1[,1][clade], biparts2[,1][isThere])
    } else if (any(orIsThere)) {
      biparts2$notSeen[orIsThere] <- 0
      matchProbs[counter,] <- c(biparts1[,1][clade], biparts2[,1][orIsThere])
    } else {
      matchProbs[counter,] <- c(biparts1[,1][[clade]], 0)
    }
    counter <- counter + 1
  }
  if(sum(biparts2$notSeen) > 0){
    matchProbs[counter:(counter+sum(biparts2$notSeen)-1),] <- cbind(0, biparts2[biparts2$notSeen == 1,1])
  }
  matchProbs <- matchProbs[rowSums(matchProbs == c(0,0)) != 2,]
  matchProbs
}

draw.contour<-function(a,alpha=0.95,plot.dens=FALSE, line.width=2, line.type=1, limits=NULL, density.res=300,spline.smooth=-1,...){
  ##a is a list or matrix of x and y coordinates (e.g., a=list("x"=rnorm(100),"y"=rnorm(100)))
  ## if a is a list or dataframe, the components must be labeled "x" and "y"
  ## if a is a matrix, the first column is assumed to be x, the second y
  ##alpha is the contour level desired
  ##if plot.dens==TRUE, then the joint density of x and y are plotted,
  ##   otherwise the contour is added to the current plot.
  ##density.res controls the resolution of the density plot
  
  ##A key assumption of this function is that very little probability mass lies outside the limits of
  ## the x and y values in "a". This is likely reasonable if the number of observations in a is large.
  
  require(MASS)
  require(ks)
  if(length(line.width)!=length(alpha)){
    line.width <- rep(line.width[1],length(alpha))
  }
  
  if(length(line.type)!=length(alpha)){
    line.type <- rep(line.type[1],length(alpha))
  }
  
  if(is.matrix(a)){
    a=list("x"=a[,1],"y"=a[,2])
  }
  ##generate approximate density values
  if(is.null(limits)){
    limits=c(range(a$x),range(a$y))
  }
  f1<-kde2d(a$x,a$y,n=density.res,lims=limits)
  
  ##plot empirical density
  if(plot.dens) image(f1,...)
  
  if(is.null(dev.list())){
    ##ensure that there is a window in which to draw the contour
    plot(a,type="n",xlim=limits[1:2],ylim=limits[3:4],...)
  }
  
  ##estimate critical contour value
  ## assume that density outside of plot is very small
  
  zdens <- rev(sort(f1$z))
  Czdens <- cumsum(zdens)
  Czdens <- (Czdens/Czdens[length(zdens)])
  for(cont.level in 1:length(alpha)){
    ##This loop allows for multiple contour levels
    crit.val <- zdens[max(which(Czdens<=alpha[cont.level]))]
    
    ##determine coordinates of critical contour
    b.full=contourLines(f1,levels=crit.val)
    for(c in 1:length(b.full)){
      ##This loop is used in case the density is multimodal or if the desired contour
      ##  extends outside the plotting region
      b=list("x"=as.vector(unlist(b.full[[c]][2])),"y"=as.vector(unlist(b.full[[c]][3])))
      
      ##plot desired contour
      line.dat<-xspline(b,shape=spline.smooth,open=TRUE,draw=FALSE)
      lines(line.dat,lty=line.type[cont.level],lwd=line.width[cont.level])
    }
  }
}

##########################################
#plotting fancier versions of the figures
##########################################

#read in data
modernHumans_M_log <- read.tree("howells_log.trees")
modernHumans_M_nonlog <- read.tree("howells_nonlog.trees")
modernHumans_F_log <- read.tree("howells_log_F.trees")
modernHumans_F_nonlog <- read.tree("howells_nonlog_F.trees")
# modernHumans_MF_log <- read.tree("howells_log_MF.trees")
# modernHumans_MF_nonlog <- read.tree("howells_nonlog_MF.trees")
# modernHumans_MwF_log <- read.tree("howells_log_MwF.trees")
# modernHumans_MwF_nonlog <- read.tree("howells_nonlog_MwF.trees")

mccTree_log_M <- maxCladeCred(modernHumans_M_log)
mccTree_nonlog_M <- maxCladeCred(modernHumans_M_nonlog)
consensus_log_M <- consensus(modernHumans_M_log, p = 0.5)
consensus_nonlog_M <- consensus(modernHumans_M_nonlog, p = 0.5)

mccTree_log_F <- maxCladeCred(modernHumans_F_log)
mccTree_nonlog_F <- maxCladeCred(modernHumans_F_nonlog)
consensus_log_F <- consensus(modernHumans_F_log, p = 0.5)
consensus_nonlog_F <- consensus(modernHumans_F_nonlog, p = 0.5)

modernHumans_log <- modernHumans_F_log
modernHumans_log <- modernHumans_M_log
modernHumans_nonlog <- modernHumans_F_nonlog
modernHumans_nonlog <- modernHumans_M_nonlog

#get stairstep CAR curve
biparts <- prop.part.df(modernHumans_M_nonlog)$postProbs
stairs <- sapply(500:1000, function(x) sum(x/1000 < biparts)/27)
plot((500:1000)/1000, stairs, type = "l", ylim = c(0,1))

#get stairstep CAR curve
biparts <- prop.part.df(modernHumans_F_nonlog)$postProbs
stairs <- sapply(500:1000, function(x) sum(x/1000 < biparts)/27)
lines((500:1000)/1000, stairs)

#MCC TREE
tree <- mccTree_nonlog_M
# tree <- root(tree, "BUSHMAN", resolve.root = T, edgelabel = T, interactive = T)
tree <- reroot(tree = tree, node.number = 15, position = 0.1) #for male tree
# tree <- reroot(tree = tree, node.number = 14, position = 0.1) #for female tree
trees <- root(modernHumans_nonlog, "BUSHMAN", resolve.root = T, edgelabel = T)
pp <- prop.clades(tree, trees, rooted = T)/length(trees)

par(mfrow = c(1,1))
# plot(tree)
# drawSupportOnEdges(round(pp, 2), col = 2)

colfunc<-colorRampPalette(c(adjustcolor( "deepskyblue", alpha.f = 0.2), adjustcolor( "darkorange", alpha.f = 0.2)))
colfunc<-colorRampPalette(c(adjustcolor( "royalblue", alpha.f = 0.2), adjustcolor( "red", alpha.f = 0.2)))
colfunc(100)

# png(filename = "howellsMCC.png", width = 1500, height = 800)
plot(tree)
title("maximum clade credibility tree, howells' data, females only, arithmetic mvBM")
edgelabels(round(pp, 2), sapply(1:tree$Nnode + length(tree$tip.label), 
                                match, tree$edge[, 2]), bg = colfunc(100)[round(pp, 2)*100], frame = "r", adj = c(0,1.5))
edgelabels(round(pp, 2), sapply(1:tree$Nnode + length(tree$tip.label), 
                                match, tree$edge[, 2]), bg = rgb(255, 255, 255, max = 255, alpha = 150, names = "white25")[[1]], frame = "r", adj = c(0,1.5))
edgelabels(round(pp, 2), sapply(1:tree$Nnode + length(tree$tip.label), 
                                match, tree$edge[, 2]), bg = colfunc(100)[round(pp, 2)*100], frame = "none", adj = c(0,1.5))
axisPhylo()
colfunc <- colorRampPalette(c("red","royalblue"))
xl <- 0.05; yb <- 5; xr <- 0.1; yt <- 20
rect(
  xl,
  head(seq(yb,yt,(yt-yb)/100),-1),
  xr,
  tail(seq(yb,yt,(yt-yb)/100),-1),
  col=colfunc(100)
)
rect(
  xl,
  head(seq(yb,yt,(yt-yb)/100),-1),
  xr,
  tail(seq(yb,yt,(yt-yb)/100),-1),
  col=rgb(255, 255, 255, max = 255, alpha = 150, names = "white25")[[1]]
)
text(labels = seq(from = 0, to = 1, length.out = 11), y = tail(seq(yb,yt,(yt-yb)/11),-1)-0.75, x = xl-(xr-xl)/2, las=2, cex=0.7)

# dev.off()

#CONSENSUS TREE
tree <- consensus(modernHumans_nonlog, p = 0.5)
tree <- root(tree, "BUSHMAN", resolve.root = T, edgelabel = T, interactive = F)
# tree <- reroot(tree = tree, node.number = 15, position = 0.1)
trees <- root(modernHumans_nonlog, "BUSHMAN", resolve.root = T, edgelabel = T)
pp <- prop.clades(tree, trees, rooted = T)/length(trees)

# png(filename = "howellsMJC.png", width = 1500, height = 800)
par(mfrow = c(1,1)); par(mar=c(2,2,2,2))
plot(tree)
title("consensus tree (p=0.5), howells' data, females only, arithmetic mvBM")
edgelabels(round(pp, 2), sapply(1:tree$Nnode + length(tree$tip.label), 
                                match, tree$edge[, 2]), bg = colfunc(100)[round(pp, 2)*100], frame = "r", adj = c(0,1.5))
edgelabels(round(pp, 2), sapply(1:tree$Nnode + length(tree$tip.label), 
                                match, tree$edge[, 2]), bg = rgb(255, 255, 255, max = 255, alpha = 150, names = "white25")[[1]], frame = "r", adj = c(0,1.5))
edgelabels(round(pp, 2), sapply(1:tree$Nnode + length(tree$tip.label), 
                                match, tree$edge[, 2]), bg = colfunc(100)[round(pp, 2)*100], frame = "none", adj = c(0,1.5))
axisPhylo()
colfunc <- colorRampPalette(c("red","royalblue"))
xl <- 2; yb <- 2; xr <-3; yt <- 17 #bottom left, top right etc. coordinates
rect(
  xl,
  head(seq(yb,yt,(yt-yb)/100),-1),
  xr,
  tail(seq(yb,yt,(yt-yb)/100),-1),
  col=colfunc(100)
)
rect(
  xl,
  head(seq(yb,yt,(yt-yb)/100),-1),
  xr,
  tail(seq(yb,yt,(yt-yb)/100),-1),
  col=rgb(255, 255, 255, max = 255, alpha = 150, names = "white25")[[1]]
)
text(labels = seq(from = 0, to = 1, length.out = 11), y = tail(seq(yb,yt,(yt-yb)/11),-1)-0.75, x = xl-(xr-xl)/2, las=2, cex=0.7)

# dev.off()

#######################################
## cophylo plots for manuscript figure
#######################################

#compute mcc trees
treeM <- reroot(tree = mccTree_nonlog_M, node.number = 15, position = 0.1) #for male tree
modernHumans_M_nonlog_R <- root(phy = modernHumans_M_nonlog, "BUSHMAN", resolve.root = T, edgelabel = T)
ppM <- prop.clades(treeM, modernHumans_M_nonlog_R, rooted = T)/length(modernHumans_M_nonlog_R)
treeM$node.label <- round(ppM,2)
treeM$node.label <- sapply(1:length(round(ppM,2)), function(x) if(round(ppM,2)[x] == 1){"1"}else{substr(round(ppM,2)[x], start = 2, stop = 4)})
plot(treeM)
edgelabels(round(ppM, 2), sapply(1:treeM$Nnode + length(treeM$tip.label), match, treeM$edge[, 2]))

treeF <- reroot(tree = mccTree_nonlog_F, node.number = 14, position = 0.1) #for female tree
modernHumans_F_nonlog_R <- root(phy = modernHumans_F_nonlog, "BUSHMAN", resolve.root = T, edgelabel = T)
ppF <- prop.clades(treeF, modernHumans_F_nonlog_R, rooted = T)/length(modernHumans_F_nonlog_R)
treeF$node.label <- round(ppF,2)
treeF$node.label <- sapply(1:length(round(ppF,2)), function(x) if(round(ppF,2)[x] == 1){"1"}else{substr(round(ppF,2)[x], start = 2, stop = 4)})
plot(treeF)
edgelabels(round(ppF, 2), sapply(1:treesF$Nnode + length(treesF$tip.label), match, treesF$edge[, 2]))

#change names of tips to preferred names
p <- c(
  "NORSE", "Northern Europe", "Norway", "Norway - Norse",
  "ZALAVAR", "Central Europe", "Hungary", "Hungary - Zalavar",
  "BERG", "Central Europe", "Austria", "Austria - Berg",
  "EGYPT", "North Africa", "Egypt", "Egypt - Egypt",
  "TEITA", "East Africa", "Kenya", "Kenya - Teita",
  "DOGON", "West Africa", "Mali", "Mali - Dogon",
  "ZULU", "South Africa", "Zulu", "South Africa - Zulu",
  "BUSHMAN", "South Africa", "San", "South Africa - San",
  "AUSTRALI", "Australia", "Australia", "Australia - Australia",
  "TASMANIA", "Australia", "Tasmania", "Australia - Tasmania",
  "TOLAI", "Melanesia", "Papua New Guinea", "Papua New Guinea - Tolai",
  "MOKAPU", "Polynesia", "Hawaii", "Hawaii - Mokapu",
  "BURIAT", "North Asia", "Siberia", "Siberia - Buriat",
  "ESKIMO", "Greenland", "Greenland", "Greenland - Eskimo",
  "PERU", "South America", "Peru", "Peru - Peru",
  "ANDAMAN", "Indonesia", "Andaman Islands", "Indonesia - Andaman",
  "EASTERI", "Polynesia", "Easter Island", "Polynesia - Easter Island",
  "ARIKARA", "North America", "USA", "USA - Arikara",
  "AINU", "East Asia", "North Japan", "Japan - Ainu",
  "NJAPAN", "East Asia", "North Japan", "Japan - North Japan",
  "SJAPAN", "East Asia", "South Japan", "Japan - South Japan",
  "HAINAN", "East Asia", "China", "China - Hainan",
  "ANYANG", "East Asia", "China", "China - Anyang",
  "ATAYAL", "East Asia", "Taiwan", "Taiwan - Atayal",
  "PHILLIPI", "Southeast Asia", "Philippines", "Phillipines - Phillipines",
  "GUAM", "Micronesia", "Guam", "Guam - Guam",
  "MORIORI", "Polynesia", "Chatham Islands", "Chatham Islands - Moriori",
  "SMAORI", "New Zealand", "New Zealand", "New Zealand - South Maori",
  "NMAORI",  "New Zealand", "New Zealand", "New Zealand - North Maori",
  "SANTACR", "North America", "USA", "USA - Santa Cruz"
)
p <- cbind(p[1:30*4-3], p[1:30*4-2], p[1:30*4-1], p[1:30*4])
treeM$tip.label <- sapply(1:length(treeM$tip.label), function(x) p[treeM$tip.label[x] == p[,1],4])
treeF$tip.label <- sapply(1:length(treeF$tip.label), function(x) p[treeF$tip.label[x] == p[,1],4])

treesMF <- cophylo(treeM, treeF, rotate = T)

png(filename = "~/mvbm_manuscript/figures/howellsTreesMF.png", width = 1000, height = 600)

plot(treesMF, scale.bar=c(.5,.5), mar=c(.1,.1,.5,.1))
title("\n\nMales                                                                                                                   Females")

#plot circles of particular shade of grey on corresponding node
colfunc<-colorRampPalette(c("white", "black"))
nodelabels.cophylo(pch=21, frame="none", bg=colfunc(100)[as.numeric(treesMF$trees[[1]]$node.label)*100], cex=1.5, which = "left")
nodelabels.cophylo(pch=21, frame="none", bg=colfunc(100)[as.numeric(treesMF$trees[[2]]$node.label)*100], cex=1.5, which= "right")
#vertical pp
# xl <- -0.5; yb <- 0; xr <- -.47; yt <- .6; 
# rect(
#   xl,
#   head(seq(yb,yt,(yt-yb)/100),-1),
#   xr,
#   tail(seq(yb,yt,(yt-yb)/100),-1),
#   col=colfunc(100)
# )
# text(labels = seq(from = 0, to = 1, length.out = 11), y = seq(yb,yt,length.out = 11), x = xl-(xr-xl)/3, las=2, cex=0.7)
# text(labels = "posterior\nprobability", x = (xl+xr)/2, y = yt+0.04)
#horizontal pp
xl <- -0.3; yb <- -.1; xr <- 0.3; yt <- -0.05; 
rect(
  head(seq(xl,xr,(xr-xl)/100),-1),
  yb,
  tail(seq(xl,xr,(xr-xl)/100),-1),
  yt,
  col=colfunc(100)
)
text(labels = seq(from = 0, to = 1, length.out = 11), y = yb-(yt-yb)/3, x = seq(xl,xr,length.out = 11), las=2, cex=1)
text(labels = "posterior probability", x = (xl+xr)/2, y = yt+(yt-yb)/3)

dev.off()

##########################################
## end cophylo plots for manuscript figure
##########################################

#plot pps on subtending edge
edgelabels.cophylo(treesMF$trees[[1]]$node.label[2:treesMF$trees[[1]]$Nnode],
                   edge=sapply(2:treesMF$trees[[1]]$Nnode+Ntip(treesMF$trees[[1]]),
                               function(n,e) which(e==n),e=treesMF$trees[[1]]$edge[,2]),
                   cex=0.7,
                   frame="none",adj=c(0.5,1.2), which = "left")
edgelabels.cophylo(treesMF$trees[[2]]$node.label[2:treesMF$trees[[2]]$Nnode],
                   edge=sapply(2:treesMF$trees[[2]]$Nnode+Ntip(treesMF$trees[[2]]),
                               function(n,e) which(e==n),e=treesMF$trees[[2]]$edge[,2]),
                   cex=0.7,
                   frame="none",adj=c(0.5,1.2),which="right")

plot(treesMF)
title("maximum clade credibility tree, howells' data, females only, arithmetic mvBM")
edgelabels.cophylo(round(treesMF$trees[[1]]$edge.length,2),cex=0.7,
                   frame="none",adj=c(0.5,1.2),which="left")
edgelabels.cophylo(round(treesMF$trees[[2]]$edge.length,2),cex=0.7,
                   frame="none",adj=c(0.5,1.2),which="right")





edgelabels.cophylo(treesMF$trees[[1]]$node.label[2:treesMF$trees[[1]]$Nnode],
                   edge=sapply(2:treesMF$trees[[1]]$Nnode+Ntip(treesMF$trees[[1]]),
                               function(n,e) which(e==n),e=treesMF$trees[[1]]$edge[,2]),
                   frame="none",adj=c(0.5,1))



plot(treesMF)
edgelabels.cophylo(treesMF$trees[[1]]$node.label[2:treesMF$trees[[1]]$Nnode],
                   edge=sapply(2:treesMF$trees[[1]]$Nnode+Ntip(treesMF$trees[[1]]),
                               function(n,e) which(e==n),e=treesMF$trees[[1]]$edge[,2]),
                   frame="none",adj=c(0.5,1))
edgelabels.cophylo(treesMF$trees[[2]]$node.label[2:treesMF$trees[[2]]$Nnode],
                   edge=sapply(2:treesMF$trees[[2]]$Nnode+Ntip(treesMF$trees[[2]]),
                               function(n,e) which(e==n),e=treesMF$trees[[2]]$edge[,2]),
                   frame="none",adj=c(0.5,1),which="right")

## show rounded edge lengths:
edgelabels.cophylo(round(treesMF$trees[[1]]$edge.length,2),cex=0.7,
                   frame="none",adj=c(0.5,1.2),which="left")
edgelabels.cophylo(round(treesMF$trees[[2]]$edge.length,2),cex=0.7,
                   frame="none",adj=c(0.5,1.2),which="right")

#consensus trees
treeM <- root(phy = consensus(modernHumans_M_nonlog, p = 0.5), "BUSHMAN", resolve.root = T, edgelabel = T) #for male tree
treeF <- root(phy = consensus(modernHumans_F_nonlog, p = 0.5), "BUSHMAN", resolve.root = T, edgelabel = T) #for female tree
plot(cophylo(treeM, treeF))

#random code snippets to delete?

layout(matrix(1:2,nrow=1),widths=c(0.8,0.2))
colfunc <- colorRampPalette(c("red","royalblue"))

par(mar=c(5.1,4.1,4.1,2.1))
plot(1:10,ann=FALSE,type="n")
grid()
points(1:10,col=colfunc(10),pch=19,cex=1.5)

xl <- 1
yb <- 1
xr <- 1.5
yt <- 2

par(mar=c(5.1,0.5,4.1,0.5))
plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(1,2),xaxt="n",yaxt="n",bty="n")
rect(
  xl,
  head(seq(yb,yt,(yt-yb)/100),-1),
  xr,
  tail(seq(yb,yt,(yt-yb)/100),-1),
  col=colfunc(100)
)
rect(
  xl,
  head(seq(yb,yt,(yt-yb)/100),-1),
  xr,
  tail(seq(yb,yt,(yt-yb)/100),-1),
  col=rgb(255, 255, 255, max = 255, alpha = 150, names = "white25")[[1]]
)
mtext(seq(from = 0, to = 1, length.out = 11),side=2,at=tail(seq(yb,yt,(yt-yb)/11),-1)-0.05,las=2,cex=0.7)

mccTree_log <- maxCladeCred(modernHumansNeandertal_log)
mccTree_nonlog <- maxCladeCred(modernHumansNeandertal_nonlog)
plot(consensus(modernHumansNeandertal_log, p = 0.5))
plot(consensus(modernHumansNeandertal_nonlog, p = 0.5))
plot(root(mccTree_log, outgroup = "NEANDERTAL"))
plot(root(mccTree_nonlog, outgroup = "NEANDERTAL"))
RF.dist(mccTree_log, mccTree_nonlog)
pp1 <- prop.clades(mccTree_log, modernHumansNeandertal_log, rooted = F) / length(modernHumansNeandertal_log)
pp2 <- prop.clades(mccTree_log, modernHumansNeandertal_nonlog, rooted = F) / length(modernHumansNeandertal_nonlog)
plot(pp1, pp2)

mccTree_log <- maxCladeCred(pan_log)
mccTree_nonlog <- maxCladeCred(pan_nonlog)
plot(root(consensus(pan_log, p = 0.5), "paniscus"))
plot(root(mccTree_log, outgroup = "paniscus"))
plot(root(consensus(pan_nonlog, p = 0.5), "paniscus"))
plot(root(mccTree_nonlog, outgroup = "paniscus"))
pp1 <- prop.clades(mccTree_log, pan_log, rooted = F) / length(modernHumansNeandertal_log)
pp2 <- prop.clades(mccTree_log, pan_nonlog, rooted = F) / length(modernHumansNeandertal_nonlog)
plot(pp1, pp2)
tree <- root(mccTree_log,  outgroup = "paniscus", resolve.root = TRUE)
# tree <- root(tree, "BUSHMAN", resolve.root = T)
plot(tree, type = "c")
pp <- prop.clades(tree, pan_nonlog, rooted = T)/length(pan_nonlog)
edgelabels(round(pp, 3), sapply(1:tree$Nnode + length(tree$tip.label), 
                                match, tree$edge[, 2]), bg = "white", frame = "r", adj = c(0,1.5))


##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

#female, male, joint, and united analyses
setwd("/Volumes/macOS/Users/nikolai/output/")

#read in data
modernHumans_M_log <- read.tree("howells_log.trees")
modernHumans_M_nonlog <- read.tree("howells_nonlog.trees")
modernHumans_F_log <- read.tree("howells_log_F.trees")
modernHumans_F_nonlog <- read.tree("howells_nonlog_F.trees")
modernHumans_MF_log <- read.tree("howells_log_MF.trees")
modernHumans_MF_nonlog <- read.tree("howells_nonlog_MF.trees")
modernHumans_MwF_log <- read.tree("howells_log_MwF.trees")
modernHumans_MwF_nonlog <- read.tree("howells_nonlog_MwF.trees")

#compute summary trees
mccTree_log_M <- maxCladeCred(modernHumans_M_log)
mccTree_nonlog_M <- maxCladeCred(modernHumans_M_nonlog)
consensus_log_M <- consensus(modernHumans_M_log, p = 0.5)
consensus_nonlog_M <- consensus(modernHumans_M_nonlog, p = 0.5)

mccTree_log_F <- maxCladeCred(modernHumans_F_log)
mccTree_nonlog_F <- maxCladeCred(modernHumans_F_nonlog)
consensus_log_F <- consensus(modernHumans_F_log, p = 0.5)
consensus_nonlog_F <- consensus(modernHumans_F_nonlog, p = 0.5)

mccTree_log_MF <- maxCladeCred(modernHumans_MF_log)
mccTree_nonlog_MF <- maxCladeCred(modernHumans_MF_nonlog)
consensus_log_MF <- consensus(modernHumans_MF_log, p = 0.5)
consensus_nonlog_MF <- consensus(modernHumans_MF_nonlog, p = 0.5)

mccTree_log_MwF <- maxCladeCred(modernHumans_MwF_log)
mccTree_nonlog_MwF <- maxCladeCred(modernHumans_MwF_nonlog)
consensus_log_MwF <- consensus(modernHumans_MwF_log, p = 0.5)
consensus_nonlog_MwF <- consensus(modernHumans_MwF_nonlog, p = 0.5)

#plot

par(mfrow = c(1,1))

plot(consensus_log_MF); title("geometric (log-transformed) mvBM, consensus tree (p > 0.5)")
pp <- prop.clades(consensus_log_MF, modernHumans_MF_log, rooted = F)/length(modernHumans_MF_log)
edgelabels(round(pp, 2), sapply(1:consensus_log_MF$Nnode + length(consensus_log_MF$tip.label), 
                                match, consensus_log_MF$edge[, 2]), bg = "white", frame = "none", adj = c(0,1.5))

plot(consensus_nonlog_MF); title("raw (not log-transformed) mvBM, consensus tree (p > 0.5)")
pp <- prop.clades(consensus_nonlog_MF, modernHumans_MF_nonlog, rooted = F)/length(modernHumans_MF_nonlog)
edgelabels(round(pp, 2), sapply(1:consensus_nonlog_MF$Nnode + length(consensus_nonlog_MF$tip.label), 
                                match, consensus_nonlog_MF$edge[, 2]), bg = "white", frame = "none", adj = c(0,1.5))

plot(mccTree_log_MF); title("geometric (log-transformed) mvBM, mcc tree (p > 0.5)")
pp <- prop.clades(mccTree_log_MF, modernHumans_MF_log, rooted = F)/length(modernHumans_MF_log)
edgelabels(round(pp, 2), sapply(1:mccTree_log_MF$Nnode + length(mccTree_log_MF$tip.label), 
                                match, mccTree_log_MF$edge[, 2]), bg = "white", frame = "none", adj = c(0,1.5))

plot(mccTree_nonlog_MF); title("raw (not log-transformed) mvBM, mcc tree (p > 0.5)")
pp <- prop.clades(mccTree_nonlog_MF, modernHumans_MF_nonlog, rooted = F)/length(modernHumans_MF_nonlog)
edgelabels(round(pp, 2), sapply(1:mccTree_nonlog_MF$Nnode + length(mccTree_nonlog_MF$tip.label), 
                                match, mccTree_nonlog_MF$edge[, 2]), bg = "white", frame = "none", adj = c(0,1.5))



#comparetrees for Ms, Fs, and combined; need to drop tips
#read in data
noFpops <- c("NMAORI", "SMAORI", "ANYANG", "PHILLIPI")
tiplabs <- c("AINU", "ARIKARA", "SANTACR", "PERU", "NJAPAN", "SJAPAN", "ATAYAL", "HAINAN", "BURIAT", "ESKIMO", "GUAM", "MORIORI", "EASTERI", "MOKAPU", "ANDAMAN", "TEITA", "BUSHMAN", "ZULU", "DOGON", "AUSTRALI", "TASMANIA", "TOLAI", "EGYPT", "NORSE", "BERG", "ZALAVAR")
modernHumans_M_log <- read.tree("howells_log.trees"); for(i in 1:length(modernHumans_M_log)){modernHumans_M_log[[i]] <- drop.tip(modernHumans_M_log[[i]], noFpops)}
modernHumans_M_nonlog <- read.tree("howells_nonlog.trees"); for(i in 1:length(modernHumans_M_nonlog)){modernHumans_M_nonlog[[i]] <- drop.tip(modernHumans_M_nonlog[[i]], noFpops)}
modernHumans_F_log <- read.tree("howells_log_F.trees")
modernHumans_F_nonlog <- read.tree("howells_nonlog_F.trees")
modernHumans_MwF_log <- read.tree("howells_log_MwF.trees"); for(i in 1:length(modernHumans_MwF_log)){modernHumans_MwF_log[[i]] <- drop.tip(modernHumans_MwF_log[[i]], noFpops)}
modernHumans_MwF_nonlog <- read.tree("howells_nonlog_MwF.trees"); for(i in 1:length(modernHumans_MwF_nonlog)){modernHumans_MwF_nonlog[[i]] <- drop.tip(modernHumans_MwF_nonlog[[i]], noFpops)}

#compute summary trees
mccTree_log_M <- maxCladeCred(modernHumans_M_log)
mccTree_nonlog_M <- maxCladeCred(modernHumans_M_nonlog)
consensus_log_M <- consensus(modernHumans_M_log, p = 0.5)
consensus_nonlog_M <- consensus(modernHumans_M_nonlog, p = 0.5)

mccTree_log_F <- maxCladeCred(modernHumans_F_log)
mccTree_nonlog_F <- maxCladeCred(modernHumans_F_nonlog)
consensus_log_F <- consensus(modernHumans_F_log, p = 0.5)
consensus_nonlog_F <- consensus(modernHumans_F_nonlog, p = 0.5)

mccTree_log_MwF <- maxCladeCred(modernHumans_MwF_log)
mccTree_nonlog_MwF <- maxCladeCred(modernHumans_MwF_nonlog)
consensus_log_MwF <- consensus(modernHumans_MwF_log, p = 0.5)
consensus_nonlog_MwF <- consensus(modernHumans_MwF_nonlog, p = 0.5)


tree.dist.matrix(c(mccTree_log_MwF, mccTree_log_F, mccTree_log_M, mccTree_nonlog_MwF, mccTree_nonlog_F, mccTree_nonlog_M))
tree.dist.matrix(c(consensus_log_MwF, consensus_log_F, consensus_log_M, consensus_nonlog_MwF, consensus_nonlog_F, consensus_nonlog_M))

groups <- list(M_log = modernHumans_M_log, M_raw = modernHumans_M_nonlog, F_log = modernHumans_F_log, F_raw = modernHumans_F_nonlog, both_log = modernHumans_MwF_log, both_raw = modernHumans_MwF_nonlog)

dev.off()
par(mar=c(2,2,2,2), mfrow = c(6,6))
for(i in 1:6){
  for(j in 1:6){
    if(j > (i-1)){
      comparison <- compareTrees(prop.part.df(groups[[i]]), prop.part.df(groups[[j]]), tiplabs)
      plot(comparison, xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1), axes = F);
      abline(b = 1,a = 0, lwd = 2, col = 2)
      axis(side = 1, at = c(0,1)); axis(side = 2, at = c(0,1))
      R2 <- round(cor(comparison)[1,2]^2, 2)
      text(labels = bquote(R^2 == .(R2)), x = .02, y = .98, srt = 0, adj = 0)
    } else {
      plot.new()
    }
    if(i == 1){
      text(labels = names(groups)[j], x = .5, y = 1.13, xpd = NA, srt = 0, adj = 0.5, font = 2, cex = 1.5)
    }
    if(j == 6){
      text(labels = names(groups)[i], x = 1.13, y = .5, xpd = NA, srt = 270, adj = 0.5, font = 2, cex = 1.5)
    }
  }
}

#compare mike's version to my version
par(mfrow = c(1,1))
plot(compareTrees(prop.part.df(modernHumans_M_log), prop.part.df(modernHumans_F_log), tiplabs), xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1), axes = F); 
abline(b = 1,a = 0, lwd = 2, col = 2)
axis(side = 1, at = c(0,1)); axis(side = 2, at = c(0,1))
treeComparisons("howells_log.trees", "howells_log_F.trees")

plot(compareTrees(prop.part.df(modernHumans_M_nonlog), prop.part.df(modernHumans_F_nonlog), tiplabs), xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1), axes = F); 
abline(b = 1,a = 0, lwd = 2, col = 2)
axis(side = 1, at = c(0,1)); axis(side = 2, at = c(0,1))
treeComparisons("howells_nonlog.trees", "howells_nonlog_F.trees")

plot(compareTrees(prop.part.df(modernHumans_M_nonlog), prop.part.df(modernHumans_MwF_nonlog), tiplabs), xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1), axes = F); 
abline(b = 1,a = 0, lwd = 2, col = 2)
axis(side = 1, at = c(0,1)); axis(side = 2, at = c(0,1))
treeComparisons("howells_nonlog.trees", "howells_nonlog_MwF.trees")

# #now look at MDS
# ### rwty stuff ###
# library(rwty)
# #write log and tree files
# modernHumans_M_log <- read.tree("howells_log.trees"); for(i in 1:length(modernHumans_M_log)){modernHumans_M_log[[i]] <- drop.tip(modernHumans_M_log[[i]], noFpops)}
# modernHumans_M_nonlog <- read.tree("howells_nonlog.trees"); for(i in 1:length(modernHumans_M_nonlog)){modernHumans_M_nonlog[[i]] <- drop.tip(modernHumans_M_nonlog[[i]], noFpops)}
# modernHumans_F_log <- read.tree("howells_log_F.trees")
# modernHumans_F_nonlog <- read.tree("howells_nonlog_F.trees")
# modernHumans_MwF_log <- read.tree("howells_log_MwF.trees"); for(i in 1:length(modernHumans_MwF_log)){modernHumans_MwF_log[[i]] <- drop.tip(modernHumans_MwF_log[[i]], noFpops)}
# modernHumans_MwF_nonlog <- read.tree("howells_nonlog_MwF.trees"); for(i in 1:length(modernHumans_MwF_nonlog)){modernHumans_MwF_nonlog[[i]] <- drop.tip(modernHumans_MwF_nonlog[[i]], noFpops)}
# 
# write.tree(modernHumans_M_log, "modernHumans_M_log_26.trees")
# write.tree(modernHumans_M_nonlog, "modernHumans_M_nonlog_26.trees")
# write.tree(modernHumans_MwF_log, "modernHumans_MwF_log_26.trees")
# write.tree(modernHumans_MwF_nonlog, "modernHumans_MwF_nonlog_26.trees")
# 
# log_modernHumans_M_log <- read.table("howells_log.log", header = T)
# log_modernHumans_M_nonlog <- read.table("howells_nonlog.log", header = T)
# log_modernHumans_F_log <- read.table("howells_log_F.log", header = T)
# log_modernHumans_F_nonlog <- read.table("howells_nonlog_F.log", header = T)
# log_modernHumans_MwF_log <- read.table("howells_log_MwF.log", header = T)
# log_modernHumans_MwF_nonlog <- read.table("howells_nonlog_MwF.log", header = T)
# 
# write.table(x = rbind(colnames(log_modernHumans_M_log[,1:4]),log_modernHumans_M_log[,1:4]), file = paste0("log_modernHumans_M_log_rwty"), row.names = F, col.names = T)
# write.table(x = rbind(colnames(log_modernHumans_M_nonlog[,1:4]),log_modernHumans_M_nonlog[,1:4]), file = paste0("log_modernHumans_M_nonlog_rwty"), row.names = F, col.names = T)
# write.table(x = rbind(colnames(log_modernHumans_F_log[,1:4]),log_modernHumans_F_log[,1:4]), file = paste0("log_modernHumans_F_log_rwty"), row.names = F, col.names = T)
# write.table(x = rbind(colnames(log_modernHumans_F_nonlog[,1:4]),log_modernHumans_F_nonlog[,1:4]), file = paste0("log_modernHumans_F_nonlog_rwty"), row.names = F, col.names = T)
# write.table(x = rbind(colnames(log_modernHumans_MwF_log[,1:4]),log_modernHumans_MwF_log[,1:4]), file = paste0("log_modernHumans_MwF_log_rwty"), row.names = F, col.names = T)
# write.table(x = rbind(colnames(log_modernHumans_MwF_nonlog[,1:4]),log_modernHumans_MwF_nonlog[,1:4]), file = paste0("log_modernHumans_MwF_nonlog_rwty"), row.names = F, col.names = T)
# 
# rwty_modernHumans_M_log <- load.trees(file = "modernHumans_M_log_26.trees", logfile = paste0("log_modernHumans_M_log_rwty"), type = "revbayes", gens.per.tree = 1)
# rwty_modernHumans_F_log <- load.trees(file = "howells_log_F.trees", logfile = paste0("log_modernHumans_F_log_rwty"), type = "revbayes", gens.per.tree = 1)
# rwty_modernHumans_MwF_log <- load.trees(file = "modernHumans_MwF_log_26.trees", logfile = paste0("log_modernHumans_MwF_log_rwty"), type = "revbayes", gens.per.tree = 1)
# rwty_modernHumans_M_nonlog <- load.trees(file = "modernHumans_M_nonlog_26.trees", logfile = paste0("log_modernHumans_M_nonlog_rwty"), type = "revbayes", gens.per.tree = 1)
# rwty_modernHumans_F_nonlog <- load.trees(file = "howells_nonlog_F.trees", logfile = paste0("log_modernHumans_F_nonlog_rwty"), type = "revbayes", gens.per.tree = 1)
# rwty_modernHumans_MwF_nonlog <- load.trees(file = "modernHumans_MwF_nonlog_26.trees", logfile = paste0("log_modernHumans_MwF_nonlog_rwty"), type = "revbayes", gens.per.tree = 1)
# 
# rwty_groups <- list(
#   rwty_modernHumans_M_log = rwty_modernHumans_M_log,
#   rwty_modernHumans_M_nonlog = rwty_modernHumans_M_nonlog,
#   rwty_modernHumans_F_log = rwty_modernHumans_F_log,
#   rwty_modernHumans_F_nonlog = rwty_modernHumans_F_nonlog,
#   rwty_modernHumans_MwF_log = rwty_modernHumans_MwF_log,
#   rwty_modernHumans_MwF_nonlog = rwty_modernHumans_MwF_nonlog)
# 
# rwty_groups <- list(
#   rwty_modernHumans_M_log = rwty_modernHumans_M_log)
# 
# # makeplot.treespace(rwty_groups, burnin = 0, n.points = 500)$treespace.heatmap
# treesp <- treespace(rwty_groups, burnin = 0, n.points = 4001)

##doing MDS without rwty
library(ape)
library(phangorn)
library(smacof)
combined_trees = c(modernHumans_M_log[1:100])
rf_matrix = RF.dist(combined_trees)
MDSpts <- cmdscale(rf_matrix) #can also use pcoa(mat)$vectors to get same
stdMDSpts <- mds(rf_matrix, init = MDSpts)$conf
plot(stdMDSpts)
plot(MDSpts)

cols <- c(brewer.pal(9, "Reds")[c(5,8)],brewer.pal(9, "Blues")[c(5,8)],brewer.pal(9, "Greens")[c(5,8)])
for(i in 1:length(unique(treesp$chain))){
  treesp$color[treesp$chain == unique(treesp$chain)[i]] <- cols[i]
}
plot(treesp[,1:2], col = treesp$color)
groupNames <- c("Male-Log", "Male-Raw", "Female-Log", "Female-Raw", "Both-Log", "Both-Raw")
legend(legend = groupNames, col = cols, x = 16.5, y = -12, pch = 1)

for(i in 1:length(unique(treesp$chain))){
  bivn <- kde2d(treesp$x[treesp$chain == unique(treesp$chain)[i]], treesp$y[treesp$chain == unique(treesp$chain)[i]], n = 100)
  contour(bivn, add = T, nlevels = 3, labcex = 1, method = "f", col = cols[i], lwd = 2, drawlabels = F)
}

dataEllipse(treesp$x, treesp$y, levels=c(0.5), add = T, col = 3, lwd = 2)

image(bivn, xlim = c(-4,4), ylim = c(-3,3), 
      main = "bivariate density of x1, x2", xlab = "x1", ylab = "x2")

contour(treesp[,1:2])
treedist(mccTree_nonlog_MF, mccTree_log_MF)
treedist(consensus_nonlog_MF, consensus_log_MF)

##### MAKE THE COMPLETE FIGURE FOR HOWELLS' DATA #####

#read in data and drop the relevant tips
groupNames <- c("Male-Log", "Male-Raw", "Female-Log", "Female-Raw", "Both-Log", "Both-Raw")
noFpops <- c("NMAORI", "SMAORI", "ANYANG", "PHILLIPI")
tiplabs <- c("AINU", "ARIKARA", "SANTACR", "PERU", "NJAPAN", "SJAPAN", "ATAYAL", "HAINAN", "BURIAT", "ESKIMO", "GUAM", "MORIORI", "EASTERI", "MOKAPU", "ANDAMAN", "TEITA", "BUSHMAN", "ZULU", "DOGON", "AUSTRALI", "TASMANIA", "TOLAI", "EGYPT", "NORSE", "BERG", "ZALAVAR")
modernHumans_M_log1 <- read.tree("howells_log.trees"); for(i in 1:length(modernHumans_M_log1)){modernHumans_M_log1[[i]] <- drop.tip(modernHumans_M_log1[[i]], noFpops)}
modernHumans_M_nonlog1 <- read.tree("howells_nonlog.trees"); for(i in 1:length(modernHumans_M_nonlog1)){modernHumans_M_nonlog1[[i]] <- drop.tip(modernHumans_M_nonlog1[[i]], noFpops)}
modernHumans_F_log1 <- read.tree("howells_log_F.trees")
modernHumans_F_nonlog1 <- read.tree("howells_nonlog_F.trees")
modernHumans_M_log2 <- read.tree("howells_log2.trees"); for(i in 1:length(modernHumans_M_log2)){modernHumans_M_log2[[i]] <- drop.tip(modernHumans_M_log2[[i]], noFpops)}
modernHumans_M_nonlog2 <- read.tree("howells_nonlog2.trees"); for(i in 1:length(modernHumans_M_nonlog2)){modernHumans_M_nonlog2[[i]] <- drop.tip(modernHumans_M_nonlog2[[i]], noFpops)}
modernHumans_F_log2 <- read.tree("howells_log_F2.trees")
modernHumans_F_nonlog2 <- read.tree("howells_nonlog_F2.trees")


groups <- list(M_log1 = modernHumans_M_log1, M_raw1 = modernHumans_M_nonlog1, F_log1 = modernHumans_F_log1, F_raw1 = modernHumans_F_nonlog1,
               M_log2 = modernHumans_M_log2, M_raw2 = modernHumans_M_nonlog2, F_log2 = modernHumans_F_log2, F_raw2 = modernHumans_F_nonlog2)
groupNames <- v <- c("log(males)", "males", "log(females)", "females")

dev.off()
par(mar=c(2,2,2,2), mfrow = c(4,4))
for(i in 1:4){
  for(j in 1:4){
    if(j > i){
      ### PLOTTING UPPER RIGHT TRI ###
      comparison <- compareTrees(prop.part.df(groups[[i]]), prop.part.df(groups[[j]]), tiplabs)
      plot(comparison, xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1), axes = F);
      abline(b = 1,a = 0, lwd = 2, col = 1)
      axis(side = 3, at = c(0,1)); axis(side = 4, at = c(0,1))
      R2 <- round(cor(comparison)[1,2]^2, 3)
      text(labels = bquote(R^2 == .(R2)), x = .02, y = .98, srt = 0, adj = 0)
      betwVwithRF <- mean(sapply(1:200, function(x) RF.dist(groups[[i]][[sample(1:length(groups[[i]]), 1)]], groups[[j]][[sample(1:length(groups[[j]]), 1)]])[1])) / 
        mean(sapply(1:200, function(x) RF.dist(groups[[i]][sample(1:length(groups[[i]]), 2)])[1])); betwVwithRF <- round(betwVwithRF, 2)
      text(labels = paste0("rf-ratio = ", betwVwithRF), x = .02, y = .90, srt = 0, adj = 0)
      mccRF <- RF.dist(maxCladeCred(groups[[i]]), maxCladeCred(groups[[j]]))[1]
      text(labels = paste0("mcc-rf = ", mccRF), x = .02, y = .82, srt = 0, adj = 0)
      if(i == 1){
        text(labels = groupNames[j], x = .5, y = 1.13, xpd = NA, srt = 0, adj = 0.5, font = 2, cex = 1.5)
      }
      if(j == 4){
        text(labels = groupNames[i], x = 1.13, y = .5, xpd = NA, srt = 270, adj = 0.5, font = 2, cex = 1.5)
      }
      
    } else if (j == i){
      ### PLOTTING DIAGONALS ###
      comparison <- compareTrees(prop.part.df(groups[[i]]), prop.part.df(groups[[i+4]]), tiplabs)
      plot(comparison, xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1), axes = F);
      abline(b = 1,a = 0, lwd = 2, col = 1)
      axis(side = 3, at = c(0,1)); axis(side = 4, at = c(0,1))
      # R2 <- round(cor(comparison)[1,2]^2, 3)
      # text(labels = bquote(R^2 == .(R2)), x = .02, y = .98, srt = 0, adj = 0)
      # betwVwithRF <- mean(sapply(1:200, function(x) RF.dist(groups[[i]][[sample(1:length(groups[[i]]), 1)]], groups[[i+4]][[sample(1:length(groups[[i+4]]), 1)]])[1])) / 
      #   mean(sapply(1:200, function(x) RF.dist(groups[[i]][sample(1:length(groups[[i]]), 2)])[1])); betwVwithRF <- round(betwVwithRF, 2)
      # text(labels = paste0("rf-ratio = ", betwVwithRF), x = .02, y = .90, srt = 0, adj = 0)
      # mccRF <- RF.dist(maxCladeCred(groups[[i]]), maxCladeCred(groups[[i+4]]))[1]
      # text(labels = paste0("mcc-rf = ", mccRF), x = .02, y = .82, srt = 0, adj = 0)
      if(i == 1){
        text(labels = groupNames[j], x = .5, y = 1.13, xpd = NA, srt = 0, adj = 0.5, font = 2, cex = 1.5)
      }
      if(j == 4){
        text(labels = groupNames[i], x = 1.13, y = .5, xpd = NA, srt = 270, adj = 0.5, font = 2, cex = 1.5)
      }
      
      sampleSize <- 500
      thinPoints <- 2
      sample.i1 <- round(seq(1, length(groups[[i]]), length.out = sampleSize))
      sample.i2 <- round(seq(1, length(groups[[i+4]]), length.out = sampleSize))
      combined_trees = c(groups[[i]][sample.i1], (groups[[i+4]][sample.i2]))
      rf_matrix = RF.dist(combined_trees)
      MDSpts <- cmdscale(rf_matrix)
      stdMDSpts <- mds(rf_matrix, init = MDSpts)$conf
      par(new = TRUE)
      plot(stdMDSpts[seq(1,(2*sampleSize), by = thinPoints),], pch = c(rep(1, sampleSize/thinPoints), rep(3, sampleSize/thinPoints)), xaxt = "n", yaxt = "n", frame.plot = F,
           col = "grey75", cex = 0.75, xlim = c(-1,1), ylim = c(-1,1))
      axis(side = 1, at = c(-1,1)); axis(side = 2, at = c(-1,1))
      bivn1 <- kde2d(stdMDSpts[1:sampleSize,1], stdMDSpts[1:sampleSize,2], n = 40, 
                     lims = c(min(stdMDSpts[1:sampleSize,1])-.1, max(stdMDSpts[1:sampleSize,1])+.1,
                              min(stdMDSpts[1:sampleSize,2])-.1, max(stdMDSpts[1:sampleSize,2])+.1))
      
      pp <- vector()
      for (k in 1:sampleSize){
        z.x <- max(which(bivn1$x < stdMDSpts[k,1]))
        z.y <- max(which(bivn1$y < stdMDSpts[k,2]))
        pp[k] <- bivn1$z[z.x, z.y]
      }
      confidencebound <- quantile(pp, probs =  0.05, na.rm = TRUE)
      contour(bivn1, add = T, levels = confidencebound, labcex = 1, col = 1, lwd = 2, drawlabels = F)
      
      bivn2 <- kde2d(stdMDSpts[(sampleSize+1):(2*sampleSize),1], stdMDSpts[(sampleSize+1):(2*sampleSize),2], n = 40, 
                     lims = c(min(stdMDSpts[(sampleSize+1):(2*sampleSize),1])-.1, max(stdMDSpts[(sampleSize+1):(2*sampleSize),1])+.1,
                              min(stdMDSpts[(sampleSize+1):(2*sampleSize),2])-.1, max(stdMDSpts[(sampleSize+1):(2*sampleSize),2])+.1))
      pp <- vector()
      for (k in (sampleSize+1):(2*sampleSize)){
        z.x <- max(which(bivn2$x < stdMDSpts[k,1]))
        z.y <- max(which(bivn2$y < stdMDSpts[k,2]))
        pp[k-sampleSize] <- bivn2$z[z.x, z.y]
      }
      confidencebound <- quantile(pp, 0.05, na.rm = TRUE)
      contour(bivn2, add = T, levels = confidencebound, labcex = 1, col = "grey50", lwd = 2, drawlabels = F)
      
      par(new = TRUE)
      plot(comparison, xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1), axes = F);
      abline(b = 1,a = 0, lwd = 2, col = 1)
      
      
    }  else if (j < i){
      ### PLOTTING LOWER LEFT TRI ###
      sampleSize <- 500
      thinPoints <- 2
      sample.i <- round(seq(1, length(groups[[i]]), length.out = sampleSize))
      sample.j <- round(seq(1, length(groups[[j]]), length.out = sampleSize))
      combined_trees = c(groups[[i]][sample.i], (groups[[j]][sample.j]))
      rf_matrix = RF.dist(combined_trees)
      MDSpts <- cmdscale(rf_matrix)
      stdMDSpts <- mds(rf_matrix, init = MDSpts)$conf
      # plot(stdMDSpts, pch = c(rep(1, sampleSize), rep(3, sampleSize)), axes = F, col = c(rep(1, sampleSize), rep("grey50", sampleSize)))
      plot(stdMDSpts[seq(1,(2*sampleSize), by = thinPoints),], pch = c(rep(1, sampleSize/thinPoints), rep(3, sampleSize/thinPoints)),
           col = c(rep(1, sampleSize/thinPoints), rep("grey50", sampleSize/thinPoints)), axes = F, cex = 0.75, xlim = c(-1,1), ylim = c(-1,1))
      axis(side = 1, at = c(-1,1)); axis(side = 2, at = c(-1,1))
      bivn1 <- kde2d(stdMDSpts[1:sampleSize,1], stdMDSpts[1:sampleSize,2], n = 40, 
                     lims = c(min(stdMDSpts[1:sampleSize,1])-.1, max(stdMDSpts[1:sampleSize,1])+.1,
                              min(stdMDSpts[1:sampleSize,2])-.1, max(stdMDSpts[1:sampleSize,2])+.1))
      
      #gets density (z) at each point in the sample
      pp <- vector()
      for (k in 1:sampleSize){
        z.x <- max(which(bivn1$x < stdMDSpts[k,1]))
        z.y <- max(which(bivn1$y < stdMDSpts[k,2]))
        pp[k] <- bivn1$z[z.x, z.y]
      }
      #find where you have to exclude probs= points to exclude those points on the countour
      confidencebound <- quantile(pp, probs =  0.05, na.rm = TRUE)
      contour(bivn1, add = T, levels = confidencebound, labcex = 1, col = 1, lwd = 2, drawlabels = F)
      
      bivn2 <- kde2d(stdMDSpts[(sampleSize+1):(2*sampleSize),1], stdMDSpts[(sampleSize+1):(2*sampleSize),2], n = 40, 
                     lims = c(min(stdMDSpts[(sampleSize+1):(2*sampleSize),1])-.1, max(stdMDSpts[(sampleSize+1):(2*sampleSize),1])+.1,
                              min(stdMDSpts[(sampleSize+1):(2*sampleSize),2])-.1, max(stdMDSpts[(sampleSize+1):(2*sampleSize),2])+.1))
      pp <- vector()
      for (k in (sampleSize+1):(2*sampleSize)){
        z.x <- max(which(bivn2$x < stdMDSpts[k,1]))
        z.y <- max(which(bivn2$y < stdMDSpts[k,2]))
        pp[k-sampleSize] <- bivn2$z[z.x, z.y]
      }
      confidencebound <- quantile(pp, 0.05, na.rm = TRUE)
      contour(bivn2, add = T, levels = confidencebound, labcex = 1, col = "grey50", lwd = 2, drawlabels = F)
      
    }
  }
}


## get the above version of the figure, just 2x2
groupNames <- v <- c("log(males)", "males", "log(females)", "females")
focus <- c("males", "females")
groups <- list(M_log1 = modernHumans_M_log1, M_raw1 = modernHumans_M_nonlog1, F_log1 = modernHumans_F_log1, F_raw1 = modernHumans_F_nonlog1,
               M_log2 = modernHumans_M_log2, M_raw2 = modernHumans_M_nonlog2, F_log2 = modernHumans_F_log2, F_raw2 = modernHumans_F_nonlog2)
groups <- groups[c(which(groupNames == focus[1]), which(groupNames == focus[2]), which(groupNames == focus[1]) + 4, which(groupNames == focus[2]) + 4)]
groupNames <- focus

dev.off()
png(filename = "~/mvbm_manuscript/figures/sensitivityHowells_2x2.png", width = 600, height = 600)
par(mar=c(2,2,2,2), mfrow = c(2,2))
for(i in 1:2){
  for(j in 1:2){
    print(c(i,j))
    if(j > i){
      ### PLOTTING UPPER RIGHT TRI ###
      comparison <- compareTrees(prop.part.df(groups[[i]]), prop.part.df(groups[[j]]), tiplabs)
      plot(comparison, xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1), axes = F);
      abline(b = 1,a = 0, lwd = 2, col = 1)
      axis(side = 3, at = c(0,1)); axis(side = 4, at = c(0,1))
      R2 <- round(cor(comparison)[1,2]^2, 3)
      text(labels = bquote(R^2 == .(R2)), x = .02, y = .98, srt = 0, adj = 0)
      betwVwithRF <- mean(sapply(1:200, function(x) RF.dist(groups[[i]][[sample(1:length(groups[[i]]), 1)]], groups[[j]][[sample(1:length(groups[[j]]), 1)]])[1])) / 
        mean(sapply(1:200, function(x) RF.dist(groups[[i]][sample(1:length(groups[[i]]), 2)])[1])); betwVwithRF <- round(betwVwithRF, 2)
      text(labels = paste0("rf-ratio = ", betwVwithRF), x = .02, y = .90, srt = 0, adj = 0)
      mccRF <- RF.dist(maxCladeCred(groups[[i]]), maxCladeCred(groups[[j]]))[1]
      text(labels = paste0("mcc-rf = ", mccRF), x = .02, y = .82, srt = 0, adj = 0)
      if(i == 1){
        text(labels = groupNames[j], x = .5, y = 1.1, xpd = NA, srt = 0, adj = 0.5, font = 2, cex = 1.5)
      }
      if(j == 2){
        text(labels = groupNames[i], x = 1.1, y = .5, xpd = NA, srt = 270, adj = 0.5, font = 2, cex = 1.5)
      }
      
    } else if (j == i){
      ### PLOTTING DIAGONALS ###
      comparison <- compareTrees(prop.part.df(groups[[i]]), prop.part.df(groups[[i+2]]), tiplabs)
      plot(comparison, xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1), axes = F);
      abline(b = 1,a = 0, lwd = 2, col = 1)
      axis(side = 3, at = c(0,1)); axis(side = 4, at = c(0,1))
      # R2 <- round(cor(comparison)[1,2]^2, 3)
      # text(labels = bquote(R^2 == .(R2)), x = .02, y = .98, srt = 0, adj = 0)
      # betwVwithRF <- mean(sapply(1:200, function(x) RF.dist(groups[[i]][[sample(1:length(groups[[i]]), 1)]], groups[[i+4]][[sample(1:length(groups[[i+4]]), 1)]])[1])) / 
      #   mean(sapply(1:200, function(x) RF.dist(groups[[i]][sample(1:length(groups[[i]]), 2)])[1])); betwVwithRF <- round(betwVwithRF, 2)
      # text(labels = paste0("rf-ratio = ", betwVwithRF), x = .02, y = .90, srt = 0, adj = 0)
      # mccRF <- RF.dist(maxCladeCred(groups[[i]]), maxCladeCred(groups[[i+4]]))[1]
      # text(labels = paste0("mcc-rf = ", mccRF), x = .02, y = .82, srt = 0, adj = 0)
      if(i == 1){
        text(labels = groupNames[j], x = .5, y = 1.1, xpd = NA, srt = 0, adj = 0.5, font = 2, cex = 1.5)
      }
      if(j == 2){
        text(labels = groupNames[i], x = 1.1, y = .5, xpd = NA, srt = 270, adj = 0.5, font = 2, cex = 1.5)
      }
      
      sampleSize <- 500
      thinPoints <- 2
      sample.i1 <- round(seq(1, length(groups[[i]]), length.out = sampleSize))
      sample.i2 <- round(seq(1, length(groups[[i+2]]), length.out = sampleSize))
      combined_trees = c(groups[[i]][sample.i1], (groups[[i+2]][sample.i2]))
      rf_matrix = RF.dist(combined_trees)
      MDSpts <- cmdscale(rf_matrix)
      stdMDSpts <- mds(rf_matrix, init = MDSpts)$conf
      par(new = TRUE)
      plot(stdMDSpts[seq(1,(2*sampleSize), by = thinPoints),], pch = c(rep(1, sampleSize/thinPoints), rep(3, sampleSize/thinPoints)), xaxt = "n", yaxt = "n", frame.plot = F,
           col = "grey75", cex = 0.75, xlim = c(-1,1), ylim = c(-1,1))
      axis(side = 1, at = c(-1,1)); axis(side = 2, at = c(-1,1))
      bivn1 <- kde2d(stdMDSpts[1:sampleSize,1], stdMDSpts[1:sampleSize,2], n = 40, 
                     lims = c(min(stdMDSpts[1:sampleSize,1])-.1, max(stdMDSpts[1:sampleSize,1])+.1,
                              min(stdMDSpts[1:sampleSize,2])-.1, max(stdMDSpts[1:sampleSize,2])+.1))
      
      pp <- vector()
      for (k in 1:sampleSize){
        z.x <- max(which(bivn1$x < stdMDSpts[k,1]))
        z.y <- max(which(bivn1$y < stdMDSpts[k,2]))
        pp[k] <- bivn1$z[z.x, z.y]
      }
      confidencebound <- quantile(pp, probs =  0.05, na.rm = TRUE)
      contour(bivn1, add = T, levels = confidencebound, labcex = 1, col = 1, lwd = 2, drawlabels = F)
      
      bivn2 <- kde2d(stdMDSpts[(sampleSize+1):(2*sampleSize),1], stdMDSpts[(sampleSize+1):(2*sampleSize),2], n = 40, 
                     lims = c(min(stdMDSpts[(sampleSize+1):(2*sampleSize),1])-.1, max(stdMDSpts[(sampleSize+1):(2*sampleSize),1])+.1,
                              min(stdMDSpts[(sampleSize+1):(2*sampleSize),2])-.1, max(stdMDSpts[(sampleSize+1):(2*sampleSize),2])+.1))
      pp <- vector()
      for (k in (sampleSize+1):(2*sampleSize)){
        z.x <- max(which(bivn2$x < stdMDSpts[k,1]))
        z.y <- max(which(bivn2$y < stdMDSpts[k,2]))
        pp[k-sampleSize] <- bivn2$z[z.x, z.y]
      }
      confidencebound <- quantile(pp, 0.05, na.rm = TRUE)
      contour(bivn2, add = T, levels = confidencebound, labcex = 1, col = "grey50", lwd = 2, drawlabels = F)
      
      par(new = TRUE)
      plot(comparison, xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1), axes = F);
      abline(b = 1,a = 0, lwd = 2, col = 1)
      
      
    }  else if (j < i){
      ### PLOTTING LOWER LEFT TRI ###
      sampleSize <- 500
      thinPoints <- 2
      sample.i <- round(seq(1, length(groups[[i]]), length.out = sampleSize))
      sample.j <- round(seq(1, length(groups[[j]]), length.out = sampleSize))
      combined_trees = c(groups[[i]][sample.i], (groups[[j]][sample.j]))
      rf_matrix = RF.dist(combined_trees)
      MDSpts <- cmdscale(rf_matrix)
      stdMDSpts <- mds(rf_matrix, init = MDSpts)$conf
      # plot(stdMDSpts, pch = c(rep(1, sampleSize), rep(3, sampleSize)), axes = F, col = c(rep(1, sampleSize), rep("grey50", sampleSize)))
      plot(stdMDSpts[seq(1,(2*sampleSize), by = thinPoints),], pch = c(rep(1, sampleSize/thinPoints), rep(3, sampleSize/thinPoints)),
           col = c(rep(1, sampleSize/thinPoints), rep("grey50", sampleSize/thinPoints)), axes = F, cex = 0.75, xlim = c(-1,1), ylim = c(-1,1))
      axis(side = 1, at = c(-1,1)); axis(side = 2, at = c(-1,1))
      bivn1 <- kde2d(stdMDSpts[1:sampleSize,1], stdMDSpts[1:sampleSize,2], n = 40, 
                     lims = c(min(stdMDSpts[1:sampleSize,1])-.1, max(stdMDSpts[1:sampleSize,1])+.1,
                              min(stdMDSpts[1:sampleSize,2])-.1, max(stdMDSpts[1:sampleSize,2])+.1))
      
      #gets density (z) at each point in the sample
      pp <- vector()
      for (k in 1:sampleSize){
        z.x <- max(which(bivn1$x < stdMDSpts[k,1]))
        z.y <- max(which(bivn1$y < stdMDSpts[k,2]))
        pp[k] <- bivn1$z[z.x, z.y]
      }
      #find where you have to exclude probs= points to exclude those points on the countour
      confidencebound <- quantile(pp, probs =  0.05, na.rm = TRUE)
      contour(bivn1, add = T, levels = confidencebound, labcex = 1, col = 1, lwd = 2, drawlabels = F)
      
      bivn2 <- kde2d(stdMDSpts[(sampleSize+1):(2*sampleSize),1], stdMDSpts[(sampleSize+1):(2*sampleSize),2], n = 40, 
                     lims = c(min(stdMDSpts[(sampleSize+1):(2*sampleSize),1])-.1, max(stdMDSpts[(sampleSize+1):(2*sampleSize),1])+.1,
                              min(stdMDSpts[(sampleSize+1):(2*sampleSize),2])-.1, max(stdMDSpts[(sampleSize+1):(2*sampleSize),2])+.1))
      pp <- vector()
      for (k in (sampleSize+1):(2*sampleSize)){
        z.x <- max(which(bivn2$x < stdMDSpts[k,1]))
        z.y <- max(which(bivn2$y < stdMDSpts[k,2]))
        pp[k-sampleSize] <- bivn2$z[z.x, z.y]
      }
      confidencebound <- quantile(pp, 0.05, na.rm = TRUE)
      contour(bivn2, add = T, levels = confidencebound, labcex = 1, col = "grey50", lwd = 2, drawlabels = F)
      
    }
  }
}
dev.off()

## other stuff
plot(mccTree_log_MF)
plot(consensus_nonlog_MF)
plot(mccTree_nonlog_MF)


plot(consensus(modernHumansNeandertal_log, p = 0.5))
plot(consensus(modernHumansNeandertal_nonlog, p = 0.5))
plot(root(mccTree_log, outgroup = "NEANDERTAL"))
plot(root(mccTree_nonlog, outgroup = "NEANDERTAL"))
RF.dist(mccTree_log, mccTree_nonlog)
pp1 <- prop.clades(mccTree_log, modernHumansNeandertal_log, rooted = F) / length(modernHumansNeandertal_log)
pp2 <- prop.clades(mccTree_log, modernHumansNeandertal_nonlog, rooted = F) / length(modernHumansNeandertal_nonlog)
plot(pp1, pp2)


#read in data for multi-Neandertal analysis
neandertal_M_nonlog <- read.tree("neandertal_nonlog_multiN.trees")
neandertal_M_nonlog_consensus <- consensus(neandertal_M_nonlog, p = .5)
plot(neandertal_M_nonlog_consensus)
pp <- prop.clades(neandertal_M_nonlog_consensus, neandertal_M_nonlog, rooted = T)/length(neandertal_M_nonlog)
drawSupportOnEdges(round(pp, 2), col = 2)
edgelabels(round(pp, 2), sapply(1:neandertal_M_nonlog_consensus$Nnode + length(neandertal_M_nonlog_consensus$tip.label), 
                                match, neandertal_M_nonlog_consensus$edge[, 2]), frame = "r", adj = c(0,1.5))

neandertal_M_nonlog_maxCladeCred <- maxCladeCred(neandertal_M_nonlog)
plot(neandertal_M_nonlog_maxCladeCred)
pp <- prop.clades(neandertal_M_nonlog_maxCladeCred, neandertal_M_nonlog, rooted = T)/length(neandertal_M_nonlog)
drawSupportOnEdges(round(pp, 2), col = 2)
edgelabels(round(pp, 2), sapply(1:neandertal_M_nonlog_maxCladeCred$Nnode + length(neandertal_M_nonlog_maxCladeCred$tip.label), 
                                match, neandertal_M_nonlog_maxCladeCred$edge[, 2]), frame = "r", adj = c(0,1.5))

# reading in trees where males and females of each population are independent tips
setwd("/Volumes/macOS/Users/nikolai/output/")
par(mfrow = c(1,1))
modernHumans_MF_log <- read.tree("howells_log_MF.trees")
modernHumans_MF_nonlog <- read.tree("howells_nonlog_MF.trees")
#nonlog
bothMFmcc_nonlog <- maxCladeCred(modernHumans_MF_nonlog); 
ppNL <- prop.clades(bothMFmcc_nonlog, modernHumans_MF_nonlog, rooted = T)/length(modernHumans_MF_nonlog)
bothMFmcc_nonlog$node.label <- round(ppNL, 2)
modernHumans_MF_nonlog_R <- root(phy = modernHumans_MF_nonlog, outgroup = c("BUSHMAN.F", 'BUSHMAN.M'))
bothMFmcc_nonlog <-  reroot(bothMFmcc_nonlog, node.number = 89, position = 0.1)
bothMFmcc_nonlog$node.label[1] <- 1 #specify pp of 1
#log
bothMFmcc_log <- maxCladeCred(modernHumans_MF_log)
ppL <- prop.clades(bothMFmcc_log, modernHumans_MF_log, rooted = T)/length(modernHumans_MF_log)
bothMFmcc_log$node.label <- round(ppL, 2)
modernHumans_MF_log_R <- root(phy = modernHumans_MF_log, outgroup = c("BUSHMAN.F", 'BUSHMAN.M'))
bothMFmcc_log <-  reroot(bothMFmcc_log, node.number = 89, position = 0.1)
bothMFmcc_log$node.label[1] <- 1 #specify pp of 1

plot(bothMFmcc_nonlog)
plot(bothMFmcc_log)

#plot mcc tree figure
logVnonlogMF <- cophylo(bothMFmcc_nonlog, bothMFmcc_log)
plot(logVnonlogMF, scale.bar=c(.5,.5), mar=c(.1,.1,.5,.1))
title("\n\nArithmetic                                                                                                                   Geometric")

#plot circles of particular shade of grey on corresponding node
colfunc<-colorRampPalette(c("white", "black"))
nodelabels.cophylo(pch=21, frame="none", bg=colfunc(100)[as.numeric(logVnonlogMF$trees[[1]]$node.label)*100], cex=1.5, which = "left")
nodelabels.cophylo(pch=21, frame="none", bg=colfunc(100)[as.numeric(logVnonlogMF$trees[[2]]$node.label)*100], cex=1.5, which= "right")

#horizontal pp
xl <- -0.3; yb <- -.1; xr <- 0.3; yt <- -0.05; 
rect(
  head(seq(xl,xr,(xr-xl)/100),-1),
  yb,
  tail(seq(xl,xr,(xr-xl)/100),-1),
  yt,
  col=colfunc(100)
)
text(labels = seq(from = 0, to = 1, length.out = 11), y = yb-(yt-yb)/3, x = seq(xl,xr,length.out = 11), las=2, cex=1)
text(labels = "posterior probability", x = (xl+xr)/2, y = yt+(yt-yb)/3)


bothMFconsensus_nonlog <- consensus(modernHumans_MF_nonlog, p = 0.5)
bothMFconsensus_nonlog <-  root(phy = bothMFconsensus_nonlog, outgroup = c("BUSHMAN.F", 'BUSHMAN.M'))
bothMFconsensus_log <- consensus(modernHumans_MF_log, p = 0.5)
bothMFconsensus_log <-  root(phy = bothMFconsensus_log, outgroup = c("BUSHMAN.F", 'BUSHMAN.M'))

plot(bothMFconsensus_nonlog)
plot(bothMFconsensus_log)
plot(cophylo(bothMFconsensus_nonlog, bothMFconsensus_log))

