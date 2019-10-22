library(ape)
library(phangorn)

convNum2Str <- function(nums, key){
  sapply(1:length(nums), function(x) key[nums[x]])
}

prop.part.df <- function(trees, cutoff = 0.001){
  numTrees <- length(trees)
  out <- prop.part(trees)
  outList <- as.list.data.frame(out)
  pps <- attributes(outList)$number/numTrees
  props <- data.frame(pps, as.matrix(outList)); colnames(props) <- c("postProbs", "clade")
  props <- props[order(-pps),]
  props <- props[props[,1] > cutoff,]
  rownames(props) <- 1:nrow(props)
  props$cladeNames <- sapply(1:length(props[,1]), function(x) sort(convNum2Str(props$clade[[x]], attributes(outList)$labels)))
  props <- props[,c(1,3)]
  props
}

badTree <- function(tree) { any(((tree$edge[,1] == 0) + (tree$edge[,2] == 0)) == 2) }

trees <- rmtree(N = 50, n = 20,rooted = F)
for(i in 1:50){
trees[[50+2*i-1]] <- rNNI(trees[[1]])
trees[[50+2*i]] <- rSPR(trees[[1]])
}
con <- consensus(trees, p = 0.5)
greedy <- con
probs <- prop.part.df(trees)
badProbs <- probs[probs[,1] < 0.5,]
tryProb <- 0
while(!is.binary(greedy)){
  tryProb <- tryProb + 1
  bipart <- badProbs[tryProb,]$cladeNames[[1]]
  tips <- sapply(1:length(bipart), function(x) which(greedy$tip.label == bipart[x]))
  nodes <- as.data.frame(table((con)$edge[,1]))
  polyNode <- as.numeric(as.character(nodes$Var1))[nodes[,2] > 2]
  
  
}
plot(con)
nodelabels()
tiplabels()
con$edge
table(multi2di(con)$edge[,1])
which(table((con)$edge[,1]) > 2)[1]

con$Nnode
greed <- rtree(n = length(con$tip.label), tip.label = con$tip.label, rooted = F)
greed$edge <- matrix(nrow = length(greed$tip.label), ncol = 2, data = c(rep(length(greed$tip.label) +1, length(greed$tip.label)), 1:length(greed$tip.label)))
greed$edge.length <- rep(1, length(greed$tip.label))
greed$Nnode <- 1
plot(greed)

greed$edge <- rbind(greed$edge, c(21, 22))
greed$Nnode <- greed$Nnode + 1
greed$edge.length <- c(greed$edge.length,1)
plot(greed)
greed$edge[sapply(1:5, function(x) which(greed$edge[,2] == x)),1] <- rep(22,5)
plot(greed)
greed$edge
rtree

multi2di <- function(phy, random = T, n)
{
  ## n: number of tips of phy
  degree <- tabulate(phy$edge[, 1])
  target <- which(degree > 2)
  if (!length(target)) return(phy)
  nb.edge <- dim(phy$edge)[1]
  nextnode <- n + phy$Nnode + 1L
  new.edge <- edge2delete <- NULL
  wbl <- FALSE
  if (!is.null(phy$edge.length)) {
    wbl <- TRUE
    new.edge.length <- NULL
  }
  
  for (node in target) {
    ind <- which(phy$edge[, 1] == node)
    N <- length(ind)
    desc <- phy$edge[ind, 2]
    if (random) {
      ## if we shuffle the descendants, we need to eventually
      ## reorder the corresponding branch lenghts (see below)
      ## so we store the result of sample()
      tmp <- sample(length(desc))
      desc <- desc[tmp]
      res <- rtree(N)$edge
    } else {
      res <- matrix(0L, 2*N - 2, 2)
      res[, 1] <- N + rep(1:(N - 1), each = 2)
      res[, 2] <- N + rep(2:N, each = 2)
      res[seq(1, by = 2, length.out = N - 1), 2] <- 1:(N - 1)
      res[length(res)] <- N
    }
    if (wbl) {
      ## keep the branch lengths coming from `node'
      el <- numeric(dim(res)[1]) # initialized with 0's
      el[res[, 2] <= N] <-
        if (random) phy$edge.length[ind][tmp] else phy$edge.length[ind]
    }
    ## now substitute the nodes in `res'
    ## `node' stays at the "root" of these new
    ## edges whereas their "tips" are `desc'
    Nodes <- c(node, nextnode:(nextnode + N - 3L))
    res[, 1] <- Nodes[res[, 1] - N]
    tmp <- res[, 2] > N
    res[tmp, 2] <- Nodes[res[tmp, 2] - N]
    res[!tmp, 2] <- desc[res[!tmp, 2]]
    new.edge <- rbind(new.edge, res)
    edge2delete <- c(edge2delete, ind)
    if (wbl) new.edge.length <- c(new.edge.length, el)
    nextnode <- nextnode + N - 2L
    phy$Nnode <- phy$Nnode + N - 2L
  }
  phy$edge <- rbind(phy$edge[-edge2delete, ], new.edge)
  if (wbl)
    phy$edge.length <- c(phy$edge.length[-edge2delete], new.edge.length)
  if (!is.null(attr(phy, "order"))) attr(phy, "order") <- NULL
  if (!is.null(phy$node.label))
    phy$node.label <-
    c(phy$node.label, rep("", phy$Nnode - length(phy$node.label)))
  phy <- reorder_ape(phy, "cladewise", FALSE, n, 1L) # fix by Klaus (2017-01-16)
  
  ## the node numbers are not in increasing order in edge[, 2]: this
  ## will confuse drop.tip and other functions (root), so renumber them
  newNb <- integer(phy$Nnode)
  newNb[1] <- n + 1L
  sndcol <- phy$edge[, 2] > n
  
  ## reorder node labels before changing edge:
  if (!is.null(phy$node.label)) {
    o <- 1 + rank(phy$edge[sndcol, 2])
    ## the root's label is not changed:
    phy$node.label <- phy$node.label[c(1, o)]
  }
  
  ## executed from right to left, so newNb is modified before phy$edge:
  phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2] - n] <-
    n + 2:phy$Nnode
  phy$edge[, 1] <- newNb[phy$edge[, 1] - n]
  phy
}

reorder_ape <- function(x, order, index.only, nb.tip, io)
{
  nb.edge <- dim(x$edge)[1]
  if (!is.null(attr(x, "order")))
    if (attr(x, "order") == order)
      if (index.only) return(1:nb.edge) else return(x)
  nb.node <- x$Nnode
  if (nb.node == 1)
    if (index.only) return(1:nb.edge) else return(x)
  
  if (io == 3) {
    x <- reorder(x)
    neworder <-
      .C(neworder_pruningwise, as.integer(nb.tip),
         as.integer(nb.node), as.integer(x$edge[, 1]),
         as.integer(x$edge[, 2]), as.integer(nb.edge),
         integer(nb.edge))[[6]]
  } else {
    neworder <- reorderRcpp(x$edge, as.integer(nb.tip),
                            as.integer(nb.tip + 1L), io)
  }
  if (index.only) return(neworder)
  x$edge <- x$edge[neworder, ]
  if (!is.null(x$edge.length))
    x$edge.length <- x$edge.length[neworder]
  attr(x, "order") <- order
  x
}

multi2di(con, T, n = length(con$tip.label))
