setwd("/Volumes/1TB/Bailey/")

data_dir <- "data_joint_thresh_means/"
data_dir <- "data/"
data_files <- list.files(data_dir)
data_files <- data_files[grep(data_files, pattern =  "Star")]
data <- lapply(1:length(data_files), function(x) "placeholder")
for(i in 1:length(data)){
  load(paste0(data_dir, data_files[i]))
  data[[i]] <- current_params
}

data <- data[grep(data_files, pattern = "NJ")]
data <- data[sapply(1:length(data), function(x) nrow(data[[x]]$means)) == 6]

mean_of_means <- data[[1]]$means
for(i in 2:length(data)){
  mean_of_means <- mean_of_means + data[[i]]$means
}
mean_of_means <- mean_of_means / length(data)

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

one.to.one.line <- function(x,y,...){
  points(x,y,...)
  abline(a = 0,b = 1, lwd = 4,...)
}

means <- sapply(1:length(data), function(run) as.vector(data[[run]]$means))
pairs(means, col = rgb(0,0,0,0.1), main = "means", lower.panel = panel.cor, upper.panel = one.to.one.line)

cors <- sapply(1:length(data), function(run) as.vector(data[[run]]$cor[upper.tri(data[[run]]$cor)]))
pairs(cors, col = rgb(0,0,0,0.3), main = "correlations", lower.panel = panel.cor, upper.panel = one.to.one.line)

threshes <- sapply(1:length(data), function(run) as.vector(data[[run]]$threshold_mat)[as.vector(data[[run]]$threshold_mat != Inf)])
pairs(threshes, col = rgb(0,0,0,0.3), main = "thresholds", lower.panel = panel.cor, upper.panel = one.to.one.line)

mean_of_cors <- data[[1]]$cor
for(i in 2:length(data)){
  mean_of_cors <- mean_of_cors + data[[i]]$cor
}
mean_of_cors <- mean_of_cors / length(data)

eigensystem <- eigen(mean_of_cors)

plot(cumsum(eigensystem$values), type = "l"); abline(v = which.max(cumsum(eigensystem$values)))


maha <- function(x1, x2, vcvm, squared = FALSE){
  squaredDist <- t(x1-x2) %*% solve(vcvm) %*% (x1-x2)
  if(!squared){
    squaredDist
  } else {
    squaredDist^0.5
  }
} 
mahaMatrix <- function(data, vcvm, squared = FALSE){
  mat <- matrix (nrow = nrow(data), ncol = nrow(data), 0)
  for(i in 1:nrow(mat)){
    for(j in i:nrow(mat)){
      mat[i,j] <- maha(data[i,], data[j,], vcvm, squared = squared)
    }
  }
  mat <- mat + t(mat)
  rownames(mat) <- colnames(mat) <-rownames(data)
  mat
}

for(i in 0:55){
  PC_thresh <- 0.5 + i/100
  PC_Scores <- t(eigensystem$vectors %*% t(mean_of_means))
  est_means_PCs <- sapply(1:which.max(cumsum(eigensystem$values)), function(ev) PC_Scores[,ev] / eigensystem$values[ev]^0.5)
  est_means_PCs <- est_means_PCs[,1:sum(cumsum(eigensystem$values / sum(eigensystem$values)) < PC_thresh)]
  mahDistMat <- mahaMatrix(est_means_PCs, diag(ncol(est_means_PCs)), squared = F)
  njTree <- ape::nj(mahDistMat)
  njTree <- phytools::reroot(tree = njTree, node.number = which(njTree$tip.label == "NEAND"), position = njTree$edge.length[which(njTree$edge[,2] == which(njTree$tip.label == "NEAND"))] / 2)
  njTree <- phytools::tipRotate(njTree, 1:length(njTree))
  plot(njTree, main = paste0(PC_thresh * 100, "% of the variance represented"))
}

#test out eigenvalue correction?
es <- eigen(mean_of_cors)
es$vectors %*% diag(es$values) %*% t(es$vectors) - mean_of_cors
o_ev <- 84
n_es <- es
n_es$values[n_es$values < 0] <- n_es$values[max(which(n_es$values > 0))]
det(n_es$vectors %*% diag(n_es$values) %*% t(n_es$vectors))
n_es$vectors %*% diag(n_es$values) %*% t(n_es$vectors) - mean_of_cors

det(mean_of_cors)
det(es$vectors[,1:n_ev] %*% diag(es$values[1:n_ev]) %*% t(es$vectors[,1:n_ev]))
es$vectors[,1:n_ev] %*% diag(es$values[1:n_ev]) %*% t(es$vectors[,1:n_ev]) - mean_of_cors


nearPD_Cor <- as.matrix(Matrix::nearPD(mean_of_cors, corr = T)$mat)
det(nearPD_Cor)
write.table(nearPD_Cor, file = "output/nearPD_Corr.txt")
write.table(mean_of_means, file = "output/mean_of_means.txt")

mahDistMat <- mahaMatrix(mean_of_means, as.matrix(Matrix::nearPD(mean_of_cors, corr = T)$mat), squared = F)
njTree <- ape::nj(mahDistMat)
njTree <- phytools::reroot(tree = njTree, node.number = which(njTree$tip.label == "NEAND"), position = njTree$edge.length[which(njTree$edge[,2] == which(njTree$tip.label == "NEAND"))] / 2)
plot(njTree)

nearPD_Cor <- as.matrix(Matrix::nearPD(mean_of_cors, corr = F)$mat)
eigensystem <- eigen(nearPD_Cor)
plot(cumsum(eigensystem$values))
# eigensystem <- eigen(mean_of_cors)
PC_thresh <- 0.99; usePCthresh = F
n_ev <- 50; useNEV = T
PC_Scores <- t(eigensystem$vectors %*% t(mean_of_means))
est_means_PCs <- sapply(1:which.max(cumsum(eigensystem$values)), function(ev) PC_Scores[,ev] / eigensystem$values[ev]^0.5)
plot(apply(est_means_PCs, 2, var)[1:n_ev])
if(usePCthresh){
  est_means_PCs <- est_means_PCs[,1:sum(cumsum(eigensystem$values / sum(eigensystem$values)) < PC_thresh)]
  mahDistMat <- mahaMatrix(est_means_PCs, diag(ncol(est_means_PCs)), squared = F)
} else if(useNEV){
  est_means_PCs <- est_means_PCs[,1:n_ev]
  mahDistMat <- mahaMatrix(est_means_PCs, diag(ncol(est_means_PCs)), squared = F)
}
njTree <- ape::nj(mahDistMat)
njTree <- phytools::reroot(tree = njTree, node.number = which(njTree$tip.label == "NEAND"), position = njTree$edge.length[which(njTree$edge[,2] == which(njTree$tip.label == "NEAND"))] / 2)
njTree <- phytools::tipRotate(njTree, 1:length(njTree))
plot(njTree, main = paste0(PC_thresh * 100, "% of the variance represented"))
plot(phangorn::upgma(mahDistMat))

meansNexus <- function (characterMatrix, outputFilePath, traitsOnColumns = T){
  if (traitsOnColumns) {ntax = dim(characterMatrix)[1]; nchar = dim(characterMatrix)[2]}
  else {ntax = dim(characterMatrix)[2]; nchar = dim(characterMatrix)[1]}
  sink(outputFilePath, append=F)
  cat("#NEXUS\n\nBegin data;\n\tDimensions ntax=", ntax, "nchar=", nchar, ";\n\tFormat datatype=Continuous missing=?;\n\tMatrix\n\n")
  for(i in 1:length(characterMatrix[,1])) {cat(rownames(characterMatrix)[i], characterMatrix[i,], "\n")}
  cat(";\nEnd;")
  sink()
}

meansNexus(est_means_PCs, outputFilePath = "/Volumes/1TB/Bailey/data/bailey_means_nearPD_STAR_nearPDnonCorr_50PCs.nex")


nV <- eigensystem$vectors
nL <- eigensystem$values
nL[nL<0] <- 0
nR <- cov2cor(nV %*% diag(nL) %*% t(nV))
n_eigensystem <- eigen(nR)
nV <- n_eigensystem$vectors
nL <- n_eigensystem$values
plot(nV[,eigensystem$values > 0], eigensystem$vectors[,eigensystem$values > 0])
plot(nL, eigensystem$values, type = "l"); abline(0,1,col = 2)
plot(as.vector(nR), as.vector(cov2cor(as.matrix(Matrix::nearPD(mean_of_cors)$mat))), xlab = "my extension to the eigenvalue extension method"); abline(0,1)
title(main = "comparing elements of outputs from the two pseudocorrelation correction methods")
legend(x = "bottomright", lwd = 1, legend = "1-to-1 line")
R2 <- cor(as.vector(nR), as.vector(cov2cor(as.matrix(Matrix::nearPD(mean_of_cors)$mat))))^2
legend(x = "topleft", legend = paste0("R^2 = ", R2))

meansNexus(est_means_PCs, outputFilePath = "~/data/bailey_means_PCA99.nex")
