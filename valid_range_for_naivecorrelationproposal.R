library(corpcor)
dim <- 10
rlkj <- function (K, eta = 1) {
  alpha <- eta + (K - 2)/2
  r12 <- 2 * rbeta(1, alpha, alpha) - 1
  R <- matrix(0, K, K)
  R[1, 1] <- 1
  R[1, 2] <- r12
  R[2, 2] <- sqrt(1 - r12^2)
  if (K > 2) 
    for (m in 2:(K - 1)) {
      alpha <- alpha - 0.5
      y <- rbeta(1, m/2, alpha)
      z <- rnorm(m, 0, 1)
      z <- z/sqrt(crossprod(z)[1])
      R[1:m, m + 1] <- sqrt(y) * z
      R[m + 1, m + 1] <- sqrt(1 - y)
    }
  return(crossprod(R))
}
plotMatrix <- function(mobject, size, location, lwd = 2, grid = T, font = 1, cex = 1, rownames = T, colnames = T, title = T, title.label = "Matrix Object"){
  lines(rbind(location, location + c(0,size[2])), lwd = lwd)
  lines(rbind(location, location + c(size[1]/8,0)), lwd = lwd)
  lines(rbind(location + c(0, size[2]), location + c(size[1]/8,size[2])), lwd = lwd)
  lines(rbind(location + c(size[1],0), location + size), lwd = lwd)
  lines(rbind(location + size, location + size - c(size[1]/8,0)), lwd = lwd)
  lines(rbind(location + c(size[1],0), location + c(size[1],0) - c(size[1]/8,0)), lwd = lwd)
  if(grid == T){
    for(i in 1:(dim(mobject)[1]-1)){
      lines(rbind(location + c(0,i*size[2]/dim(mobject)[1]), location + c(size[1], i*size[2]/dim(mobject)[1])))
    }
    for(j in 1:(dim(mobject)[2]-1)){
      lines(rbind(location + c(j*size[1]/dim(mobject)[2],0), location + c(j*size[1]/dim(mobject)[2], size[2])))
    }
  }
  if(class(mobject[1,1]) != "expression" & class(mobject[1,1]) != "character"){mobject <- matrix(as.character(mobject), nrow = dim(mobject)[1], ncol = dim(mobject)[2])}
  for(i in 1:(dim(mobject)[1])){
    for(j in 1:dim(mobject)[2]){
      text(labels = mobject[i,j], x = location[1] + (j-1/2)*size[1]/dim(mobject)[2], y = location[2] + size[2] - (i-1/2)*size[2]/dim(mobject)[1], font = font, cex = cex)
    }
  }
  if(title){
    text(title.label, x = location[1] + size[1]/2, y = location[2] + size[2] + strheight(title.label, font = 2, cex = 1.5)/1.5, cex = 1.5, font = 2)
  }
  if(rownames){
    for(i in 1:dim(mobject)[1]){
      text(rownames(mobject)[i], x = location[1] - strwidth(rownames(mobject)[i])/2 - size[1]/(ncol(mobject)*6), y = location[2] + size[2] - (i-1/2)*size[2]/dim(mobject)[2])
    }
  }
  if(colnames){
    for(i in 1:dim(mobject)[1]){
      text(colnames(mobject)[i], x = location[1] + (i-1/2)*size[1]/dim(mobject)[1], y = location[2] - strheight(colnames(mobject)[i])/2- size[2]/(nrow(mobject)*6))
    }
  }
}

R <- rlkj(K = dim) 
valid_range <- function(R, inds = c(1,2), incr = 0.01){
  R_test <- R
  up <- 0
  while(matrixcalc::is.positive.definite(R_test)){
    up <- up + incr
    R_test[inds[1], inds[2]] <- R_test[inds[2], inds[1]] <- R[inds[1], inds[2]] + up
  }
  R_test <- R
  down <- 0
  while(matrixcalc::is.positive.definite(R_test)){
    down <- down - incr
    R_test[inds[1], inds[2]] <- R_test[inds[2], inds[1]] <- R[inds[1], inds[2]] + down
  }
  return((up - down) - 2*incr)
}


dims <- round(1.15^(1:45))[-(1:5)]

# save(ranges, file = "ranges_of_expected_wiggleroom.txt")
ranges <- sapply(dims, function(dim) replicate(n = 2000, valid_range(rlkj(K = dim), incr = 0.001)))
load("ranges_of_expected_wiggleroom.txt")
par(mar = c(5,5,0,2))
png(filename = "/Users/nikolai/Documents/Harvati_Reanalysis_Manuscript/figures/expected_wiggle_room.png", width = 800, height = 600)
plot(dims, log10(apply(ranges, 2, mean)), ylim = c(-2.1,0.5),
     type = "l", xlab = "matrix dimensionality", 
     ylab = "expected interval width", cex.lab = 1.25, lwd = 3, bty = "n", col = "darkred", xaxt = "n", yaxt = "n")
axis(side = 2, labels = c(2*2^(-(0:3)), sapply(-(3:9), function(i) as.expression(bquote(2^ .(i))))), at = log10(2*2^(-(0:10))), las = 2)
axis(side = 1, labels = -1:11*50, at = -1:11*50)
dev.off()

#partial correlation matrix test
makeSymmetric <- function(R){
  return((R + t(R)) / 2)
}
valid_range_pc <- function(R, inds = c(1,2), incr = 0.01){
  R_test <- corpcor::cor2pcor(R)
  up <- 0
  yPSD <- T
  while(yPSD){
    up <- up + incr
    R_test[inds[1], inds[2]] <- R_test[inds[2], inds[1]] <- R[inds[1], inds[2]] + up
    m = R_test
    m = -m
    diag(m) = -diag(m)
    m = solve(m)
    if(any(diag(m) < 0)){
      yPSD <- F
    } else{
      yPSD <- matrixcalc::is.positive.definite(makeSymmetric(corpcor::pcor2cor(R_test)))
    }
  }
  R_test <- corpcor::cor2pcor(R)
  down <- 0
  yPSD <- T
  while(yPSD){
    down <- down - incr
    R_test[inds[1], inds[2]] <- R_test[inds[2], inds[1]] <- R[inds[1], inds[2]] + down
    m = R_test
    m = -m
    diag(m) = -diag(m)
    m = solve(m)
    if(any(diag(m) < 0)){
      yPSD <- F
    } else{
      yPSD <- matrixcalc::is.positive.definite(makeSymmetric(corpcor::pcor2cor(R_test)))
    }
  }
  return((up - down) - 2*incr)
}

#not all partial correlations are valid

#generate partial correlation matrix
set.seed(1)
dim <- 5
pcm <- matrix(runif(dim^2, -1, 1), dim, dim)
diag(pcm) <- 1

#transform to correlation matrix
cm = -pcm
diag(cm) <- -diag(cm)
cm <- solve(cm)
is.positive.definite(cov2cor(cm))


R <- rlkj(K = 10)
valid_range_pc(R, incr = 0.001)
ranges <- sapply(dims, function(dim) replicate(n = 2000, valid_range_pc(rlkj(K = dim), incr = 0.001)))
