library(mvMORPH)
library(adephylo)
library(rethinking)
library(parallel)
library(foreach)
library(doParallel)

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

n_traits <- 100
n_tips <- 20

noise_vars <- c(0,0.001,0.01,0.1,0.3,0.5,1,2,5)
exp_cor_evs <- matrix(0, nrow = length(noise_vars), ncol = n_traits)
from_RM_dist <- T

cl <- makeCluster(5, outfile="")
registerDoParallel(cl)
getDoParWorkers()
# stopCluster(cl)

foreach(nv=1:length(noise_vars), .packages = c("mvMORPH", "adephylo", "rethinking")) %dopar% {
# for(nv in 1:length(noise_vars)){
  nrep = 1E6
  cs_ev <- matrix(0, nrow = nrep, ncol = n_traits)
  cor_ev <- matrix(0, nrow = nrep, ncol = n_traits)
  for(i in 1:nrep){
    if(i %% 10 == 0){cat(paste0(i, " "))}
    
    R <- rlkjcorr(1, n_traits)
    tree <- pbtree(n = n_tips)
    traits <- mvSIM(tree, model="BM1", param=list(sigma=R, theta=rep(0, n_traits))) 
    patdists <- as.matrix(distTips(tree))
    
    noise_var <- noise_vars[nv]
    if(!from_RM_dist){
      noise <- matrix(rnorm(n = n_tips * n_traits, mean = 0, sd = sqrt(noise_var)), nrow = n_tips, ncol = n_traits)
    } else {
      noise <- rmvnorm(n = n_tips, mean = rep(0, n_traits), sigma = R * noise_var)
      
    }
    traits <- traits + noise
    
    V <- eigen(R)
    L <- V$values
    V <- V$vectors
    # VarExpl <- round(L / sum(L), 4) * 100
    
    #plot(cumsum(VarExpl), type = "l")
    PC_Scores <- t((t(V) %*% t(traits)))
    pcdists <- sapply(1:n_traits, function(pc) as.matrix(dist(PC_Scores[,pc]))[upper.tri(diag(n_tips))])
    cors_with_patdists <- sapply(1:n_traits, function(pc) cor(pcdists[,pc], patdists[upper.tri(patdists)]))
    #plot(cors_with_patdists)
    
    # alltrmahadists <- mahaMatrix(traits, vcvm = R)
    # alltrEucPCdists <- as.matrix(dist(PC_Scores))
    # cor(alltrmahadists[upper.tri(alltrmahadists)], patdists[upper.tri(patdists)])
    # cor(alltrEucPCdists[upper.tri(alltrEucPCdists)], alltrmahadists[upper.tri(alltrmahadists)])
    # cor(alltrEucPCdists[upper.tri(alltrEucPCdists)], patdists[upper.tri(patdists)])
    
    # cs_ev[i,] <- VarExpl
    cor_ev[i,] <- cors_with_patdists
    
  }
  exp_cor_ev <- apply(cor_ev, 2, mean)
  write.table(exp_cor_ev, file = paste0("~/exp_cor_ev_", ifelse(from_RM_dist, "RMdist_", ""), noise_var, ".txt"))
  # exp_cor_evs[nv,] <- exp_cor_ev
}

exp_cor_evs <- matrix(0, nrow = length(noise_vars), ncol = n_traits)
for(nv in 1:length(noise_vars)){
  noise_var <- noise_vars[nv]
  exp_cor_ev <- c(read.table(file = paste0("~/exp_cor_ev_", noise_var, ".txt")))$x
  exp_cor_evs[nv,] <- exp_cor_ev
}

exp_cor_evs_rm <- matrix(0, nrow = length(noise_vars), ncol = n_traits)
for(nv in 1:length(noise_vars)){
  noise_var <- noise_vars[nv]
  exp_cor_ev <- c(read.table(file = paste0("~/exp_cor_ev_RMdist_", noise_var, ".txt")))$x
  exp_cor_evs_rm[nv,] <- exp_cor_ev
}

egval <- replicate(1E3, eigen(rlkjcorr(1, n_traits))$values)
varexpl <- function(x){x / sum(x) * 100}
egval <- apply(egval, 2, varexpl)
expected_treeheight <- mean(replicate(1E4, max(distTips(pbtree(n = n_tips))) / 2))
sd_treeheight <- sd(replicate(1E4, max(distTips(pbtree(n = n_tips))) / 2))
library(latex2exp)
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


grDevices::cairo_pdf("~/dissertation/figures/noisydecay_of_association.pdf", width = 12, height = 4)
layout(matrix(c(1,1,1,2,2,2,3,3), nrow = 1))

par(mar = c(4.5,5.25,4,1))
cols = viridis::viridis(n = length(noise_vars), end = 0.95)
plot(exp_cor_evs[1,], type = "l", ylim = c(-0.0075,0.3), col = cols[1], lwd = 3, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
box(lwd=2)
box(which = "figure", lwd = 2, lty = 2)
axis(2, at = 0:3 / 10, labels = rep("", 4), lwd = 2, cex.axis = 3, tck = -0.025)
mtext(text = 0:3 / 10, side = 2, at = 0:3 / 10, cex = 1, line = 0.75, las = 2)
axis(1, at = 0:5 * 20, labels =  rep("", 6), lwd = 2, cex.axis = 3, tck = -0.025)
mtext(text = 0:5 * 20, side = 1, at = 0:5 * 20, cex = 1, line = 0.75, las = 2)
title(xlab = "Ordered Eigenvector Index", cex.lab = 1.65)
title(ylab = TeX(paste0("E(r) between d_{patristic} and d_{Euclidean}")), cex.lab = 1.65)
for(i in 2:nrow(exp_cor_evs)){
  lines(exp_cor_evs[i,], lwd = 3, col = cols[i])
}
legend(x = "bottomleft", col = cols, lwd = 3, legend = noise_vars, title = TeX(paste0("Gaussian Noise ($\\sigma^2$)")), cex = 1.15, box.lty = 2, box.lwd = 2, ncol = 2)
title(main = "Expected Correlation Between Distances Along \nSuccessive Eigenvectors with Patristic Distance", cex.main  = 1.5)
fig_label(text = "a)", region = "figure", cex = 2, shrinkY = 0.975)
box(lwd=2)
text(x = 30.5, y = -0.005, TeX("r_{ij} = 0"), cex  = 1.15)
text(x = 77, y = -0.01, TeX(paste0("tree ~ pure-birth($\\lambda = 1$, n_{tips} = ", n_tips, ")")), cex  = 1.15)

###
plot(exp_cor_evs_rm[1,], type = "l", ylim = c(-0.0075,0.3), col = cols[1], lwd = 3, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
box(lwd=2)
box(which = "figure", lwd = 2, lty = 2)
axis(2, at = 0:3 / 10, labels = rep("", 4), lwd = 2, cex.axis = 3, tck = -0.025)
mtext(text = 0:3 / 10, side = 2, at = 0:3 / 10, cex = 1, line = 0.75, las = 2)
axis(1, at = 0:5 * 20, labels =  rep("", 6), lwd = 2, cex.axis = 3, tck = -0.025)
mtext(text = 0:5 * 20, side = 1, at = 0:5 * 20, cex = 1, line = 0.75, las = 2)
title(xlab = "Ordered Eigenvector Index", cex.lab = 1.65)
title(ylab = TeX(paste0("E(r) between d_{patristic} and d_{Euclidean}")), cex.lab = 1.65)
for(i in 2:nrow(exp_cor_evs_rm)){
  lines(exp_cor_evs_rm[i,], lwd = 3, col = cols[i])
}
legend(x = "bottomleft", col = cols, lwd = 3, legend = noise_vars, title = TeX(paste0("Gaussian Noise ($\\sigma^2$)")), cex = 1.15, box.lty = 2, box.lwd = 2, ncol = 2)
title(main = "Expected Correlation Between Distances Along \nSuccessive Eigenvectors with Patristic Distance", cex.main  = 1.5)
fig_label(text = "b)", region = "figure", cex = 2, shrinkY = 0.975)
box(lwd=2)
text(x = 30.5, y = -0.005, TeX("r_{ij} = R_{ij}"), cex  = 1.15)
text(x = 77, y = -0.01, TeX(paste0("tree ~ pure-birth($\\lambda = 1$, n_{tips} = ", n_tips, ")")), cex  = 1.15)


###
par(mar = c(4.5,3.75,4,1))
plot(egval[,1], type = "l", ylim = c(-0,4.5), xlim = c(0,100), col = rgb(0,0,0,0.05), lwd = 3, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
box(lwd=2)
box(which = "figure", lwd = 2, lty = 2)
axis(2, at = 0:4, labels = rep("", 5), lwd = 2, cex.axis = 3, tck = -0.025)
mtext(text = 0:4, side = 2, at = 0:4, cex = 1, line = 0.75, las = 2)
axis(1, at = 0:5 * 20, labels =  rep("", 6), lwd = 2, cex.axis = 3, tck = -0.025)
mtext(text = 0:5 * 20, side = 1, at = 0:5 * 20, cex = 1, line = 0.75, las = 2)
title(xlab = "Ordered Eigenvector Index", cex.lab = 1.65)
title(ylab = TeX(paste0("Percent of Variance in LKJ($\\eta = 1$)")), cex.lab = 1.65, line = 1.75)
for(i in 2:nrow(egval)){
  lines(egval[,i], lwd = 3, col = rgb(0,0,0,0.05))
}
title(main = (paste0("Percent of Variance in Each\n Eigenvector in a Flat LKJ")), cex.main  = 1.5)
fig_label(text = "c)", region = "figure", cex = 2, shrinkY = 0.975)
box(lwd=2)

dev.off()
