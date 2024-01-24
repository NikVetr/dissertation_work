library(magick)
library(purrr)

setwd("~/Documents/mv_threshold_pics/")
file.remove(list.files(pattern=".png"))
### generate single plot ###


library(mvtnorm)
mu <- c(.75,.85)                         
Sigma <- matrix(c(1, .4, .4, 1), 2)  

incr <- 0.01
x <- seq(from = -3, to = 5, by = incr)
y <- seq(from = -3, to = 5, by = incr)
z <- sapply(1:length(x), function(h) dmvnorm(cbind(rep(x[h], length(y)),y), mean = mu, sigma = Sigma))
bivn <- list(x = x, y = y, z = z)

nts_x <- 3 
ts_x <- runif(nts_x-1, 0.6,1.5) #locations of thresholds
ts_x <- cumsum(c(0,ts_x))

nts_y <- 3 
ts_y <- runif(nts_y-1, 0.6,1.5) #locations of thresholds
ts_y <- cumsum(c(0,ts_y))

png(filename = "bivN.png", width = 525, height = 575)

plot(1,1,xlim = c(-3,5), ylim = c(-3,5), col = "white", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
axis(1, -3:5)
axis(2, -3:5)
abline(h = ts_y, lty = 2)
abline(v = ts_x, lty = 2)

text(paste0("state ", 0:(nts_x)), x = sapply(2:5, function(t) (c(-3, ts_x, 5)[t] + c(-3, ts_x, 5)[t-1])/2), y = 5.5, xpd = T, cex = 1.2)
text(paste0("state ", 0:(nts_y)), x = 5.5, y = sapply(2:5, function(t) (c(-3, ts_y, 5)[t] + c(-3, ts_y, 5)[t-1])/2), xpd = T, cex = 1.2, srt=270)
text("state (0,0)", x = -2.5, y = -3)
text(paste0("state (", nts_x, ",", nts_y, ")"), x = 4.5, y = 5)

title(xlab = "TRAIT 1", ylab = "TRAIT 2", line = 2.5, cex.lab = 1.25, font.lab = 2)


clevels <- 2^(-10:0)
contour(bivn, add = T, levels = clevels, labels = paste0("2^", -10:0), col = "black", lwd = 2)

dev.off()

png(filename = "heat.png", width = 525, height = 575)
par(mfrow = c(1,1), mar = c(1,1,8,5), xpd = F)

probs <- matrix(0, nrow = nts_x+1, ncol = nts_y+1)
total_prob <- sum(bivn$z)
tsx <- c(-Inf, ts_x, Inf)
tsy <- c(-Inf, ts_y, Inf)
for(i in 1:(nts_x+1)){
  for(j in 1:(nts_y+1)){
    probs[i,j] <- sum(z[x > tsx[i] & x < tsx[i+1], y > tsy[j] & y < tsy[j+1]]) / total_prob
  }
}

row_i <- as.vector(matrix(1:nrow(probs), nrow = nrow(probs), ncol = ncol(probs)))
col_j <- as.vector(matrix(1:ncol(probs), nrow = nrow(probs), ncol = ncol(probs), byrow = T))
maxProb <- 0.75
cols <- viridisLite::plasma(maxProb * 100)
plot(x = 0:(dim(probs)[1]+1), y = 0:(dim(probs)[2]+1), col = "white", xlab = "", ylab = "", bty = "n", axes = F)
points(x = row_i, y = col_j, col = cols[ceiling(as.vector((probs)) * 100)], cex = 14, pch = 15)
text(x = row_i, y = col_j, labels = round(as.vector((probs)) * 100), col = "white", cex = 1.3)


axis(1, 1:dim(probs)[1], labels = paste0("state ", 0: (dim(probs)[1]-1)), tick = F, line = -5, cex.axis = 1.3)
axis(2, 1:dim(probs)[1], labels = paste0("state ", 0: (dim(probs)[1]-1)), tick = F, line = -5, cex.axis = 1.3)
title(xlab = "TRAIT 1", ylab = "TRAIT 2", line = -2.5, cex.lab = 1.25, font.lab = 2)
rect(xleft = dim(probs)[1] + 0.675, xright = dim(probs)[1] + 0.675 + apply(probs, 2, sum), 
     ybottom = 1:(dim(probs)[2]) - 0.375 , ytop = 1:(dim(probs)[2]) + 0.375,
     col = cols[ceiling(apply(probs, 2, sum)*100)], density = 20, border = "black", lwd = 1.5, xpd = T)
rect(xleft = 1:dim(probs)[1] - 0.375, xright = 1:dim(probs)[1] + 0.375, 
     ybottom = (dim(probs)[2]) + 0.65 , ytop = (dim(probs)[2]) + 0.65 + apply(probs, 1, sum),
     col = cols[ceiling(apply(probs, 2, sum)*100)], density = 20, border = "black", lwd = 1.5, xpd = T)
# barplot(height = apply(probs, 1, sum), horiz = T, axes = F, 
#         add = T, xlim = c(1,(dim(probs)[1])), ylim = c(2,(dim(probs)[2])),
#         width = 0.75, xpd = T, space = 0.25, offset = 5)

plotrix::color.legend( xl = 0.55, xr = dim(probs)[1] + 0.45, yb = (dim(probs)[2]) * 1.15, yt = (dim(probs)[2]) * 1.125 , 
                       gradient="x", rect.col=cols, legend = "")
plotrix::color.legend( xl = (dim(probs)[1]) * 1.155, xr = (dim(probs)[1]) * 1.1255, yb = 0.55, yt = dim(probs)[2] + 0.45 , 
                       gradient="y", rect.col=cols, legend = "")
text(x = 0.45, y = (dim(probs)[2]) * 1.14, labels = 0)
text(x = (dim(probs)[1]) * 1.14, y = 0.45, labels = 0)
text(x = dim(probs)[1] + 0.575, y = (dim(probs)[2]) * 1.14, labels = maxProb * 100)

dev.off()

file.remove(list.files(pattern=".png"))

###########################
### animate moving mean ###
###########################

t_maxProb <- 0.50
prop_to_stretch <- 0.2
cols <- viridisLite::plasma(t_maxProb * 100)
cols <- c(cols[1:(t_maxProb*100-prop_to_stretch*t_maxProb*100)], 
          rep(cols[(t_maxProb*100-prop_to_stretch*t_maxProb*100+1):(t_maxProb*100)], each = (100-t_maxProb*100) / (prop_to_stretch*t_maxProb*100)+1))
maxProb <- 1

t_loop <- 1
rad <- 1
unit <- 0.025
xmu <- seq(-rad, rad, by = unit)
ymu <- sqrt(rad^2 - xmu^2)
xmu <- c(xmu, rev(xmu))
ymu <- c(ymu, -rev(ymu))
mu_disp <- cbind(xmu, ymu)
mu_disp <- rbind(mu_disp[ceiling(length(xmu)*0.875):length(xmu),], mu_disp[1:(ceiling(length(xmu)*0.875)-1),])
mu_disp <- t(t(mu_disp) - mu_disp[1,])
mus <- t(t(mu_disp) + mu)

nts_x <- 3 
ts_x <- runif(nts_x-1, 0.6,1.5) #locations of thresholds
ts_x <- cumsum(c(0,ts_x))

nts_y <- 3 
ts_y <- runif(nts_y-1, 0.6,1.5) #locations of thresholds
ts_y <- cumsum(c(0,ts_y))

Sigma <- matrix(c(1, .4, .4, 1), 2)  
incr <- 0.01


for(i in 1:length(xmu)){
  cat(paste0(i, " "))
  mu <- mus[i,]                      

  
  x <- seq(from = -3, to = 5, by = incr)
  y <- seq(from = -3, to = 5, by = incr)
  z <- sapply(1:length(x), function(h) dmvnorm(cbind(rep(x[h], length(y)),y), mean = mu, sigma = Sigma))
  bivn <- list(x = x, y = y, z = z)
  
  png(filename = paste0("bivN_", i, ".png"), width = 525, height = 575)
  
  plot(1,1,xlim = c(-3,5), ylim = c(-3,5), col = "white", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  axis(1, -3:5)
  axis(2, -3:5)
  abline(h = ts_y, lty = 2)
  abline(v = ts_x, lty = 2)
  
  text(paste0("state ", 0:(nts_x)), x = sapply(2:5, function(t) (c(-3, ts_x, 5)[t] + c(-3, ts_x, 5)[t-1])/2), y = 5.5, xpd = T, cex = 1.2)
  text(paste0("state ", 0:(nts_y)), x = 5.5, y = sapply(2:5, function(t) (c(-3, ts_y, 5)[t] + c(-3, ts_y, 5)[t-1])/2), xpd = T, cex = 1.2, srt=270)
  text("state (0,0)", x = -2.5, y = -3)
  text(paste0("state (", nts_x, ",", nts_y, ")"), x = 4.5, y = 5)
  
  title(xlab = "TRAIT 1", ylab = "TRAIT 2", line = 2.5, cex.lab = 1.25, font.lab = 2)
  
  
  clevels <- 2^(-10:0)
  # contour(bivn, add = T, levels = clevels, labels = paste0("2^", -10:0), col = "black", lwd = 2)
  contour(bivn, add = T, levels = clevels, labels = paste0("2^", -10:0), col = "black", lwd = 2, method = "edge")
  
  dev.off()
}


files <- list.files(pattern = "*.png", full.names = T)
files <- files[order(as.integer(sapply(1:length(files), function(x) strsplit(strsplit(files[x], split = "_")[[1]][2], split = ".png")[[1]][1])))]

files %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=25, loop = t_loop) %>% # animates, can opt for number of loops
  image_write("mv_threshold_mu_smooth_bivN.gif") # write to current dir
# file.rename(from = files[1], to = "hide_image")
file.remove(list.files(pattern=".png"))


for(i in 1:length(xmu)){

  cat(paste0(i, " "))
  mu <- mus[i,]                      
  
  x <- seq(from = -3, to = 5, by = incr)
  y <- seq(from = -3, to = 5, by = incr)
  z <- sapply(1:length(x), function(h) dmvnorm(cbind(rep(x[h], length(y)),y), mean = mu, sigma = Sigma))
  bivn <- list(x = x, y = y, z = z)
    
  png(filename = paste0("heat_", i, ".png"), width = 525, height = 575)
  par(mfrow = c(1,1), mar = c(1,1,8,5), xpd = F)

  probs <- matrix(0, nrow = nts_x+1, ncol = nts_y+1)
  total_prob <- sum(bivn$z)
  tsx <- c(-Inf, ts_x, Inf)
  tsy <- c(-Inf, ts_y, Inf)
  for(i in 1:(nts_x+1)){
    for(j in 1:(nts_y+1)){
      probs[i,j] <- sum(z[x > tsx[i] & x < tsx[i+1], y > tsy[j] & y < tsy[j+1]]) / total_prob
    }
  }
  
  row_i <- as.vector(matrix(1:nrow(probs), nrow = nrow(probs), ncol = ncol(probs)))
  col_j <- as.vector(matrix(1:ncol(probs), nrow = nrow(probs), ncol = ncol(probs), byrow = T))
  plot(x = 0:(dim(probs)[1]+1), y = 0:(dim(probs)[2]+1), col = "white", xlab = "", ylab = "", bty = "n", axes = F)
  points(x = row_i, y = col_j, col = cols[ceiling(as.vector((probs)) * 100)], cex = 14, pch = 15)
  text(x = row_i, y = col_j, labels = round(as.vector((probs)) * 100), col = "white", cex = 1.3)
  
  
  axis(1, 1:dim(probs)[1], labels = paste0("state ", 0: (dim(probs)[1]-1)), tick = F, line = -5, cex.axis = 1.3)
  axis(2, 1:dim(probs)[1], labels = paste0("state ", 0: (dim(probs)[1]-1)), tick = F, line = -5, cex.axis = 1.3)
  title(xlab = "TRAIT 1", ylab = "TRAIT 2", line = -2.5, cex.lab = 1.25, font.lab = 2)
  rect(xleft = dim(probs)[1] + 0.675, xright = dim(probs)[1] + 0.675 + apply(probs, 2, sum), 
       ybottom = 1:(dim(probs)[2]) - 0.375 , ytop = 1:(dim(probs)[2]) + 0.375,
       col = cols[ceiling(apply(probs, 2, sum)*100)], density = 20, border = "black", lwd = 1.5, xpd = T)
  rect(xleft = 1:dim(probs)[1] - 0.375, xright = 1:dim(probs)[1] + 0.375, 
       ybottom = (dim(probs)[2]) + 0.65 , ytop = (dim(probs)[2]) + 0.65 + apply(probs, 1, sum),
       col = cols[ceiling(apply(probs, 2, sum)*100)], density = 20, border = "black", lwd = 1.5, xpd = T)
  # barplot(height = apply(probs, 1, sum), horiz = T, axes = F, 
  #         add = T, xlim = c(1,(dim(probs)[1])), ylim = c(2,(dim(probs)[2])),
  #         width = 0.75, xpd = T, space = 0.25, offset = 5)
  
  plotrix::color.legend( xl = 0.55, xr = dim(probs)[1] + 0.45, yb = (dim(probs)[2]) * 1.15, yt = (dim(probs)[2]) * 1.125 , 
                         gradient="x", rect.col=cols, legend = "")
  plotrix::color.legend( xl = (dim(probs)[1]) * 1.155, xr = (dim(probs)[1]) * 1.1255, yb = 0.55, yt = dim(probs)[2] + 0.45 , 
                         gradient="y", rect.col=cols, legend = "")
  text(x = 0.45, y = (dim(probs)[2]) * 1.14, labels = 0)
  text(x = (dim(probs)[1]) * 1.14, y = 0.45, labels = 0)
  text(x = dim(probs)[1] + 0.585, y = (dim(probs)[2]) * 1.14, labels = maxProb * 100)
  
  dev.off()
  
}

files <- list.files(pattern = "*.png", full.names = T)
files <- files[order(as.integer(sapply(1:length(files), function(x) strsplit(strsplit(files[x], split = "_")[[1]][2], split = ".png")[[1]][1])))]

files %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=25, loop = t_loop) %>% # animates, can opt for number of loops
  image_write("mv_threshold_mu_smooth_heat.gif") # write to current dir
# file.rename(from = files[1], to = "hide_image")
file.remove(list.files(pattern=".png"))

###########################
### animate moving thresholds ###
###########################

t_loop <- 1

ts_x
thresh_to_vary_x <- ceiling(median(1:nts_x))
amount_per_step_x <- 0.025
ts_xs <- c(seq(from = ts_x[thresh_to_vary_x], to = ts_x[thresh_to_vary_x-1]+amount_per_step_x, by = -amount_per_step_x), 
             seq(from = ts_x[thresh_to_vary_x-1]+amount_per_step_x, to = ts_x[thresh_to_vary_x+1]-amount_per_step_x, by = amount_per_step_x),
             seq(from = ts_x[thresh_to_vary_x+1]-amount_per_step_x, to = ts_x[thresh_to_vary_x], by = -amount_per_step_x))

ts_y
thresh_to_vary_y <- ceiling(median(1:nts_y))
amount_per_step_y <- 2 * (ts_y[thresh_to_vary_y+1] - ts_y[thresh_to_vary_y-1] - 4 * amount_per_step_x) / length(ts_xs)
ts_ys <- c(seq(from = ts_y[thresh_to_vary_y], to = ts_y[thresh_to_vary_y-1]+amount_per_step_y, by = -amount_per_step_y), 
           seq(from = ts_y[thresh_to_vary_y-1]+amount_per_step_y, to = ts_y[thresh_to_vary_y+1]-amount_per_step_y, by = amount_per_step_y),
           seq(from = ts_y[thresh_to_vary_y+1]-amount_per_step_y, to = ts_y[thresh_to_vary_y], by = -amount_per_step_y))
ts_ys <- ts_ys[1:length(ts_xs)]

for(i in 1:length(ts_xs)){
  cat(paste0(i, " "))
  
  ts_x[thresh_to_vary_x] <- ts_xs[i]
  ts_y[thresh_to_vary_y] <- ts_ys[i]
  
  x <- seq(from = -3, to = 5, by = incr)
  y <- seq(from = -3, to = 5, by = incr)
  z <- sapply(1:length(x), function(h) dmvnorm(cbind(rep(x[h], length(y)),y), mean = mu, sigma = Sigma))
  bivn <- list(x = x, y = y, z = z)
  
  png(filename = paste0("bivN_", i, ".png"), width = 525, height = 575)
  
  plot(1,1,xlim = c(-3,5), ylim = c(-3,5), col = "white", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  axis(1, -3:5)
  axis(2, -3:5)
  abline(h = ts_y, lty = 2)
  abline(v = ts_x, lty = 2)
  
  text(paste0("state ", 0:(nts_x)), x = sapply(2:5, function(t) (c(-3, ts_x, 5)[t] + c(-3, ts_x, 5)[t-1])/2), y = 5.5, xpd = T, cex = 1.2)
  text(paste0("state ", 0:(nts_y)), x = 5.5, y = sapply(2:5, function(t) (c(-3, ts_y, 5)[t] + c(-3, ts_y, 5)[t-1])/2), xpd = T, cex = 1.2, srt=270)
  text("state (0,0)", x = -2.5, y = -3)
  text(paste0("state (", nts_x, ",", nts_y, ")"), x = 4.5, y = 5)
  
  title(xlab = "TRAIT 1", ylab = "TRAIT 2", line = 2.5, cex.lab = 1.25, font.lab = 2)
  
  
  clevels <- 2^(-10:0)
  # contour(bivn, add = T, levels = clevels, labels = paste0("2^", -10:0), col = "black", lwd = 2)
  contour(bivn, add = T, levels = clevels, labels = paste0("2^", -10:0), col = "black", lwd = 2, method = "edge")
  
  dev.off()
}


files <- list.files(pattern = "*.png", full.names = T)
files <- files[order(as.integer(sapply(1:length(files), function(x) strsplit(strsplit(files[x], split = "_")[[1]][2], split = ".png")[[1]][1])))]

files %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=25, loop = t_loop) %>% # animates, can opt for number of loops
  image_write("mv_threshold_thresh_smooth_bivN.gif") # write to current dir
# file.rename(from = files[1], to = "hide_image")
file.remove(list.files(pattern=".png"))

for(i in 1:length(ts_xs)){
  cat(paste0(i, " "))
  
  ts_x[thresh_to_vary_x] <- ts_xs[i]
  ts_y[thresh_to_vary_y] <- ts_ys[i]           
  
  x <- seq(from = -3, to = 5, by = incr)
  y <- seq(from = -3, to = 5, by = incr)
  z <- sapply(1:length(x), function(h) dmvnorm(cbind(rep(x[h], length(y)),y), mean = mu, sigma = Sigma))
  bivn <- list(x = x, y = y, z = z)
  
  png(filename = paste0("heat_", i, ".png"), width = 525, height = 575)
  par(mfrow = c(1,1), mar = c(1,1,8,5), xpd = F)
  
  probs <- matrix(0, nrow = nts_x+1, ncol = nts_y+1)
  total_prob <- sum(bivn$z)
  tsx <- c(-Inf, ts_x, Inf)
  tsy <- c(-Inf, ts_y, Inf)
  for(i in 1:(nts_x+1)){
    for(j in 1:(nts_y+1)){
      probs[i,j] <- sum(z[x > tsx[i] & x < tsx[i+1], y > tsy[j] & y < tsy[j+1]]) / total_prob
    }
  }
  
  row_i <- as.vector(matrix(1:nrow(probs), nrow = nrow(probs), ncol = ncol(probs)))
  col_j <- as.vector(matrix(1:ncol(probs), nrow = nrow(probs), ncol = ncol(probs), byrow = T))
  plot(x = 0:(dim(probs)[1]+1), y = 0:(dim(probs)[2]+1), col = "white", xlab = "", ylab = "", bty = "n", axes = F)
  points(x = row_i, y = col_j, col = cols[ceiling(as.vector((probs)) * 100)], cex = 14, pch = 15)
  text(x = row_i, y = col_j, labels = round(as.vector((probs)) * 100), col = "white", cex = 1.3)
  
  
  axis(1, 1:dim(probs)[1], labels = paste0("state ", 0: (dim(probs)[1]-1)), tick = F, line = -5, cex.axis = 1.3)
  axis(2, 1:dim(probs)[1], labels = paste0("state ", 0: (dim(probs)[1]-1)), tick = F, line = -5, cex.axis = 1.3)
  title(xlab = "TRAIT 1", ylab = "TRAIT 2", line = -2.5, cex.lab = 1.25, font.lab = 2)
  rect(xleft = dim(probs)[1] + 0.675, xright = dim(probs)[1] + 0.675 + apply(probs, 2, sum), 
       ybottom = 1:(dim(probs)[2]) - 0.375 , ytop = 1:(dim(probs)[2]) + 0.375,
       col = cols[ceiling(apply(probs, 2, sum)*100)], density = 20, border = "black", lwd = 1.5, xpd = T)
  rect(xleft = 1:dim(probs)[1] - 0.375, xright = 1:dim(probs)[1] + 0.375, 
       ybottom = (dim(probs)[2]) + 0.65 , ytop = (dim(probs)[2]) + 0.65 + apply(probs, 1, sum),
       col = cols[ceiling(apply(probs, 2, sum)*100)], density = 20, border = "black", lwd = 1.5, xpd = T)
  # barplot(height = apply(probs, 1, sum), horiz = T, axes = F, 
  #         add = T, xlim = c(1,(dim(probs)[1])), ylim = c(2,(dim(probs)[2])),
  #         width = 0.75, xpd = T, space = 0.25, offset = 5)
  
  plotrix::color.legend( xl = 0.55, xr = dim(probs)[1] + 0.45, yb = (dim(probs)[2]) * 1.15, yt = (dim(probs)[2]) * 1.125 , 
                         gradient="x", rect.col=cols, legend = "")
  plotrix::color.legend( xl = (dim(probs)[1]) * 1.155, xr = (dim(probs)[1]) * 1.1255, yb = 0.55, yt = dim(probs)[2] + 0.45 , 
                         gradient="y", rect.col=cols, legend = "")
  text(x = 0.45, y = (dim(probs)[2]) * 1.14, labels = 0)
  text(x = (dim(probs)[1]) * 1.14, y = 0.45, labels = 0)
  text(x = dim(probs)[1] + 0.585, y = (dim(probs)[2]) * 1.14, labels = maxProb * 100)
  
  dev.off()
  
}

files <- list.files(pattern = "*.png", full.names = T)
files <- files[order(as.integer(sapply(1:length(files), function(x) strsplit(strsplit(files[x], split = "_")[[1]][2], split = ".png")[[1]][1])))]

files %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=25, loop = t_loop) %>% # animates, can opt for number of loops
  image_write("mv_threshold_thresh_smooth_heat.gif") # write to current dir
# file.rename(from = files[1], to = "hide_image")
file.remove(list.files(pattern=".png"))

###########################
### animate correlation ###
###########################

t_loop <- 1

incr_r <- 0.025
rs <- Sigma[1,2]
rs <- c(seq(rs, 0.9, by = incr_r), seq(0.9, -0.9, by = -incr_r), seq(-0.9, rs, by = incr_r))


for(i in 1:length(rs)){
  cat(paste0(i, " "))
  
  Sigma[1,2] <- Sigma[2,1] <- rs[i]
  
  x <- seq(from = -3, to = 5, by = incr)
  y <- seq(from = -3, to = 5, by = incr)
  z <- sapply(1:length(x), function(h) dmvnorm(cbind(rep(x[h], length(y)),y), mean = mu, sigma = Sigma))
  bivn <- list(x = x, y = y, z = z)
  
  png(filename = paste0("bivN_", i, ".png"), width = 525, height = 575)
  
  plot(1,1,xlim = c(-3,5), ylim = c(-3,5), col = "white", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  axis(1, -3:5)
  axis(2, -3:5)
  abline(h = ts_y, lty = 2)
  abline(v = ts_x, lty = 2)
  
  text(paste0("state ", 0:(nts_x)), x = sapply(2:5, function(t) (c(-3, ts_x, 5)[t] + c(-3, ts_x, 5)[t-1])/2), y = 5.5, xpd = T, cex = 1.2)
  text(paste0("state ", 0:(nts_y)), x = 5.5, y = sapply(2:5, function(t) (c(-3, ts_y, 5)[t] + c(-3, ts_y, 5)[t-1])/2), xpd = T, cex = 1.2, srt=270)
  text("state (0,0)", x = -2.5, y = -3)
  text(paste0("state (", nts_x, ",", nts_y, ")"), x = 4.5, y = 5)
  
  title(xlab = "TRAIT 1", ylab = "TRAIT 2", line = 2.5, cex.lab = 1.25, font.lab = 2)
  
  
  clevels <- 2^(-10:0)
  # contour(bivn, add = T, levels = clevels, labels = paste0("2^", -10:0), col = "black", lwd = 2)
  contour(bivn, add = T, levels = clevels, labels = paste0("2^", -10:0), col = "black", lwd = 2, method = "edge")
  
  dev.off()
}


files <- list.files(pattern = "*.png", full.names = T)
files <- files[order(as.integer(sapply(1:length(files), function(x) strsplit(strsplit(files[x], split = "_")[[1]][2], split = ".png")[[1]][1])))]

files %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=25, loop = t_loop) %>% # animates, can opt for number of loops
  image_write("mv_threshold_corr_smooth_bivN.gif") # write to current dir
# file.rename(from = files[1], to = "hide_image")
file.remove(list.files(pattern=".png"))


for(i in 1:length(rs)){
  cat(paste0(i, " "))
  
  Sigma[1,2] <- Sigma[2,1] <- rs[i]
  
  x <- seq(from = -3, to = 5, by = incr)
  y <- seq(from = -3, to = 5, by = incr)
  z <- sapply(1:length(x), function(h) dmvnorm(cbind(rep(x[h], length(y)),y), mean = mu, sigma = Sigma))
  bivn <- list(x = x, y = y, z = z)
  
  png(filename = paste0("heat_", i, ".png"), width = 525, height = 575)
  par(mfrow = c(1,1), mar = c(1,1,8,5), xpd = F)
  
  probs <- matrix(0, nrow = nts_x+1, ncol = nts_y+1)
  total_prob <- sum(bivn$z)
  tsx <- c(-Inf, ts_x, Inf)
  tsy <- c(-Inf, ts_y, Inf)
  for(i in 1:(nts_x+1)){
    for(j in 1:(nts_y+1)){
      probs[i,j] <- sum(z[x > tsx[i] & x < tsx[i+1], y > tsy[j] & y < tsy[j+1]]) / total_prob
    }
  }
  
  row_i <- as.vector(matrix(1:nrow(probs), nrow = nrow(probs), ncol = ncol(probs)))
  col_j <- as.vector(matrix(1:ncol(probs), nrow = nrow(probs), ncol = ncol(probs), byrow = T))
  plot(x = 0:(dim(probs)[1]+1), y = 0:(dim(probs)[2]+1), col = "white", xlab = "", ylab = "", bty = "n", axes = F)
  points(x = row_i, y = col_j, col = cols[ceiling(as.vector((probs)) * 100)], cex = 14, pch = 15)
  text(x = row_i, y = col_j, labels = round(as.vector((probs)) * 100), col = "white", cex = 1.3)
  
  
  axis(1, 1:dim(probs)[1], labels = paste0("state ", 0: (dim(probs)[1]-1)), tick = F, line = -5, cex.axis = 1.3)
  axis(2, 1:dim(probs)[1], labels = paste0("state ", 0: (dim(probs)[1]-1)), tick = F, line = -5, cex.axis = 1.3)
  title(xlab = "TRAIT 1", ylab = "TRAIT 2", line = -2.5, cex.lab = 1.25, font.lab = 2)
  rect(xleft = dim(probs)[1] + 0.675, xright = dim(probs)[1] + 0.675 + apply(probs, 2, sum), 
       ybottom = 1:(dim(probs)[2]) - 0.375 , ytop = 1:(dim(probs)[2]) + 0.375,
       col = cols[ceiling(apply(probs, 2, sum)*100)], density = 20, border = "black", lwd = 1.5, xpd = T)
  rect(xleft = 1:dim(probs)[1] - 0.375, xright = 1:dim(probs)[1] + 0.375, 
       ybottom = (dim(probs)[2]) + 0.65 , ytop = (dim(probs)[2]) + 0.65 + apply(probs, 1, sum),
       col = cols[ceiling(apply(probs, 2, sum)*100)], density = 20, border = "black", lwd = 1.5, xpd = T)
  # barplot(height = apply(probs, 1, sum), horiz = T, axes = F, 
  #         add = T, xlim = c(1,(dim(probs)[1])), ylim = c(2,(dim(probs)[2])),
  #         width = 0.75, xpd = T, space = 0.25, offset = 5)
  
  plotrix::color.legend( xl = 0.55, xr = dim(probs)[1] + 0.45, yb = (dim(probs)[2]) * 1.15, yt = (dim(probs)[2]) * 1.125 , 
                         gradient="x", rect.col=cols, legend = "")
  plotrix::color.legend( xl = (dim(probs)[1]) * 1.155, xr = (dim(probs)[1]) * 1.1255, yb = 0.55, yt = dim(probs)[2] + 0.45 , 
                         gradient="y", rect.col=cols, legend = "")
  text(x = 0.45, y = (dim(probs)[2]) * 1.14, labels = 0)
  text(x = (dim(probs)[1]) * 1.14, y = 0.45, labels = 0)
  text(x = dim(probs)[1] + 0.585, y = (dim(probs)[2]) * 1.14, labels = maxProb * 100)
  
  dev.off()
  
}

files <- list.files(pattern = "*.png", full.names = T)
files <- files[order(as.integer(sapply(1:length(files), function(x) strsplit(strsplit(files[x], split = "_")[[1]][2], split = ".png")[[1]][1])))]

files %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=25, loop = t_loop) %>% # animates, can opt for number of loops
  image_write("mv_threshold_corr_smooth_heat.gif") # write to current dir
# file.rename(from = files[1], to = "hide_image")
file.remove(list.files(pattern=".png"))

###########################
### animate wiggling mean ###
###########################

t_loop <- 1

n_steps <- 100
mu_disp <- rmvnorm(n_steps-1, sigma = Sigma/n_steps, mean = c(0,0))
mu_disp <- rbind(c(0,0), cbind(cumsum(mu_disp[,1]), cumsum(mu_disp[,2])))
mus <- t(t(mu_disp) + mu)


for(i in 1:n_steps){
  cat(paste0(i, " "))
  mu <- mus[i,]                      
  
  
  x <- seq(from = -3, to = 5, by = incr)
  y <- seq(from = -3, to = 5, by = incr)
  z <- sapply(1:length(x), function(h) dmvnorm(cbind(rep(x[h], length(y)),y), mean = mu, sigma = Sigma))
  bivn <- list(x = x, y = y, z = z)
  
  png(filename = paste0("bivN_", i, ".png"), width = 525, height = 575)
  
  plot(1,1,xlim = c(-3,5), ylim = c(-3,5), col = "white", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  axis(1, -3:5)
  axis(2, -3:5)
  abline(h = ts_y, lty = 2)
  abline(v = ts_x, lty = 2)
  
  text(paste0("state ", 0:(nts_x)), x = sapply(2:5, function(t) (c(-3, ts_x, 5)[t] + c(-3, ts_x, 5)[t-1])/2), y = 5.5, xpd = T, cex = 1.2)
  text(paste0("state ", 0:(nts_y)), x = 5.5, y = sapply(2:5, function(t) (c(-3, ts_y, 5)[t] + c(-3, ts_y, 5)[t-1])/2), xpd = T, cex = 1.2, srt=270)
  text("state (0,0)", x = -2.5, y = -3)
  text(paste0("state (", nts_x, ",", nts_y, ")"), x = 4.5, y = 5)
  
  title(xlab = "TRAIT 1", ylab = "TRAIT 2", line = 2.5, cex.lab = 1.25, font.lab = 2)
  
  
  clevels <- 2^(-10:0)
  # contour(bivn, add = T, levels = clevels, labels = paste0("2^", -10:0), col = "black", lwd = 2)
  contour(bivn, add = T, levels = clevels, labels = paste0("2^", -10:0), col = "black", lwd = 2, method = "edge")
  
  dev.off()
}


files <- list.files(pattern = "*.png", full.names = T)
files <- files[order(as.integer(sapply(1:length(files), function(x) strsplit(strsplit(files[x], split = "_")[[1]][2], split = ".png")[[1]][1])))]

files %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=25, loop = t_loop) %>% # animates, can opt for number of loops
  image_write("mv_threshold_mu_bm_bivN.gif") # write to current dir
# file.rename(from = files[1], to = "hide_image")
file.remove(list.files(pattern=".png"))


for(i in 1:n_steps){
  
  cat(paste0(i, " "))
  mu <- mus[i,]                      
  
  x <- seq(from = -3, to = 5, by = incr)
  y <- seq(from = -3, to = 5, by = incr)
  z <- sapply(1:length(x), function(h) dmvnorm(cbind(rep(x[h], length(y)),y), mean = mu, sigma = Sigma))
  bivn <- list(x = x, y = y, z = z)
  
  png(filename = paste0("heat_", i, ".png"), width = 525, height = 575)
  par(mfrow = c(1,1), mar = c(1,1,8,5), xpd = F)
  
  probs <- matrix(0, nrow = nts_x+1, ncol = nts_y+1)
  total_prob <- sum(bivn$z)
  tsx <- c(-Inf, ts_x, Inf)
  tsy <- c(-Inf, ts_y, Inf)
  for(i in 1:(nts_x+1)){
    for(j in 1:(nts_y+1)){
      probs[i,j] <- sum(z[x > tsx[i] & x < tsx[i+1], y > tsy[j] & y < tsy[j+1]]) / total_prob
    }
  }
  
  row_i <- as.vector(matrix(1:nrow(probs), nrow = nrow(probs), ncol = ncol(probs)))
  col_j <- as.vector(matrix(1:ncol(probs), nrow = nrow(probs), ncol = ncol(probs), byrow = T))
  plot(x = 0:(dim(probs)[1]+1), y = 0:(dim(probs)[2]+1), col = "white", xlab = "", ylab = "", bty = "n", axes = F)
  points(x = row_i, y = col_j, col = cols[ceiling(as.vector((probs)) * 100)], cex = 14, pch = 15)
  text(x = row_i, y = col_j, labels = round(as.vector((probs)) * 100), col = "white", cex = 1.3)
  
  
  axis(1, 1:dim(probs)[1], labels = paste0("state ", 0: (dim(probs)[1]-1)), tick = F, line = -5, cex.axis = 1.3)
  axis(2, 1:dim(probs)[1], labels = paste0("state ", 0: (dim(probs)[1]-1)), tick = F, line = -5, cex.axis = 1.3)
  title(xlab = "TRAIT 1", ylab = "TRAIT 2", line = -2.5, cex.lab = 1.25, font.lab = 2)
  rect(xleft = dim(probs)[1] + 0.675, xright = dim(probs)[1] + 0.675 + apply(probs, 2, sum), 
       ybottom = 1:(dim(probs)[2]) - 0.375 , ytop = 1:(dim(probs)[2]) + 0.375,
       col = cols[ceiling(apply(probs, 2, sum)*100)], density = 20, border = "black", lwd = 1.5, xpd = T)
  rect(xleft = 1:dim(probs)[1] - 0.375, xright = 1:dim(probs)[1] + 0.375, 
       ybottom = (dim(probs)[2]) + 0.65 , ytop = (dim(probs)[2]) + 0.65 + apply(probs, 1, sum),
       col = cols[ceiling(apply(probs, 2, sum)*100)], density = 20, border = "black", lwd = 1.5, xpd = T)
  # barplot(height = apply(probs, 1, sum), horiz = T, axes = F, 
  #         add = T, xlim = c(1,(dim(probs)[1])), ylim = c(2,(dim(probs)[2])),
  #         width = 0.75, xpd = T, space = 0.25, offset = 5)
  
  plotrix::color.legend( xl = 0.55, xr = dim(probs)[1] + 0.45, yb = (dim(probs)[2]) * 1.15, yt = (dim(probs)[2]) * 1.125 , 
                         gradient="x", rect.col=cols, legend = "")
  plotrix::color.legend( xl = (dim(probs)[1]) * 1.155, xr = (dim(probs)[1]) * 1.1255, yb = 0.55, yt = dim(probs)[2] + 0.45 , 
                         gradient="y", rect.col=cols, legend = "")
  text(x = 0.45, y = (dim(probs)[2]) * 1.14, labels = 0)
  text(x = (dim(probs)[1]) * 1.14, y = 0.45, labels = 0)
  text(x = dim(probs)[1] + 0.585, y = (dim(probs)[2]) * 1.14, labels = maxProb * 100)
  
  dev.off()
  
}

files <- list.files(pattern = "*.png", full.names = T)
files <- files[order(as.integer(sapply(1:length(files), function(x) strsplit(strsplit(files[x], split = "_")[[1]][2], split = ".png")[[1]][1])))]

files %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=25, loop = t_loop) %>% # animates, can opt for number of loops
  image_write("mv_threshold_mu_bm_heat.gif") # write to current dir
# file.rename(from = files[1], to = "hide_image")
file.remove(list.files(pattern=".png"))

cols_wiggles <- viridisLite::cividis(n_steps)
for(i in 1:n_steps){
  
  cat(paste0(i, " "))
  
  png(filename = paste0("wiggle_", i, ".png"), width = 525, height = 575)
  
  plot(mus[1:i,], type = "l", xlim = c(min(mus[,1]), max(mus[,1])), ylim = c(min(mus[,2]), max(mus[,2])), 
       lwd = 2, xlab = "TRAIT 1", ylab = "TRAIT 2", col = "white", cex.lab = 1.425, font.lab = 2, xpd = T)
  for(j in 2:i){
    lines(c(mus[j-1,1], mus[j,1]), c(mus[j-1,2], mus[j,2]), col = cols_wiggles[j], lwd = 2)
  }
  
  
  plotrix::color.legend( xl = min(mus[,1]), xr = max(mus[,1]), yb = max(mus[,2]) + 0.1, yt = max(mus[,2]) + 0.2 , 
                         gradient="x", rect.col=cols_wiggles, legend = "")
  title(main = "TIME", font = 2, cex.main = 2, line = 0.9, xpd = T, col.main = "white")
  
  text(x = min(mus[,1])-0.04, max(mus[,2]) + 0.15, labels = 0, xpd = T, cex = 1.55)
  text(x = max(mus[,1])+0.04, max(mus[,2]) + 0.15, labels = 1, xpd = T, cex = 1.55)
  
  dev.off()
  
}


files <- list.files(pattern = "*.png", full.names = T)
files <- files[order(as.integer(sapply(1:length(files), function(x) strsplit(strsplit(files[x], split = "_")[[1]][2], split = ".png")[[1]][1])))]

files %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=25, loop = t_loop) %>% # animates, can opt for number of loops
  image_write("mv_threshold_mu_bm_wiggle.gif") # write to current dir
# file.rename(from = files[1], to = "hide_image")
file.remove(list.files(pattern=".png"))
