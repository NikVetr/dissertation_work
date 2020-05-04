setwd("~/Documents/1D_BrownianMotion")
file.remove(list.files(pattern=".png"))
library(magick)
library(purrr)

ntips <- 30
cols <- sample(viridisLite::magma(ntips))
lw = 3 #line width
wid <- 1500
hei <- 750
decr_loc <- 0.8
tl <- 10
sigma2 <- 1
breakpoint <- 0.25
increments <- 150

rOU <- function(n, x0 = 0, t = 1, sigma2 = 1, alpha = 0.1, theta = 0){
  mean <- theta + (x0-theta)*exp(-1*alpha*t)
  var <- sigma2 / (2*alpha) * (1-exp(-2*alpha*t))
  rnorm(n = n, mean = mean, sd = sqrt(var))
}


nlines <- 1 #number of realizations of the process to plot
x0 <- 0 #starting state of the process
theta <- 0 #location of the optimum
alpha = 2.75 #strength of the rubber band / returning force
t <- (0:(increments*breakpoint))*tl/increments #discrete time steps over which to simulate the process
x <- matrix(0,nrow = nlines, ncol = length(t)); x[,1] <- x0 #initialize the matrix of realizations
for(j in 1:nlines){
  for(i in 2:length(t)){
    x[j,i] <- rOU(n=1, x0=x[j,i-1], t = t[i]-t[i-1], sigma2=sigma2, alpha=alpha, theta=theta)
  }
}
for(i in 1:(breakpoint*increments)){
  png(filename = paste0("wiggle_", i, ".png"), width = wid, height = hei)
  plot(t[1:i],x[1,1:i], type = "l", ylim = c(-5,5), xlim = c(0,tl), col=rgb(0,0,0,1), lwd = lw, xlab = "", ylab = "", yaxt="n", xaxt="n", bty="n")
  dev.off()
}
anc_t <- t
anc_x <- x
nlines <- ntips #number of realizations of the process to plot
x0 <- x[1,breakpoint*increments+1] #starting state of the process
theta <- 0 #location of the optimum
alpha = 0.15 #strength of the rubber band / returning force
t <- ((increments*breakpoint):increments)*tl/increments #discrete time steps over which to simulate the process
x <- matrix(0,nrow = nlines, ncol = length(t)); x[,1] <- x0 #initialize the matrix of realizations
for(j in 1:nlines){
  for(i in 2:length(t)){
    x[j,i] <- rOU(n=1, x0=x[j,i-1], t = t[i]-t[i-1], sigma2=sigma2, alpha=alpha, theta=theta)
  }
}
for(i in 1:(length(t)-1)){
  png(filename = paste0("wiggle_", increments*breakpoint+i, ".png"), width = wid, height = hei)
  plot(anc_t,anc_x, type = "l", ylim = c(-5,5), xlim = c(0,tl), col=rgb(0,0,0,1), lwd = lw, xlab = "", ylab = "", yaxt="n", xaxt="n", bty="n")
  for(j in 1:(ntips)){
    # wiggler <- ifelse(i == length(t),0,sample(-1:1,1))
    wiggler <- 0
    colc <- col2rgb(cols[j])
    opacity <- 1
    colc <- rgb(colc[1], colc[2], colc[3], opacity*256, maxColorValue = 256)
    if(i < (increments*decr_loc - increments*breakpoint)){    
      lines(t[1:(i+wiggler)],x[j,1:(i+wiggler)], type = "l", ylim = c(-5,5), xlim = c(0,tl), col=colc, lwd = lw, xlab = "", ylab = "", yaxt="n", bty="n")
    } else {
      solid_i <- (increments*decr_loc - increments*breakpoint)
      lines(t[1:(solid_i+wiggler)],x[j,1:(solid_i+wiggler)], type = "l", ylim = c(-5,5), xlim = c(0,tl), col=colc, lwd = lw, xlab = "", ylab = "", yaxt="n", bty="n")
      for(i_fade in solid_i:(i)){
        opacity <- round(1 - (i_fade - (increments*decr_loc - increments*breakpoint))/(increments*(1-decr_loc) + 1), digits = 5)
        colc <- col2rgb(cols[j])
        colc <- rgb(colc[1], colc[2], colc[3], opacity*256, maxColorValue = 256)
        lines(t[i_fade:(i_fade+1+wiggler)],x[j,i_fade:(i_fade+1+wiggler)], type = "l", ylim = c(-5,5), xlim = c(0,tl), col=colc, lwd = lw, xlab = "", ylab = "", yaxt="n", bty="n")
      }
    }
  }
  dev.off()
}
files <- list.files(pattern = "*.png", full.names = T)
files <- files[order(as.integer(sapply(1:length(files), function(x) strsplit(strsplit(files[x], split = "_")[[1]][2], split = ".png")[[1]][1])))]
files <- c(files, rev(files))
files %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=25) %>% # animates, can opt for number of loops
  image_write("pretty_uBM.gif") # write to current dir
file.remove(list.files(pattern=".png"))

png(filename = paste0("tree.png"), width = 400, height = 1000)
plot(0:11, 0:11/11*tl, col = "white", axes = F, xlab = "", ylab = "")
segments(x0 = 5, x1 = 5, y0 = 1, y1 = tl*breakpoint, col = "black", lwd = lw)
for(i in 1:ntips){
  segments(x0 = 5, x1 = i/ntips*10, y0 = tl*breakpoint, y1 = tl, col = cols[i], lwd = lw)
}
segments(x0 = 5, x1 = 5, y0 = 1, y1 = tl*breakpoint, col = "black", lwd = lw)
dev.off()