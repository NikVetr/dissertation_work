library(magick)
library(purrr)

setwd("~/Documents/univ_threshold_pics/")
file.remove(list.files(pattern=".png"))
### generate single plot ###

par(mfrow = c(1,2), mar = c(4,3,3,2), xpd = F)
x <- (-500:700)/100

mu <- 1 #location of std norm mean
y <- dnorm(x, mean = mu)
plot(x, y, type = "l", lwd = 2, ylim = c(0,0.45), yaxt = "n", ylab = "", xlab = "")
title(xlab = "latent liability", line=2.25)

nts <- 3 #number of thresholds
ts <- runif(nts-1, 0.6,1.5) #locations of thresholds
ts <- cumsum(c(0,ts))

#specify colors
cols <- viridisLite::viridis(nts+1, alpha = 0.8)

#shade underneath curve
for(i in 1:(nts+1)){
  
  if(i == 1){
    xs <- c(x[x <= ts[i]], rev(x[x <= ts[i]]))
    ys <- c(y[x <= ts[i]], rep(0, length(y[x <= ts[i]])))
  }
  
  if(1 < i & i <= nts) {
    xs <- c(x[ts[i-1] <= x & x <= ts[i]], rev(x[ts[i-1] <= x & x <= ts[i]]))
    ys <- c(y[ts[i-1] <= x & x <= ts[i]], rep(0, length(y[ts[i-1] <= x & x <= ts[i]])))
  }
  
  if(i == (nts+1)){
    xs <- c(x[x >= ts[i-1]], rev(x[x >= ts[i-1]]))
    ys <- c(y[x >= ts[i-1]], rep(0, length(y[x >= ts[i-1]])))
  }
  
  polygon(x = xs, y = ys, col = cols[i])
}

abline(v = ts, lwd = 3, lty = 2)


#get frequencies of discrete traits
discr_props <- c(0, pnorm(q = ts, mean = mu, sd = 1), 1)
discr_props <- diff(discr_props)
discr_props
barplot(discr_props, col = cols, ylim = c(0,0.6), xlab = "")
title(xlab = "discrete trait frequency", line=1)

################
### vary mus ###
################

mus <- (8:40)/20; mus <- c(mus, rev(mus)[-1])
nts <- 3 #number of thresholds
ts <- runif(nts-1, 0.6,1.5) #locations of thresholds
ts <- cumsum(c(0,ts))
x <- (-500:700)/100
t_loop <- 1

#specify colors
cols <- viridisLite::viridis(nts+1, alpha = 0.8)

for(i in 1:length(mus)){
  cat(paste0(i, " "))
  mu <- mus[i] #location of std norm mean
  y <- dnorm(x, mean = mu)
  
  png(filename = paste0("univ_threshold_plot_", i, ".png"), width = 1000, height = 350)
  par(mfrow = c(1,2), mar = c(4,4,3,2), xpd = F)
  
  plot(x, y, type = "l", lwd = 2, ylim = c(0,0.45), yaxt = "n", ylab = "", xlab = "")
  title(xlab = "latent liability", line=2.25, cex.lab = 1.5)
  
  #shade underneath curve
  for(i in 1:(nts+1)){
    
    if(i == 1){
      xs <- c(x[x <= ts[i]], rev(x[x <= ts[i]]))
      ys <- c(y[x <= ts[i]], rep(0, length(y[x <= ts[i]])))
    }
    
    if(1 < i & i <= nts) {
      xs <- c(x[ts[i-1] <= x & x <= ts[i]], rev(x[ts[i-1] <= x & x <= ts[i]]))
      ys <- c(y[ts[i-1] <= x & x <= ts[i]], rep(0, length(y[ts[i-1] <= x & x <= ts[i]])))
    }
    
    if(i == (nts+1)){
      xs <- c(x[x >= ts[i-1]], rev(x[x >= ts[i-1]]))
      ys <- c(y[x >= ts[i-1]], rep(0, length(y[x >= ts[i-1]])))
    }
    
    polygon(x = xs, y = ys, col = cols[i])
  }
  
  abline(v = ts, lwd = 3, lty = 2)
  
  
  #get frequencies of discrete traits
  discr_props <- c(0, pnorm(q = ts, mean = mu, sd = 1), 1)
  discr_props <- diff(discr_props)
  discr_props
  barplot(discr_props, col = cols, ylim = c(0,0.6), xlab = "")
  title(xlab = "discrete trait frequency", line=1, cex.lab = 1.5)
  
  dev.off()
}

files <- list.files(pattern = "*.png", full.names = T)
files <- files[order(as.integer(sapply(1:length(files), function(x) strsplit(strsplit(files[x], split = "_")[[1]][4], split = ".png")[[1]][1])))]

files %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=25, loop = t_loop) %>% # animates, can opt for number of loops
  image_write("univ_threshold_mu_smooth.gif") # write to current dir
file.rename(from = files[1], to = "hide_image")
file.remove(list.files(pattern=".png"))


#######################
### vary thresholds ###
#######################

thresh_to_vary <- ceiling(median(1:nts))
amount_per_step <- 0.05
ts_locs <- c(seq(from = ts[thresh_to_vary], to = ts[thresh_to_vary-1], by = -amount_per_step), 
             seq(from = ts[thresh_to_vary-1], to = ts[thresh_to_vary+1], by = amount_per_step),
             seq(from = ts[thresh_to_vary+1], to = ts[thresh_to_vary], by = -amount_per_step))

for(i in 1:length(ts_locs)){
  cat(paste0(i, " "))
  y <- dnorm(x, mean = mu)
  ts[thresh_to_vary] <- ts_locs[i]
  
  png(filename = paste0("univ_threshold_plot_", i, ".png"), width = 1000, height = 350)
  par(mfrow = c(1,2), mar = c(4,4,3,2), xpd = F)
  
  plot(x, y, type = "l", lwd = 2, ylim = c(0,0.45), yaxt = "n", ylab = "", xlab = "")
  title(xlab = "latent liability", line=2.25, cex.lab = 1.5)
  
  #shade underneath curve
  for(i in 1:(nts+1)){
    
    if(i == 1){
      xs <- c(x[x <= ts[i]], rev(x[x <= ts[i]]))
      ys <- c(y[x <= ts[i]], rep(0, length(y[x <= ts[i]])))
    }
    
    if(1 < i & i <= nts) {
      xs <- c(x[ts[i-1] <= x & x <= ts[i]], rev(x[ts[i-1] <= x & x <= ts[i]]))
      ys <- c(y[ts[i-1] <= x & x <= ts[i]], rep(0, length(y[ts[i-1] <= x & x <= ts[i]])))
    }
    
    if(i == (nts+1)){
      xs <- c(x[x >= ts[i-1]], rev(x[x >= ts[i-1]]))
      ys <- c(y[x >= ts[i-1]], rep(0, length(y[x >= ts[i-1]])))
    }
    
    polygon(x = xs, y = ys, col = cols[i])
  }
  
  abline(v = ts, lwd = 3, lty = 2)
  
  
  #get frequencies of discrete traits
  discr_props <- c(0, pnorm(q = ts, mean = mu, sd = 1), 1)
  discr_props <- diff(discr_props)
  discr_props
  barplot(discr_props, col = cols, ylim = c(0,0.6), xlab = "")
  title(xlab = "discrete trait frequency", line=1, cex.lab = 1.5)
  
  dev.off()
}

files <- list.files(pattern = "*.png", full.names = T)
files <- files[order(as.integer(sapply(1:length(files), function(x) strsplit(strsplit(files[x], split = "_")[[1]][4], split = ".png")[[1]][1])))]
files %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=25, loop = t_loop) %>% # animates, can opt for number of loops
  image_write("univ_threshold_thresh_smooth.gif") # write to current dir
file.remove(list.files(pattern=".png"))


########################
### evolve mus by bm ###
########################

mus <- cumsum(c(mus[length(mus)], rnorm(n = 75, mean = 0, sd = 0.05)))

for(i in 1:length(mus)){
  cat(paste0(i, " "))
  mu <- mus[i] #location of std norm mean
  y <- dnorm(x, mean = mu)
  
  png(filename = paste0("univ_threshold_plot_", i, ".png"), width = 1000, height = 350)
  par(mfrow = c(1,2), mar = c(4,4,3,2), xpd = F)
  
  plot(x, y, type = "l", lwd = 2, ylim = c(0,0.45), yaxt = "n", ylab = "", xlab = "")
  title(xlab = "latent liability", line=2.25, cex.lab = 1.5)
  
  #shade underneath curve
  for(i in 1:(nts+1)){
    
    if(i == 1){
      xs <- c(x[x <= ts[i]], rev(x[x <= ts[i]]))
      ys <- c(y[x <= ts[i]], rep(0, length(y[x <= ts[i]])))
    }
    
    if(1 < i & i <= nts) {
      xs <- c(x[ts[i-1] <= x & x <= ts[i]], rev(x[ts[i-1] <= x & x <= ts[i]]))
      ys <- c(y[ts[i-1] <= x & x <= ts[i]], rep(0, length(y[ts[i-1] <= x & x <= ts[i]])))
    }
    
    if(i == (nts+1)){
      xs <- c(x[x >= ts[i-1]], rev(x[x >= ts[i-1]]))
      ys <- c(y[x >= ts[i-1]], rep(0, length(y[x >= ts[i-1]])))
    }
    
    polygon(x = xs, y = ys, col = cols[i])
  }
  
  abline(v = ts, lwd = 3, lty = 2)
  
  
  #get frequencies of discrete traits
  discr_props <- c(0, pnorm(q = ts, mean = mu, sd = 1), 1)
  discr_props <- diff(discr_props)
  discr_props
  barplot(discr_props, col = cols, ylim = c(0,0.6), xlab = "")
  title(xlab = "discrete trait frequency", line=1, cex.lab = 1.5)
  
  dev.off()
  
}

files <- list.files(pattern = "*.png", full.names = T)
files <- files[order(as.integer(sapply(1:length(files), function(x) strsplit(strsplit(files[x], split = "_")[[1]][4], split = ".png")[[1]][1])))]
files %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=25, loop = t_loop) %>% # animates, can opt for number of loops
  image_write("univ_threshold_thresh_bm.gif") # write to current dir
file.remove(list.files(pattern=".png"))


#plot the brownian motion too
for(i in 1:length(mus)){
  
  png(filename = paste0("univ_threshold_plot_", i, ".png"), width = 500, height = 350)
  
  plot(seq(from = 0, to = 1, length.out = length(mus))[1:i], mus[1:i], 
       type = "l", xlim = c(0,1), ylim = c(min(mus), max(mus)),
       xlab = "", ylab = "")
  title(xlab = "time", ylab = "mean liability", line=2.25, cex.lab = 1.5)
  
  dev.off()
  
}

files <- list.files(pattern = "*.png", full.names = T)
files <- files[order(as.integer(sapply(1:length(files), function(x) strsplit(strsplit(files[x], split = "_")[[1]][4], split = ".png")[[1]][1])))]
files %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=25, loop = t_loop) %>% # animates, can opt for number of loops
  image_write("bm.gif") # write to current dir
file.remove(list.files(pattern=".png"))
file.rename(from = "hide_image", to = files[1])
