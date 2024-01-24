library(plot3D)
library(TeachingDemos)
library(mvtnorm)
library(magick)
library(purrr)
library(ape)
library(RColorBrewer)
setwd("/Users/nikolai/3DBM_visual")


## testing stuff out

nsteps <- 100
t <- seq(0, 1, by = 0.01)
h <- c(4, rnorm(n = nsteps, mean = 0, sd = 0.1))
h <- cumsum(h)
w <- c(4, rnorm(n = nsteps, mean = 0, sd = 0.1))
w <- cumsum(w)

scatter3D(x = h, y = w, z = t, phi = 0, type = "l", 
          ticktype = "simple", lwd = 1, xlab = "height", ylab = "width", zlab = zlabel, add = T, col = colorChoices[1])

scatter3D(x = h, y = w, z = t, type = "l", xlim = c(0,8), ylim = c(0,8), 
          ticktype = "detailed", lwd = 1, xlab = "height", ylab = "width", zlab = zlabel, add = F, col = colorChoices[1], phi = 45, theta = 45)


setwd("/Users/nikolai/3DBM_visual")
x <- 1
y <- 1
for(i in 1:length(h)){
  png(width = pngWidth, height = pngHeight,  paste0("cyl", paste0(rep(0, (6 - nchar(i))), collapse = ""), i, ".png"))
  plot(x,y, xlim=c(-.5,2.5), ylim=c(-.5,2.5), xlab = i, type='n') 
  my.symbols(x,y,ms.cylinder, height=h[i], width=w[i], angle=.5)
  dev.off()
}

list.files(pattern = "*.png", full.names = T) %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=5) %>% # animates, can opt for number of loops
  image_write("test.gif") # write to current dir

system(command = "convert -delay 10 *.png cyl.gif")
file.remove(list.files(pattern=".png"))

####################################################################################################
####################################################################################################
####################################################################################################
#try animating evolution on a simple 3-tip tree
####################################################################################################
####################################################################################################
####################################################################################################

#SIMULATE DATA

#subtending branch to root -- 0.25 BLs
stepsPerBL <- 100
rm <- matrix(data = c(1,0.5,0.5,2), nrow = 2)
rootMean <- c(4,4)
BL1 <- 0.25
nsteps <- BL1*stepsPerBL
t1 <- seq(0, BL1, by = 1 / stepsPerBL)
hw1 <- rbind(rootMean, rmvnorm(n = nsteps, mean = c(0,0), sigma = rm / stepsPerBL))
hw1 <- cbind(cumsum(hw1[,1]), cumsum(hw1[,2]))

#split 1, 0.5 BL on branch 2 to split 2, 1 BL on branch 3 to tip 1
BL2 <- 0.5
nsteps <- BL2*stepsPerBL
t2 <- seq(tail(t1, 1), tail(t1, 1) + BL2, by = 1 / stepsPerBL)
hw2 <- rbind(hw1[dim(hw1)[1],], rmvnorm(n = nsteps, mean = c(0,0), sigma = rm / stepsPerBL))
hw2 <- cbind(cumsum(hw2[,1]), cumsum(hw2[,2]))

BL3 <- 1
nsteps <- BL3*stepsPerBL
t3 <- seq(tail(t1, 1), tail(t1, 1) + BL3, by = 1 / stepsPerBL)
hw3 <- rbind(hw1[dim(hw1)[1],], rmvnorm(n = nsteps, mean = c(0,0), sigma = rm / stepsPerBL))
hw3 <- cbind(cumsum(hw3[,1]), cumsum(hw3[,2]))

# split 2 to tips 2 and 3, each branch with 0.5 BLs
BL4 <- 0.5
nsteps <- BL4*stepsPerBL
t4 <- seq(tail(t2, 1), tail(t2, 1) + BL4, by = 1 / stepsPerBL)
hw4 <- rbind(hw2[dim(hw2)[1],], rmvnorm(n = nsteps, mean = c(0,0), sigma = rm / stepsPerBL))
hw4 <- cbind(cumsum(hw4[,1]), cumsum(hw4[,2]))

BL5 <- 0.5
nsteps <- BL5*stepsPerBL
t5 <- seq(tail(t2, 1), tail(t2, 1) + BL5, by = 1 / stepsPerBL)
hw5 <- rbind(hw2[dim(hw2)[1],], rmvnorm(n = nsteps, mean = c(0,0), sigma = rm / stepsPerBL))
hw5 <- cbind(cumsum(hw5[,1]), cumsum(hw5[,2]))

tree <- read.tree(text = "((A:0.5,B:0.5):0.5,C:1):0.25;")
plot(tree, type = "cladogram", root.edge = T, edge.width = 4, font = 2, cex = 3, direction = "up", tip.color = colorChoices[c(4,5,3)], edge.color = colorChoices[c(2,4,5,3)], srt = -90, label.offset = .06)
axisPhylo(2,backward = F)

####################################################################################################
####################################################################################################
####################################################################################################

#ANIMATE
colorChoices <- c(1,2,3,4,5)
colorChoices <- brewer.pal(5, "Dark2")
# colorChoices <- brewer.pal(5, "Set1")
pngWidth <- 960
pngHeight <- 960
zlabel <- ""

#find animation parameters
project = c("xy0", "yz0", "xz0", "xy1", "yz1", "xz1") #all planes projected upon
project = c("xz1", "yz0")
# project = ""
wigglyLineThickness <- 2
all <- rbind(hw1, hw2, hw3, hw4, hw5)
xl <- c(min(all[,1]), max(all[,1]))
yl <- c(min(all[,2]), max(all[,2]))
rotate <- TRUE
bounce <- TRUE
if(rotate){
  if (bounce){
    rotateDegrees <- seq(0, 22.5, length.out = (length(t1) + length(t2) + length(t4) - 3)/2)
    rotateDegrees <- c(rotateDegrees, rev(rotateDegrees))
  } else {
    rotateDegrees <- seq(0, 90, length.out = (length(t1) + length(t2) + length(t4) - 3))
  }
} else {
  rotateDegrees <- seq(0, 0, length.out = (length(t1) + length(t2) + length(t4) - 3))
}

counter <- 0
for(i in 2:length(t1)){
  counter <- counter + 1
  png(width = pngWidth, height = pngHeight,  paste0("mvBM", paste0(rep(0, (6 - nchar(counter))), collapse = ""), counter, ".png"))
  #main squigglies
  scatter3D(x = hw1[1:i,1], y = hw1[1:i,2], z = t1[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t1[i]), 
            ticktype = "detailed", lwd = wigglyLineThickness, xlab = "height", ylab = "width", zlab = zlabel, add = F, col = colorChoices[1], phi = 45, theta = 45 + rotateDegrees[counter])
  #projections
  projLWD <- 0.4
  if(any(project == "xy0")){
    scatter3D(x = hw1[1:i,1], y = hw1[1:i,2], z = rep(0, length(hw1[1:i,1])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
  }
  if(any(project == "xy1")){
    scatter3D(x = hw1[1:i,1], y = hw1[1:i,2], z = rep(t1[i], length(hw1[1:i,1])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
   }
  if(any(project == "xz0")){
    scatter3D(x = hw1[1:i,1], y =rep(yl[1], length(hw1[1:i,1])), z = t1[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
  }
  if(any(project == "xz1")){
    scatter3D(x = hw1[1:i,1], y =rep(yl[2], length(hw1[1:i,1])), z = t1[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
  }
  if(any(project == "yz0")){
    scatter3D(x = rep(xl[1], length(hw1[1:i,1])), y = hw1[1:i,2], z = t1[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
  }
  if(any(project == "yz1")){
    scatter3D(x = rep(xl[2], length(hw1[1:i,1])), y = hw1[1:i,2], z = t1[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
  }  
  dev.off()
}

for(i in 2:length(t2)){
  counter <- counter + 1
  png(width = pngWidth, height = pngHeight,  paste0("mvBM", paste0(rep(0, (6 - nchar(counter))), collapse = ""), counter, ".png"))
  #main squigglies
  scatter3D(x = hw1[,1], y = hw1[,2], z = t1, type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
            ticktype = "detailed", lwd = 2, xlab = "height", ylab = "width", zlab = zlabel, add = F, col = colorChoices[1], phi = 45, theta = 45 + rotateDegrees[counter])
  scatter3D(x = hw2[1:i,1], y = hw2[1:i,2], z = t2[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
            ticktype = "detailed", lwd = 2, xlab = "height", ylab = "width", zlab = zlabel, col = colorChoices[2], phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  scatter3D(x = hw3[1:i,1], y = hw3[1:i,2], z = t3[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
            ticktype = "detailed", lwd = 2, xlab = "height", ylab = "width", zlab = zlabel, col = colorChoices[3], phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  #projections
  projLWD <- 0.3
  if(any(project == "xy0")){
    scatter3D(x = hw1[,1], y = hw1[,2], z = rep(0, length(hw1[,1])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
    scatter3D(x = hw2[1:i,1], y = hw2[1:i,2], z = rep(0, length(hw2[1:i,1])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[2], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw3[1:i,1], y = hw3[1:i,2], z = rep(0, length(hw3[1:i,1])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[3], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  }
  if(any(project == "xy1")){
    scatter3D(x = hw1[,1], y = hw1[,2], z = rep(t2[i], length(hw1[,1])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
    scatter3D(x = hw2[1:i,1], y = hw2[1:i,2], z = rep(t2[i], length(hw2[1:i,1])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[2], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw3[1:i,1], y = hw3[1:i,2], z = rep(t2[i], length(hw3[1:i,1])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[3], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  }
  if(any(project == "xz0")){
    scatter3D(x = hw1[,1], y =rep(yl[1], length(hw1[,1])), z = t1, type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
    scatter3D(x = hw2[1:i,1], y = rep(yl[1], length(hw2[1:i,1])), z = t2[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[2], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw3[1:i,1], y = rep(yl[1], length(hw3[1:i,1])), z = t3[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[3], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  }
  if(any(project == "xz1")){
    scatter3D(x = hw1[,1], y =rep(yl[2], length(hw1[,1])), z = t1, type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
    scatter3D(x = hw2[1:i,1], y = rep(yl[2], length(hw2[1:i,1])), z = t2[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[2], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw3[1:i,1], y = rep(yl[2], length(hw3[1:i,1])), z = t3[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[3], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  }
  if(any(project == "yz0")){
    scatter3D(x = rep(xl[1], length(hw1[,1])), y = hw1[,2], z = t1, type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
    scatter3D(x = rep(xl[1], length(hw2[1:i,1])), y = hw2[1:i,2], z = t2[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[2], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = rep(xl[1], length(hw3[1:i,1])), y = hw3[1:i,2], z = t3[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[3], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  }
  if(any(project == "yz1")){
    scatter3D(x = rep(xl[2], length(hw1[,1])), y = hw1[,2], z = t1, type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
    scatter3D(x = rep(xl[2], length(hw2[1:i,1])), y = hw2[1:i,2], z = t2[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[2], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = rep(xl[2], length(hw3[1:i,1])), y = hw3[1:i,2], z = t3[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[3], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  }  
  dev.off()
}

for(i in 2:length(t4)){
  counter <- counter + 1
  png(width = pngWidth, height = pngHeight,  paste0("mvBM", paste0(rep(0, (6 - nchar(counter))), collapse = ""), counter, ".png"))
  #main squigglies
  scatter3D(x = hw1[,1], y = hw1[,2], z = t1, type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
            ticktype = "detailed", lwd = wigglyLineThickness, xlab = "height", ylab = "width", zlab = zlabel, add = F, col = colorChoices[1], phi = 45, theta = 45 + rotateDegrees[counter])
  scatter3D(x = hw2[,1], y = hw2[,2], z = t2, type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
            ticktype = "detailed", lwd = wigglyLineThickness, xlab = "height", ylab = "width", zlab = zlabel, col = colorChoices[2], phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  scatter3D(x = hw3[1:(length(t2)+i-1),1], y = hw3[1:(length(t2)+i-1),2], z = t3[1:(length(t2)+i-1)], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
            ticktype = "detailed", lwd = wigglyLineThickness, xlab = "height", ylab = "width", zlab = zlabel, col = colorChoices[3], phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  scatter3D(x = hw4[1:i,1], y = hw4[1:i,2], z = t4[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
            ticktype = "detailed", lwd = wigglyLineThickness, xlab = "height", ylab = "width", zlab = zlabel, col = colorChoices[4], phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  scatter3D(x = hw5[1:i,1], y = hw5[1:i,2], z = t5[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
            ticktype = "detailed", lwd = wigglyLineThickness, xlab = "height", ylab = "width", zlab = zlabel, col = colorChoices[5], phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  #projections
  projLWD <- 0.2
  if(any(project == "xy0")){
    scatter3D(x = hw1[,1], y = hw1[,2], z = rep(0,length(t1)), type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
    scatter3D(x = hw2[,1], y = hw2[,2], z = rep(0,length(t2)), type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[2], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw3[1:(length(t2)+i-1),1], y = hw3[1:(length(t2)+i-1),2], z = rep(0,length(t3[1:(length(t2)+i-1)])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[3], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw4[1:i,1], y = hw4[1:i,2], z = rep(0,length(t4[1:i])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[4], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw5[1:i,1], y = hw5[1:i,2], z = rep(0,length(t5[1:i])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[5], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  }
  if(any(project == "xy1")){
    scatter3D(x = hw1[,1], y = hw1[,2], z = rep(t4[i],length(t1)), type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
    scatter3D(x = hw2[,1], y = hw2[,2], z = rep(t4[i],length(t2)), type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[2], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw3[1:(length(t2)+i-1),1], y = hw3[1:(length(t2)+i-1),2], z = rep(t4[i],length(t3[1:(length(t2)+i-1)])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[3], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw4[1:i,1], y = hw4[1:i,2], z = rep(t4[i], length(hw4[1:i,1])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[4], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw5[1:i,1], y = hw5[1:i,2], z = rep(t5[i], length(hw5[1:i,1])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = projLWD, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[5], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  }
  if(any(project == "xz0")){
    scatter3D(x = hw1[,1], y = rep(yl[1], length(hw1[,2])), z = t1, type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, xlab = "height", ylab = "width", zlab = zlabel, add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
    scatter3D(x = hw2[,1], y = rep(yl[1], length(hw2[,2])), z = t2, type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[2], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw3[1:(length(t2)+i-1),1], y = rep(yl[1], length(hw3[1:(length(t2)+i-1),2])), z = t3[1:(length(t2)+i-1)], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[3], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw4[1:i,1], y = rep(yl[1], length(hw4[1:i,2])), z = t4[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[4], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw5[1:i,1], y = rep(yl[1], length(hw5[1:i,2])), z = t5[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[5], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  }
  if(any(project == "xz1")){
    scatter3D(x = hw1[,1], y = rep(yl[2], length(hw1[,2])), z = t1, type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, xlab = "height", ylab = "width", zlab = zlabel, add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
    scatter3D(x = hw2[,1], y = rep(yl[2], length(hw2[,2])), z = t2, type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[2], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw3[1:(length(t2)+i-1),1], y = rep(yl[2], length(hw3[1:(length(t2)+i-1),2])), z = t3[1:(length(t2)+i-1)], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[3], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw4[1:i,1], y = rep(yl[2], length(hw4[1:i,2])), z = t4[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[4], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw5[1:i,1], y = rep(yl[2], length(hw5[1:i,2])), z = t5[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[5], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  }
  if(any(project == "yz0")){
    scatter3D(x = rep(xl[1], length(hw1[,1])), y = hw1[,2], z = t1, type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, xlab = "height", ylab = "width", zlab = zlabel, add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
    scatter3D(x = rep(xl[1], length(hw2[,1])), y = hw2[,2], z = t2, type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[2], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = rep(xl[1], length(hw3[1:(length(t2)+i-1),1])), y = hw3[1:(length(t2)+i-1),2], z = t3[1:(length(t2)+i-1)], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[3], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = rep(xl[1], length(hw4[1:i,1])), y = hw4[1:i,2], z = t4[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[4], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = rep(xl[1], length(hw5[1:i,1])), y = hw5[1:i,2], z = t5[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[5], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  }
  if(any(project == "yz1")){
    scatter3D(x = rep(xl[2], length(hw1[,1])), y = hw1[,2], z = t1, type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, xlab = "height", ylab = "width", zlab = zlabel, add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
    scatter3D(x = rep(xl[2], length(hw2[,1])), y = hw2[,2], z = t2, type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[2], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = rep(xl[2], length(hw3[1:(length(t2)+i-1),1])), y = hw3[1:(length(t2)+i-1),2], z = t3[1:(length(t2)+i-1)], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[3], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = rep(xl[2], length(hw4[1:i,1])), y = hw4[1:i,2], z = t4[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[4], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = rep(xl[2], length(hw5[1:i,1])), y = hw5[1:i,2], z = t5[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, xlab = "height", ylab = "width", zlab = zlabel, col = adjustcolor(colorChoices[5], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  }  
  dev.off()
}

list.files(pattern = "*.png", full.names = T) %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=10) %>% # animates, can opt for number of loops
  image_write("treeGrow.gif") # write to current dir
file.remove(list.files(pattern=".png"))

####################################################################################################
#draw rectangles
####################################################################################################

gradSegment <- function(x0, y0, x1, y1, numFade = 100, col = colorChoices[1], exponent = 1, lwd = 1){
  for(i in 1:numFade){
    prop <- (i/numFade)^exponent
    xp <- (1-prop)*x0 + prop*x1
    yp <- (1-prop)*y0 + prop*y1
    segments(x0,y0,xp,yp, lwd = lwd, col = adjustcolor(col, 1/numFade))
  }
}

counter <- 0
for(i in 1:length(t1)){
  counter <- counter + 1
  png(width = pngWidth, height = pngHeight,  paste0("rec", paste0(rep(0, (6 - nchar(counter))), collapse = ""), counter, ".png"))
  plot(0,0,xlim=c(0,30), ylim=c(0,30), type='n', axes = F, xlab = "", ylab = "") 
  h=hw1[i,1]; w=hw1[i,2]
  bottomLeft <- c(0,0)
  rect(bottomLeft[1],bottomLeft[2],bottomLeft[1]+w,bottomLeft[2]+h, col = adjustcolor(colorChoices[1], .8))
  dev.off()
}

for(j in 1:length(t2)){
  counter <- counter + 1
  png(width = pngWidth, height = pngHeight,  paste0("rec", paste0(rep(0, (6 - nchar(counter))), collapse = ""), counter, ".png"))
  plot(0,0,xlim=c(0,30), ylim=c(0,30), type='n', axes = F, xlab = "", ylab = "") 
  
  #1 lineage
  h=hw1[i,1]; w=hw1[i,2]
  bottomLeft <- c(0,0)
  rect(bottomLeft[1],bottomLeft[2],bottomLeft[1]+w,bottomLeft[2]+h, col = adjustcolor(colorChoices[1], .8))
  
  #2 lineages
  h=hw2[j,1]; w=hw2[j,2]
  bottomLeft <- c(0,10)
  rect(bottomLeft[1],bottomLeft[2],bottomLeft[1]+w,bottomLeft[2]+h, col = adjustcolor(colorChoices[2], .8))
  arrowFrom <- c((0 + hw1[i,2])/2, (0 + hw1[i,1])+.3)
  arrowTo <- c((bottomLeft[1]+w)/2, (bottomLeft[2])-.3)
  arrows(arrowFrom[1],arrowFrom[2],arrowTo[1],arrowTo[2], col = colorChoices[2], lwd = 2, length = 0.1)
  gradSegment(arrowFrom[1],arrowFrom[2],arrowTo[1],arrowTo[2], col = colorChoices[1], exponent = 1, lwd = 3)
  
  h=hw3[j,1]; w=hw3[j,2]
  bottomLeft <- c(10,10)
  rect(bottomLeft[1],bottomLeft[2],bottomLeft[1]+w,bottomLeft[2]+h, col = adjustcolor(colorChoices[3], .8))
  arrowFrom <- c((0 + hw1[i,2])/2, (0 + hw1[i,1])+.3)
  arrowTo <- c((2*bottomLeft[1]+w)/2, (bottomLeft[2])-.3)
  arrows(arrowFrom[1],arrowFrom[2],arrowTo[1],arrowTo[2], col = colorChoices[3], lwd = 2, length = 0.1)
  gradSegment(arrowFrom[1],arrowFrom[2],arrowTo[1],arrowTo[2], col = colorChoices[1], exponent = 1, lwd = 3)
  
  dev.off()
  
}

for(k in 1:length(t4)){
  counter <- counter + 1
  png(width = pngWidth, height = pngHeight,  paste0("rec", paste0(rep(0, (6 - nchar(counter))), collapse = ""), counter, ".png"))
  plot(0,0,xlim=c(0,30), ylim=c(0,30), type='n', axes = F, xlab = "", ylab = "") 
  
  #1 lineage
  h=hw1[i,1]; w=hw1[i,2]
  bottomLeft <- c(0,0)
  rect(bottomLeft[1],bottomLeft[2],bottomLeft[1]+w,bottomLeft[2]+h, col = adjustcolor(colorChoices[1], .8))
  
  #2 lineages
  h=hw2[j,1]; w=hw2[j,2]
  bottomLeft <- c(0,10)
  rect(bottomLeft[1],bottomLeft[2],bottomLeft[1]+w,bottomLeft[2]+h, col = adjustcolor(colorChoices[2], .8))
  arrowFrom <- c((0 + hw1[i,2])/2, (0 + hw1[i,1])+.3)
  arrowTo <- c((bottomLeft[1]+w)/2, (bottomLeft[2])-.3)
  arrows(arrowFrom[1],arrowFrom[2],arrowTo[1],arrowTo[2], col = colorChoices[2], lwd = 2, length = 0.1)
  gradSegment(arrowFrom[1],arrowFrom[2],arrowTo[1],arrowTo[2], col = colorChoices[1], exponent = 1, lwd = 3)
  
  h=hw3[j,1]; w=hw3[j,2]
  bottomLeft <- c(10,10)
  rect(bottomLeft[1],bottomLeft[2],bottomLeft[1]+w,bottomLeft[2]+h, col = adjustcolor(colorChoices[3], .8))
  arrowFrom <- c((0 + hw1[i,2])/2, (0 + hw1[i,1])+.3)
  arrowTo <- c((2*bottomLeft[1]+w)/2, (bottomLeft[2])-.3)
  arrows(arrowFrom[1],arrowFrom[2],arrowTo[1],arrowTo[2], col = colorChoices[3], lwd = 2, length = 0.1)
  gradSegment(arrowFrom[1],arrowFrom[2],arrowTo[1],arrowTo[2], col = colorChoices[1], exponent = 1, lwd = 3)
  
  
  #3 lineages
  h=hw3[(k + length(t2) - 1),1]; w=hw3[(k + length(t2) - 1),2]
  bottomLeft <- c(20,20)
  rect(bottomLeft[1],bottomLeft[2],bottomLeft[1]+w,bottomLeft[2]+h, col = adjustcolor(colorChoices[3], .8))
  arrowFrom <- c((10*2 + hw3[j,2])/2, (10 + hw3[j,1])+.3)
  arrowTo <- c((bottomLeft[1]*2+w)/2, (bottomLeft[2])-.3)
  arrows(arrowFrom[1],arrowFrom[2],arrowTo[1],arrowTo[2], col = colorChoices[3], lwd = 2, length = 0.1)
  gradSegment(arrowFrom[1],arrowFrom[2],arrowTo[1],arrowTo[2], col = colorChoices[3], exponent = 1, lwd = 3)
  
  h=hw4[k,1]; w=hw4[k,2]
  bottomLeft <- c(0,20)
  rect(bottomLeft[1],bottomLeft[2],bottomLeft[1]+w,bottomLeft[2]+h, col = adjustcolor(colorChoices[4], .8))
  arrowFrom <- c((0*2 + hw2[j,2])/2, (10 + hw2[j,1])+.3)
  arrowTo <- c((bottomLeft[1]*2+w)/2, (bottomLeft[2])-.3)
  arrows(arrowFrom[1],arrowFrom[2],arrowTo[1],arrowTo[2], col = colorChoices[4], lwd = 2, length = 0.1)
  gradSegment(arrowFrom[1],arrowFrom[2],arrowTo[1],arrowTo[2], col = colorChoices[2], exponent = 1, lwd = 3)
  
  
  h=hw5[k,1]; w=hw5[k,2]
  bottomLeft <- c(10,20)
  rect(bottomLeft[1],bottomLeft[2],bottomLeft[1]+w,bottomLeft[2]+h, col = adjustcolor(colorChoices[5], .8))
  arrowFrom <- c((0*2 + hw2[j,2])/2, (10 + hw2[j,1])+.3)
  arrowTo <- c((bottomLeft[1]*2+w)/2, (bottomLeft[2])-.3)
  arrows(arrowFrom[1],arrowFrom[2],arrowTo[1],arrowTo[2], col = colorChoices[5], lwd = 2, length = 0.1)
  gradSegment(arrowFrom[1],arrowFrom[2],arrowTo[1],arrowTo[2], col = colorChoices[2], exponent = 1, lwd = 3)
  
  dev.off()
  
}


list.files(pattern = "*.png", full.names = T) %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=10) %>% # animates, can opt for number of loops
  image_write("recMorph.gif") # write to current dir
file.remove(list.files(pattern=".png"))

####################################################################################################
####################################################################################################
####################################################################################################





#draw cylinders??
ms.cylinder <- function(height, width, angle) { 
  theta <- seq(2*pi, 0, length=200) 
  x1 <- cos(theta) * width 
  y1 <- sin(theta) * width * angle 
  
  x <- c(x1, x1[1:100], x1[100]) 
  y <- c(y1 + height/2, y1[1:100]-height/2, y1[100]+height/2) 
  
  return(cbind(x,y)) 
}

x <- 1
y <- 1
for(i in 1:length(t1)){
  png(width = pngWidth, height = pngHeight,  paste0("cyl", paste0(rep(0, (6 - nchar(i))), collapse = ""), i, ".png"))
  plot(x,y, xlim=c(-5,25), ylim=c(-5,25), xlab = i, type='n') 
  my.symbols(x,y,ms.cylinder, height=hw1[i,2], width=hw1[i,2], angle=.5, col =2)
  dev.off()
}
list.files(pattern = "*.png", full.names = T) %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=10) %>% # animates, can opt for number of loops
  image_write("cylMorph.gif") # write to current dir
file.remove(list.files(pattern=".png"))

