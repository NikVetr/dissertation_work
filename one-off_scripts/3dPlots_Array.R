library(plot3D)
library(TeachingDemos)
library(mvtnorm)
library(magick)
library(purrr)
library(ape)
library(RColorBrewer)
setwd("/Users/nikolai/3DBM_visual")


#SIMULATE DATA

for(run in 1:17){

print(run)

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
colorChoices <- sample(c(brewer.pal(9, "Set1"),brewer.pal(8, "Set2"),brewer.pal(12, "Set3"), brewer.pal(8, "Dark2")), size = 5)
pngWidth <- 360
pngHeight <- 360
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
    ticktype = "detailed", lwd = wigglyLineThickness, zlab = zlabel, xlab = "", ylab = "",  add = F, col = colorChoices[1], phi = 45, theta = 45 + rotateDegrees[counter])
  #projections
  projLWD <- 0.4
  if(any(project == "xy0")){
    scatter3D(x = hw1[1:i,1], y = hw1[1:i,2], z = rep(0, length(hw1[1:i,1])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
  }
  if(any(project == "xy1")){
    scatter3D(x = hw1[1:i,1], y = hw1[1:i,2], z = rep(t1[i], length(hw1[1:i,1])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
  }
  if(any(project == "xz0")){
    scatter3D(x = hw1[1:i,1], y =rep(yl[1], length(hw1[1:i,1])), z = t1[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
  }
  if(any(project == "xz1")){
    scatter3D(x = hw1[1:i,1], y =rep(yl[2], length(hw1[1:i,1])), z = t1[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
  }
  if(any(project == "yz0")){
    scatter3D(x = rep(xl[1], length(hw1[1:i,1])), y = hw1[1:i,2], z = t1[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
  }
  if(any(project == "yz1")){
    scatter3D(x = rep(xl[2], length(hw1[1:i,1])), y = hw1[1:i,2], z = t1[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
  }  
  dev.off()
}

for(i in 2:length(t2)){
  counter <- counter + 1
  png(width = pngWidth, height = pngHeight,  paste0("mvBM", paste0(rep(0, (6 - nchar(counter))), collapse = ""), counter, ".png"))
  #main squigglies
  scatter3D(x = hw1[,1], y = hw1[,2], z = t1, type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
            ticktype = "detailed", lwd = 2, zlab = zlabel, xlab = "", ylab = "",  add = F, col = colorChoices[1], phi = 45, theta = 45 + rotateDegrees[counter])
  scatter3D(x = hw2[1:i,1], y = hw2[1:i,2], z = t2[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
            ticktype = "detailed", lwd = 2, zlab = zlabel, xlab = "", ylab = "",  col = colorChoices[2], phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  scatter3D(x = hw3[1:i,1], y = hw3[1:i,2], z = t3[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
            ticktype = "detailed", lwd = 2, zlab = zlabel, xlab = "", ylab = "",  col = colorChoices[3], phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  #projections
  projLWD <- 0.3
  if(any(project == "xy0")){
    scatter3D(x = hw1[,1], y = hw1[,2], z = rep(0, length(hw1[,1])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
    scatter3D(x = hw2[1:i,1], y = hw2[1:i,2], z = rep(0, length(hw2[1:i,1])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[2], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw3[1:i,1], y = hw3[1:i,2], z = rep(0, length(hw3[1:i,1])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[3], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  }
  if(any(project == "xy1")){
    scatter3D(x = hw1[,1], y = hw1[,2], z = rep(t2[i], length(hw1[,1])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
    scatter3D(x = hw2[1:i,1], y = hw2[1:i,2], z = rep(t2[i], length(hw2[1:i,1])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[2], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw3[1:i,1], y = hw3[1:i,2], z = rep(t2[i], length(hw3[1:i,1])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[3], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  }
  if(any(project == "xz0")){
    scatter3D(x = hw1[,1], y =rep(yl[1], length(hw1[,1])), z = t1, type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
    scatter3D(x = hw2[1:i,1], y = rep(yl[1], length(hw2[1:i,1])), z = t2[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[2], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw3[1:i,1], y = rep(yl[1], length(hw3[1:i,1])), z = t3[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[3], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  }
  if(any(project == "xz1")){
    scatter3D(x = hw1[,1], y =rep(yl[2], length(hw1[,1])), z = t1, type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
    scatter3D(x = hw2[1:i,1], y = rep(yl[2], length(hw2[1:i,1])), z = t2[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[2], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw3[1:i,1], y = rep(yl[2], length(hw3[1:i,1])), z = t3[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[3], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  }
  if(any(project == "yz0")){
    scatter3D(x = rep(xl[1], length(hw1[,1])), y = hw1[,2], z = t1, type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
    scatter3D(x = rep(xl[1], length(hw2[1:i,1])), y = hw2[1:i,2], z = t2[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[2], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = rep(xl[1], length(hw3[1:i,1])), y = hw3[1:i,2], z = t3[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[3], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  }
  if(any(project == "yz1")){
    scatter3D(x = rep(xl[2], length(hw1[,1])), y = hw1[,2], z = t1, type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
    scatter3D(x = rep(xl[2], length(hw2[1:i,1])), y = hw2[1:i,2], z = t2[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[2], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = rep(xl[2], length(hw3[1:i,1])), y = hw3[1:i,2], z = t3[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t2[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[3], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  }  
  dev.off()
}

for(i in 2:length(t4)){
  counter <- counter + 1
  png(width = pngWidth, height = pngHeight,  paste0("mvBM", paste0(rep(0, (6 - nchar(counter))), collapse = ""), counter, ".png"))
  #main squigglies
  scatter3D(x = hw1[,1], y = hw1[,2], z = t1, type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
            ticktype = "detailed", lwd = wigglyLineThickness, zlab = zlabel, xlab = "", ylab = "",  add = F, col = colorChoices[1], phi = 45, theta = 45 + rotateDegrees[counter])
  scatter3D(x = hw2[,1], y = hw2[,2], z = t2, type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
            ticktype = "detailed", lwd = wigglyLineThickness, zlab = zlabel, xlab = "", ylab = "",  col = colorChoices[2], phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  scatter3D(x = hw3[1:(length(t2)+i-1),1], y = hw3[1:(length(t2)+i-1),2], z = t3[1:(length(t2)+i-1)], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
            ticktype = "detailed", lwd = wigglyLineThickness, zlab = zlabel, xlab = "", ylab = "",  col = colorChoices[3], phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  scatter3D(x = hw4[1:i,1], y = hw4[1:i,2], z = t4[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
            ticktype = "detailed", lwd = wigglyLineThickness, zlab = zlabel, xlab = "", ylab = "",  col = colorChoices[4], phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  scatter3D(x = hw5[1:i,1], y = hw5[1:i,2], z = t5[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
            ticktype = "detailed", lwd = wigglyLineThickness, zlab = zlabel, xlab = "", ylab = "",  col = colorChoices[5], phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  #projections
  projLWD <- 0.2
  if(any(project == "xy0")){
    scatter3D(x = hw1[,1], y = hw1[,2], z = rep(0,length(t1)), type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
    scatter3D(x = hw2[,1], y = hw2[,2], z = rep(0,length(t2)), type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[2], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw3[1:(length(t2)+i-1),1], y = hw3[1:(length(t2)+i-1),2], z = rep(0,length(t3[1:(length(t2)+i-1)])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[3], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw4[1:i,1], y = hw4[1:i,2], z = rep(0,length(t4[1:i])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[4], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw5[1:i,1], y = hw5[1:i,2], z = rep(0,length(t5[1:i])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[5], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  }
  if(any(project == "xy1")){
    scatter3D(x = hw1[,1], y = hw1[,2], z = rep(t4[i],length(t1)), type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
    scatter3D(x = hw2[,1], y = hw2[,2], z = rep(t4[i],length(t2)), type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[2], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw3[1:(length(t2)+i-1),1], y = hw3[1:(length(t2)+i-1),2], z = rep(t4[i],length(t3[1:(length(t2)+i-1)])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[3], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw4[1:i,1], y = hw4[1:i,2], z = rep(t4[i], length(hw4[1:i,1])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[4], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw5[1:i,1], y = hw5[1:i,2], z = rep(t5[i], length(hw5[1:i,1])), type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = projLWD, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[5], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  }
  if(any(project == "xz0")){
    scatter3D(x = hw1[,1], y = rep(yl[1], length(hw1[,2])), z = t1, type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, zlab = zlabel, xlab = "", ylab = "",  add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
    scatter3D(x = hw2[,1], y = rep(yl[1], length(hw2[,2])), z = t2, type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[2], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw3[1:(length(t2)+i-1),1], y = rep(yl[1], length(hw3[1:(length(t2)+i-1),2])), z = t3[1:(length(t2)+i-1)], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[3], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw4[1:i,1], y = rep(yl[1], length(hw4[1:i,2])), z = t4[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[4], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw5[1:i,1], y = rep(yl[1], length(hw5[1:i,2])), z = t5[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[5], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  }
  if(any(project == "xz1")){
    scatter3D(x = hw1[,1], y = rep(yl[2], length(hw1[,2])), z = t1, type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, zlab = zlabel, xlab = "", ylab = "",  add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
    scatter3D(x = hw2[,1], y = rep(yl[2], length(hw2[,2])), z = t2, type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[2], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw3[1:(length(t2)+i-1),1], y = rep(yl[2], length(hw3[1:(length(t2)+i-1),2])), z = t3[1:(length(t2)+i-1)], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[3], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw4[1:i,1], y = rep(yl[2], length(hw4[1:i,2])), z = t4[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[4], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = hw5[1:i,1], y = rep(yl[2], length(hw5[1:i,2])), z = t5[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[5], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  }
  if(any(project == "yz0")){
    scatter3D(x = rep(xl[1], length(hw1[,1])), y = hw1[,2], z = t1, type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, zlab = zlabel, xlab = "", ylab = "",  add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
    scatter3D(x = rep(xl[1], length(hw2[,1])), y = hw2[,2], z = t2, type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[2], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = rep(xl[1], length(hw3[1:(length(t2)+i-1),1])), y = hw3[1:(length(t2)+i-1),2], z = t3[1:(length(t2)+i-1)], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[3], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = rep(xl[1], length(hw4[1:i,1])), y = hw4[1:i,2], z = t4[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[4], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = rep(xl[1], length(hw5[1:i,1])), y = hw5[1:i,2], z = t5[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[5], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  }
  if(any(project == "yz1")){
    scatter3D(x = rep(xl[2], length(hw1[,1])), y = hw1[,2], z = t1, type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, zlab = zlabel, xlab = "", ylab = "",  add = T, col = adjustcolor(colorChoices[1], .5), phi = 45, theta = 45 + rotateDegrees[counter])
    scatter3D(x = rep(xl[2], length(hw2[,1])), y = hw2[,2], z = t2, type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[2], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = rep(xl[2], length(hw3[1:(length(t2)+i-1),1])), y = hw3[1:(length(t2)+i-1),2], z = t3[1:(length(t2)+i-1)], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[3], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = rep(xl[2], length(hw4[1:i,1])), y = hw4[1:i,2], z = t4[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[4], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
    scatter3D(x = rep(xl[2], length(hw5[1:i,1])), y = hw5[1:i,2], z = t5[1:i], type = "l", xlim = xl, ylim = yl, zlim = c(0, t4[i]),
              ticktype = "detailed", lwd = 0.2, zlab = zlabel, xlab = "", ylab = "",  col = adjustcolor(colorChoices[5], .5), phi = 45, theta = 45 + rotateDegrees[counter], add = T)
  }  
  dev.off()
}

list.files(pattern = "*.png", full.names = T) %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=10) %>% # animates, can opt for number of loops
  image_write(paste0("manyTree", run, ".gif")) # write to current dir
file.remove(list.files(pattern=".png"))

}
