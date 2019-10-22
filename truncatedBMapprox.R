# truncated BM numerical approximation

library(mvtnorm)

#function for computing conditional multivariate normals
condMVN <- function(means, cov, obs, inds){
  A <- cov[-inds, -inds]
  C <- cov[inds, inds]
  B <- cov[-inds, inds]
  condCov <- A - B%*%solve(C)%*%t(B)
  condMean <- as.vector(means[-inds] + B%*%solve(C)%*%(obs - means[inds]))
  return(list(means = condMean, cov = condCov))
}


#bivariate BM on single branch

#in a circle
nhists <- 5e2
for(k in 1:nhists){
  print(k)
  
  #specify bounds of truncation, in this case a circle
  circleRadius <- 10
  circleCenter <- c(0,0)
  
  #specify starting and ending positions
  start <- c(-5,-4)
  end <- c(4,4.5)
  
  #specify brownian motion rate matrix
  R <- matrix(c(2,1,1,3), ncol = 2)
  
  #specify branch length separating start and end points
  BL <- 50
  
  #specify the number of intermediate timepoints over which to sample history
  n_IntermediateTimePoints <- 2^4-1
  
  #initialize matrix that will contain the bivariate history 
  positions <- matrix(0, ncol = 2, nrow = n_IntermediateTimePoints + 2)
  #specify the number of locations in total
  n_locations <- dim(positions)[1]
  #and specify start and end points
  positions[1,] <- start; positions[n_locations,] <- end
  
  #specify the position of each timepoint, where the starting branchlength is 0
  BLs <- 0:(n_IntermediateTimePoints+1)/(n_IntermediateTimePoints+1)*BL
  
  
  
  #keep track of which locations have observations already
  sampledInds <- rep(F, n_locations); sampledInds[1] <- TRUE; sampledInds[n_locations] <- TRUE
  
  #start sampling locations in log_2 stages, since under mvBM we can use already-sampled intermediate locations
  #to sample further intermediate locations
  for(i in 1:(log2(n_IntermediateTimePoints+1))){
    #which indices do we want to sample in this stage?
    indsToSample <- seq(from = 1, to = n_locations, length.out = 2^i+1)[-c(1,2^i+1)]
    if(i != 1){
      #make sure we haven't already sampled them!
      indsToSample <- setdiff(indsToSample, which(sampledInds))
      #and then iterate over unsampled indices
      for(j in 1:length(indsToSample)){
        indBeingSampled <- indsToSample[j]
        #find the two closest *already sampled* indices -- since mvBM is Markovian and evolution on branch i
        #independent of evolution on branch j, one can find conditionals only based on the two closest points 
        dists <- abs(which(sampledInds) - indBeingSampled)
        startEndInds <- which(sampledInds)[which(dists == min(dists))]
        strt <- positions[startEndInds[1],]
        ennd <- positions[startEndInds[2],]
        condNorm <- condMVN(means = rep(strt, 2), 
                            cov = kronecker(matrix(c(rep(BLs[indBeingSampled] - BLs[startEndInds[1]], 3), 
                            BLs[startEndInds[2]] - BLs[startEndInds[1]]), ncol = 2), R), obs = ennd, inds = c(3,4))
        samp <- rmvnorm(n = 1, mean = condNorm$means, sigma = condNorm$cov)
        #use rejection sampling to make sure we're in the allowable space
        while(norm(samp-circleCenter, type = "F") > circleRadius){
          samp <- rmvnorm(n = 1, mean = condNorm$means, sigma = condNorm$cov)
        }
        positions[indBeingSampled,] <- samp
      }
    } else {
      #special case when sampling first value
      condNorm <- condMVN(means = rep(start, 2), cov = kronecker(matrix(c(rep(BLs[indsToSample], 3), 
                          BLs[n_locations]), ncol = 2), R), obs = end, inds = c(3,4))
      samp <- rmvnorm(n = 1, mean = condNorm$means, sigma = condNorm$cov)
      while(norm(samp-circleCenter, type = "F") > circleRadius){
        samp <- rmvnorm(n = 1, mean = condNorm$means, sigma = condNorm$cov)
      }
      positions[indsToSample,] <- samp
    }
    sampledInds[indsToSample] <- T
  }
  
  #draw circle
  if(k==1){
    plot(x = circleCenter[1], y = circleCenter[2], pch = 3, xlim = c(-circleRadius, circleRadius), 
         ylim = c(-circleRadius, circleRadius), xlab = "", ylab = "")
    circX <- (-100:100/100*circleRadius)
    circY <- sqrt(circleRadius^2 - circX^2)
    lines(circX, circY, lwd = 2); lines(circX, -circY, lwd = 2)
  }
  
  #draw sampled path
  points(rbind(start, end), col = c(3,2), lwd = 2, pch = 19) #green marks start, red marks end, 
  lines(positions, col = rgb(0,0,0,mean(c(rep(1,10),rep(30/nhists, nhists)))))

}


#in two circles, connected by a rectangle
nhists <- 1e2
for(k in 1:nhists){
  print(k)
  
  #specify bounds of truncation, in this case: two circles and a rectangle
  circle1Radius <- 8
  circle1Center <- c(-11,0)
  
  circle2Radius <- 8
  circle2Center <- c(11,0)
  
  rectangleCorners <- rbind(c(-3.3,-1), c(3.3,1))
  
  #function to check if point is within these bounds
  withinBoundaries <- function(point, circle1Radius, circle1Center, circle2Radius, circle2Center, rectangleCorners){
    (norm(point-circle1Center, type = "F") < circle1Radius) |
    (norm(point-circle2Center, type = "F") < circle2Radius) |
    ((point[1] > rectangleCorners[1,1]) & 
     (point[1] < rectangleCorners[2,1]) &
     (point[2] > rectangleCorners[1,2]) & 
     (point[2] < rectangleCorners[2,2]))   
  }

  
  #specify starting and ending positions
  start <- c(-12,-6)
  end <- c(14,4)
  
  #specify brownian motion rate matrix
  R <- matrix(c(2,1,1,3), ncol = 2)
  
  #specify branch length separating start and end points
  BL <- 40
  
  #specify the number of intermediate timepoints over which to sample history
  n_IntermediateTimePoints <- 2^9-1
  
  #initialize matrix that will contain the bivariate history 
  positions <- matrix(0, ncol = 2, nrow = n_IntermediateTimePoints + 2)
  #specify the number of locations in total
  n_locations <- dim(positions)[1]
  #and specify start and end points
  positions[1,] <- start; positions[n_locations,] <- end
  
  #specify the position of each timepoint, where the starting branchlength is 0
  BLs <- 0:(n_IntermediateTimePoints+1)/(n_IntermediateTimePoints+1)*BL
  
  
  
  #keep track of which locations have observations already
  sampledInds <- rep(F, n_locations); sampledInds[1] <- TRUE; sampledInds[n_locations] <- TRUE
  
  #start sampling locations in log_2 stages, since under mvBM we can use already-sampled intermediate locations
  #to sample further intermediate locations
  for(i in 1:(log2(n_IntermediateTimePoints+1))){
    #which indices do we want to sample in this stage?
    indsToSample <- seq(from = 1, to = n_locations, length.out = 2^i+1)[-c(1,2^i+1)]
    if(i != 1){
      #make sure we haven't already sampled them!
      indsToSample <- setdiff(indsToSample, which(sampledInds))
      #and then iterate over unsampled indices
      for(j in 1:length(indsToSample)){
        indBeingSampled <- indsToSample[j]
        #find the two closest *already sampled* indices -- since mvBM is Markovian and evolution on branch i
        #independent of evolution on branch j, one can find conditionals only based on the two closest points 
        dists <- abs(which(sampledInds) - indBeingSampled)
        startEndInds <- which(sampledInds)[which(dists == min(dists))]
        strt <- positions[startEndInds[1],]
        ennd <- positions[startEndInds[2],]
        condNorm <- condMVN(means = rep(strt, 2), 
                            cov = kronecker(matrix(c(rep(BLs[indBeingSampled] - BLs[startEndInds[1]], 3), 
                                                     BLs[startEndInds[2]] - BLs[startEndInds[1]]), ncol = 2), R), obs = ennd, inds = c(3,4))
        samp <- rmvnorm(n = 1, mean = condNorm$means, sigma = condNorm$cov)
        #use rejection sampling to make sure we're in the allowable space
        while(!withinBoundaries(point = samp, circle1Radius, circle1Center, circle2Radius, circle2Center, rectangleCorners)){
          samp <- rmvnorm(n = 1, mean = condNorm$means, sigma = condNorm$cov)
        }
        positions[indBeingSampled,] <- samp
      }
    } else {
      #special case when sampling first value
      condNorm <- condMVN(means = rep(start, 2), cov = kronecker(matrix(c(rep(BLs[indsToSample], 3), 
                                                                          BLs[n_locations]), ncol = 2), R), obs = end, inds = c(3,4))
      samp <- rmvnorm(n = 1, mean = condNorm$means, sigma = condNorm$cov)
      while(!withinBoundaries(point = samp, circle1Radius, circle1Center, circle2Radius, circle2Center, rectangleCorners)){
        samp <- rmvnorm(n = 1, mean = condNorm$means, sigma = condNorm$cov)
      }
      positions[indsToSample,] <- samp
    }
    sampledInds[indsToSample] <- T
  }
  
  #draw circle
  if(k==1){
    plot(x = circle1Center[1], y = circle1Center[2], pch = 3, xlim = c(circle1Center[1]-circle1Radius, circle2Center[1]+circle2Radius), 
         ylim = c(circle1Center[2]-circle1Radius, circle1Center[2]+circle1Radius), xlab = "", ylab = "")
    points(x = circle2Center[1], y = circle2Center[2], pch = 3)
    
    #circle1
    circ1X <- (-100:100/100*circle1Radius)
    circ1Y <- sqrt(circle1Radius^2 - circ1X^2)
    lines(circ1X + circle1Center[1], circ1Y + circle1Center[2], lwd = 2); lines(circ1X + circle1Center[1], -circ1Y + circle1Center[2], lwd = 2)
    
    #circle2
    circ1X <- (-100:100/100*circle2Radius)
    circ1Y <- sqrt(circle2Radius^2 - circ1X^2)
    lines(circ1X + circle2Center[1], circ1Y + circle2Center[2], lwd = 2); lines(circ1X + circle2Center[1], -circ1Y + circle2Center[2], lwd = 2)
    
    #rectangle
    rect(xleft = rectangleCorners[1,1], ybottom = rectangleCorners[1,2], xright = rectangleCorners[2,1], ytop = rectangleCorners[2,2], lwd = 2)
    rect(border = NA, col = "white", xleft = rectangleCorners[1,1]-0.1, ybottom = rectangleCorners[1,2]+0.1, xright = rectangleCorners[2,1]+0.1, ytop = rectangleCorners[2,2]-0.1, lwd = 2)
  }
  
  #draw sampled path
  points(rbind(start, end), col = c(3,2), lwd = 2, pch = 19) #green marks start, red marks end, 
  lines(positions, col = rgb(0,0,0,mean(c(rep(1,5),rep(10/nhists, nhists)))))
  
}
