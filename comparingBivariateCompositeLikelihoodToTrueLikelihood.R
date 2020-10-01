setwd("/Volumes/1TB/Bailey/")

library(mvtnorm)
library(pbivnorm)

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

detbivcor <- function(rhos){
  sapply(rhos, function(rho) det(diag(2)*(1-rho) + rho))
}

pbivnorm_infs <- function(upper1, upper2, rhos){
  nInf <- (upper1 == -Inf | upper2 == -Inf)
  Inf1 <- upper1 == Inf
  Inf2 <- upper2 == Inf
  n2_1 <- Inf1 & !Inf2
  n1_2 <- Inf2 & !Inf1
  include <- !(nInf | Inf1 | Inf2)
  
  results <- rep(0, length(include))
  results[Inf1 & Inf2] <- 1
  results[n2_1] <- pnorm(q = upper2[n2_1])
  results[n1_2] <- pnorm(q = upper1[n1_2])
  results[include] <- pbivnorm(upper1[include], upper2[include], rhos[include])
  return(results)
} #checks for buggy edge cases and corrects

compareCompVsFull <- function(dimr, cor = NA, meanl = -1, meanu = 1, lower = NA, upper = NA, uniformWindow = F, unifL = c(-1,0), unifU = c(0,1), repError = F, applyDetCorr = F){
  if(all(is.na(cor))){
    x <- rlkj(dimr)
  } else {
    x <- cor
    dimr <- dim(x)[1]
  }
  
  if(!uniformWindow){
    if(all(is.na(lower))){
      lower <- rnorm(dimr, mean = meanl)
    }
    if(all(is.na(upper))){
      upper <- rnorm(dimr, mean = meanu)
    }
    swapi <- lower > upper
    swap <- upper[swapi]
    upper[swapi] <- lower[swapi]
    lower[swapi] <- swap
  }
  
  if(uniformWindow){
    lower <- runif(n = dimr, unifL[1], unifL[2])
    upper <- lower + runif(n = dimr, unifU[1], unifU[2])
  }
  
  trueprob <- (pmvnorm(lower = lower, upper = upper, corr = x, algorithm = GenzBretz()))
  absError <- attr(trueprob, "error")
  errorMag <- absError / trueprob[1]
  trueprob <- log(trueprob[1])
  
  uppers <- matrix(0, choose(dimr, 2), 2)
  lowers <- matrix(0, choose(dimr, 2), 2)
  rhos <- rep(0, choose(dimr, 2))
  counter <- 0
  for(i in 1:(dimr-1)){
    for(j in (i+1):dimr){
      counter <- counter + 1
      uppers[counter,] <- c(upper[i], upper[j])
      lowers[counter,] <- c(lower[i], lower[j])
      rhos[counter] <- x[i,j]
    }
  }
  
  probs <- pbivnorm_infs(uppers[,1], uppers[,2], rhos) - 
    pbivnorm_infs(lowers[,1], uppers[,2], rhos) - 
    pbivnorm_infs(uppers[,1], lowers[,2], rhos) + 
    pbivnorm_infs(lowers[,1], lowers[,2], rhos)
  
  probs[probs < 0] <- 0 #get rid of floating point errors
  if(applyDetCorr){
    probs <- probs / sqrt(detbivcor(x[upper.tri(x)]))
  }
  prob1 <- mean(log(probs)) #geometric mean
  # prob2 <- log(1 / (mean(1 / (probs)))) #harmonic mean
  # prob3 <- log(mean(probs)) #arithmetic mean
  prob <- prob1
  prob <- prob * (dimr) / 2
  if(applyDetCorr){
    prob <- prob + 0.5 * log(det(x))
  }
  if(!repError){
    return(c(trueprob, prob))
  } else if(repError){
    return(c(trueprob, prob, errorMag, absError))
  }
}

r2_1to1 <- function(var1, var2, uplow = F){
  if(!uplow){return(1 - (var(var1 - var2) / var(var1)))}
  if(uplow){
    tot <- 1 - (var(var1 - var2) / var(var1))
    low <- var1 > var2
    up <- 1 - (var(var1[!low] - var2[!low]) / var(var1[!low]))
    low <- 1 - (var(var1[low] - var2[low]) / var(var1[low]))
    return(c(tot, up, low))
  }
}


# R <- rlkj(10)
# (det(R))  / prod(detbivcor(R[upper.tri(R)]))

dimr <- 118
nrep <- 1000
error_threshold <- 10000
comparingTrueVsComp <- t(replicate(nrep, compareCompVsFull(dimr, cor = NA, lower = runif(dimr, -1, 1), upper = runif(dimr, -1, 1), repError = T, applyDetCorr = F)))
# comparingTrueVsComp <- t(replicate(nrep, compareCompVsFull(dimr, cor = diag(dimr) + 0.9 - diag(dimr)*0.9, lower = runif(dimr, -1, 1), upper = runif(dimr, -1, 1), repError = T, applyDetCorr = F)))
comparingTrueVsComp_orig <- comparingTrueVsComp
if(any(comparingTrueVsComp[,3] > error_threshold)){
  comparingTrueVsComp <- comparingTrueVsComp[-which(comparingTrueVsComp[,3] > error_threshold),]
}
if(any(is.na(comparingTrueVsComp[,3]))){
  comparingTrueVsComp <- comparingTrueVsComp[-which(is.na(comparingTrueVsComp[,3])),]
}
# comparingTrueVsComp[,2] <- comparingTrueVsComp[,2] / dimr
# lm(comparingTrueVsComp[,1] ~ comparingTrueVsComp[,2])
cor(comparingTrueVsComp[,1], comparingTrueVsComp[,2])^2
r2_1to1(comparingTrueVsComp[,1], comparingTrueVsComp[,2], uplow = T)
# lm(exp(comparingTrueVsComp[,1]) ~ exp(comparingTrueVsComp[,2]))

par(mar = c(5,5,3,3))
plot(comparingTrueVsComp[,1:2], xlab = "\"full\" multivariate normal integral", ylab = "composite integral obtained from \nall pairwise bivariate integrals")
title(paste0("dimension = ", dimr, ", r = ", round(cor(comparingTrueVsComp)[1,2], 3)))
abline(0, 1, lwd = 2, col = 2)
legend(x = "topleft", legend = "1-to-1 line", col = 2, lwd = 2, lty = 1)
for(i in 1:nrow(comparingTrueVsComp)){
  segments(y0 = comparingTrueVsComp[i,2], y1 = comparingTrueVsComp[i,2], 
           x0 = log(exp(comparingTrueVsComp[i,1] - comparingTrueVsComp[i,4])), x1 = log(exp(comparingTrueVsComp[i,1] + comparingTrueVsComp[i,4])))
}


##############################
#### making a nice figure ####
##############################
load(file = "comparingTrueVsComp_118d0.9r ")
comparingTrueVsComp09r <- comparingTrueVsComp  
load(file = "comparingTrueVsComp_118dLKJr")
comparingTrueVsCompLKJr <- comparingTrueVsComp  

fig_label <- function(text, region="figure", pos="topleft", cex=NULL, shrinkX = 1, shrinkY = 1, ...) {
  
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

#figure 4.2
grDevices::cairo_pdf(filename = "/Volumes/macOS/Users/nikolai/dissertation/figures/chpt4_figure2.pdf", width = 1500 / 72, height = 750 / 72)
par(mfrow = c(1,2), mar = c(5,7,3,1))
plot(comparingTrueVsCompLKJr[,1:2], xlab = "\"full\" multivariate normal integral", ylab = "composite integral obtained from \nall pairwise bivariate integrals",
     xlim = c(-350, -150), cex.axis = 1.5, pch = 16, col = rgb(0,0,0,0.5), cex = 2.5, cex.lab = 2)
title(main =  latex2exp::TeX("Approximate Integral Comparison: LKJ($\\eta$ = 1) Correlations"), cex.main = 2)
# text(paste0("dimension = 118, r = ", round(cor(comparingTrueVsCompLKJr)[1,2], 3)), x = -319, y = -177, cex = 1.5)
fig_label(paste0("dimension = 118, r = ", round(cor(comparingTrueVsCompLKJr)[1,2], 3)), cex = 1.5, region = "plot", shrinkX = 1.01, shrinkY = 1.005)
abline(0, 1, lwd = 3, col = 2)
legend(x = "bottomright", legend = "1-to-1 line", col = 2, lwd = 3, lty = 1, box.lty = 2, box.lwd = 2, cex = 2)
fig_label("a)", shrinkX = 0.992, shrinkY = 1.005, cex = 2.25)
box(which = "figure", lty = 2, col = rgb(0,0,0,0.5))

plot(comparingTrueVsComp09r[,1:2], xlab = "\"full\" multivariate normal integral", ylab = "composite integral obtained from \nall pairwise bivariate integrals",
     xlim = c(-240, -120), cex.axis = 1.5, pch = 16, col = rgb(0,0,0,0.5), cex = 2.5, cex.lab = 2)
title(main =  latex2exp::TeX("Approximate Integral Comparison: r$_{ij}$ = 0.9 Correlations"), cex.main = 2)
# text(paste0("dimension = 118, r = ", round(cor(comparingTrueVsComp09r)[1,2], 3)), x = -217, y = -167, cex = 1.5)
fig_label(paste0("dimension = 118, r = ", round(cor(comparingTrueVsComp09r)[1,2], 3)), cex = 1.5, region = "plot", shrinkX = 1.01, shrinkY = 1.0075)
abline(0, 1, lwd = 3, col = 2)
legend(x = "bottomright", legend = "1-to-1 line", col = 2, lwd = 3, lty = 1, box.lty = 2, box.lwd = 2, cex = 2)
box(which = "figure", lty = 2, col = rgb(0,0,0,0.5))
fig_label("b)", shrinkX = 0.9925, shrinkY = 1.0075, cex = 2.25)
dev.off()

#approximation will underestimate how tiny very low full mvn probs are
#but this might also be running up against the limits of GenzBretz' abilities

dims <- 2:150
r2s <- matrix(0, nrow = c(length(dims)), ncol = 4)
nrep <- 2000
error_magnitude_threshold <- 0.25

for(dimr in dims){
  cat(paste0(dimr, " "))
  comparingTrueVsComp <- t(replicate(nrep, compareCompVsFull(dimr, cor = NA, uniformWindow = T, unifL = c(-1,0), unifU = c(0,1), repError = T)))
  if(length(which(comparingTrueVsComp == -Inf, arr.ind = T)[,1]) > 0){
    comparingTrueVsComp <- comparingTrueVsComp[-which(comparingTrueVsComp == -Inf, arr.ind = T)[,1],] #remove -Inf errors
  }
  if(length(which(comparingTrueVsComp[,3] > error_magnitude_threshold)) > 0){
    comparingTrueVsComp <- comparingTrueVsComp[-which(comparingTrueVsComp[,3] > error_magnitude_threshold),] #remove high error estimates for full multiv
  }
  if(length(which(comparingTrueVsComp[,3] == error_magnitude_threshold)) > 0){
    comparingTrueVsComp <- comparingTrueVsComp[-which(comparingTrueVsComp[,3] == error_magnitude_threshold),] #remove high error estimates for full multiv
  }
  
  r2s[dimr-1,] <- c(cor(comparingTrueVsComp[,1], comparingTrueVsComp[,2])^2, r2_1to1(comparingTrueVsComp[,1], comparingTrueVsComp[,2], uplow = T))
}

plot(r2s[,3])

# bigdiff <- F
# while(!bigdiff){
#   r <- rlkj(dimr)
#   lower <- runif(n = dimr, unifL[1], unifL[2])
#   upper <- lower + runif(n = dimr, unifU[1], unifU[2])
#   vals <- compareCompVsFull(cor = r, lower = lower, upper = upper) #oh wow the full one is super unstable!
#   if(diff(vals) > 50){
#     bigdiff = T
#   }
# }

# dev.off()


dimrange <- 1:10
intercepts <- rep(0, length(dimrange))
slopes <- rep(0, length(dimrange))
for(dimr in dimrange){
  print(dimr)
  nrep <- 10000
  comparingTrueVsComp <- t(replicate(nrep, compareCompVsFull(dimr*10, cor = NA)))
  lmx <- lm(comparingTrueVsComp[,1] ~ comparingTrueVsComp[,2])
  intercepts[dimr] <- lmx[1]$coefficients[1]
  slopes[dimr] <- lmx[1]$coefficients[2]
}
par(mfrow = c(1,2))
plot(dimrange*10, slopes, xlab = "dimension", ylab = "least squares estimate of slope for \n\"full mvn integral ~ mean pairwise bivn integral\"", main = "slope")
tdmr <- dimrange*10
lm(slopes ~ tdmr)
mean(slopes / (dimrange*10))
plot(dimrange*10, intercepts, xlab = "dimension", ylab = "least squares estimate of intercept for \n\"full mvn integral ~ mean pairwise bivn integral\"", main = "intercept")
lm(intercepts ~ tdmr)
(intercepts / (dimrange*2))
diff(intercepts)

# 
# par(mfrow = c(1,1))
# dimr <- 2:20
# corcof <- 0.2
# vals <- sapply(dimr, function(x) compareCompVsFull(dimr = x, cor = (diag(x)*(1-corcof) + corcof), lower = rep(-0.5, x), upper = rep(0.5, x)))
# plot(dimr, vals[1,]); abline()
# dimr2 <- dimr - 2
# lm(vals[1,] ~ (dimr2))
# 
# 
# plot(log(sapply(2:50, function(x) det(diag(x)*(1-corcof) + corcof))))
# lm(log(sapply(dimr, function(x) det(diag(x)*(1-corcof) + corcof)))* 0.5 ~ dimr)
# 
# plot(t(sapply(1:9/10, function(x) compareCompVsFull(dimr = 20, cor = (diag(20)*(1-x) + x), lower = rep(-0.5, 20), upper = rep(0.5, 20)))))
# 
# f <- function(){
#   R <- rlkj(20)
#   c(log(det(R)), sum(log(detbivcor(R[upper.tri(R)]))))
# }
# 
# plot(t(replicate(100, f())))

#does the bivariate probability mirror the univariate?
nrep = 100
vals <- matrix(data = 0, nrow = nrep, ncol = 3)
low_extra <- runif(5, -0.5, 0.5)
high_extra <- runif(5, 0.5, 1.5)
for(i in 1:nrep){
  corr_param <- 0.75
  low = c(runif(1, -0.5, 0.5), 0.5, low_extra)
  high = c(runif(1, 0.5, 1.5), 1.1, high_extra)
  dimr <- length(low)
  R <- diag(c(rep(1-corr_param, dimr))) + corr_param
  R <- rlkj(dimr)
  vals[i,1] <- compareCompVsFull(dimr= 2, cor = R, lower = low, upper = high)[1]
  vals[i,2] <- compareCompVsFull(dimr= 2, cor = R, lower = low, upper = high)[2]
  vals[i,3] <- log(pnorm(q = high[1]) - pnorm(q = low[1]))
}
pairs(vals)
