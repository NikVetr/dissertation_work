#load libraries
library(foreach)
library(doParallel)
library(parallel)

#specify functions
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

moments_from_grid <- function(x, d, m = c(1:2), ll = F){
  
  #prevent underflow
  if(ll){
    d <- d - max(d) + 10
    d <- exp(d)
  } 
  
  m1 <- sum(x*d) / sum(d)
  m_eval <- setdiff(m, 1)
  ms <- sapply(m_eval, function(m_i) sum((x-m1)^m_i*d) / sum(d))
  if(1 %in% m){return(unlist(c(m1, ms)))}else{return(ms)}  
  
}

bvn_posterior_dist <- function(n, r, corr_val_interval = 2E-3){

  #set sim params
  p <- 2
  x <- matrix(rnorm(n*p), n, p)
  R <- diag(p)*(1-r) + r
  
  #transform data to specified sample moments
  L_old <- t(chol(cov(x)))
  L_new <- t(chol(R))
  x <- t(L_new %*% solve(L_old) %*% t(x))
  x <- x - matrix(apply(x, 2, mean), n, p, T)
  
  #compute likelihood surface over correlation from -1 to 1 & normalize
  corr_vals <- seq(-1+corr_val_interval,1-corr_val_interval,by=corr_val_interval)
  lls <- sapply(corr_vals, function(ri) sum(mvtnorm::dmvnorm(x = x, mean = rep(0,p), sigma = diag(p)*(1-ri) + ri, log = T)))
  # lls <- lls - log(sum(exp(lls)*rep(corr_val_interval, length(lls))))
  lls <- lls - (log(sum(exp(lls - max(lls) + 10)*rep(corr_val_interval, length(lls)))) - 10 + max(lls))
  
  
  #get some useful quantities
  posterior_mode <- corr_vals[which.max(lls)]
  posterior_moments <- moments_from_grid(corr_vals, lls, ll = T)
  posterior_mean <- posterior_moments[1]
  posterior_var <- posterior_moments[2]
  posterior_mean_dens <- exp(lls[which.min(abs(corr_vals - posterior_mean))])
  
  return(list(corr_vals = corr_vals, 
              lls = lls, 
              posterior_mode = posterior_mode,
              posterior_mean = posterior_mean,
              posterior_var = posterior_var,
              posterior_mean_dens = posterior_mean_dens))

}

cvi = 2E-4
dens_info <- bvn_posterior_dist(n = 3, r = 0.4, corr_val_interval = cvi)
inds_to_plot <- which(dens_info$corr_vals > -1 & dens_info$corr_vals < 1)
plot(dens_info$corr_vals[inds_to_plot], exp(dens_info$lls[inds_to_plot]), type = "l", xlab = "correlation", ylab = "prior density")
polygon(x = c(dens_info$corr_vals[inds_to_plot], rev(dens_info$corr_vals[inds_to_plot])), 
        y = c(exp(dens_info$lls[inds_to_plot]), rep(0,length(exp(dens_info$lls[inds_to_plot])))),
        col = adjustcolor("#091F92", 0.25))
segments(x0 = dens_info$posterior_mode, x1 = dens_info$posterior_mode, y0 = 0, 
         y1 = exp(max(dens_info$lls)), lwd = 2, col = "#170a47")
segments(x0 = dens_info$posterior_mean, x1 = dens_info$posterior_mean, y0 = 0, 
         y1 = dens_info$posterior_mean_dens, lwd = 2, col = "#3964C3")

#make a quick animation
n_vals <- 3:250
r = 0.4
dens_infos <- mclapply(n_vals, function(n) bvn_posterior_dist(n = n, r = r, corr_val_interval = cvi), mc.cores = 12)
lls_mat <- t(sapply(n_vals, function(n) dens_infos[[n-min(n_vals)+1]]$lls))
posterior_modes <- sapply(n_vals, function(n) dens_infos[[n-min(n_vals)+1]]$posterior_mode)
posterior_vars <- sapply(n_vals, function(n) dens_infos[[n-min(n_vals)+1]]$posterior_var)
posterior_means <- sapply(n_vals, function(n) dens_infos[[n-min(n_vals)+1]]$posterior_mean)
posterior_mean_densities <- sapply(n_vals, function(n) dens_infos[[n-min(n_vals)+1]]$posterior_mean_dens)

#### do the plotting ####

#test out basic plot
n = 250
layout(matrix(c(1,1,2,1,1,3,1,1,4), ncol = 3, nrow = 3, byrow = T))
par(mar = c(4.5,5,2,2))
plot(dens_infos[[n-min(n_vals)+1]]$corr_vals[inds_to_plot], exp(dens_infos[[n-min(n_vals)+1]]$lls[inds_to_plot]), 
     type = "l", xlab = "correlation", ylab = "prior density",
     xlim = c(-1,1), ylim = c(0, max(exp(lls_mat))))
polygon(x = c(dens_infos[[n-min(n_vals)+1]]$corr_vals[inds_to_plot], rev(dens_infos[[n-min(n_vals)+1]]$corr_vals[inds_to_plot])), 
        y = c(exp(dens_infos[[n-min(n_vals)+1]]$lls[inds_to_plot]), rep(0,length(exp(dens_infos[[n-min(n_vals)+1]]$lls[inds_to_plot])))),
        col = adjustcolor("#091F92", 0.25))
segments(x0 = dens_infos[[n-min(n_vals)+1]]$posterior_mode, x1 = dens_infos[[n-min(n_vals)+1]]$posterior_mode, y0 = 0, 
         y1 = exp(max(dens_infos[[n-min(n_vals)+1]]$lls)), lwd = 2, col = "#170a47")
segments(x0 = dens_infos[[n-min(n_vals)+1]]$posterior_mean, x1 = dens_infos[[n-min(n_vals)+1]]$posterior_mean, y0 = 0, 
         y1 = dens_infos[[n-min(n_vals)+1]]$posterior_mean_dens, lwd = 2, col = "#3964C3")
text(x = -1.075, y = max(exp(lls_mat)), col = "#b36200", cex = 1.4, pos = 4,
     label = paste0("Sample Size = ", n, ", Correlation = ", r))

par(mar = c(3.5,4,2,2))
plot(n_vals, posterior_modes, type = "l", lwd = 2, col = "#170a47",
     ylab = "prior mode", xlab = "sample size", xpd = NA, cex.axis = 0.75)
points(x = n, y = dens_infos[[n-min(n_vals)+1]]$posterior_mode, pch = 21, bg = "#dc7800", cex = 1.5)

par(mar = c(4,4,1.5,2))
plot(n_vals, posterior_means, type = "l", lwd = 2, col = "#3964C3",
     ylab = "prior mean", xlab = "sample size", xpd = NA, cex.axis = 0.75)
points(x = n, y = dens_infos[[n-min(n_vals)+1]]$posterior_mean, pch = 21, bg = "#dc7800", cex = 1.5)

par(mar = c(4.5,4,1,2))
plot(n_vals, posterior_vars, type = "l", lwd = 2, col = "#3964C3",
     ylab = "prior mean", xlab = "sample size", xpd = NA, log = 'y', cex.axis = 0.75)
points(x = n, y = dens_infos[[n-min(n_vals)+1]]$posterior_var, pch = 21, bg = "#dc7800", cex = 1.5)



#####
#now iterate through all these plots
dir.create(path = "~/Documents/Documents - nikolai/bivnorm_posterior_animation")
file.remove(paste0("~/Documents/Documents - nikolai/bivnorm_posterior_animation/frames/", 
  list.files("~/Documents/Documents - nikolai/bivnorm_posterior_animation/frames", pattern = "*.png")))
dir.create(path = "~/Documents/Documents - nikolai/bivnorm_posterior_animation/frames")

n_frames <- 600
rate = 4
ns <- cumsum(c(0, rev(diff(pexp(q = seq(0,0.75,length.out = n_frames), rate = rate))))) / pexp(0.75, rate = rate) * diff(range(n_vals)) + min(n_vals)
diff(head(ns))

if(!exists("cl")){
  cl <- makeCluster(12, outfile="")
  registerDoParallel(cl)
}
getDoParWorkers()

foreach(ni=(1:n_frames)) %dopar% {
# for(ni in (1:n_frames)){
  cat(paste0(ni, " "))
  
  #retrieve parameters
  if(abs(ns[ni] %% 1) < 1E-6){
    n <- round(ns[ni])
    
    corr_vals <- dens_infos[[n-min(n_vals)+1]]$corr_vals[inds_to_plot]
    lls <- dens_infos[[n-min(n_vals)+1]]$lls[inds_to_plot]
    posterior_mode <- dens_infos[[n-min(n_vals)+1]]$posterior_mode
    posterior_mean <- dens_infos[[n-min(n_vals)+1]]$posterior_mean
    posterior_mean_dens <- dens_infos[[n-min(n_vals)+1]]$posterior_mean_dens
    posterior_var <- dens_infos[[n-min(n_vals)+1]]$posterior_var
  
  } else {
    n <- ns[ni] 
    
    prop_ceiling <- ns[ni]%%1
    prop_floor <- 1-prop_ceiling
    corr_vals <- dens_infos[[floor(ns[ni])-min(n_vals)+1]]$corr_vals[inds_to_plot] * prop_floor + 
      dens_infos[[ceiling(ns[ni])-min(n_vals)+1]]$corr_vals[inds_to_plot] * prop_ceiling
    lls <- dens_infos[[floor(ns[ni])-min(n_vals)+1]]$lls[inds_to_plot] * prop_floor + 
      dens_infos[[ceiling(ns[ni])-min(n_vals)+1]]$lls[inds_to_plot] * prop_ceiling
    posterior_mode <- dens_infos[[floor(ns[ni])-min(n_vals)+1]]$posterior_mode * prop_floor + 
      dens_infos[[ceiling(ns[ni])-min(n_vals)+1]]$posterior_mode * prop_ceiling
    posterior_mean <- dens_infos[[floor(ns[ni])-min(n_vals)+1]]$posterior_mean * prop_floor + 
      dens_infos[[ceiling(ns[ni])-min(n_vals)+1]]$posterior_mean * prop_ceiling
    posterior_mean_dens <- dens_infos[[floor(ns[ni])-min(n_vals)+1]]$posterior_mean_dens * prop_floor + 
      dens_infos[[ceiling(ns[ni])-min(n_vals)+1]]$posterior_mean_dens * prop_ceiling
    posterior_var <- dens_infos[[floor(ns[ni])-min(n_vals)+1]]$posterior_var * prop_floor + 
      dens_infos[[ceiling(ns[ni])-min(n_vals)+1]]$posterior_var * prop_ceiling
  }
  
  png(filename = paste0("~/Documents/Documents - nikolai/bivnorm_posterior_animation/frames/", 
                        paste0(rep(0, 5-nchar(((ni - 1)) + 1)), collapse = ""), ((ni - 1)) + 1,".png"),
      width = 1500, height = 900, res = 200)
  
  ##PLOTTING CODE
  layout(matrix(c(1,1,2,1,1,3,1,1,4), ncol = 3, nrow = 3, byrow = T))
  par(mar = c(4.5,5,2,2))
  plot(corr_vals, exp(lls), 
       type = "l", xlab = "correlation", ylab = "marginal prior density",
       xlim = c(-1,1), ylim = c(0, max(exp(lls_mat))))
  polygon(x = c(corr_vals, rev(corr_vals)), 
          y = c(exp(lls), rep(0,length(lls))),
          col = adjustcolor("#091F92", 0.25))
  segments(x0 = posterior_mode, x1 = posterior_mode, y0 = 0, 
           y1 = exp(max(lls)), lwd = 2, col = "#170a47")
  segments(x0 = posterior_mean, x1 = posterior_mean, y0 = 0, 
           y1 = posterior_mean_dens, lwd = 2, col = "#3964C3")
  text(x = -1.075, y = max(exp(lls_mat)), col = "#b36200", cex = 1.4, pos = 4,
       label = paste0("Sample Size = ", round(ns[ni]), ", Correlation = ", r))
  
  par(mar = c(3.5,4,2,2))
  plot(n_vals, posterior_modes, type = "l", lwd = 2, col = "#170a47",
       ylab = "prior mode", xlab = "sample size", xpd = NA)
  points(x = ns[ni], y = posterior_mode, pch = 21, bg = "#dc7800", cex = 1.5)
  
  par(mar = c(4,4,1.5,2))
  plot(n_vals, posterior_means, type = "l", lwd = 2, col = "#3964C3",
       ylab = "prior mean", xlab = "sample size", xpd = NA)
  points(x = ns[ni], y = posterior_mean, pch = 21, bg = "#dc7800", cex = 1.5)
  
  par(mar = c(4.5,4,1,2))
  plot(n_vals, posterior_vars, type = "l", lwd = 2, col = "#530000",
       ylab = "prior variance", xlab = "sample size", xpd = NA)
  points(x = ns[ni], y = posterior_var, pch = 21, bg = "#dc7800", cex = 1.5)
  
  dev.off()
    
} 


#stich together animation
file.remove(paste0("~/Documents/Documents - nikolai/bivnorm_posterior_animation/", 
                   list.files("~/Documents/Documents - nikolai/bivnorm_posterior_animation", pattern = "*.mp4")))
system(paste0("cd \'Documents/Documents - nikolai/bivnorm_posterior_animation\'; ffmpeg -r ", 
              60," -f image2 -s 1500x900 -i frames/%05d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p forw_bivnorm_animation.mp4"))

#reverse and append and loop
system(paste0("cd \'Documents/Documents - nikolai/bivnorm_posterior_animation\'; ", 
             "ffmpeg -i forw_bivnorm_animation.mp4 -vf reverse rev_bivnorm_animation.mp4; ",
             "touch input.txt;",
             "echo \"file forw_bivnorm_animation.mp4\nfile rev_bivnorm_animation.mp4\" > input.txt;",
             "ffmpeg -f concat -i input.txt -codec copy once_bivnorm_animation.mp4; ",
             "ffmpeg -stream_loop 1 -i once_bivnorm_animation.mp4 -c copy bivnorm_animation.mp4"))


#### now let's try it in 3D! ####
#first confirm I know how to do this lol w/ a trivariate normal
p <- 3
r <- 0.4
R <- diag(p)*(1-r) + r
x <- seq(-5,5,length.out = 100)
xc <- expand.grid(list(x,x,x))
dens <- mvtnorm::dmvnorm(x = xc, sigma = R, log = T)
m_dens <- sapply(x, function(xi) {
  vals <- dens[xc[,1] == xi]
  log(mean(exp(vals - max(vals) + 10))) - 10 + max(vals) #underflow hack
})
m_dens <- m_dens - log(sum(exp(m_dens)*rep(diff(x[1:2]), length(m_dens))))
par(mfrow = c(1,2))
plot(x, exp(m_dens), type = "l")
lines(x, dnorm(x, log = F), type = "l", col = "green")
plot(m_dens, dnorm(x, log = T)); abline(0,1)
exp(m_dens) - dnorm(x, log = F) #discrete approximation error
#ok seems that I do, 

cvg_to_mat_3D <- function(cvg){
  cmat <- diag(3)
  cmat[1,2] <- cvg[1]
  cmat[1,3] <- cvg[2]
  cmat[2,3] <- cvg[3]
  cmat + t(cmat) - diag(3)
  
}

invlogit <- function(x) exp(x)/(1+exp(x))

logit <- function(p) log(p / (1-p))

tvn_posterior_dist <- function(n, r, corr_val_interval = 1E-1){
  
  #set sim params
  p <- 3
  x <- matrix(rnorm(n*p), n, p)
  R <- diag(p)*(1-r) + r
  
  #transform data to specified sample moments
  L_old <- t(chol(cov(x)))
  L_new <- t(chol(R))
  x <- t(L_new %*% solve(L_old) %*% t(x))
  x <- x - matrix(apply(x, 2, mean), n, p, T)
  
  #compute likelihood surface over correlation from -1 to 1 & normalize
  corr_vals <- seq(-1+corr_val_interval,1-corr_val_interval,by=corr_val_interval)
  corr_vals_grid <- expand.grid(list(corr_vals,corr_vals,corr_vals))
  lls <- sapply(1:nrow(corr_vals_grid), function(ri) 
    sum(mvtnorm::dmvnorm(x = x, mean = rep(0,p), sigma = cvg_to_mat_3D(unlist(corr_vals_grid[ri,])), log = T)))
  lls <- sapply(corr_vals, function(cv) {
    vals <- lls[corr_vals_grid[,1] == cv]
    log(mean(exp(vals - max(vals) + 10))) - 10 + max(vals) #underflow hack
  })
  
  # lls <- lls - log(sum(exp(lls)*rep(corr_val_interval, length(lls))))
  lls <- lls - (log(sum(exp(lls - max(lls) + 10)*rep(corr_val_interval, length(lls)))) - 10 + max(lls))
  
  #get some useful quantities
  posterior_mode <- corr_vals[which.max(lls)]
  posterior_moments <- moments_from_grid(corr_vals, lls, ll = T)
  posterior_mean <- posterior_moments[1]
  posterior_var <- posterior_moments[2]
  posterior_mean_dens <- exp(lls[which.min(abs(corr_vals - posterior_mean))])
  
  return(list(corr_vals = corr_vals, 
              lls = lls, 
              posterior_mode = posterior_mode,
              posterior_mean = posterior_mean,
              posterior_var = posterior_var,
              posterior_mean_dens = posterior_mean_dens))
  
}

plot(corr_vals, exp(tvn_posterior_dist(10, 0.4, 1E-1)$lls), type = "l")

#collect output for 
n_vals <- 4:10
r = 0.4
cvi = 1E-1
dens_infos <- mclapply(n_vals, function(n) bvn_posterior_dist(n = n, r = r, corr_val_interval = cvi), mc.cores = 12)
lls_mat <- t(sapply(n_vals, function(n) dens_infos[[n-min(n_vals)+1]]$lls))
posterior_modes <- sapply(n_vals, function(n) dens_infos[[n-min(n_vals)+1]]$posterior_mode)
posterior_vars <- sapply(n_vals, function(n) dens_infos[[n-min(n_vals)+1]]$posterior_var)
posterior_means <- sapply(n_vals, function(n) dens_infos[[n-min(n_vals)+1]]$posterior_mean)
posterior_mean_densities <- sapply(n_vals, function(n) dens_infos[[n-min(n_vals)+1]]$posterior_mean_dens)

dens_infos_tvn <- mclapply(n_vals, function(n) tvn_posterior_dist(n = n, r = r, corr_val_interval = cvi), mc.cores = 12)
lls_mat_tvn <- t(sapply(n_vals, function(n) dens_infos_tvn[[n-min(n_vals)+1]]$lls))
posterior_modes_tvn <- sapply(n_vals, function(n) dens_infos_tvn[[n-min(n_vals)+1]]$posterior_mode)
posterior_vars_tvn <- sapply(n_vals, function(n) dens_infos_tvn[[n-min(n_vals)+1]]$posterior_var)
posterior_means_tvn <- sapply(n_vals, function(n) dens_infos_tvn[[n-min(n_vals)+1]]$posterior_mean)
posterior_mean_densities_tvn <- sapply(n_vals, function(n) dens_infos_tvn[[n-min(n_vals)+1]]$posterior_mean_dens)

#### now do the plotting ####
#test out basic plot
n = 8
inds_to_plot <- which(dens_infos[[1]]$corr_vals > -1 & dens_infos[[1]]$corr_vals < 1)

layout(matrix(c(1,1,2,1,1,3,1,1,4), ncol = 3, nrow = 3, byrow = T))
par(mar = c(4.5,5,2,2))
plot(dens_infos[[n-min(n_vals)+1]]$corr_vals[inds_to_plot], exp(dens_infos[[n-min(n_vals)+1]]$lls[inds_to_plot]), 
     type = "l", xlab = "correlation", ylab = "prior density",
     xlim = c(-1,1), ylim = c(0, max(exp(c(lls_mat, lls_mat_tvn)))))

#lines
segments(x0 = dens_infos[[n-min(n_vals)+1]]$posterior_mode, x1 = dens_infos[[n-min(n_vals)+1]]$posterior_mode, y0 = 0, 
         y1 = exp(max(dens_infos[[n-min(n_vals)+1]]$lls)), lwd = 2, col = "#170a47")
segments(x0 = dens_infos[[n-min(n_vals)+1]]$posterior_mean, x1 = dens_infos[[n-min(n_vals)+1]]$posterior_mean, y0 = 0, 
         y1 = dens_infos[[n-min(n_vals)+1]]$posterior_mean_dens, lwd = 2, col = "#3964C3")
segments(x0 = dens_infos_tvn[[n-min(n_vals)+1]]$posterior_mode, x1 = dens_infos_tvn[[n-min(n_vals)+1]]$posterior_mode, y0 = 0, 
         y1 = exp(max(dens_infos_tvn[[n-min(n_vals)+1]]$lls)), lwd = 2, col = "#C83200")
segments(x0 = dens_infos_tvn[[n-min(n_vals)+1]]$posterior_mean, x1 = dens_infos_tvn[[n-min(n_vals)+1]]$posterior_mean, y0 = 0, 
         y1 = dens_infos_tvn[[n-min(n_vals)+1]]$posterior_mean_dens, lwd = 2, col = "#FF9912")

#polygons
polygon(x = c(dens_infos[[n-min(n_vals)+1]]$corr_vals[inds_to_plot], rev(dens_infos[[n-min(n_vals)+1]]$corr_vals[inds_to_plot])), 
        y = c(exp(dens_infos[[n-min(n_vals)+1]]$lls[inds_to_plot]), rep(0,length(exp(dens_infos[[n-min(n_vals)+1]]$lls[inds_to_plot])))),
        col = adjustcolor("#091F92", 0.25))
polygon(x = c(dens_infos_tvn[[n-min(n_vals)+1]]$corr_vals[inds_to_plot], rev(dens_infos_tvn[[n-min(n_vals)+1]]$corr_vals[inds_to_plot])), 
        y = c(exp(dens_infos_tvn[[n-min(n_vals)+1]]$lls[inds_to_plot]), rep(0,length(exp(dens_infos_tvn[[n-min(n_vals)+1]]$lls[inds_to_plot])))),
        col = adjustcolor("#C83200", 0.25))

#beta lkj marginal regularization
d=4
beta_ldensities_lkj1 <- dbeta(x = (dens_infos[[1]]$corr_vals + 1)/2, 1+(d-2)/2, 1+(d-2)/2, log = T) - log(2)
new_corr_dens <- dens_infos[[n-min(n_vals)+1]]$lls[inds_to_plot] + beta_ldensities_lkj1
new_corr_dens <- new_corr_dens - (log(sum(exp(new_corr_dens - max(new_corr_dens) + 10)*rep(corr_val_interval, length(new_corr_dens)))) - 10 + max(new_corr_dens))

polygon(x = c(dens_infos[[n-min(n_vals)+1]]$corr_vals[inds_to_plot], rev(dens_infos[[n-min(n_vals)+1]]$corr_vals[inds_to_plot])), 
        y = c(exp(new_corr_dens), rep(0,length(exp(dens_infos[[n-min(n_vals)+1]]$lls[inds_to_plot])))),
        col = adjustcolor("green", 0.25))

text(x = -1.075, y = max(exp(lls_mat)), col = "#b36200", cex = 1.4, pos = 4,
     label = paste0("Sample Size = ", n, ", Correlation = ", r))

par(mar = c(3.5,4,2,2))
plot(n_vals, posterior_modes, type = "l", lwd = 2, col = "#170a47",
     ylab = "prior mode", xlab = "sample size", xpd = NA, cex.axis = 0.75)
points(x = n, y = dens_infos[[n-min(n_vals)+1]]$posterior_mode, pch = 21, bg = "#dc7800", cex = 1.5)

par(mar = c(4,4,1.5,2))
plot(n_vals, posterior_means, type = "l", lwd = 2, col = "#3964C3",
     ylab = "prior mean", xlab = "sample size", xpd = NA, cex.axis = 0.75)
points(x = n, y = dens_infos[[n-min(n_vals)+1]]$posterior_mean, pch = 21, bg = "#dc7800", cex = 1.5)

par(mar = c(4.5,4,1,2))
plot(n_vals, posterior_vars, type = "l", lwd = 2, col = "#3964C3",
     ylab = "prior mean", xlab = "sample size", xpd = NA, log = 'y', cex.axis = 0.75)
points(x = n, y = dens_infos[[n-min(n_vals)+1]]$posterior_var, pch = 21, bg = "#dc7800", cex = 1.5)


#####

#set sim params
foo <- function(n){
  p <- 3
  x <- matrix(rnorm(n*p), n, p)
  R <- diag(p)*(1-r) + r
  
  #transform data to specified sample moments
  L_old <- t(chol(cov(x)))
  L_new <- t(chol(R))
  x <- t(L_new %*% solve(L_old) %*% t(x))
  x <- x - matrix(apply(x, 2, mean), n, p, T)
  
  R2 <- diag(p)*(1-0.7) + 0.7
  
  
  dens1 <- sum(mvtnorm::dmvnorm(x, sigma = R2, log = T))
  dens1
}
#check data duplication
r = 0.4
n = 10
p <- 5
x <- matrix(rnorm(n*p), n, p)
R <- diag(p)*(1-r) + r

#transform data to specified sample moments
L_old <- t(chol(cov(x)))
L_new <- t(chol(R))
x <- t(L_new %*% solve(L_old) %*% t(x))
x <- x - matrix(apply(x, 2, mean), n, p, T)

cov(x)
nreps <- 1:100
var_reduct <- sapply(nreps, function(nrep) cov(do.call(rbind, lapply(1:nrep,function(foo)x)))[1,1])
MASS::fractions(var_reduct)
plot(nreps, var_reduct, type = "l")

a = 0.1
R2 <- diag(p)*(1-0.7) + 0.7
# plot(sapply(1:20, function(i) dnorm(x = a*i, log = T)), type = "l")
plot((1:20)^2 - 1, sapply(1:20, function(i) sum(mvtnorm::dmvnorm(x*i, sigma = R2, log = T))), type = "l")
(lm(ll ~ tdim, 
   data = data.frame(list(tdim = (1:20)^2 - 1, ll = sapply(1:20, function(i) sum(mvtnorm::dmvnorm(x*i, sigma = R2, log = T)))))))
b = sum(mvtnorm::dmvnorm(x, sigma = R2, log = T))+1/2*log(det(2*pi*R2))*n
a = sum(mvtnorm::dmvnorm(x, sigma = R2, log = T))
i = 17
a
b

sum(mvtnorm::dmvnorm(x*i, sigma = R2, log = T))
a + (i^2 -1)* b



i = 0.5
sum(mvtnorm::dmvnorm(x*i, sigma = R2, log = T))
ldens <- sum(mvtnorm::dmvnorm(x, sigma = R2, log = T))
ldens + (ldens+1/2*log(det(2*pi*R2))*n) * (i^2 - 1)

cov_scale <- function(c,n) c*(n-1) / (n*c-1)

#set sample parameters
r = 0.4
p <- 5
R <- diag(p)*(1-r) + r
L_new <- t(chol(R))

#simulate low sample data
n1 = 7
x1 <- matrix(rnorm(n1*p), n1, p)
L_old_1 <- t(chol(cov(x1)))
x1 <- t(L_new %*% solve(L_old_1) %*% t(x1))
x1 <- x1 - matrix(apply(x1, 2, mean), n1, p, T)

#simulate high sample data
n2 = 111
x2 <- matrix(rnorm(n2*p), n2, p)
L_old_2 <- t(chol(cov(x2)))
x2 <- t(L_new %*% solve(L_old_2) %*% t(x2))
x2 <- x2 - matrix(apply(x2, 2, mean), n2, p, T)

#important
c <- n2 / n1

#evaluate log-density
R2 <- rlkj(p, 1)
#high-sample
sum(mvtnorm::dmvnorm(x2, sigma = R2, log = T))
#low-sample data duplication
sum(mvtnorm::dmvnorm(x1 /sqrt(c*(n1-1) / (n1*c-1)), sigma = R2, log = T)) * c
#alternatively
a = sum(mvtnorm::dmvnorm(x1, sigma = R2, log = T)) * c
log_det_bit <- log(det(2*pi*R2))
c <- n2 / n1
i = 1 / (c*(n1-1) / (n1*c-1))
b = a+1/2*log_det_bit*n2
a + (i-1)* b

#alternatively
# a = sum(mvtnorm::dmvnorm(x2, sigma = R2, log = T))
# b = a+1/2*log(det(2*pi*R2))*n2
# i = (c*(n1-1) / (n1*c-1))
# a
# b
# sum(mvtnorm::dmvnorm(x2*i, sigma = R2, log = T))
# a + (i^2-1)* b


#benchmark
microbenchmark::microbenchmark(
  sum(mvtnorm::dmvnorm(x2, sigma = R2, log = T)),
  sum(mvtnorm::dmvnorm(x1 /sqrt(c*(n1-1) / (n1*c-1)), sigma = R2, log = T)) * c,
  {c <- n2 / n1
  i = 1 / (c*(n1-1) / (n1*c-1))
  b = a+1/2*log_det_bit*n2
  a + (i-1)* b}
)

#possible to calculate just from sample cov?
x <- x1
n <- n1
R <- R2
sum(mvtnorm::dmvnorm(x, sigma = R, log = T))
-1/2*log(det(2*pi*R))*n - 1/2*sum(sapply(1:n, function(i) t(x[i,])%*%solve(R)%*%(x[i,])))
-1/2*log(det(2*pi*R))*n - 1/2*sum(diag((x)%*%solve(R)%*%t(x)))
-1/2*log(det(2*pi*R))*n - 1/2*sum(diag((n-1)*cov(x)%*%solve(R)))

sum(diag(1/2*(x)%*%solve(R)%*%t(x))) #trace
sum(diag(1/2*t(x)%*%x%*%solve(R))) #Tr(ABC) = Tr(CAB) = Tr(BCA)
sum(diag(1/2*(n-1)*cov(x)%*%solve(R))) #Tr(ABC) = Tr(CAB) = Tr(BCA)

# cov(x)
# t(x)%*%x/(n-1)

microbenchmark::microbenchmark(sum(mvtnorm::dmvnorm(x, sigma = R, log = T)),
-1/2*log(det(2*pi*R))*n - 1/2*sum(diag((n-1)*cov(x)%*%solve(R))))


# ldens + (ldens+1/2*log(det(2*pi*R2))*n2) * (i^2 - 1)
# 
# a <- ldens_low_scaled + 1/2*log(det(2*pi*R2))*n2
# ddf <- n2 / n1 #data duplication factor
# 
# 
# ddf * ldens_low + (1 - cov_scale(ddf,n1)^2) * a * ddf
# 
# sum(mvtnorm::dmvnorm(do.call(rbind, lapply(1:c,function(foo)x1)), sigma = R2, log = T))
#         
# 
# cov(do.call(rbind, lapply(1:ddf,function(foo)x1)))


##### see if we can find any telling pattern with the integral?
#set sample parameters
r = 0.4
p <- 3
n <- 10
R <- diag(p)*(1-r) + r

#new corrmat
R2 <- rethinking::rlkjcorr(1, p)

unnormalized_log_density <- function(corrmat, sample_corrmat, n) {
  -1/2*log(det(2*pi*corrmat))*n - 1/2*sum(diag((n-1)*sample_corrmat%*%solve(corrmat)))  
}

unnormalized_log_density(R2, R, n)

#use monte carlo sim to compute integral over all R2 

take_integral_and_determinant <- function(p, n, nR2 = 1E3){
  R <- rethinking::rlkjcorr(1, p)
  R2s <- rethinking::rlkjcorr(nR2, p)
  densities <- sapply(1:nR2, function(i) unnormalized_log_density(R2s[i,,], sample_corrmat = R, n))
  MCintegral_over_densities <- log(mean(exp(densities - max(densities) + 10))) - 10 + max(densities)
  ldetR <- determinant(R, logarithm= T)$modulus[1]
  ldetRtransf <- log(1 / sqrt(det(R * 2 * pi)))
  return(c(ldetR = ldetR, MCintegral_over_densities = MCintegral_over_densities, ldetRtransf = ldetRtransf))
}

p <- 3
n <- 10
logdet_VS_mcIntegral <- t(replicate(50, take_integral_and_determinant(p = p, n = n, nR2 = 5E4)))
# logdet_VS_mcIntegral <- logdet_VS_mcIntegral[-which.min(logdet_VS_mcIntegral[,1]),]
plot(logdet_VS_mcIntegral, xlab = "log determinant of sample correlation matrix",
     ylab = "log monte carlo integral of unnormalized distribution")
lmfit <- MASS::rlm(logdet_VS_mcIntegral[,2] ~ logdet_VS_mcIntegral[,1])
abline(lmfit, col = 2, lwd = 2)

par(mfrow = c(1,1))
plot(logdet_VS_mcIntegral[,3], logdet_VS_mcIntegral[,2], xlab = "log determinant of sample correlation matrix",
     ylab = "log monte carlo integral of unnormalized distribution")
plot(logdet_VS_mcIntegral[,"ldetRtransf"], logdet_VS_mcIntegral[,"MCintegral_over_densities"],
     xlab = latex2exp::TeX("$log(|2\\pi R|^{-1/2})$"),
     ylab = "log(average density)", pch = 19, col = adjustcolor(1, 0.4), cex = 1.5)
lmfit <- MASS::rlm(logdet_VS_mcIntegral[,"MCintegral_over_densities"] ~ logdet_VS_mcIntegral[,"ldetRtransf"])
abline(lmfit, col = adjustcolor(2, 0.5), lwd = 2)
title(main = paste0("Dimension = ", p, "; Sample Size = ", n))
title(main = paste0("Intercept = ", round(coefficients(lmfit)["(Intercept)"], 3), 
                    "; Slope = ", round(coefficients(lmfit)[2], 3)), cex.main = 1, line = 0.5, col.main = "darkred")

#####

fit_model_to_integral_and_determinant <- function(n, p, nrep = 50, nR2 = 1E3){
  logdet_VS_mcIntegral <- t(replicate(nrep, take_integral_and_determinant(n = n, p = p, nR2 = nR2)))
  # lmfit <- MASS::rlm(logdet_VS_mcIntegral[,2] ~ logdet_VS_mcIntegral[,1])
  lmfit <- MASS::rlm(logdet_VS_mcIntegral[,2] ~ logdet_VS_mcIntegral[,3])
  return(lmfit$coefficients)
}

p <- 4
sample_sizes <- 1:100
coeff_ests <- do.call(rbind, mclapply(sample_sizes, function(n){
  system(sprintf('echo "%s"', paste0(n, collapse="")));
  fit_model_to_integral_and_determinant(n = n, p = p, nrep = 100, nR2 = 5E3)}, mc.cores = 14))

par(mfrow = c(1,2))
plot(sample_sizes, coeff_ests[,"(Intercept)"], pch = 19, col = adjustcolor(1, 0.4), cex = 1,
     xlab = "Sample Size", ylab = "Estimated Intercept")
lmfit1 <- MASS::rlm(coeff_ests[,"(Intercept)"] ~ sample_sizes)
abline(lmfit1, col = adjustcolor(2, 0.5), lwd = 2)
title(main = paste0("Intercept = ", round(coefficients(lmfit1)["(Intercept)"], 3), 
                    "; Slope = ", round(coefficients(lmfit1)["sample_sizes"], 3)), 
      cex.main = 1, line = 0.5, col.main = "darkred")
text(par("usr")[2] + diff(par("usr")[c(1,2)]) * 0.2, 
     par("usr")[4] + diff(par("usr")[c(3,4)]) * 0.125, 
     labels = paste0("Dimension = ", p), pos = 3, cex = 2, xpd = NA)

plot(sample_sizes, coeff_ests[,"logdet_VS_mcIntegral[, 3]"], pch = 19, col = adjustcolor(1, 0.4), cex = 1,
     xlab = "Sample Size", ylab = "Estimated Slope")
lmfit2 <- MASS::rlm(coeff_ests[,"logdet_VS_mcIntegral[, 3]"] ~ sample_sizes)
abline(lmfit2, col = adjustcolor(2, 0.5), lwd = 2)
title(main = paste0("Intercept = ", round(coefficients(lmfit2)["(Intercept)"], 3), 
                    "; Slope = ", round(coefficients(lmfit2)["sample_sizes"], 3)), 
      cex.main = 1, line = 0.5, col.main = "darkred")

#####
# R <- rethinking::rlkjcorr(1, p)
# n_R2s <- 1E4
# R2s <- rethinking::rlkjcorr(n_R2s, p, eta = 1)
# densities <- sapply(1:n_R2s, function(i) rethinking::dlkjcorr(R2s[i,,], eta = 1.1, log = T))
# MCintegral_over_densities <- log(sum(exp(densities) / n_R2s))
# MCintegral_over_densities

# log(sum(exp(densities) / n_R2s))
# 
# log(sum(exp(densities - log(1/n_R2s))))
# 
# 


#####
nrunif <- 1E4
bound <- 0.34
sum(dnorm(runif(nrunif, min = -bound, max = bound)) / nrunif) * 2 * bound
mean(dnorm(runif(nrunif, min = -bound, max = bound))) * 2 * bound
pnorm(bound) - pnorm(-bound)

nrunif <- 5E5
bound <- 10
x <- matrix(runif(nrunif*p, min = -bound, max = bound), nrunif, p)
mean(mvtnorm::dmvnorm(x, sigma = R, mean = rep(0, 3))) * (2 * bound)^p


##### now try to implement this in Stan

#check #1 -- confirm that the average density is the same for all distributions, since the volume is the same
#check #2 -- implement in Stan and sample from the prior, confirm that we get the prior on 'n' back

## mvgamma functions ##
mvgamma <- function(a,p) prod(sapply(1:p, function(j) gamma(a + (1-j)/2))) * pi^(p*(p-1)/4)
lmvgamma <- function(a,p) sum(sapply(1:p, function(j) lgamma(a + (1-j)/2))) + (p*(p-1)/4)*log(pi)
#check
a = 7.1
p = 5
mvgamma(a, p)
exp(lmvgamma(a, p))

#density functions
log_unnormalized_density <- function(R, C, N) {
  -1/2*N*log(det(2*pi*R)) - 1/2*sum(diag((N-1)*C%*%solve(R)))  
}

normalizing_constant <- function(C, N){
  p <- dim(C)[1]
  2^(p*(1-(N-p-1)/2)) * det((N-1)*C)^((N-p-1)/2) / mvgamma(a = (N-p-1)/2, p = p)
}

log_normalizing_constant <- function(C, N){
  p <- dim(C)[1]
  p*(1-(N-p-1)/2)*log(2) + 
    (N-p-1)/2*determinant((N-1)*C, logarithm= T)$modulus[1] - 
    lmvgamma(a = (N-p-1)/2, p = p)
}

log_normalized_density <- function(R, C, N){
  log_unnormalized_density(R, C, N) + log_normalizing_constant(C, N)
}

#check
r = 0.4
p <- 3
R <- diag(p)*(1-r) + r
N <- 15
r_c <- 0.6
C <- diag(p)*(1-r_c) + r_c
normalizing_constant(R, N)
exp(log_normalizing_constant(R, N))

## confirm equality of mean densities ##
rlkj <- function (d, eta = 1) {
  alpha <- eta + (d - 2)/2
  r12 <- 2 * rbeta(1, alpha, alpha) - 1
  R <- matrix(0, d, d)
  R[1, 1] <- 1
  R[1, 2] <- r12
  R[2, 2] <- sqrt(1 - r12^2)
  if (d > 2) 
    for (m in 2:(d - 1)) {
      alpha <- alpha - 0.5
      y <- rbeta(1, m/2, alpha)
      z <- rnorm(m, 0, 1)
      z <- z/sqrt(crossprod(z)[1])
      R[1:m, m + 1] <- sqrt(y) * z
      R[m + 1, m + 1] <- sqrt(1 - y)
    }
  return(crossprod(R))
}

Ns <- c(30, 35)
nrunif <- 1E4
runif_corrmats <- replicate(nrunif, rlkj(p, eta = 1))
sapply(Ns, function(N) {
  log_densities <- sapply(1:nrunif, function(i) log_normalized_density(runif_corrmats[,,i], C, N))
  mean_density <- log(mean(exp(log_densities - max(log_densities) + 10))) - 10 + max(log_densities)
  mean_density
})

