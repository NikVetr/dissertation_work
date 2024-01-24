library(coda)
library(rethinking)
setwd("/Volumes/1TB/Bailey")
nreps <- 1:500

n_monomorphic <- nreps
n_ok_traits <- nreps
for(i in nreps){
  load(file = paste0("traits_indiv_discr/traits_indiv_discr_", i))
  all_traits <- do.call(rbind, traits_indiv_discr)
  n_monomorphic[i] <- sum(sapply(1:118, function(colj) length(table(all_traits[,colj])) == 1))
  ldis <- which(sapply(1:118, function(colj) sum(table(all_traits[,colj]) > 10)) > 1)
  n_ok_traits[i] <- length(ldis)
  print(paste0(i, ": ", n_monomorphic[i], ", ", n_ok_traits[i]))
}
quantile(n_monomorphic, 0:20/20)
quantile(118 - n_ok_traits - n_monomorphic, 0:20/20)
for(i in nreps){
  means <- read.table(file = paste0("means/true_means_", i))
  load(file = paste0("data_sim/param_estimates_replicate_", i, "_optimrep_1"))
  print(cor(c(current_params$means), c(as.matrix(means))))
  plot(c(current_params$means), c(as.matrix(means)))
  abline(0, 1)
}

extraFilename <- "_noPanAsian"
corrMat <- as.matrix(read.table(paste0("output/nearPD_Corr", extraFilename, ".txt")))
n_traits <- dim(corrMat)[1]
weight <- 0.02
corrMat <- (corrMat + weight*diag(n_traits)) / (1 + weight)
for(i in nreps){
  load(file = paste0("traits_indiv_discr/traits_indiv_discr_", i))
  all_traits <- do.call(rbind, traits_indiv_discr)
  sapply(1:118, function(colj) sum(table(all_traits[,colj]) > 10)) > 3
  
  load(file = paste0("data_sim/param_estimates_replicate_", i, "_optimrep_1"))
  cor1 <- current_params$cor
  load(file = paste0("data_sim/param_estimates_replicate_", i, "_optimrep_2"))
  cor2 <- current_params$cor
  load(file = paste0("data_sim/param_estimates_replicate_", i, "_optimrep_3"))
  cor3 <- current_params$cor
  load(file = paste0("data_sim/param_estimates_replicate_", i, "_optimrep_4"))
  cor4 <- current_params$cor
  
  # print(cor(cor1[upper.tri(cor1)], cor2[upper.tri(cor2)]))
  
  cor_est <- as.matrix(Matrix::nearPD((cor1 + cor2 + cor3 + cor4)/4, corr = T)$mat)
  cor_est <- (cor_est + weight*diag(n_traits)) / (1 + weight)
  print(cor(cor_est[upper.tri(cor_est)], corrMat[upper.tri(corrMat)]))
  
  #lotsa data inds
  ldis <- which(sapply(1:118, function(colj) sum(table(all_traits[,colj]) > 5)) > 3)
  cor_est_ld <- cor_est[ldis, ldis]
  corrMat_ld <- corrMat[ldis, ldis]
  print(cor(cor_est_ld[upper.tri(cor_est_ld)], corrMat_ld[upper.tri(corrMat_ld)]))
  plot(cor_est_ld[upper.tri(cor_est_ld)], corrMat_ld[upper.tri(corrMat_ld)])
  
  plot(cor_est[upper.tri(cor_est)], corrMat[upper.tri(corrMat)])
  abline(0, 1)
}


thresholds <- as.matrix(read.table(paste0("output/mean_of_thresholds", extraFilename, ".txt")))
for(i in nreps){
  load(file = paste0("data_sim/param_estimates_replicate_", i, "_optimrep_1"))
  
  print(cor(current_params$threshold_mat[,-1][current_params$threshold_mat[,-1] != Inf], thresholds[,-1][thresholds[,-1] != Inf]))
  plot(current_params$threshold_mat[,-1][current_params$threshold_mat[,-1] != Inf], thresholds[,-1][thresholds[,-1] != Inf])
  abline(0, 1)
}

setwd("/Volumes/1TB/Bailey_LOTSADATA/")

for(i in 501:600){
  means <- read.table(file = paste0("means/true_means_", i))
  load(file = paste0("data_sim/param_estimates_replicate_", i, "_optimrep_2"))
  print(cor(c(current_params$means), c(as.matrix(means))))
  plot(c(current_params$means), c(as.matrix(means)))
  abline(0, 1)
}
#357, 417 are bad

extraFilename <- "_noPanAsian"
corrMat <- as.matrix(read.table(paste0("output/nearPD_Corr", extraFilename, ".txt")))
n_traits <- dim(corrMat)[1]
weight <- 0.02
corrMat <- (corrMat + weight*diag(n_traits)) / (1 + weight)
for(i in 501:600){
  load(file = paste0("traits_indiv_discr/traits_indiv_discr_", i))
  all_traits <- do.call(rbind, traits_indiv_discr)
  sapply(1:118, function(colj) sum(table(all_traits[,colj]) > 10)) > 3
  
  load(file = paste0("data_sim/param_estimates_replicate_", i, "_optimrep_1"))
  cor1 <- current_params$cor
  load(file = paste0("data_sim/param_estimates_replicate_", i, "_optimrep_2"))
  cor2 <- current_params$cor
  load(file = paste0("data_sim/param_estimates_replicate_", i, "_optimrep_3"))
  cor3 <- current_params$cor
  load(file = paste0("data_sim/param_estimates_replicate_", i, "_optimrep_4"))
  cor4 <- current_params$cor
  
  # print(cor(cor1[upper.tri(cor1)], cor2[upper.tri(cor2)]))
  
  cor_est <- as.matrix(Matrix::nearPD((cor1 + cor2 + cor3 + cor4)/4, corr = T)$mat)
  cor_est <- as.matrix(Matrix::nearPD((cor1 + cor2)/2, corr = T)$mat)
  cor_est <- (cor_est + weight*diag(n_traits)) / (1 + weight)
  print(cor(cor_est[upper.tri(cor_est)], corrMat[upper.tri(corrMat)]))
  
  #lotsa data inds
  ldis <- which(sapply(1:118, function(colj) sum(table(all_traits[,colj]) > 5)) > 3)
  cor_est_ld <- cor_est[ldis, ldis]
  corrMat_ld <- corrMat[ldis, ldis]
  print(cor(cor_est_ld[upper.tri(cor_est_ld)], corrMat_ld[upper.tri(corrMat_ld)]))
  plot(cor_est_ld[upper.tri(cor_est_ld)], corrMat_ld[upper.tri(corrMat_ld)])
  
  plot(cor_est[upper.tri(cor_est)], corrMat[upper.tri(corrMat)])
  abline(0, 1)
}


thresholds <- as.matrix(read.table(paste0("output/mean_of_thresholds", extraFilename, ".txt")))
for(i in 1:100){
  load(file = paste0("data_sim/param_estimates_replicate_", i, "_optimrep_1"))
  
  print(cor(current_params$threshold_mat[,-1][current_params$threshold_mat[,-1] != Inf], thresholds[,-1][thresholds[,-1] != Inf]))
  plot(current_params$threshold_mat[,-1][current_params$threshold_mat[,-1] != Inf], thresholds[,-1][thresholds[,-1] != Inf])
  abline(0, 1)
}

setwd("/Volumes/1TB/Bailey/")
library(phangorn)
for(i in 501:600){
  true_tree <- read.tree(file = paste0("trees/tree_", i))
  t1 <- read.tree(file = paste0("output_sim/fixCorrs_infRates_allTraits_trees_regRates_replicate", i, "_c1.trees"))
  t2 <- read.tree(file = paste0("output_sim/fixCorrs_infRates_allTraits_trees_regRates_replicate", i, "_c2.trees"))
  t3 <- read.tree(file = paste0("output_sim/fixCorrs_infRates_allTraits_trees_regRates_replicate", i, "_c3.trees"))
  t4 <- read.tree(file = paste0("output_sim/fixCorrs_infRates_allTraits_trees_regRates_replicate", i, "_c4.trees"))
  est_tree <- maxCladeCred(c(t1,t2,t3,t4))
  est_tree <- maxCladeCred(t4)
  print(phangorn::RF.dist(est_tree, unroot(true_tree), normalize = T))
  # RF.dist(maxCladeCred(t1), maxCladeCred(t3), normalize = T)
  
  # true_rates <- read.table(file = paste0("rates/rates_", i))
  # r1 <- read.table(file = paste0("output_sim/fixCorrs_infRates_allTraits_rates_regRates_replicate", i, "_c1.txt"))
  # plot((as.numeric(true_rates)), as.vector(apply(r1, 2, mean)))
  
}

for(i in 1:100){
  true_tree <- read.tree(file = paste0("trees/tree_", i))
  t1 <- read.tree(file = paste0("output_sim/fixCorrs_infRates_allTraits_trees_regRates_replicate_TRUEVALS_", i, "_c1.trees"))
  t2 <- read.tree(file = paste0("output_sim/fixCorrs_infRates_allTraits_trees_regRates_replicate_TRUEVALS_", i, "_c2.trees"))
  t3 <- read.tree(file = paste0("output_sim/fixCorrs_infRates_allTraits_trees_regRates_replicate_TRUEVALS_", i, "_c3.trees"))
  t4 <- read.tree(file = paste0("output_sim/fixCorrs_infRates_allTraits_trees_regRates_replicate_TRUEVALS_", i, "_c4.trees"))
  est_tree <- maxCladeCred(c(t1,t2,t3,t4))
  est_tree <- maxCladeCred(t4)
  print(phangorn::RF.dist(est_tree, unroot(true_tree), normalize = T))
  # RF.dist(maxCladeCred(t1), maxCladeCred(t3), normalize = T)
  
  # true_rates <- read.table(file = paste0("rates/rates_", i))
  # r1 <- read.table(file = paste0("output_sim/fixCorrs_infRates_allTraits_rates_regRates_replicate", i, "_c1.txt"))
  # plot((as.numeric(true_rates)), as.vector(apply(r1, 2, mean)))
  
}

for(i in 1:100){
  t1 <- read.tree(file = paste0("output_sim/fixCorrs_infRates_allTraits_trees_regRates_replicate", sample(1:500, 1), "_c1.trees"))
  t2 <- read.tree(file = paste0("output_sim/fixCorrs_infRates_allTraits_trees_regRates_replicate", sample(1:500, 1), "_c2.trees"))
  print(RF.dist(maxCladeCred(t2), maxCladeCred(t1), normalize = T))
}

########################
#### making figures ####
########################

#### figure 4.6: comparing true vs inferred values
extraFilename <- "_noPanAsian"
corrMat <- as.matrix(read.table(paste0("output/nearPD_Corr", extraFilename, ".txt")))
n_traits <- dim(corrMat)[1]
nTips <- 8
weight <- 0.02
corrMat <- (corrMat + weight*diag(n_traits)) / (1 + weight)

thresholds <- as.matrix(read.table(paste0("output/mean_of_thresholds", extraFilename, ".txt")))

reps <- c(1:1000)
cor_cors <- rep(0, length(reps))
cor_means <- rep(0, length(reps))
cor_threshes <- rep(0, length(reps))
all_cors <- lapply(1:length(reps), function(x) NA)
all_means <- lapply(1:length(reps), function(x) NA)
all_threshes <- lapply(1:length(reps), function(x) NA)

for(i_ind in 1:length(reps)){
  cat(paste0(" ", i_ind))
  i <- reps[i_ind]
  means <- as.matrix(read.table(file = paste0("means/true_means_", i)))
  if(file.exists(paste0("data_sim/param_estimates_replicate_", i, "_optimrep_1"))){
    load(file = paste0("data_sim/param_estimates_replicate_", i, "_optimrep_1"))
    cor1 <- current_params$cor
    m1 <- current_params$means
    t1 <- current_params$threshold_mat
    h1 <- T
  } else {
    cor1 <- 0
    m1 <- 0
    t1 <- 0
    h1 <- F
  }
  if(file.exists(paste0("data_sim/param_estimates_replicate_", i, "_optimrep_2"))){
    load(file = paste0("data_sim/param_estimates_replicate_", i, "_optimrep_2"))
    cor2 <- current_params$cor
    m2 <- current_params$means
    t2 <- current_params$threshold_mat
    h2 <- T
  } else {
    cor2 <- 0
    m1 <- 0
    t1 <- 0
    h2 <- F
  }
  if(file.exists(paste0("data_sim/param_estimates_replicate_", i, "_optimrep_3"))){
    load(file = paste0("data_sim/param_estimates_replicate_", i, "_optimrep_3"))
    cor3 <- current_params$cor
    m3 <- current_params$means
    t3 <- current_params$threshold_mat
    h3 <- T
  } else {
    cor3 <- 0
    m1 <- 0
    t1 <- 0
    h3 <- F
  }
  if(file.exists(paste0("data_sim/param_estimates_replicate_", i, "_optimrep_4"))){
    load(file = paste0("data_sim/param_estimates_replicate_", i, "_optimrep_4"))
    cor4 <- current_params$cor
    m4 <- current_params$means
    t4 <- current_params$threshold_mat
    h4 <- T
  } else {
    cor4 <- 0
    m1 <- 0
    t1 <- 0
    h4 <- F
  }

  cor_est <- as.matrix(Matrix::nearPD((cor1 + cor2 + cor3 + cor4) / sum(h1, h2, h3, h4), corr = T)$mat)
  cor_est <- (cor_est + weight*diag(n_traits)) / (1 + weight)
  cor_cors[i_ind] <- cor(cor_est[upper.tri(cor_est)], corrMat[upper.tri(corrMat)])
  all_cors[[i_ind]] <- cbind(cor_est[upper.tri(cor_est)], corrMat[upper.tri(corrMat)])
  
  cat("c")
  
  means_est <- (m1 + m2 + m3 + m4) / sum(h1, h2, h3, h4)
  cor_means[i_ind] <- cor(c(means_est), c(means))
  all_means[[i_ind]] <- cbind(c(means_est), c(means))
  
  cat("m")
  
  thresh_est <- (t1 + t2 + t3 + t4) / sum(h1, h2, h3, h4)
  cor_threshes[i_ind] <- cor(thresh_est[,-1][thresh_est[,-1] != Inf], thresholds[,-1][thresholds[,-1] != Inf])
  all_threshes[[i_ind]] <- cbind(thresh_est[,-1][thresh_est[,-1] != Inf], thresholds[,-1][thresholds[,-1] != Inf])
  
  cat("t")
  
}

fins <- rep(0, length(reps))
for(i in reps){
  fins[i] <- (sum(file.exists(paste0("data_sim/param_estimates_replicate_", i, "_optimrep_1")),
      file.exists(paste0("data_sim/param_estimates_replicate_", i, "_optimrep_2")),
      file.exists(paste0("data_sim/param_estimates_replicate_", i, "_optimrep_3")),
      file.exists(paste0("data_sim/param_estimates_replicate_", i, "_optimrep_4"))))
}
table(fins)
fails <- which(fins != 4)
plot(cor_threshes[-which(fins != 4)])
plot(cor_means[-which(fins != 4)])
plot(cor_cors[-which(fins != 4)])

all_threshes_1 <- all_threshes[1:500]
all_cors_1 <- all_cors[1:500]
all_means_1 <- all_means[1:500]
cor_means_1 <- cor_means[1:500]
cor_cors_1 <- cor_cors[1:500]
cor_threshes_1 <- cor_threshes[1:500]

if(any(fails < 500)){
  all_threshes_1 <- all_threshes_1[-fails[fails < 500]]
  all_cors_1 <- all_cors_1[-fails[fails < 500]]
  all_means_1 <- all_means_1[-fails[fails < 500]]
  cor_means_1 <- cor_means_1[-fails[fails < 500]]
  cor_cors_1 <- cor_cors_1[-fails[fails < 500]]
  cor_threshes_1 <- cor_threshes_1[-fails[fails < 500]]
}

all_threshes_1 <- do.call(all_threshes_1, what = rbind)
all_cors_1 <- do.call(all_cors_1, what = rbind)
all_means_1 <- do.call(all_means_1, what = rbind)

#

all_threshes_2 <- all_threshes[501:1000]
all_cors_2 <- all_cors[501:1000]
all_means_2 <- all_means[501:1000]
cor_means_2 <- cor_means[501:1000]
cor_cors_2 <- cor_cors[501:1000]
cor_threshes_2 <- cor_threshes[501:1000]


if(any(fails > 500)){
  all_threshes_2 <- all_threshes_2[-(fails[fails > 500]-500)]
  all_cors_2 <- all_cors_2[-(fails[fails > 500]-500)]
  all_means_2 <- all_means_2[-(fails[fails > 500]-500)]
  cor_means_2 <- cor_means_2[-(fails[fails > 500]-500)]
  cor_cors_2 <- cor_cors_2[-(fails[fails > 500]-500)]
  cor_threshes_2 <- cor_threshes_2[-(fails[fails > 500]-500)]
}

all_threshes_2 <- do.call(all_threshes_2, what = rbind)
all_cors_2 <- do.call(all_cors_2, what = rbind)
all_means_2 <- do.call(all_means_2, what = rbind)

all_data_pts <- list(all_means_1, all_cors_1, all_threshes_1, all_means_2, all_cors_2, all_threshes_2)
all_corcofs <- list(cor_means_1, cor_means_2, cor_cors_1, cor_cors_2, cor_threshes_1, cor_threshes_2)
xaxis_labels <- list(-1:2*3, -1:1/4, 0:5, -1:2*3, -1:1/2, 0:5)
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
dens2 <- function (x, adj = 0.5, norm.comp = FALSE, main = "", show.HPDI = FALSE, 
          show.zero = FALSE, rm.na = TRUE, add = FALSE, shade.col = "black", ...) {
  if (inherits(x, "data.frame")) {
    n <- ncol(x)
    cnames <- colnames(x)
    set_nice_margins()
    par(mfrow = make.grid(n))
    for (i in 1:n) {
      dens(x[, i], adj = adj, norm.comp = norm.comp, show.HPDI = show.HPDI, 
           show.zero = TRUE, xlab = cnames[i], ...)
    }
  }
  else {
    if (rm.na == TRUE) 
      x <- x[!is.na(x)]
    thed <- density(x, adjust = adj)
    if (add == FALSE) {
      set_nice_margins()
      plot(thed, main = main, ...)
    }
    else lines(thed$x, thed$y, ...)
    if (show.HPDI != FALSE) {
      hpd <- HPDI(x, prob = show.HPDI)
      shade(thed, hpd, col = shade.col)
    }
    if (norm.comp == TRUE) {
      mu <- mean(x)
      sigma <- sd(x)
      curve(dnorm(x, mu, sigma), col = "white", lwd = 2, 
            add = TRUE)
      curve(dnorm(x, mu, sigma), add = TRUE)
    }
    if (show.zero == TRUE) {
      lines(c(0, 0), c(0, max(thed$y) * 2), lty = 2)
    }
  }
}





grDevices::cairo_pdf(filename = "/Volumes/macOS/Users/nikolai/dissertation/figures/chpt4_figure6.pdf", width = 2000 / 72, height = 2000 / 72)
par(mfrow = c(3,3))
for(i in 1:6){
  par(mar = c(ifelse(any(i == c(4,5,6)), 7, 4),
              ifelse(any(i == c(1,4)), 7, 6),
              ifelse(any(i == c(1,2,3)), 6, 2),
              1))
  par(xpd = F)
  subdata <- all_data_pts[[i]][sample(1:nrow(all_data_pts[[i]]), replace = F, size = 1E4),]
  plot(subdata, pch = 16, 
       cex = 5, xlab = "", ylab = "", cex.axis = 3, col = rgb(0,0,0,0.025), xaxt = "n", yaxt = "n", frame.plot = F)
  axis(side = 1, tick = -2, lwd = 3, lwd.ticks = 3, cex.axis = 3, labels = rep("", length(xaxis_labels[[i]])), at = xaxis_labels[[i]])
  mtext(text = xaxis_labels[[i]], side = 1, line = 2, at = xaxis_labels[[i]], cex = 2)
  axis(side = 2, tick = -2, lwd = 3, lwd.ticks = 3, cex.axis = 3)
  box("plot", lwd = 3, lty = 1)
  if(i < 4){
    title(main = c("Means", "Correlations", "Thresholds")[i], cex.main = 6)
  }
  if(i > 3){
    title(xlab = "Estimated Values", cex.lab = 4, line = 5)
  }
  if(any(i == c(1,4))){
    title(ylab = "True Values", cex.lab = 4, line = 4)
  }
  if(i == 1){
    legend(x = "topleft", col = 2, lwd = 4, lty = 2, legend = "1-to-1 line", box.lty = 2, box.lwd = 4, cex = 3)
  }
  fig_label(paste0(letters[i], ")"), shrinkX = ifelse(any(i == c(3,6)), 0.95, 0.98), shrinkY = ifelse(any(i == c(3,6)), 0.99, 0.98), cex = 4.25)
  fig_label(text = latex2exp::TeX(paste0("R$^2$ = ", round(cor(subdata)[1,2]^2, 3))), region = "plot", pos = "top",
            shrinkY = ifelse(any(i == c(3,6)), 0.99, 0.98), cex = 2.5)
  fig_label(text = ifelse(i < 4, "Empirical Parameterization", "High Sample-Size Parameterization"), region = "plot", pos = "bottomright", 
            cex = 2.75, shrinkX = 1, shrinkY = ifelse(any(i == c(3,6)), 0.625, 0.98))#, shrinkX = , )
  box(which = "figure", lty = 1, col = rgb(0,0,0,0.75), lwd = 2)
  abline(0, 1, col = 2, lwd = 4, lty = 2)
}
for(i in 1:6){
  if(i < 3){
    par(mar = c(7,7,2,1))
  } else {
    par(mar = c(7,6,2,1))
  }
  densobj <- density(x = all_corcofs[[i]]^2, adjust = 2)
  if(i %% 2 != 0){
    plot(densobj, xlim = c(0,1), ylim = c(0,c(16,13,34)[(i+1)/2]), lwd = 3, main = "", xaxt = "n", yaxt = "n", frame.plot = F, xlab = "", ylab = "")
    axis(side = 1, tick = -2, lwd = 3, lwd.ticks = 3, cex.axis = 3, labels = rep("", length(xaxis_labels[[i]])), at = xaxis_labels[[i]])
    mtext(text = 0:5*2/10, side = 1, line = 2, at = 0:5*2/10, cex = 2)
    axis(side = 2, tick = -2, lwd = 3, lwd.ticks = 3, cex.axis = 3)
    box("plot", lwd = 3, lty = 1)
    title(xlab = latex2exp::TeX(paste0("R$^2$ Across Replicates")), cex.lab = 4, line = 6)
    shade(densobj, col = rgb(0,0,0,0.4), lim = c(0, 2))
    fig_label(paste0(c("g", "h", "i")[(i+1)/2], ")"), shrinkX = 0.9, shrinkY = 0.9875, cex = 4.25)
  } else{
    lines(densobj, lwd = 3)
    shade(densobj, col = rgb(0,0,0,0.8), lim = c(0, 2))
    # segments(x0 = -0.02, x1= 1.02, y0 = 0, y1 = 0, col = "white", lwd = 4)
  }
  if(i == 1){
    title(ylab = "Density", cex.lab = 4, line = 4)
    legend(x = "topleft", fill = c(rgb(0,0,0,0.4), rgb(0,0,0,0.8)), 
           legend = c("Empirical Parameterization", "High Sample-Size Parameterization"), box.lty = 2, box.lwd = 4, cex = 3)    
  }
  box(which = "figure", lty = 1, col = rgb(0,0,0,0.75), lwd = 3)
}
dev.off()


## now on to generate figure 7!

missing_probs <- read.table(file = "/Volumes/1TB/Bailey/estimated_missing_probs/estimated_missing_probs_replicate_1_optimrep_1")
reps <- 1:500
missing_probs <- matrix(0, nrow = length(reps), ncol = nrow(missing_probs))

for(i_ind in 1:length(reps)){
  cat(paste0(" ", i_ind))
  i <- reps[i_ind]
  means <- as.matrix(read.table(file = paste0("means/true_means_", i)))
  if(file.exists(paste0("estimated_missing_probs/estimated_missing_probs_replicate_", i, "_optimrep_1"))){
    mp1 <- read.table(file = paste0("estimated_missing_probs/estimated_missing_probs_replicate_", i, "_optimrep_1"))
    h1 <- T
  } else {
    mp1 <- 0
    h1 <- F
  }
  if(file.exists(paste0("estimated_missing_probs/estimated_missing_probs_replicate_", i, "_optimrep_2"))){
    mp2 <- read.table(file = paste0("estimated_missing_probs/estimated_missing_probs_replicate_", i, "_optimrep_2"))
    h2 <- T
  } else {
    mp2 <- 0
    h2 <- F
  }
  if(file.exists(paste0("estimated_missing_probs/estimated_missing_probs_replicate_", i, "_optimrep_3"))){
    mp3 <- read.table(file = paste0("estimated_missing_probs/estimated_missing_probs_replicate_", i, "_optimrep_3"))
    h3 <- T
  } else {
    mp3 <- 0
    h3 <- F
  }
  if(file.exists(paste0("estimated_missing_probs/estimated_missing_probs_replicate_", i, "_optimrep_4"))){
    mp4 <- read.table(file = paste0("estimated_missing_probs/estimated_missing_probs_replicate_", i, "_optimrep_4"))
    h4 <- T
  } else {
    mp4 <- 0
    h4 <- F
  }
  mp <- (mp1 + mp2 + mp3 + mp4) / sum(c(h1, h2, h3, h4))
  missing_probs[i_ind,] <- rev(c(as.vector(mp))$x)
  
}

true_missing_probs <- c(0.7858551, 0.6900000, 0.5744717, 0.4501955, 0.3318387, 0.2314964, 0.1544808)[1:7]
1:7 /
library(vioplot)

grDevices::cairo_pdf(filename = "/Volumes/macOS/Users/nikolai/dissertation/figures/chpt4_figure7.pdf", width = 750 / 72, height = 750 / 72)
par(mar = c(6,6,6,4))
vioplot(missing_probs, at = rev(true_missing_probs) * 10, xaxt = "n", yaxt = "n", frame.plot = F, rectCol = NA, 
        wex = 1.25, cex  = 0.00001, cex.lab = "n", lwd = 3, col = rgb(0,0,0,0.35))
points(rev(true_missing_probs) * 10, apply(missing_probs, 2, median), pch = 16, col = "black", cex = 0.5)
box(which = "plot", lwd = 3)
axis(side = 1, tick = -2, lwd = 3, lwd.ticks = 3, cex.axis = 3, labels = rep("", 11), at = 0:10)
mtext(text = 2:8/10, side = 1, line = 2, at = 2:8, cex = 2.5)
axis(side = 2, tick = -2, lwd = 3, lwd.ticks = 3, cex.axis = 2.5)
title(main = "Pooled State-Dependent Missingness", cex.main = 2.75)
title(xlab = "True Pr(Missing | State)", cex.lab = 3, line = 4.5)
title(ylab = "Estimated Pr(Missing | State)", cex.lab = 3, line = 3.75)
abline(a = 0, b = 1/10, lwd = 3, col = 2, lty = 2)
legend(x = "topleft", col = 2, lwd = 3, lty = 2, legend = "1-to-1 line", box.lty = 2, box.lwd = 4, cex = 2)
dev.off()


###################################
## and now on to generate figure 8!
###################################
shorter <- function(x1, x2){
  if(length(x1) > length(x2)){
    return(x2)
  } else{
    return(x1)
  }
}
startsWith2 <- function(x, pre){apply(sapply(1:length(pre), function(n) startsWith(x, pre[n])), 1, any)}
greedyCT <- function(trees){
  tipNames <- trees[[1]]$tip.label
  tipNums <- paste0("t", 1:length(tipNames))
  for(i in 1:length(trees)){
    trees[[i]]$tip.label <- tipNums[match(trees[[i]]$tip.label, tipNames)]
  }
  ape::write.tree(trees, file = "gct_trees.txt")
  writeLines("gct_trees.txt\nY", con = "consense_script.txt")
  if(file.exists("outfile")){file.remove("outfile")}
  if(file.exists("outtree")){file.remove("outtree")}
  system(paste0("/Applications/phylip-3.695/exe/consense < consense_script.txt"), intern = T)
  gct <- ape::read.tree("outtree")
  if(file.exists("outfile")){file.remove("outfile")}
  if(file.exists("outtree")){file.remove("outtree")}
  if(file.exists("gct_trees.txt")){file.remove("gct_trees.txt")}
  if(file.exists("consense_script.txt")){file.remove("consense_script.txt")}
  gct$tip.label <- tipNames[match(gct$tip.label, tipNums)]
  return(gct)
}
convNum2Str <- function(nums, key){
  sapply(1:length(nums), function(x) key[nums[x]])
}
prop.part.df <- function(trees, cutoff = 0.01, bs = T){
  if(class(trees) == "multiPhylo"){
    if(trees[[1]]$Nnode == (length(trees[[1]]$tip.label) - 1)){
      trees <- unroot(trees) #unroot rooted trees
    }
    tipLabs <- trees[[1]]$tip.label
    numTrees <- length(trees)
    if(bs) {
      out <- as.prop.part(bitsplits(trees))
    } else {
      out <- prop.part(trees)
    }
    outList <- as.list.data.frame(out)
    pps <- attributes(outList)$number/numTrees
    props <- data.frame(pps, as.matrix(outList)); colnames(props) <- c("postProbs", "clade")
    props <- props[order(-pps),]
    props <- props[props[,1] > cutoff,]
    rownames(props) <- 1:nrow(props)
    props$cladeNames <- sapply(1:length(props[,1]), function(x) sort(convNum2Str(props$clade[[x]], attributes(outList)$labels)))
    props <- props[,c(1,3)]
    if(!bs) {
      props <- props[-1,]
    }
    allClades <- c(props$cladeNames, lapply(1:length(props$postProbs), function(x) sort(setdiff(tipLabs, props$cladeNames[[x]]))))
    if(any(duplicated(allClades))){
      # print("duplicate clades found")
      dupClades <- allClades[duplicated(allClades)]
      dupCladesComp <- lapply(1:length(dupClades), function(x) sort(setdiff(tipLabs, dupClades[[x]])))
      matchedDupes <- cbind(1:length(dupClades), sapply(1:length(dupClades), function(x) which(sapply(1:length(dupCladesComp), function(y) setequal(dupClades[[x]], dupCladesComp[[y]])))))
      dupes <- matrix(data = 0, nrow = length(matchedDupes[,1])/2, ncol = 2)
      for(i in 1:length(matchedDupes[,1])){
        pair <- sort(matchedDupes[i,])
        if(!any(sapply(1:length(dupes[,1]), function(x) pair == dupes[x,]))){
          dupes[min(which(apply(dupes, 1, sum) == 0)),] <- pair
        }
      }
      for(i in 1:length(dupes[,1])){
        clade1 <- dupClades[dupes[i,1]][[1]]
        clade2 <- dupClades[dupes[i,2]][[1]]
        propsInd1 <- which(sapply(1:length(props[,1]), function(x) setequal(clade1, props$cladeNames[[x]])))
        propsInd2 <- which(sapply(1:length(props[,1]), function(x) setequal(clade2, props$cladeNames[[x]])))
        props[propsInd1,1] <- props[propsInd1,1] + props[propsInd2,1]
        props <- props[-propsInd2,]
      }
      props <- props[order(props[,1], decreasing = T),]
    }
    props
  } else if (class(trees) == "phylo"){
    if(trees$Nnode == (length(trees$tip.label) - 1)){
      trees <- unroot(trees) #unroot rooted trees
    }
    tipLabs <- trees$tip.label
    numTrees <- 1
    if(bs) {
      out <- as.prop.part(bitsplits(trees))
    } else {
      out <- prop.part(trees)
    }
    outList <- as.list.data.frame(out)
    pps <- attributes(outList)$number/numTrees
    props <- data.frame(pps, as.matrix(outList)); colnames(props) <- c("postProbs", "clade")
    props <- props[order(-pps),]
    props <- props[props[,1] > cutoff,]
    rownames(props) <- 1:nrow(props)
    props$cladeNames <- sapply(1:length(props[,1]), function(x) sort(convNum2Str(props$clade[[x]], attributes(outList)$labels)))
    props <- props[,c(1,3)]
    if(!bs) {
      props <- props[-1,]
    }
    props
  }
}
compareTrees <- function(biparts1, biparts2, tipLabs, returnCladeNames = F){
  matchProbs <- matrix(0, nrow = sum(length(biparts1[,1]), length(biparts2[,1])), ncol = ifelse(returnCladeNames, 3, 2))
  counter <- 1
  biparts2$notSeen <- 1
  for(clade in 1:length(biparts1[,2])){
    cladeName <- biparts1[,2][[clade]]
    isThere <- sapply(1:length(biparts2[,1]), function(x) identical(cladeName, biparts2[x,2][[1]]))
    altCladeName <- sort(setdiff(tipLabs, cladeName))
    shortCladeName <- paste0(sort(shorter(cladeName, altCladeName)), collapse = ", ")
    orIsThere <- sapply(1:length(biparts2[,1]), function(x) identical(altCladeName, biparts2[x,2][[1]]))
    if(any(isThere)){
      biparts2$notSeen[isThere] <- 0
      if(returnCladeNames){
        matchProbs[counter,] <- c(biparts1[,1][clade], biparts2[,1][isThere], shortCladeName)
      } else {
        matchProbs[counter,] <- c(biparts1[,1][clade], biparts2[,1][isThere])
      }
    } else if (any(orIsThere)) {
      biparts2$notSeen[orIsThere] <- 0
      if(returnCladeNames){
        matchProbs[counter,] <- c(biparts1[,1][clade], biparts2[,1][orIsThere], shortCladeName)
      } else {
        matchProbs[counter,] <- c(biparts1[,1][clade], biparts2[,1][orIsThere])
      }
    } else {
      if(returnCladeNames){
        matchProbs[counter,] <- c(biparts1[,1][[clade]], 0, shortCladeName)
      } else {
        matchProbs[counter,] <- c(biparts1[,1][[clade]], 0)
      }
    }
    counter <- counter + 1
  }
  if(sum(biparts2$notSeen) > 0){
    if(returnCladeNames){
      cladeNames <- sapply(1:length(biparts2[biparts2$notSeen == 1,2]), function(bip) 
        paste0(sort(shorter(biparts2[biparts2$notSeen == 1,2][[bip]], sort(setdiff(tipLabs, biparts2[biparts2$notSeen == 1,2][[bip]])))), collapse = ", "))
      matchProbs[counter:(counter+sum(biparts2$notSeen)-1),] <- cbind(0, cbind(biparts2[biparts2$notSeen == 1,1], cladeNames)) #needs to be an sapply
    } else {
      matchProbs[counter:(counter+sum(biparts2$notSeen)-1),] <- cbind(0, biparts2[biparts2$notSeen == 1,1])
    }
  }
  if(returnCladeNames){
    matchProbs <- matchProbs[rowSums(matchProbs == c(0,0)) != ifelse(returnCladeNames,3,2),]
  }
  matchProbs
}
clade_prob <- function(tipNames, trees, allTipNames = NA, partfreqs = NA){
  if(is.na(partfreqs)[1]){
    clade_monophyly <- prop.part.df(trees, cutoff = 0.00001)
    if(is.na(allTipNames)[1]){
      allTipNames <- trees[[1]]$tip.label
    }
  } else {
    clade_monophyly <- partfreqs
    if(is.na(allTipNames)[1]){
      allTipNames <- unique(unlist(bipart_probs$cladeNames))
    }
  }
  cladeName <- tipNames
  clade <- list(allTipNames[startsWith2(allTipNames, cladeName)], allTipNames[!startsWith2(allTipNames, cladeName)])
  clade_ind <- (sapply(1:length(clade_monophyly$cladeNames), function(x) setequal(clade_monophyly$cladeNames[[x]], clade[[1]]) |
                         setequal(clade_monophyly$cladeNames[[x]], clade[[2]])))
  if(all(!clade_ind)){return(0)}
  # clade_monophyly$cladeNames[clade_ind]
  return(clade_monophyly$postProbs[clade_ind])
}
internal_nodes_probs <- function(tree, trees_to_use){
  bipart_probs <- prop.part.df(trees_to_use, cutoff = 0.00001)
  ti <- 1:length(tree$tip.label)
  tips <- sapply(1:Nnode(tree), function(node) as.numeric(na.omit(ti[getDescendants(tree, length(tree$tip.label) + node)])))
  tips <- sapply(1:length(tips), function(node) tree$tip.label[tips[[node]]])
  node_pps <- sapply(1:length(tips), function(node) clade_prob(tipNames = tips[[node]], allTipNames = tree$tip.label, partfreqs = bipart_probs))
  return(node_pps)
}
combine <- function(t1, t2, t3, t4){
  t1h <- all(!is.na(t1))
  t2h <- all(!is.na(t2))
  t3h <- all(!is.na(t3))
  t4h <- all(!is.na(t4))
  altr <- list(t1, t2, t3, t4)
  altr <- altr[which(c(t1h, t2h, t3h, t4h))]
  return(do.call(c, altr))
}
cpts <- function(t1, t2, t3, t4, tipLabs){
  t1h <- all(!is.na(t1))
  t2h <- all(!is.na(t2))
  t3h <- all(!is.na(t3))
  t4h <- all(!is.na(t4))
  altr <- list(t1, t2, t3, t4)
  altr <- altr[which(c(t1h, t2h, t3h, t4h))]
  pps_altr <- lapply(altr, function(trs) prop.part.df(trs))
  for(i in 1:length(pps_altr)){
    pps_altr[[i]]$cladeNames <- sapply(1:length(pps_altr[[i]]$cladeNames), function(clade)
      sapply(1:length(strsplit(pps_altr[[i]]$cladeNames[[clade]], split = "-")), function(name)
        paste0(strsplit(pps_altr[[i]]$cladeNames[[clade]], split = "-")[[name]], collapse = "")))
  }
  return(unlist(sapply(1:(length(altr)-1), function(rowi) sapply((rowi+1):length(altr), function(colj) cor(compareTrees(pps_altr[[rowi]], pps_altr[[colj]], tipLabs)[,1:2])[1,2]))))
}

reps <- 1:1500
post_prob_comp <- lapply(1:length(reps), function(X) NA)
compTreesChains <- lapply(1:length(reps), function(X) NA)
for(i_ind in 1:length(reps)){
  cat(paste0(" ", i_ind))
  i <- reps[i_ind]
  true_tree <- read.tree(file = paste0("trees/tree_", ifelse(i < 1001, i, i-1000)))
  if(file.exists(paste0("output_sim/fixCorrs_infRates_allTraits_trees_regRates_replicate", ifelse(i < 1001, i, paste0("_TRUEVALS_", i-1000)), "_c1.trees"))){
    t1 <- read.tree(file = paste0("output_sim/fixCorrs_infRates_allTraits_trees_regRates_replicate", ifelse(i < 1001, i, paste0("_TRUEVALS_", i-1000)), "_c1.trees"))
  } else {
    t1 <- NA
  }
  if(file.exists(paste0("output_sim/fixCorrs_infRates_allTraits_trees_regRates_replicate", ifelse(i < 1001, i, paste0("_TRUEVALS_", i-1000)), "_c2.trees"))){
    t2 <- read.tree(file = paste0("output_sim/fixCorrs_infRates_allTraits_trees_regRates_replicate", ifelse(i < 1001, i, paste0("_TRUEVALS_", i-1000)), "_c2.trees"))
  } else {
    t2 <- NA
  }
  if(file.exists(paste0("output_sim/fixCorrs_infRates_allTraits_trees_regRates_replicate", ifelse(i < 1001, i, paste0("_TRUEVALS_", i-1000)), "_c3.trees"))){
    t3 <- read.tree(file = paste0("output_sim/fixCorrs_infRates_allTraits_trees_regRates_replicate", ifelse(i < 1001, i, paste0("_TRUEVALS_", i-1000)), "_c3.trees"))
  } else {
    t3 <- NA
  }
  if(file.exists(paste0("output_sim/fixCorrs_infRates_allTraits_trees_regRates_replicate", ifelse(i < 1001, i, paste0("_TRUEVALS_", i-1000)), "_c4.trees"))){
    t4 <- read.tree(file = paste0("output_sim/fixCorrs_infRates_allTraits_trees_regRates_replicate", ifelse(i < 1001, i, paste0("_TRUEVALS_", i-1000)), "_c4.trees"))
  } else {
    t4 <- NA
  }
  if(any(!is.na(c(t1,t2,t3,t4)))){
    tipLabs <- c("NEAND",  "AUSNG", "WAS", "NEAS",  "AMER",   "SSAF",   "EUR", "SAS")
    trees <- combine(t1, t2, t3, t4)
    if(sum(is.na(c(t1,t2,t3,t4))) < 3){
      compTreesChains[[i_ind]] <- cpts(t1, t2, t3, t4, tipLabs = tipLabs)
    } else {
      compTreesChains[[i_ind]] <- NA
    }
    pp1 <- prop.part.df(trees, cutoff = 0.001)
    pp1$cladeNames <- sapply(1:length(pp1$cladeNames), function(clade)
      sapply(1:length(strsplit(pp1$cladeNames[[clade]], split = "-")), function(name)
        paste0(strsplit(pp1$cladeNames[[clade]], split = "-")[[name]], collapse = "")))
    pp2 <- prop.part.df(true_tree)
    pp2$cladeNames <- sapply(1:length(pp2$cladeNames), function(clade)
      sapply(1:length(strsplit(pp2$cladeNames[[clade]], split = "-")), function(name)
        paste0(strsplit(pp2$cladeNames[[clade]], split = "-")[[name]], collapse = "")))
    post_prob_comp[[i_ind]] <- compareTrees(pp1, pp2, tipLabs = tipLabs)
  } else {
    post_prob_comp[[i_ind]] <- NA
    compTreesChains[[i_ind]] <- NA
  }
}

reps <- 1:1500
ratePercentile <- lapply(1:length(reps), function(X) NA)
compRatesChains <- lapply(1:length(reps), function(X) NA)
rateCorrs <- rep(NA, length(reps))
for(i_ind in 1:length(reps)){
  cat(paste0(" ", i_ind))
  i <- reps[i_ind]
  true_rates <- as.numeric(read.table(file = paste0("rates/rates_", ifelse(i < 1001, i, i-1000))))
  if(file.exists(paste0("output_sim/fixCorrs_infRates_allTraits_rates_regRates_replicate", ifelse(i < 1001, i, paste0("_TRUEVALS_", i-1000)), "_c1.txt"))){
    r1 <- read.table(file = paste0("output_sim/fixCorrs_infRates_allTraits_rates_regRates_replicate", ifelse(i < 1001, i, paste0("_TRUEVALS_", i-1000)), "_c1.txt"))
  } else {
    r1 <- NA
  }
  if(file.exists(paste0("output_sim/fixCorrs_infRates_allTraits_rates_regRates_replicate", ifelse(i < 1001, i, paste0("_TRUEVALS_", i-1000)), "_c2.txt"))){
    r2 <- read.table(file = paste0("output_sim/fixCorrs_infRates_allTraits_rates_regRates_replicate", ifelse(i < 1001, i, paste0("_TRUEVALS_", i-1000)), "_c2.txt"))
  } else {
    r2 <- NA
  }
  if(file.exists(paste0("output_sim/fixCorrs_infRates_allTraits_rates_regRates_replicate", ifelse(i < 1001, i, paste0("_TRUEVALS_", i-1000)), "_c3.txt"))){
    r3 <- read.table(file = paste0("output_sim/fixCorrs_infRates_allTraits_rates_regRates_replicate", ifelse(i < 1001, i, paste0("_TRUEVALS_", i-1000)), "_c3.txt"))
  } else {
    r3 <- NA
  }
  if(file.exists(paste0("output_sim/fixCorrs_infRates_allTraits_rates_regRates_replicate", ifelse(i < 1001, i, paste0("_TRUEVALS_", i-1000)), "_c4.txt"))){
    r4 <- read.table(file = paste0("output_sim/fixCorrs_infRates_allTraits_rates_regRates_replicate", ifelse(i < 1001, i, paste0("_TRUEVALS_", i-1000)), "_c4.txt"))
  } else {
    r4 <- NA
  }
  
  if(any(!is.na(c(r1,r2,r3,r4)))){
    rates <- list(r1,r2,r3,r4)
    if(any(is.na(rates))){
      rates <- rates[!is.na(rates)]
    }
    rates <- do.call(rbind, rates)
    compRatesChains[[i_ind]] <- effectiveSize(rates)
    ratePercentile[[i_ind]] <- sapply(1:length(true_rates), function(trait) sum(true_rates[trait] > rates[,trait]) / length(rates[,trait]))
    rateCorrs[i_ind] <- cor(apply(rates, 2, mean), true_rates)
  } else {
    compRatesChains[[i_ind]] <- NA
    ratePercentile[[i_ind]] <- NA
    rateCorrs[i_ind] <- NA
  }

}

source("~/scripts/bailey_processingSimStudy_estMuTrueCorr.R")
source("~/scripts/bailey_processingSimStudy_truMuEstCorr.R")

sum(sapply(1:500, function(x) all(compTreesChains_truMu[[x]] > 0.99)))

ratePercentile_truMu
compRatesChains_truMu
rateCorrs_truMu
post_prob_comp_truMu
compTreesChains_truMu
fails_truMu <- which(!sapply(1:length(compTreesChains_truMu), function(comper) all(compTreesChains_truMu[[comper]]^2 > 0.99)))
fails_truMu <- unique(c(fails_truMu, which(sapply(1:length(compRatesChains_truMu), function(chain) !all(compRatesChains_truMu[[chain]] > 1000) | any(is.na(compRatesChains_truMu[[chain]]))))))

ratePercentile_truCorr
compRatesChains_truCorr
rateCorrs_truCorr
post_prob_comp_truCorr
compTreesChains_truCorr
fails_truCorr <- which(!sapply(1:length(compTreesChains_truCorr), function(comper) all(compTreesChains_truCorr[[comper]]^2 > 0.99)))
fails_truCorr <- unique(c(fails_truCorr, which(sapply(1:length(compRatesChains_truCorr), function(chain) !all(compRatesChains_truCorr[[chain]] > 1000) | any(is.na(compRatesChains_truCorr[[chain]]))))))


fails <- which(sapply(1:length(compRatesChains), function(chain) !all(compRatesChains[[chain]] > 1000) | any(is.na(compRatesChains[[chain]]))))
rate_percs_1 <- ratePercentile[1:500]
rate_percs_1 <- unlist(rate_percs_1[-(fails[fails > 0 & fails < 501])])
rate_percs_2 <- ratePercentile[501:1000]
rate_percs_2 <- unlist(rate_percs_2[-(fails[fails > 500 & fails < 1001] - 500)])
rate_percs_3 <- ratePercentile[1001:1500]
rate_percs_3 <- unlist(rate_percs_3[-(fails[fails > 1000 & fails < 1501] - 1500)])

rate_corrs_1 <- rateCorrs[1:500]
rate_corrs_1 <- (rate_corrs_1[-fails[fails > 0 & fails < 501]])
rate_corrs_2 <- rateCorrs[501:1000]
rate_corrs_2 <- (rate_corrs_2[-(fails[fails > 500 & fails < 1001] - 500)])
rate_corrs_3 <- rateCorrs[1001:1500]
rate_corrs_3 <- (rate_corrs_3[-(fails[fails > 1000 & fails < 1501] - 1500)])

hist(rate_percs_1)
hist(rate_percs_2)
hist(rate_percs_3)

hist(rate_corrs_1)
hist(rate_corrs_2)
hist(rate_corrs_3)


fails <- which(!sapply(1:length(compTreesChains), function(comper) all(compTreesChains[[comper]]^2 > 0.99)))
post_probs_1 <- post_prob_comp[1:500]
post_probs_1 <- post_probs_1[-(fails[fails < 501])]
post_probs_2 <- post_prob_comp[501:1000]
post_probs_2 <- post_probs_2[-(fails[fails > 500 & fails < 1001] - 500)]
post_probs_3 <- post_prob_comp[1001:1500]
post_probs_3 <- post_probs_3[-(fails[fails > 1000 & fails < 1501] - 500)]

post_probs_1 <- post_prob_comp[1:500]
post_probs_1 <- post_probs_1[!is.na(post_probs_1)]
post_probs_2 <- post_prob_comp[501:1000]
post_probs_2 <- post_probs_2[!is.na(post_probs_2)]
post_probs_3 <- post_prob_comp[1001:1500]

post_probs_all <- list(post_probs_1, post_probs_2, post_probs_3)
rate_corrs_all <- rev(list(rate_corrs_1, rate_corrs_2, rate_corrs_3))
rate_percs_all <- rev(list(rate_percs_1, rate_percs_2, rate_percs_3))

if(F){
  post_probs_all <- list(post_probs_1, post_probs_2, post_probs_3, post_prob_comp_truMu[-fails_truMu], post_prob_comp_truCorr[-fails_truCorr])
  rate_corrs_all <- rev(list(rate_corrs_1, rate_corrs_2, rate_corrs_3, rateCorrs_truMu[-fails_truMu], rateCorrs_truCorr[-fails_truCorr]))
  rate_percs_all <- rev(list(rate_percs_1, rate_percs_2, rate_percs_3, ratePercentile_truMu[-fails_truMu], ratePercentile_truCorr[-fails_truCorr]))
}

grDevices::cairo_pdf(filename = "/Volumes/macOS/Users/nikolai/dissertation/figures/chpt4_figure8_2.pdf", width = 2250 / 72, height = 750 / 72)
layout(mat = matrix(c(1,2,5,1,3,5,1,4,5), nrow = 3, ncol = 3, byrow = T))
par(mar = c(8.5,8.5,6,2))
# par(mar = c(6,6,6,2), mfrow = c(1,2))

for(i in 1:length(post_probs_all)){
  pp_cat <- do.call(rbind, post_probs_all[[i]])
  pp_cat <- pp_cat[order(pp_cat[,1], decreasing = F),]
  interval_width = 0.1 + 2/30
  intervals <- seq(0, 1, by = interval_width)
  iv_means <- sapply(1:(length(intervals) - 1), function(iv) mean(pp_cat[pp_cat[,1] > intervals[iv] & pp_cat[,1] < intervals[iv+1],1]))
  iv_freqs_true <- sapply(1:(length(intervals) - 1), function(iv) mean(pp_cat[pp_cat[,1] > intervals[iv] & pp_cat[,1] < intervals[iv+1],2]))
  if(i == 1){
    plot(iv_means, iv_freqs_true, xlim = c(0,1.11), ylim = c(0,1), type = "l", lwd = 3, xaxt = "n", yaxt = "n", frame = F, xlab = "", ylab = "")
    box(which = "plot", lwd = 3)
    box(which = "figure", lwd = 2, lty = 2)
    axis(side = 1, tck = -0.02, lwd = 3, lwd.ticks = 3, cex.axis = 3, labels = rep("", 11), at = 0:10/10)
    mtext(text = 0:10/10, side = 1, line = 3, at = 0:10/10, cex = 2)
    axis(side = 2, tck = -0.02, lwd = 3, lwd.ticks = 3, cex.axis = 3, labels = rep("", 11), at = 0:10/10)
    mtext(text = 0:10/10, side = 2, line = 1.75, at = 0:10/10, cex = 2)
    title(main = "Calibration Curves", cex.main = 5.75)
    title(xlab = "Estimated Bipartition Probability", cex.lab = 4, line = 6.5)
    title(ylab = "True Bipartition Frequency", cex.lab = 4, line = 5)
    legend(x = "topleft", col = c(2,1), lwd = 3, lty = c(2,1), legend = c("1-to-1 line", "calibration curve"), box.lty = 2, box.lwd = 4, cex = 2.75)
    fig_label(text = "a)", cex = 5, shrinkX = 0.9, shrinkY = 0.98)
  } else {
    lines(iv_means, iv_freqs_true, lwd = 3)
  }
  segments(x0 = 0, x1 = 1, y0  = 0, y1 = 1, lwd = 3, col = 2, lty = 2)
  if(i == 3){
    text(x = iv_means[length(iv_means)] + 0.083, y = iv_freqs_true[length(iv_freqs_true)] + 0.025, labels = "known means,\ncorrelations", cex = 2.1)
  } else if(i == 2){
    text(x = iv_means[length(iv_means)] + 0.07, y = iv_freqs_true[length(iv_freqs_true)] + 0.025, labels = "high sample-size \nestimates", cex = 2.1)
  } else if(i == 1){
    text(x = iv_means[length(iv_means)] + 0.07, y = iv_freqs_true[length(iv_freqs_true)] - 0.015, labels = "empirical \nestimates", cex = 2.1)
  }
}
par(mar = c(8.5,9.5,6,1))
for(i in 1:length(rate_percs_all)){
  par(lwd = 3)
  hist(rate_percs_all[[i]], main = "", xaxt = "n", yaxt = "n", xlab = "", ylab = "", breaks = 10, freq = F)
  box(which = "figure", lwd = 2, lty = 2)
  fig_label(text = paste0(letters[i+1], ")"), cex = 5, shrinkX = 0.9, shrinkY = 0.9)
  axis(side = 1, tck = -0.1, lwd = 3, lwd.ticks = 3, cex.axis = 3, labels = rep("", 11), at = 0:10/10)
  mtext(text = 0:10/10, side = 1, line = 3, at = 0:10/10, cex = 2)
  title(xlab = "True Rate Percentile in Marginal Posterior Distribution", cex.lab = 2.5, line = 6.5)
  title(ylab = "Probability Mass   ", cex.lab = 2.5, line = 2)
  if(i == 1){
    title(main = "Known Means, Correlations", cex.main = 3)
    axis(side = 2, tck = -0.1, lwd = 3, lwd.ticks = 3, cex.axis = 3, labels = rep("", 4), at = 0:3/2.5, line = -2)
    mtext(text = 0:3/2.5, side = 2, line = -0.5, at = 0:3/2.5, cex = 1.3)
  } 
  if(i == 2){
    title(main = "High Sample-Size Estimates", cex.main = 3)
    axis(side = 2, tck = -0.1, lwd = 3, lwd.ticks = 3, cex.axis = 3, labels = rep("", 3), at = 0:2, line = -2)
    mtext(text = 0:2, side = 2, line = -0.5, at = 0:2, cex = 1.3)
  }
  if(i == 3){
    title(main = "Empirical Estimates", cex.main = 3)
    axis(side = 2, tck = -0.1, lwd = 3, lwd.ticks = 3, cex.axis = 3, labels = rep("", 3), at = 0:2, line = -2)
    mtext(text = 0:2, side = 2, line = -0.5, at = 0:2, cex = 1.3)
  }
}
par(mar = c(8.5,8.5,6,1))
for(i in 1:3){
  densobj <- density(x = rate_corrs_all[[i]]^2, adjust = 2, from = 0)
  if(i == 1){
    plot(densobj, xlim = c(0,1), ylim = c(0,12), lwd = 3, main = "", xaxt = "n", yaxt = "n", frame.plot = F, xlab = "", ylab = "")
    axis(side = 1, tck = -0.02, lwd = 3, lwd.ticks = 3, cex.axis = 3, labels = rep("", 11), at = 0:10/10)
    mtext(text = 0:10/10, side = 1, line = 3, at = 0:10/10, cex = 2)
    axis(side = 2, tck = -0.02, lwd = 3, lwd.ticks = 3, cex.axis = 3, labels = rep("", 7), at = 0:6 * 2)
    mtext(text = 0:6 * 2, side = 2, line = 1.75, at = 0:6 * 2, cex = 2)

    box("plot", lwd = 3, lty = 1)
    title(xlab = latex2exp::TeX(paste0("R$^2$ Across Replicates")), cex.lab = 4, line = 7.5)
    shade(densobj, col = rgb(0,0,0,0.2), lim = c(0, 2))
    fig_label(text = paste0("e", ")"), cex = 5, shrinkX = 0.925, shrinkY = 0.975)
    title(main = latex2exp::TeX(paste0("R$^2$ between True Rates & Posterior Mean Rates")), cex.main = 3.75)
  } 
  if(i == 2){
    lines(densobj, lwd = 3)
    shade(densobj, col = rgb(0,0,0,0.5), lim = c(0, 2))
  }
  if(i == 3){
    lines(densobj, lwd = 3)
    shade(densobj, col = rgb(0,0,0,0.8), lim = c(0, 2))
  }
  if(i == 1){
    title(ylab = "Density", cex.lab = 4, line = 5)
    legend(x = "topleft", fill = c(rgb(0,0,0,0.2), rgb(0,0,0,0.5), rgb(0,0,0,0.8)), 
           legend = c("known means, correlations", "high sample-size estimates", "empirical estimates"), box.lty = 2, box.lwd = 4, cex = 3)    
  }
  box(which = "figure", lwd = 2, lty = 2)
}
dev.off()