library(mvtnorm)
library(pbivnorm)

rho <- runif(1,0,1)
upper <- runif(2,1,3)
lower <- runif(2,-2,1)

rho <- rhos[5]
upper <- uppers[79,]
lower <- lowers[79,]

# pbivnorm(x = upper[1], y = upper[1], rho = rho)
pbivnorm_infs(uppers[5,1], uppers[5,2], 0.94)
pbivnorm_infs(upper1 = uppers[,1], upper2 = uppers[,2], rhos = rep(0.99, length(uppers[,1])))



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
}

a <- pbivnorm_infs(uppers[,1], uppers[,2], rhos) - 
  pbivnorm_infs(lowers[,1], uppers[,2], rhos) - 
  pbivnorm_infs(uppers[,1], lowers[,2], rhos) + 
  pbivnorm_infs(lowers[,1], lowers[,2], rhos)

pbivnorm_infs(uppers[15,1], uppers[15,2], rhos[15]) - 
  pbivnorm_infs(lowers[15,1], uppers[15,2], rhos[15]) - 
  pbivnorm_infs(uppers[15,1], lowers[15,2], rhos[15]) + 
  pbivnorm_infs(lowers[15,1], lowers[15,2], rhos[15])

log(pbivnorm_infs(uppers[,1], uppers[,2], rhos) - 
      pbivnorm_infs(lowers[,1], uppers[,2], rhos) - 
      pbivnorm_infs(uppers[,1], lowers[,2], rhos) + 
      pbivnorm_infs(lowers[,1], lowers[,2], rhos))

rhos <- -1
pbivnorm(upper1, upper2, rhos)
pbivnorm_infs(upper1, upper2, rhos)
pmvnorm(upper = c(100, upper2[5]), lower = c(-Inf, -Inf), corr = matrix(c(1,rho[5],rho[5],1),2,2), mean = c(0,0))[1]

rho <- -100:100/100
for(i in 1:length(rho)){
  print(pmvnorm(upper = c(100, upper2), lower = c(-Inf, -Inf), corr = matrix(c(1,rho[i],rho[i],1),2,2), mean = c(0,0))[1])
}
pmv <- pmvnorm(upper = upper, lower = lower, corr = matrix(c(1,rho,rho,1),2,2), mean = c(0,0))[1]
pbiv1 <- pbivnorm(upper[1], upper[2], rho) - pbivnorm(lower[1], upper[2], rho) - pbivnorm(upper[1], lower[2], rho) + pbivnorm(lower[1], lower[2], rho)
pbiv1 <- pbivnorm_infs(upper[1], upper[2], rho) - pbivnorm_infs(lower[1], upper[2], rho) - pbivnorm_infs(upper[1], lower[2], rho) + pbivnorm_infs(lower[1], lower[2], rho)
pbiv2 <- pbivnorm(upper[1], upper[2], rho) - pbivnorm(lower[1], Inf, rho) - pbivnorm(Inf, lower[2], rho) + pbivnorm(lower[1], lower[2], rho)

pmv - pbiv1
pmv - pbiv2

nrep <- 150 * 10 * 10

mbm <- microbenchmark(pbivnorm(rep(upper[1], nrep), rep(upper[2], nrep), rho) - pbivnorm(rep(lower[1], nrep), rep(upper[2], nrep), rho) - pbivnorm(rep(upper[1], nrep), rep(lower[2], nrep), rho) + pbivnorm(rep(lower[1], nrep), rep(lower[2], nrep), rho), 
                      replicate(nrep, pmvnorm(upper = upper, lower = lower, corr = matrix(c(1,rho,rho,1),2,2), mean = c(0,0))[1]))

print(mbm)

.Fortran("PBIVNORM", prob, lower = as.double(-100,-100), uppera = as.double(0,0), upperb = as.double(0,0), infin, rho = 0.2, lt, PACKAGE = "pbivnorm")[[1]]

pbivnorm2 <- function (x, y, rho = 0, recycle = TRUE) 
{
  if (length(dim(x))) {
    if (ncol(x) != 2) 
      stop("'x' must have two columns if specified as a matrix")
    if (!missing(y) && !is.null(y)) 
      warning("'x' was specified as a matrix, so 'y' will be ignored")
    y <- x[, 2]
    x <- x[, 1]
  }
  if (any(abs(rho) > 1)) 
    stop("'rho' must be a valid correlation (-1 <= rho <= 1)")
  x <- replace(x, x == Inf, .Machine$double.xmax)
  x <- replace(x, x == -Inf, -.Machine$double.xmax)
  y <- replace(y, y == Inf, .Machine$double.xmax)
  y <- replace(y, y == -Inf, -.Machine$double.xmax)
  lengths <- sapply(list(x, y, rho), length)
  if (recycle) {
    x <- rep(x, length.out = max(lengths))
    y <- rep(y, length.out = max(lengths))
    rho <- rep(rho, length.out = max(lengths))
  } else if (diff(range(lengths)) > 0) {
    stop("'x', 'y', and 'rho' must be the same length when recycling is disallowed")
  }
  lower <- as.double(c(0, 0))
  infin <- as.integer(c(0, 0))
  uppera <- as.double(x)
  upperb <- as.double(y)
  lt <- as.integer(length(x))
  prob <- double(lt)
  rho <- as.double(rho)
  ans <- .Fortran("PBIVNORM", prob, lower, uppera, upperb, infin, rho, lt, PACKAGE = "pbivnorm")[[1]]
  return(ans)
}

f1 <- function() pmvnorm(upper = up, lower = low, corr = matrix(c(1,rho,rho,1),2,2), mean = c(0,0), algorithm = "GenzBretz")
f2 <- function() pmvnorm(upper = up, lower = low, corr = matrix(c(1,rho,rho,1),2,2), mean = c(0,0), algorithm = "Miwa")

f1()
f2()
microbenchmark(f1, f2)

x <- rnorm(10)
y <- rnorm(10)
rho <- runif(-1,1,10)

plot(-10:10/10, pbivnorm(0, 3,-10:10/10))

f1 <- function() {
  uniqueSitePatterns <- sapply(1:(d_trait-1), function(trait) uniqueSP_redux(dat[dat$tip == mean_tip, c(mean_trait, (1:d_trait)[-mean_trait][trait])], dat[dat$tip == mean_tip,d_trait+1]))
  nOrds <- sapply(1:(d_trait-1), function(trait) nrow(uniqueSitePatterns[1,trait][[1]]))
  nOrds_traits <- rep(1:(d_trait-1), nOrds)
  new_std_thresholds <- cbind(rep(-Inf, d_trait), sapply(1:ncol(threshold_mat), function(col) threshold_mat[,col] - means[mean_tip,]), rep(Inf, d_trait))
  SPs <- do.call(rbind, uniqueSitePatterns[1,])
  uppers <- t(sapply(1:length(nOrds_traits), function(pat) c(new_std_thresholds[mean_trait,SPs[pat,1]+2], new_std_thresholds[-mean_trait,][nOrds_traits[pat],SPs[pat,2]+2])))
  lowers <- t(sapply(1:length(nOrds_traits), function(pat) c(new_std_thresholds[mean_trait,SPs[pat,1]+1], new_std_thresholds[-mean_trait,][nOrds_traits[pat],SPs[pat,2]+1])))
  rhos = rep(corrs_of_trait, nOrds)
  counts = unlist(uniqueSitePatterns[2,])
  logLL <- -sum(log(pbivnorm_infs(uppers[,1], uppers[,2], rhos) - pbivnorm_infs(lowers[,1], uppers[,2], rhos) - pbivnorm_infs(uppers[,1], lowers[,2], rhos) + pbivnorm_infs(lowers[,1], lowers[,2], rhos)) * counts)
  return(logLL)
}

f2 <- function() {-sum(unlist(sapply(1:(d_trait-1), function(trait) log(multi_pmvnorm_optim_multithresh(mean = c(par_to_opt, cont_traits_of_tip[trait]),
                                                                                                        sigma = matrix(c(1,corrs_of_trait[trait],corrs_of_trait[trait],1),2,2),
                                                                                                        ordinals = uniqueSP_redux(dat[dat$tip == mean_tip, c(mean_trait, (1:d_trait)[-mean_trait][trait])], dat[dat$tip == mean_tip,d_trait+1])$SPs,
                                                                                                        thresholds = rbind(threshold_mat[mean_trait,], threshs_of_other_traits[trait,]))) * 
                                       uniqueSP_redux(dat[dat$tip == mean_tip, c(mean_trait, (1:d_trait)[-mean_trait][trait])], dat[dat$tip == mean_tip,d_trait+1])$counts)))
}

f3 <- function() {-sum(log(pbivnorm_infs(uppers[,1], uppers[,2], rhos) - pbivnorm_infs(lowers[,1], uppers[,2], rhos) - pbivnorm_infs(uppers[,1], lowers[,2], rhos) + pbivnorm_infs(lowers[,1], lowers[,2], rhos)) * counts)}

microbenchmark(f1(), f2(), f3(), times = 10)
f1()
f2()

microbenchmark(  uniqueSitePatterns <- sapply(1:(d_trait-1), function(trait) uniqueSP_fast(dat[dat$tip == mean_tip, c(mean_trait, (1:d_trait)[-mean_trait][trait])], dat[dat$tip == mean_tip,d_trait+1])),
                 nOrds <- sapply(1:(d_trait-1), function(trait) nrow(uniqueSitePatterns[1,trait][[1]])),
                 nOrds_traits <- rep(1:(d_trait-1), nOrds),
                 new_std_thresholds <- cbind(rep(-Inf, d_trait), sapply(1:ncol(threshold_mat), function(col) threshold_mat[,col] - means[mean_tip,]), rep(Inf, d_trait)),
                 SPs <- do.call(rbind, uniqueSitePatterns[1,]),
                 uppers <- t(sapply(1:length(nOrds_traits), function(pat) c(new_std_thresholds[mean_trait,SPs[pat,1]+2], new_std_thresholds[-mean_trait,][nOrds_traits[pat],SPs[pat,2]+2]))),
                 lowers <- t(sapply(1:length(nOrds_traits), function(pat) c(new_std_thresholds[mean_trait,SPs[pat,1]+1], new_std_thresholds[-mean_trait,][nOrds_traits[pat],SPs[pat,2]+1]))),
                 rhos = rep(corrs_of_trait, nOrds),
                 counts = unlist(uniqueSitePatterns[2,]),
                 logLL <- -sum(log(pbivnorm_infs(uppers[,1], uppers[,2], rhos) - pbivnorm_infs(lowers[,1], uppers[,2], rhos) - pbivnorm_infs(uppers[,1], lowers[,2], rhos) + pbivnorm_infs(lowers[,1], lowers[,2], rhos)) * counts),
                 times = 10)

#uniqueSP_redux most time consuming -- maybe preprocess?
#otherwise computing upper and lowe bounds more time consuming

microbenchmark(uniqueSitePatterns <- sapply(1:(d_trait-1), function(trait) uniqueSP_redux(dat[dat$tip == mean_tip, c(mean_trait, (1:d_trait)[-mean_trait][trait])], dat[dat$tip == mean_tip,d_trait+1])), times = 10)

sapply(1:(d_trait-1), function(trait) uniqueSP_redux(dat[dat$tip == mean_tip, c(mean_trait, (1:d_trait)[-mean_trait][trait])], dat[dat$tip == mean_tip,d_trait+1]))

uniqueSP_redux <- function(indivs, freqs){
  indivs_str <- sapply(1:length(indivs[,1]), function(indiv) paste0(indivs[indiv,], collapse = "")) 
  uniq_indivs_str <- unique(indivs_str)
  counts <- sapply(1:length(uniq_indivs_str), function(string) sum(freqs[indivs_str == uniq_indivs_str[string]]))
  SPs <- t(apply((do.call(rbind, (strsplit(uniq_indivs_str, "")))), 1, as.numeric))
  return(list(SPs = SPs, counts = counts))
}

uniqueSP_fast <- function(indivs, freqs){
  newFreqs <- as.data.frame(wtd.table(indivs[,1], indivs[,2], weights = freqs))
  newFreqs <- (newFreqs[newFreqs$Freq > 0,])
  newFreqs$Var1 <- as.numeric(levels(newFreqs$Var1))[newFreqs$Var1]
  newFreqs$Var2 <- as.numeric(levels(newFreqs$Var2))[newFreqs$Var2]
  return(list(SPs = as.matrix(newFreqs[,1:2]), counts = newFreqs[,3]))
}


microbenchmark(uniqueSP_redux(dat[dat$tip == mean_tip, c(mean_trait, (1:d_trait)[-mean_trait][trait])], dat[dat$tip == mean_tip,d_trait+1]),
               uniqueSP_fast(dat[dat$tip == mean_tip, c(mean_trait, (1:d_trait)[-mean_trait][trait])], dat[dat$tip == mean_tip,d_trait+1]))

microbenchmark(uniqueSitePatterns <- sapply(1:(d_trait-1), function(trait) uniqueSP_redux(dat[dat$tip == mean_tip, c(mean_trait, (1:d_trait)[-mean_trait][trait])], dat[dat$tip == mean_tip,d_trait+1])), 
               uniqueSitePatterns2 <- sapply(1:(d_trait-1), function(trait) uniqueSP_fast(dat[dat$tip == mean_tip, c(mean_trait, (1:d_trait)[-mean_trait][trait])], dat[dat$tip == mean_tip,d_trait+1])),
               times = 10)

microbenchmark(uniqueSitePatterns <- lapply(1:nTaxa, function(tip) sapply(1:(d_trait-1), function(trait) uniqueSP_fast(dat[dat$tip == tip, c(thresh_trait, (1:d_trait)[-thresh_trait][trait])], dat[dat$tip == tip,d_trait+1]))),
               nOrds <- lapply(1:nTaxa, function(tip) sapply(1:(d_trait-1), function(trait) nrow(uniqueSitePatterns[[tip]][1,trait][[1]]))),
               nOrds_traits <- lapply(1:nTaxa, function(tip) rep(1:(d_trait-1), nOrds[[tip]])),
               new_std_thresholds <- lapply(1:nTaxa, function(tip) cbind(rep(-Inf, d_trait), sapply(1:ncol(threshold_mat), function(col) threshold_mat[,col] - means[tip,]), rep(Inf, d_trait))),
               SPs <- lapply(1:nTaxa, function(tip) do.call(rbind, uniqueSitePatterns[[tip]][1,])),
               uppers_lowers <- lapply(1:nTaxa, function(tip) t(sapply(1:length(nOrds_traits[[tip]]), 
                                                                       function(pat) c(new_std_thresholds[[tip]][thresh_trait, SPs[[tip]][pat,1]+c(2,1)], 
                                                                                       new_std_thresholds[[tip]][nOrds_traits[[tip]][pat] + ifelse(nOrds_traits[[tip]][pat] < thresh_trait, 0, 1), SPs[[tip]][pat,2]+c(2,1)])))),
               uppers <- do.call(rbind, lapply(1:nTaxa, function(tip) uppers_lowers[[tip]][,c(1,3)])),
               lowers <- do.call(rbind, lapply(1:nTaxa, function(tip) uppers_lowers[[tip]][,c(2,4)])),
               rhos = corrs_of_thresh[unlist(nOrds_traits)],
               counts = unlist(lapply(1:nTaxa, function(tip) unlist(uniqueSitePatterns[[tip]][2,]))),
               logLL <- -sum(log(pbivnorm_infs(uppers[,1], uppers[,2], rhos) - 
                                   pbivnorm_infs(lowers[,1], uppers[,2], rhos) - 
                                   pbivnorm_infs(uppers[,1], lowers[,2], rhos) + 
                                   pbivnorm_infs(lowers[,1], lowers[,2], rhos)) * counts),
               times = 10)


microbenchmark(dat <- list(dat_orig, k_par, prunes_init, rownames(traits), n_thresh),
               n_thresh <- dat[[5]],
               tipNames <- dat[[4]],
               prunes <- dat[[3]],
               k_par <- dat[[2]],
               dat2 <- dat[[1]])


pbivnorm_infs(uppers[,1], uppers[,2], rhos = rep(0.94, length(uppers[,1])))
pbivnorm_infs(uppers[5,1], uppers[5,2], rhos[5])
