# CREATING REVBAYES REPLICATE SCRIPTS
setwd("/Volumes/2TB/mvBM_sims_Howells_8tips")

library(ape)
library(coda)
library(rwty)
library(phangorn)
library(lattice)
library(latticeExtra)
library(MCMCpack)
library(mvMORPH)
library(phytools)


#################################################
###### Getting Modern Human Var-Cov Matrix ######
#################################################

d.orig <- read.csv("C:\\Users\\Nikolai\\Documents\\data\\Howell.csv")
d <- d.orig
# str(d)
# d[d == 0] <- NA
# d <- d[complete.cases(d),]
men <- d[d$Sex == "M",]
#men <- d
#head(as.matrix(men))
#unique(men$Population)[1]
pops <- sapply(1:length(unique(men$Population)), function (x) as.matrix((t(men[men$Population == unique(men$Population)[x],-(1:4)]))))

#get traits
traits <- list()
for(i in 1:82){ 
  print(i)
  trait <- list()
  trait[[1]] <- pops[[1]][i,]
  for (j in 2:30){
    trait[[j]] <- pops[[j]][i,]
  } 
  traits[[i]] <- trait
}

#make cov matrix
cov <- matrix(nrow = 82, ncol = 82)
for (i in 1:82){
  print(i)
  for(j in 1:82){
    trait1 <- traits[[i]]
    trait2 <- traits[[j]]
    trait1means <- sapply(1:length(trait1), function(x) mean(trait1[[x]]))
    trait2means <- sapply(1:length(trait2), function(x) mean(trait2[[x]]))
    trait1deviations <- sapply(1:length(trait1), function(x) trait1[[x]] - trait1means[x])
    trait2deviations <- sapply(1:length(trait2), function(x) trait2[[x]] - trait2means[x])
    deviationsProducts <- sapply(1:length(trait1), function(x) trait1deviations[[x]] * trait2deviations[[x]])
    sumALL <- 0
    for (k in (1:length(deviationsProducts))){
      sumALL <- sumALL + sum(deviationsProducts[[k]])
    }
    indiv <- sum(sapply(1:length(trait1), function(x) length(trait1[[x]])))
    covariance <- sumALL/(indiv - length(deviationsProducts))
    cov[i,j] <- covariance
  }
}

npop <- length(pops)
#compute matrix of population means
traitsHOWELL <- matrix(nrow=npop, ncol=82)
for(i in 1:npop){
  for(j in 1:82){
    traitsHOWELL[i,j] <- mean(pops[[i]][j,])
  }
}
rownames(traitsHOWELL) <- unique(men$Population)
colnames(traitsHOWELL) <- colnames(men)[-(1:4)]

rownames(cov) <- colnames(cov) <- colnames(men)[-(1:4)]
covHuman <- cov

#taken from roseman and weaver 2004, only BPC -> BPL (due to a typo? it's the basion-prosthion length )
linMeasNames <- c("GOL", "NOL", "BNL", "BBH", "XCB", "XFB", "STB", "ZYB", "AUB", "WCB", "ASB", "BPL", "NPH", "NLH", "OBH", "OBB", "JUB", "NLB", "MAB", "MDH", "MDB", "ZMB", "SSS", "FMB", "NAS", "EKB", "DKB", "NDS", "WNB", "SIS", "IML", "XML", "MLS", "WMH", "FOL", "FRC", "FRS", "FRF", "PAC", "PAS", "PAF", "OCC", "OCS", "OCF", "VRR", "NAR", "SSR", "PRR", "DKR", "ZOR", "FMR", "EKR", "ZMR", "AVR", "DKS", "SOS", "GLS")
linMeasTraitsHuman <- traitsHOWELL[,linMeasNames]
covFULL <- linMeasCovHuman <- covHuman[linMeasNames, linMeasNames]

##########################################
###### Getting Chimp Var-Cov Matrix ######
##########################################


############# DATA PREPARATION #############

d.orig <- read.csv("C:\\Users\\Nikolai\\Documents\\data\\PanSorted.csv")
d <- d.orig
men <- d[d$Sex == "M",]
d <- men
str(d)


head(as.matrix(d))
unique(d$Taxon)[1]
sp <- sapply(1:length(unique(d$Taxon)), function (x) as.matrix((t(d[d$Taxon == unique(d$Taxon)[x],-(1:2)]))))
ntraits <- 27
nspecies <- 4

#get traits
traits <- list()
for(i in 1:ntraits){ 
  print(i)
  trait <- list()
  trait[[1]] <- sp[[1]][i,]
  for (j in 2:nspecies){
    trait[[j]] <- sp[[j]][i,]
  } 
  traits[[i]] <- trait
}

#make cov matrix
cov <- matrix(nrow = ntraits, ncol = ntraits)
for (i in 1:ntraits){
  print(i)
  for(j in 1:ntraits){
    trait1 <- traits[[i]]
    trait2 <- traits[[j]]
    trait1means <- sapply(1:length(trait1), function(x) mean(trait1[[x]]))
    trait2means <- sapply(1:length(trait2), function(x) mean(trait2[[x]]))
    trait1deviations <- sapply(1:length(trait1), function(x) trait1[[x]] - trait1means[x])
    trait2deviations <- sapply(1:length(trait2), function(x) trait2[[x]] - trait2means[x])
    deviationsProducts <- sapply(1:length(trait1), function(x) trait1deviations[[x]] * trait2deviations[[x]])
    sumALL <- 0
    for (k in (1:length(deviationsProducts))){
      sumALL <- sumALL + sum(deviationsProducts[[k]])
    }
    indiv <- sum(sapply(1:length(trait1), function(x) length(trait1[[x]])))
    covariance <- sumALL/(indiv - length(deviationsProducts))
    cov[i,j] <- covariance
  }
}

offDiagCorr <- c(cov2cor(cov))
offDiagCorr <- offDiagCorr[offDiagCorr < 1] 
par(mfrow = c(1,2))
hist(offDiagCorr, main = "")
title("Correlations Between Traits")
hist(diag(cov), main = "")
title("Variances of Traits")

npop <- length(sp)
#compute matrix of Species means
traitsCHIMP <- matrix(nrow=npop, ncol=ntraits)
for(i in 1:npop){
  for(j in 1:ntraits){
    traitsCHIMP[i,j] <- mean(sp[[i]][j,])
  }
}
rownames(traitsCHIMP) <- unique(d$Taxon)
colnames(traitsCHIMP) <- colnames(d)[-(1:2)]

rownames(cov) <- colnames(cov) <- colnames(d)[-(1:2)]
covChimp <- cov



################################################################
# Mitteroecker and Bookstein Distance for Regional Differences #
################################################################

pops <- levels(men$Population)
#DEFINE REGIONS
ASIA <- c("AINU", "ATAYAL", "HAINAN", "ANYANG", "PHILLIPI", "S JAPAN", "N JAPAN", "BURIAT", "ANDAMAN")
AMERICAS <- c("ESKIMO", "ARIKARA", "SANTA CR", "PERU")
AUSTRALASIA <- c("TOLAI", "TASMANIA", "AUSTRALI")
EUROPE_nAFRICA <- c("NORSE", "EGYPT", "ZALAVAR", "BERG")
SUBSAHARAN_AFRICA <- c("TEITA", "BUSHMAN", "ZULU", "DOGON")
POLYNESIA_MICRONESIA <- c("GUAM", "N MAORI", "S MAORI", "MORIORI", "MOKAPU", "EASTER I")
regions <- list(ASIA, AMERICAS, AUSTRALASIA, EUROPE_nAFRICA, SUBSAHARAN_AFRICA, POLYNESIA_MICRONESIA)

linMeasNames <- c("GOL", "NOL", "BNL", "BBH", "XCB", "XFB", "STB", "ZYB", "AUB", "WCB", "ASB", "BPL", "NPH", "NLH", "OBH", "OBB", "JUB", "NLB", "MAB", "MDH", "MDB", "ZMB", "SSS", "FMB", "NAS", "EKB", "DKB", "NDS", "WNB", "SIS", "IML", "XML", "MLS", "WMH", "FOL", "FRC", "FRS", "FRF", "PAC", "PAS", "PAF", "OCC", "OCS", "OCF", "VRR", "NAR", "SSR", "PRR", "DKR", "ZOR", "FMR", "EKR", "ZMR", "AVR", "DKS", "SOS", "GLS")

d <- d.orig <- read.csv("C:\\Users\\Nikolai\\Documents\\data\\Howell.csv")
regionCovariances <- list()
for (x in 1:6){
  print(x)
  
  region1 <- regions[[x]]
  
  men <- d[d$Sex == "M",]
  
  region1men <- data.frame()
  for (k in 1:length(region1)){
    region1men <- rbind(region1men, men[men$Population == region1[k],])
  }
  
  pops1 <- sapply(1:length(unique(region1men$Population)), function (x) as.matrix((t(region1men[region1men$Population == unique(region1men$Population)[x],-(1:4)]))))
  
  #get traits
  traits1 <- list()
  for(i in 1:82){ 
    trait <- list()
    trait[[1]] <- pops1[[1]][i,]
    for (j in 2:length(region1)){
      trait[[j]] <- pops1[[j]][i,]
    } 
    traits1[[i]] <- trait
  }
  
  #make cov matrix
  cov1 <- matrix(nrow = 82, ncol = 82)
  for (i in 1:82){
    for(j in 1:82){
      trait1 <- traits1[[i]]
      trait2 <- traits1[[j]]
      trait1means <- sapply(1:length(trait1), function(x) mean(trait1[[x]]))
      trait2means <- sapply(1:length(trait2), function(x) mean(trait2[[x]]))
      trait1deviations <- sapply(1:length(trait1), function(x) trait1[[x]] - trait1means[x])
      trait2deviations <- sapply(1:length(trait2), function(x) trait2[[x]] - trait2means[x])
      deviationsProducts <- sapply(1:length(trait1), function(x) trait1deviations[[x]] * trait2deviations[[x]])
      sumALL <- 0
      for (k in (1:length(deviationsProducts))){
        sumALL <- sumALL + sum(deviationsProducts[[k]])
      }
      indiv <- sum(sapply(1:length(trait1), function(x) length(trait1[[x]])))
      covariance <- sumALL/(indiv - length(deviationsProducts))
      cov1[i,j] <- covariance
    }
  }
  
  rownames(cov1) <- colnames(cov1) <- colnames(region1men)[-(1:4)]
  #taken from roseman and weaver 2004, only BPC -> BPL (due to a typo? it's the basion-prosthion length )
  linMeasNames <- c("GOL", "NOL", "BNL", "BBH", "XCB", "XFB", "STB", "ZYB", "AUB", "WCB", "ASB", "BPL", "NPH", "NLH", "OBH", "OBB", "JUB", "NLB", "MAB", "MDH", "MDB", "ZMB", "SSS", "FMB", "NAS", "EKB", "DKB", "NDS", "WNB", "SIS", "IML", "XML", "MLS", "WMH", "FOL", "FRC", "FRS", "FRF", "PAC", "PAS", "PAF", "OCC", "OCS", "OCF", "VRR", "NAR", "SSR", "PRR", "DKR", "ZOR", "FMR", "EKR", "ZMR", "AVR", "DKS", "SOS", "GLS")
  cov1 <- cov1[linMeasNames, linMeasNames]
  
  regionCovariances[[x]] <- cov1
}


avgDistRegAll <- vector()
for(i in 57:57){
  print(i)
  dist <- vector()
  for(j in 1:30000){
    subsample <- sample(x = rownames(regionCovariances[[1]]), size = i, replace = F)
    regions <- sample(x = 1:6, size = 2, replace = F)
    reduCov1 <- regionCovariances[[regions[1]]][subsample, subsample]
    reduCov2 <- regionCovariances[[regions[2]]][subsample, subsample]
    dist[j] <- sum(log(eigen(solve(reduCov1) %*% reduCov2)$values)^2)^0.5
  }
  avgDistRegAll[i] <- mean(dist)
}
plot(y= avgDistRegAll, pch = 1, x =  1:57, ylab = "Average Mitteroecker and Bookstein Distance", xlab = "Number of Traits in Subsample", main = "Regional Comparison")

numTraitsLM <- 1:57
regionalLM <- lm(formula = avgDistRegAll ~ numTraitsLM + I(numTraitsLM^2))
lines(numTraitsLM, fitted(regionalLM), col='red', type='b')
new <- data.frame(numTraitsLM = c(2,4,8,16,32,64))
modelReg_MBD <- predict(regionalLM, new, se.fit = F)
empiricalReg_MBD <- avgDistRegAll[c(2,4,8,16,32,64)]
targetRegionalMBD <- sapply(1:6, function (x) if(!is.na(empiricalReg_MBD[x])) { empiricalReg_MBD[x] } else { modelReg_MBD[x] })
targetRegionalMBD <- c(0.3016502, 0.5492471, 1.0165738, 1.9375714, 3.8145961, 7.8122877)
targetRegionalMBD <- c(0.3016502, 0.5492471, 1.0165738, 1.9375714, 3.8145961, 6.939882, 7.8122877) #incl 57


#####################################################################
## Mitteroecker and Bookstein Distance for Chimp-Human Differences ##
#####################################################################

sharedHumanCov <- covHuman[rownames(covChimp), rownames(covChimp)]
sharedChimpCov <- covChimp

#Mitteroecker and Bookstein Distance
avgDist <- vector()
for(i in 1:27){
  print(i)
  dist <- vector()
  for(j in 1:30000){
    subsample <- sample(x = rownames(sharedHumanCov), size = i, replace = F)
    reduSHC <- sharedHumanCov[subsample, subsample]
    reduSCC <- sharedChimpCov[subsample, subsample]
    dist[j] <- sum(log(eigen(solve(reduSHC) %*% reduSCC)$values)^2)^0.5
  }
  avgDist[i] <- mean(dist)
}
plot(y= avgDist, x =  1:27, ylab = "Average Mitteroecker and Bookstein Distance", xlab = "Number of Traits in Subsample", main = "Chimp-Human Comparison")

numTraitsLM <- 1:27
chimpHumanLM <- lm(formula = avgDist ~ numTraitsLM + I(numTraitsLM^2))
lines(numTraitsLM, fitted(chimpHumanLM), col='red', type='b')
new <- data.frame(numTraitsLM = c(2,4,8,16,32,64))
modelCH_MBD <- predict(chimpHumanLM, new, se.fit = F)
empiricalCH_MBD <- avgDist[c(2,4,8,16,32,64)]
targetChimpHumanMBD <- sapply(1:6, function (x) if(!is.na(empiricalCH_MBD[x])) { empiricalCH_MBD[x] } else { modelCH_MBD[x] })
targetChimpHumanMBD <- c(0.7828651, 1.2172830, 1.9157832, 3.1369000, 5.3419662, 8.0803367)

#incl 57 traits
numTraitsLM <- 1:27
plot(y= avgDist, x =  1:27, xlim = c(0,64), ylim = c(0,8), ylab = "Average Mitteroecker and Bookstein Distance", xlab = "Number of Traits in Subsample", main = "Chimp-Human Comparison")
chimpHumanLM <- lm(formula = avgDist ~ numTraitsLM + I(numTraitsLM^2))
traitCHnums <- data.frame(numTraitsLM = 1:64)
lines(1:64, predict(chimpHumanLM, traitCHnums, se.fit = F), col='red', type='b')
new <- data.frame(numTraitsLM = c(2,4,8,16,32,57,64))
modelCH_MBD <- predict(chimpHumanLM, new, se.fit = F)
empiricalCH_MBD <- avgDist[c(2,4,8,16,32,64)]
targetChimpHumanMBD <- sapply(1:7, function (x) if(!is.na(empiricalCH_MBD[x])) { empiricalCH_MBD[x] } else { modelCH_MBD[x] })
targetChimpHumanMBD <- c(0.7828651, 1.2172830, 1.9157832, 3.1369000, 5.3419662, 7.6690310, 8.0803367)


##########################################
############ USEFUL FUNCTIONS ############
##########################################

#write means to file function
meansNexus <- function (characterMatrix, outputFilePath, traitsOnColumns = T){
  if (traitsOnColumns) {ntax = dim(characterMatrix)[1]; nchar = dim(characterMatrix)[2]}
  else {ntax = dim(characterMatrix)[2]; nchar = dim(characterMatrix)[1]}
  sink(outputFilePath, append=F)
  cat("#NEXUS\n\nBegin data;\n\tDimensions ntax=", ntax, "nchar=", nchar, ";\n\tFormat datatype=Continuous missing=?;\n\tMatrix\n\n")
  for(i in 1:length(characterMatrix[,1])) {cat(rownames(characterMatrix)[i], characterMatrix[i,], "\n")}
  cat(";\nEnd;")
  sink()
}

# simulate traits
# library(phytools)
# tree <- rtree(n = 30, rooted = T)
# write.nexus(tree, file="A:\\Google Drive\\data\\testFuncTree.nex")
# library(mvMORPH)
# traits <- mvSIM(tree = tree, nsim = 1, model = "BM1", 
#                 param = list(ntraits = dim(cov)[1], sigma=round(cov, 4), mu= rep(0, times = dim(cov)[1])))

#write to file
# meansNexus(traits, "A:\\Google Drive\\data\\testFuncMeans.nex")

ratesPhy <- function (rateMatrix, outputFilePath) {
  if(!isSymmetric(rateMatrix)) {stop("Error: This matrix is not symmetric.")}
  sink(outputFilePath, append=F)
  cat(dim(rateMatrix)[1],"\n")
  for(i in 1:length(rateMatrix[,1])) {
    cat(rownames(rateMatrix)[i], paste0(as.vector(rateMatrix[i, (1:(length(rateMatrix[i,])-1))])))
    cat(" ")
    cat(rateMatrix[i, (length(rateMatrix[i,]))], "\n", sep = "")
  }
  sink()
}

ratesTSV <- function (rateMatrix, outputFilePath) {
  if(!isSymmetric(rateMatrix)) {stop("Error: This matrix is not symmetric.")}
  sink(outputFilePath, append=F)
  for(i in 1:length(rateMatrix[,1])) {
    cat(paste0(as.vector(rateMatrix[i, (1:(length(rateMatrix[i,])-1))])), sep = "\t")
    cat("\t")
    cat(rateMatrix[i, (length(rateMatrix[i,]))], "\n\n", sep = "")
  }
  sink()
}

# ratesPhy(round(linMeasCovHuman, 4), "A:\\Google Drive\\data\\testFuncRates.phy")


extendedOnionExpansion <- function(cov, d, eta){
  sds <- diag(cov)^0.5
  cor <- cov2cor(cov)
  a <- dim(cor)[1]
  beta <- eta + (d-2)/2
  beta <- beta - (1/2)*(a-2)
  for(k in (a+1):d){
    beta <- beta - 1/2
    y = rbeta(n = 1, shape1 = (k-1)/2, shape2 = beta)
    r <- y^0.5
    theta <- rnorm(n = (k-1), mean = 0.5, sd = 1)
    theta <- theta / (sum(theta^2)^0.5)
    w <- r * theta
    L <- t(chol(cor))
    q <- as.vector(L %*% w)
    traitNames <- c(rownames(cor), paste0("trait", k)) 
    cor <- rbind(cor, q)
    cor <- cbind(cor, c(q,1))
    rownames(cor) <- colnames(cor) <- traitNames
  }
  newSDS <- c(sds, sample(x = sds, size = (d-a), replace = T))
  ncov <- diag(newSDS) %*% cor %*% diag(newSDS)
  ncov
}

#########################################################################################
## figure out how many degrees of freedom are needed to draw from Wishart distribution ##
#########################################################################################

traitNumIter <- c(2,4,8,16,32,64)

wishartSamplingIntensity <- 5000

traitNumIter <- c(2,4,8,16,32,64)

dfMatrix <- matrix(nrow = length(traitNumIter), ncol = length(degreeMisspecification))

degreeMisspecification <- c("perfect", "interregional", "chimpHuman")

# k iterates over degree of misspecification
for (k in 2:length(degreeMisspecification)) {
  
  # j iterates over trait numbers
  for (j in 1:length(traitNumIter)) {
    
    print(c(k,j))
    
    
    rateMatrixSpecification <- degreeMisspecification[k]
    ntraits <- traitNumIter[j]
    
    if (ntraits < 58){
      traitsToUse <- sample(linMeasNames, size = ntraits, replace = F)
    } else {
      traitsToUse <- linMeasNames
    }
    
    cov <- covFULL[traitsToUse, traitsToUse]
    
    if (ntraits > 57){
      cov <- extendedOnionExpansion(cov = cov, d = ntraits, eta = 1)
    }
    
    if (rateMatrixSpecification == "interregional"){
      targetMBD <- targetRegionalMBD[j]
    } else if (rateMatrixSpecification == "chimpHuman"){
      targetMBD <- targetChimpHumanMBD[j]
    } 
    
    
    
    
    avgDist <- vector()
    dist <- vector()
    dfWish <- ntraits
    for(r in 1:wishartSamplingIntensity){
      sampCov <- rwish(dfWish, cov)/dfWish
      dist[r] <- sum(log(eigen(solve(sampCov) %*% cov)$values)^2)^0.5
    }
    avgDist[1] <- mean(dist)
    
    while(avgDist[length(avgDist)] > targetMBD){
      dfWish <- dfWish + 20
      dist <- vector()  
      for(r in 1:wishartSamplingIntensity){
        sampCov <- rwish(dfWish, cov)/dfWish
        dist[r] <- sum(log(eigen(solve(sampCov) %*% cov)$values)^2)^0.5
      }
      # print(mean(dist))
      # print(sd(dist)/sqrt(length(dist)))
      avgDist <- c(avgDist, mean(dist))
    }
    
    dfWish <- dfWish - 20
    avgDist <- avgDist[-length(avgDist)]
    
    while(avgDist[length(avgDist)] > targetMBD){
      dfWish <- dfWish + 10
      dist <- vector()  
      for(r in 1:wishartSamplingIntensity){
        sampCov <- rwish(dfWish, cov)/dfWish
        dist[r] <- sum(log(eigen(solve(sampCov) %*% cov)$values)^2)^0.5
      }
      #print(mean(dist))
      avgDist <- c(avgDist, mean(dist))
    }
    
    dfWish <- dfWish - 10
    avgDist <- avgDist[-length(avgDist)]
    
    while(avgDist[length(avgDist)] > targetMBD){
      dfWish <- dfWish + 5
      dist <- vector()  
      for(r in 1:wishartSamplingIntensity){
        sampCov <- rwish(dfWish, cov)/dfWish
        dist[r] <- sum(log(eigen(solve(sampCov) %*% cov)$values)^2)^0.5
      }
      #print(mean(dist))
      avgDist <- c(avgDist, mean(dist))
    }
    
    dfWish <- dfWish - 5
    avgDist <- avgDist[-length(avgDist)]
    
    while(avgDist[length(avgDist)] > targetMBD){
      dfWish <- dfWish + 1
      dist <- vector()  
      for(r in 1:wishartSamplingIntensity){
        sampCov <- rwish(dfWish, cov)/dfWish
        dist[r] <- sum(log(eigen(solve(sampCov) %*% cov)$values)^2)^0.5
      }
      #print(mean(dist))
      avgDist <- c(avgDist, mean(dist))
    }
    
    if (which.min(abs(avgDist - targetMBD)) != length(avgDist)) {
      dfWish <- dfWish - 1
    }
    
    dfMatrix[j,k] <- dfWish
    
  }
  
}

##########################################
## figure out df when traitnum includes 57
##########################################

traitNumIter <- c(2,4,8,16,32,57,64)

wishartSamplingIntensity <- 10000

dfMatrix <- matrix(nrow = length(traitNumIter), ncol = length(degreeMisspecification))

degreeMisspecification <- c("perfect", "interregional", "chimpHuman")

# k iterates over degree of misspecification
for (k in 2:length(degreeMisspecification)) {
  
  # j iterates over trait numbers
  for (j in 6:6) {
    
    print(c(k,j))
    
    
    rateMatrixSpecification <- degreeMisspecification[k]
    ntraits <- traitNumIter[j]
    
    if (ntraits < 58){
      traitsToUse <- sample(linMeasNames, size = ntraits, replace = F)
    } else {
      traitsToUse <- linMeasNames
    }
    
    cov <- covFULL[traitsToUse, traitsToUse]
    
    if (ntraits > 57){
      cov <- extendedOnionExpansion(cov = cov, d = ntraits, eta = 1)
    }
    
    if (rateMatrixSpecification == "interregional"){
      targetMBD <- targetRegionalMBD[j]
    } else if (rateMatrixSpecification == "chimpHuman"){
      targetMBD <- targetChimpHumanMBD[j]
    } 
    
    
    
    
    avgDist <- vector()
    dist <- vector()
    dfWish <- ntraits
    for(r in 1:wishartSamplingIntensity){
      sampCov <- rwish(dfWish, cov)/dfWish
      dist[r] <- sum(log(eigen(solve(sampCov) %*% cov)$values)^2)^0.5
    }
    avgDist[1] <- mean(dist)
    
    while(avgDist[length(avgDist)] > targetMBD){
      dfWish <- dfWish + 20
      dist <- vector()  
      for(r in 1:wishartSamplingIntensity){
        sampCov <- rwish(dfWish, cov)/dfWish
        dist[r] <- sum(log(eigen(solve(sampCov) %*% cov)$values)^2)^0.5
      }
      # print(mean(dist))
      # print(sd(dist)/sqrt(length(dist)))
      avgDist <- c(avgDist, mean(dist))
    }
    
    dfWish <- dfWish - 20
    avgDist <- avgDist[-length(avgDist)]
    
    while(avgDist[length(avgDist)] > targetMBD){
      dfWish <- dfWish + 10
      dist <- vector()  
      for(r in 1:wishartSamplingIntensity){
        sampCov <- rwish(dfWish, cov)/dfWish
        dist[r] <- sum(log(eigen(solve(sampCov) %*% cov)$values)^2)^0.5
      }
      #print(mean(dist))
      avgDist <- c(avgDist, mean(dist))
    }
    
    dfWish <- dfWish - 10
    avgDist <- avgDist[-length(avgDist)]
    
    while(avgDist[length(avgDist)] > targetMBD){
      dfWish <- dfWish + 5
      dist <- vector()  
      for(r in 1:wishartSamplingIntensity){
        sampCov <- rwish(dfWish, cov)/dfWish
        dist[r] <- sum(log(eigen(solve(sampCov) %*% cov)$values)^2)^0.5
      }
      #print(mean(dist))
      avgDist <- c(avgDist, mean(dist))
    }
    
    dfWish <- dfWish - 5
    avgDist <- avgDist[-length(avgDist)]
    
    while(avgDist[length(avgDist)] > targetMBD){
      dfWish <- dfWish + 1
      dist <- vector()  
      for(r in 1:wishartSamplingIntensity){
        sampCov <- rwish(dfWish, cov)/dfWish
        dist[r] <- sum(log(eigen(solve(sampCov) %*% cov)$values)^2)^0.5
      }
      #print(mean(dist))
      avgDist <- c(avgDist, mean(dist))
    }
    
    if (which.min(abs(avgDist - targetMBD)) != length(avgDist)) {
      dfWish <- dfWish - 1
    }
    
    dfMatrix[j,k] <- dfWish
    
  }
  
}

#w/out 57
dfMatrix <- matrix(c(0,0,0,0,0,0,58,67,76,86,98,118,10,16,26,41,62,114), nrow = 7, ncol = 3)

#incl 57
dfMatrix <- matrix(c(0,0,0,0,0,0,0,58,67,76,86,98,113,118,10,16,26,41,62,101,114), nrow = 7, ncol = 3)


#getting the trees
howellsMC2trees1 <- read.tree(file = "A:\\Google Drive\\outputGD\\howellsLogUnifMCMC_1.trees")
howellsMC2trees1 <- howellsMC2trees1[2001:10000] #discard burn-in

howellsMC2trees2 <- read.tree(file = "A:\\Google Drive\\outputGD\\howellsLogUnifMCMC_2.trees")
howellsMC2trees2 <- howellsMC2trees2[2001:10000] #discard burn-in

trees <- c(howellsMC2trees1, howellsMC2trees2)


## actually preparing the files ##

traitNumIter <- c(2,4,8,16,32,57,64)

degreeMisspecification <- c("perfect", "interregional", "chimpHuman")

nrun <- 2

# i iterates over replicates
for(i in 1:100){

# k iterates over degree of misspecification
for (k in 1:length(degreeMisspecification)) {


# j iterates over trait numbers
for (j in 1:length(traitNumIter)) {

    
print(c(j, k, i))

ntraits <- traitNumIter[j]
if (ntraits < 58){
  traitsToUse <- sample(linMeasNames, size = ntraits, replace = F)
} else {
  traitsToUse <- linMeasNames
}
cov <- covFULL[traitsToUse, traitsToUse]
rateMatrixSpecification <- degreeMisspecification[k]
replicateNum <- i

# generate more traits if needed
if (traitNumIter[j] > 57){
    cov <- extendedOnionExpansion(cov = cov, d = ntraits, eta = 1)
}

covOrig <- cov

# wiggle traits if needed

if (rateMatrixSpecification != "perfect"){
  dfWish <- dfMatrix[j,k]
  cov <- rwish(dfWish, cov)/dfWish
}

fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_", "replicate_", replicateNum)

# write the covariance matrices
ratesPath <- paste0("A:\\Google Drive\\mvBM_sims\\rates\\", fileName, "_rates.Phy")
ratesPhy(round(covOrig, 4), ratesPath)

ratesPathTSV <- paste0("A:\\Google Drive\\mvBM_sims\\rates\\", fileName, "_rates.tsv")
ratesTSV(round(covOrig, 4), ratesPathTSV)

# simulate the traits and write to file
traitsPath <- paste0("A:\\Google Drive\\mvBM_sims\\data\\", fileName, "_traits.nex")

#use howells trees
tree <- trees[[sample(size = 1, x = 1:length(trees))]]
tree <- root(tree, outgroup = "BUSHMAN", resolve.root = T)

#use pbtrees
# tree <- pbtree(n=7)

traits <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = dim(cov)[1], sigma=round(cov, 4), mu= rep(0, times = dim(cov)[1])))
meansNexus(traits, traitsPath)
treePath <- paste0("A:\\Google Drive\\mvBM_sims\\trees\\", fileName, "_tree.nex")
write.nexus(tree, file = treePath, translate = T)
treePathF <- paste0("A:\\Google Drive\\mvBM_sims\\trees\\translateEqualsFalse_", fileName, "_tree.nex")
write.nexus(tree, file = treePathF, translate = F)

# write the master scripts
for (l in 1:nrun) {
scriptPath <- paste0("A:\\Google Drive\\mvBM_sims\\scripts\\", fileName, "_script_run_", l, ".Rev")
sink(scriptPath, append = F)

  cat("## Analysis of Simulated Data:", ntraits, "traits,", rateMatrixSpecification, "rate matrix specification, replication number", replicateNum, "##\n\n")
  cat("replicateNum <-", replicateNum, "\n")
  cat("runNum <-", l, "\n")
  cat(paste0("rateMatrixSpecification <- \"", rateMatrixSpecification, "\"\n\n"))
  cat("mi = 0\nseed(runNum)\n\n")
  
  cat("# read in the data\n")
  cat(paste0("means <- readContinuousCharacterData(\"data/", fileName, "_traits.nex\")\n\n"))

  cat("# read in the rate matrix\n")
  cat(paste0("rate = readDistanceMatrix(\"rates/", fileName, "_rates.Phy\")\n"))
  cat("rateMatrix = rate.symmetricMatrix()\n\n")
  
  cat("# get some useful variables from the data\nnumTips = means.ntaxa()\nnumBranches = 2 * numTips - 3\nntraits <- means.nchar()\ntaxa = means.taxa()\nnumNodes = numTips * 2 - 1\n\n")
  
  cat("# specify tree topology prior and moves\nsource(\"simulationsTreeTopology.Rev\")\n\n")

  cat("# specify branch length priors and moves\nsource(\"simulationsBranchLengths.Rev\")\n\n")
  
  cat("# assemble the tree\nphylogeny := treeAssembly(topology, br_lens)\n\n")
  
  cat("# put it all together\nmvBNdata ~ dnPhyloMultivariateBrownianREML(phylogeny, rep(1.0, numBranches), rateMatrix, nSites=ntraits)\nmvBNdata.clamp(means)\n\n")
  
  cat("# specify monitors\nsource(\"simulationsMonitors.Rev\")\n\n")
  
  cat("# define model\nmymodel = model(phylogeny)\n\n")
  
  cat("# run MCMC\nsource(\"simulationsMCMC.Rev\")\n\n")
  
  cat("# analyze mcmc output\nsource(\"simulationsAnalysis.Rev\")\n\n")
  
  cat("q()")
  
sink()

# add master script name to list

scriptPath <- paste0("scripts/", fileName, "_script_run_", l, ".Rev")
sink("A:\\Google Drive\\mvBM_sims\\job.txt", append = T)
cat(scriptPath, "\n")
sink()

}

}

}
  
}
