###### Getting Chimp Var-Cov Matrix ######
getwd()
setwd("/Users/nikolai/data")

meansNexus <- function (characterMatrix, outputFilePath, traitsOnColumns = T){
  if (traitsOnColumns) {ntax = dim(characterMatrix)[1]; nchar = dim(characterMatrix)[2]}
  else {ntax = dim(characterMatrix)[2]; nchar = dim(characterMatrix)[1]}
  sink(outputFilePath, append=F)
  cat("#NEXUS\n\nBegin data;\n\tDimensions ntax=", ntax, "nchar=", nchar, ";\n\tFormat datatype=Continuous missing=?;\n\tMatrix\n\n")
  for(i in 1:length(characterMatrix[,1])) {cat(rownames(characterMatrix)[i], characterMatrix[i,], "\n")}
  cat(";\nEnd;")
  sink()
}

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

############# DATA PREPARATION #############

#non-log analyses!

d.orig <- read.csv("chimp.csv")
d.orig <- read.csv("PanSorted.csv")
d <- d.orig
women <- d[d$Sex == "F",]
d <- women
sum(d == 0)

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

#### write to file in a way compatible with revbayes ####

ratesTSV(covChimp, "nonlogChimpRate_F.tsv")
meansNexus(traitsCHIMP, "nonlogChimpTraits_F.nex")

#################################################
###### Getting Modern Human Var-Cov Matrix ######
#################################################


d.orig <- read.csv("Howell.csv")
d <- d.orig
women <- d[d$Sex == "F",]
linMeasNames <- c("GOL", "NOL", "BNL", "BBH", "XCB", "XFB", "STB", "ZYB", "AUB", "WCB", "ASB", "BPL", "NPH", "NLH", 
                  "OBH", "OBB", "JUB", "NLB", "MAB", "MDH", "MDB", "ZMB", "SSS", "FMB", "NAS", "EKB", "DKB", "NDS", 
                  "WNB", "SIS", "IML", "XML", "MLS", "WMH", "FOL", "FRC", "FRS", "FRF", "PAC", "PAS", "PAF", "OCC", 
                  "OCS", "OCF", "VRR", "NAR", "SSR", "PRR", "DKR", "ZOR", "FMR", "EKR", "ZMR", "AVR", "DKS", "SOS", "GLS")
women <- cbind(women[,c(1,2,3,4)], women[,linMeasNames])
women
women[women == 0] <- NA
women <- women[complete.cases(women),]
women
pops <- sapply(1:length(unique(women$Population)), function (x) as.matrix((t(women[women$Population == unique(women$Population)[x],-(1:4)]))))

#get traits
traits <- list()
for(i in 1:length(linMeasNames)){ 
  print(i)
  trait <- list()
  trait[[1]] <- pops[[1]][i,]
  for (j in 2:26){
    trait[[j]] <- pops[[j]][i,]
  } 
  traits[[i]] <- trait
}

#make cov matrix
cov <- matrix(nrow = 57, ncol = 57)
for (i in 1:57){
  print(i)
  for(j in 1:57){
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
traitsHOWELL <- matrix(nrow=npop, ncol=57)
for(i in 1:npop){
  for(j in 1:57){
    traitsHOWELL[i,j] <- mean(pops[[i]][j,])
  }
}
rownames(traitsHOWELL) <- unique(women$Population)
colnames(traitsHOWELL) <- colnames(women)[-(1:4)]

rownames(cov) <- colnames(cov) <- colnames(women)[-(1:4)]
covHuman <- cov



# populations <- c("N JAPAN", "S JAPAN", "S MAORI", "N MAORI", "BUSHMAN", "AUSTRALI", "TASMANIA", "ARIKARA", "SANTA CR", "PERU", "BERG", "ZALAVAR", "AINU", "NORSE", "EASTER I", "MOKAPU")
#populations <- c("N JAPAN", "S JAPAN", "S MAORI", "N MAORI", "BUSHMAN", "AUSTRALI", "TASMANIA", "ARIKARA", "SANTA CR", "PERU")
#taken from roseman and weaver 2004, only BPC -> BPL (due to a typo? it's the basion-prosthion length )
linMeasNames <- c("GOL", "NOL", "BNL", "BBH", "XCB", "XFB", "STB", "ZYB", "AUB", "WCB", "ASB", "BPL", "NPH", "NLH", "OBH", "OBB", "JUB", "NLB", "MAB", "MDH", "MDB", "ZMB", "SSS", "FMB", "NAS", "EKB", "DKB", "NDS", "WNB", "SIS", "IML", "XML", "MLS", "WMH", "FOL", "FRC", "FRS", "FRF", "PAC", "PAS", "PAF", "OCC", "OCS", "OCF", "VRR", "NAR", "SSR", "PRR", "DKR", "ZOR", "FMR", "EKR", "ZMR", "AVR", "DKS", "SOS", "GLS")
linMeasTraitsHuman <- traitsHOWELL[,linMeasNames]
rownames(linMeasTraitsHuman) <- gsub(" ", "", rownames(linMeasTraitsHuman))

linMeasCovHuman <- covHuman[linMeasNames, linMeasNames]

#### write to file in a way compatible with revbayes ####

ratesTSV(linMeasCovHuman, "nonlogHowellsRate_F.tsv")
meansNexus(linMeasTraitsHuman, "nonlogHowellsTraits_F.nex")


##########################################
###### Get Neandertals in there Too ######
##########################################

# d.orig <- read.csv("neandertal.csv")
# d <- d.orig
# d <- as.matrix(d[1:20, 3:29])
# covNean <- cov(d, use = "pairwise.complete.obs")
# colmeans <- sapply(1:length(d[1,]), function (x) mean(d[,x], na.rm = T))
# d.imp <- matrix(nrow = 20, ncol = 27, dimnames = list(rownames(d), colnames(d)))
# for (i in 1:20){ 
#   for (j in 1:27) { 
#     print(j)
#     if (is.na(d[i,j])) {d.imp[i,j] <- colmeans[j]} 
#     else {d.imp[i,j] <- d[i,j]}  
#   }
# }
# covNean <- cov(d.imp)
# 
# NeanHumanTraits <- intersect(rownames(covNean), rownames(covHuman))
# sharedNeanCov <- covNean[NeanHumanTraits, NeanHumanTraits]
# sharedHumanCov <- covHuman[NeanHumanTraits, NeanHumanTraits]
# 
# ## write neandertal means in with humans and the appropriate rate matrix
# 
# humNeanTraits <- linMeasTraitsHuman[,NeanHumanTraits]
# NEANDERTAL <- colmeans
# humNeanTraits <- rbind(humNeanTraits, NEANDERTAL)
# all(colnames(d) == colnames(humNeanTraits))
# rownames(humNeanTraits) <- gsub(" ", "", rownames(humNeanTraits))
# 
# sharedHumanCovr <- round(sharedHumanCov, digits = 4)
# isSymmetric(sharedHumanCovr)
# 
# #adding in the additional neandertal traits
# 
# neanTraits <- "
# + XX X 111.1
# + X X 93.4
# + X X X X X 98.9
# + X X X X X 112.5
# + X X X X X 119.4
# + XX X 116.9
# + X X X X X X X X 19.9
# + X X X X X 85.3
# + X X X X X 201.8
# + X X X X X X 197.0
# + X X 126.6
# + X X X X X X 148.7
# + X X X X X X X X 121.6
# + X X X X X X X X 134.0
# + X X X X 121.1
# + X X X X X 62.1
# + X X X X X X 36.4
# + X X X X X 44.1
# + X X X X X X 32.5
# + X X X X X X X 74.3
# + X X X X X X X 22.1
# + X X XXXX 111.5
# + X X XXXX 36.4
# + X X X X X X X X XXX 112.7
# + X X X X X X X 22.5
# + X X X X X 109.5
# + X X X X X X 28.1
# + X X X X X X X X 111.0
# + X X X X X X X X 57.3
# + X X X X X X X 108.5
# + X X X X X X X 18.6
# + X X X X X X X 54.7
# + X X 33.2
# + X X 40.7
# + X X X X X 118.4
# + X X X X X X X 107.3
# + X X X X 77.2"
# neanTraits <- gsub("X", "", neanTraits); neanTraits <- gsub("\n", "", neanTraits)
# neanTraits <- strsplit(neanTraits, split = " ")
# neanTraits <- neanTraits[[1]][which(neanTraits[[1]] != "")][-1]
# neanTraits[length(neanTraits)] <- paste0(neanTraits[length(neanTraits)], "+")
# neanTraits <- gsub('.{1}$', '', neanTraits)
# as.numeric(neanTraits)
# 
# 
# ntnames <- "BNL Basion-nasion length X X X X
# OCC Lambda-opisthion chord X X X X
# AVR Molar alveolus radius X X X
# SSR Subspinale radius X X X X X
# PRR Prosthion radius X X X X X
# BPL Basion-prosthion length X X X X
# FRS Nasion-bregma subtense X X X X X X X
# NPH Nasion-prosthion height X X X X X
# GOL Glabella-occipital length X X X X X X
# NOL Nasio-occipital length X X X X X X
# BBH Basion-bregma height X X X X
# XCB Maximum cranial breadth X X X X X X X
# XFB Maximum frontal breadth X X X X X X X X
# AUB Biauricular breadth X X X X X X
# ASB Biasterionic breadth X X X X X X
# NLH Nasal height X X X X X
# OBH Orbital height X X X X X X
# OBB Orbit breadth X X X X X X
# NLB Nasal breadth X X X X X
# MAB Palate breadth X X X X X X
# MDH Mastoid height X X X X X
# ZMB Bimaxillary breadth X X X X X X X
# SSS Subspinale radius X X X X X
# FMB Bifrontal breadth X X X X X X X X
# NAS Nasio-frontal subtense X X X X X X X X
# EKB Biorbital breadth X X X X X X X
# DKB Interorbital breadth X X X X X X
# FRC Nasion-bregma chord X X X X X X X
# FRF Nasion-subtense fraction X X X X X X X
# PAC Bregma-lambda chord X X X X X X
# PAS Bregma-lambda subtense X X X X X X
# PAF Bregma-subtense fraction X X X X X X
# OCS Lambda-opisthion subtense X X X X
# OCF Lambda-subtense fraction X X X X
# VRR Vertex radius X X X X X X
# NAR Nasion radius X X X X X X X
# ZMR"
# ntnames <- strsplit(gsub("X\n", "", ntnames), " ")
# ntnames <- ntnames[[1]][sapply(1:length(ntnames[[1]]), function (x) nchar(ntnames[[1]][x])) == 3]
# 
# 
# NeanHumanTraits <- intersect(ntnames, rownames(covHuman))
# sharedHumanCov <- covHuman[ntnames, ntnames]
# 
# humNeanTraits <- linMeasTraitsHuman[,ntnames]
# NEANDERTAL <- as.numeric(neanTraits)
# humNeanTraits <- rbind(humNeanTraits, NEANDERTAL)
# rownames(humNeanTraits) <- gsub(" ", "", rownames(humNeanTraits))
# humNeanTraits <- round(humNeanTraits, 2)
# 
# sharedHumanCovr <- round(sharedHumanCov, digits = 4)
# isSymmetric(sharedHumanCovr)


## excluding female and indeterminate neandertals from the sample

d.orig <- read.csv("Neandertal_meas37.csv")
d <- d.orig
m <- d[d$SEX == "F",]
m <- m[,-(1:2)]
NEANDERTAL <- sapply(1:length(m[1,]), function (x) mean(m[,x], na.rm = T))
traitnames <- colnames(m)
sharedHumanCov <- covHuman[traitnames, traitnames]

humNeanTraits <- linMeasTraitsHuman[,traitnames]
humNeanTraits <- rbind(humNeanTraits, NEANDERTAL)
rownames(humNeanTraits) <- gsub(" ", "", rownames(humNeanTraits))
humNeanTraits <- round(humNeanTraits, 4)

sharedHumanCovr <- round(sharedHumanCov, digits = 4)
isSymmetric(sharedHumanCovr)

#### write to file in a way compatible with revbayes ####

ratesTSV(sharedHumanCovr, "nonlogMaleNeandertalModHumRate_F.tsv")
meansNexus(humNeanTraits, "nonlogMaleNeandertalModHumTraits_F.nex")










############### ############### ############### 
############### log analyses! ################# 
############### ############### ###############

d.orig <- read.csv("chimp.csv")
d.orig <- read.csv("PanSorted.csv")
d <- d.orig
men <- d[d$Sex == "F",]
d <- men
sum(d == 0)
d[,-c(1,2)] <- log(d[,-c(1,2)])
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

npop <- length(sp)
#compute matrix of Species means
traitsCHIMPlog <- matrix(nrow=npop, ncol=ntraits)
for(i in 1:npop){
  for(j in 1:ntraits){
    traitsCHIMPlog[i,j] <- mean(sp[[i]][j,])
  }
}
rownames(traitsCHIMPlog) <- unique(d$Taxon)
colnames(traitsCHIMPlog) <- colnames(d)[-(1:2)]
rownames(cov) <- colnames(cov) <- colnames(d)[-(1:2)]
covChimpLog <- cov

#### write to file in a way compatible with revbayes ####

ratesTSV(covChimpLog, "logChimpRate_F.tsv")
meansNexus(traitsCHIMPlog, "logChimpTraits_F.nex")

#################################################
###### Getting Modern Human Var-Cov Matrix ######
#################################################


d.orig <- read.csv("Howell.csv")
d <- d.orig
men <- d[d$Sex == "F",]
linMeasNames <- c("GOL", "NOL", "BNL", "BBH", "XCB", "XFB", "STB", "ZYB", "AUB", "WCB", "ASB", "BPL", "NPH", "NLH", 
                  "OBH", "OBB", "JUB", "NLB", "MAB", "MDH", "MDB", "ZMB", "SSS", "FMB", "NAS", "EKB", "DKB", "NDS", 
                  "WNB", "SIS", "IML", "XML", "MLS", "WMH", "FOL", "FRC", "FRS", "FRF", "PAC", "PAS", "PAF", "OCC", 
                  "OCS", "OCF", "VRR", "NAR", "SSR", "PRR", "DKR", "ZOR", "FMR", "EKR", "ZMR", "AVR", "DKS", "SOS", "GLS")
men <- cbind(men[,c(1,2,3,4)], men[,linMeasNames])
men
men[men == 0] <- NA
men <- men[complete.cases(men),]
men
men[,-c(1:4)] <- log(men[,-c(1:4)])

pops <- sapply(1:length(unique(men$Population)), function (x) as.matrix((t(men[men$Population == unique(men$Population)[x],-(1:4)]))))


#get traits
traits <- list()
for(i in 1:57){ 
  print(i)
  trait <- list()
  trait[[1]] <- pops[[1]][i,]
  for (j in 2:26){
    trait[[j]] <- pops[[j]][i,]
  } 
  traits[[i]] <- trait
}

#make cov matrix
cov <- matrix(nrow = 57, ncol = 57)
for (i in 1:57){
  print(i)
  for(j in 1:57){
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
traitsHOWELLlog <- matrix(nrow=npop, ncol=57)
for(i in 1:npop){
  for(j in 1:57){
    traitsHOWELLlog[i,j] <- mean(pops[[i]][j,])
  }
}
rownames(traitsHOWELLlog) <- unique(men$Population)
colnames(traitsHOWELLlog) <- colnames(men)[-(1:4)]

rownames(cov) <- colnames(cov) <- colnames(men)[-(1:4)]
covHumanLog <- cov



# populations <- c("N JAPAN", "S JAPAN", "S MAORI", "N MAORI", "BUSHMAN", "AUSTRALI", "TASMANIA", "ARIKARA", "SANTA CR", "PERU", "BERG", "ZALAVAR", "AINU", "NORSE", "EASTER I", "MOKAPU")
#populations <- c("N JAPAN", "S JAPAN", "S MAORI", "N MAORI", "BUSHMAN", "AUSTRALI", "TASMANIA", "ARIKARA", "SANTA CR", "PERU")
#taken from roseman and weaver 2004, only BPC -> BPL (due to a typo? it's the basion-prosthion length )
linMeasNames <- c("GOL", "NOL", "BNL", "BBH", "XCB", "XFB", "STB", "ZYB", "AUB", "WCB", "ASB", "BPL", "NPH", "NLH", "OBH", "OBB", "JUB", "NLB", "MAB", "MDH", "MDB", "ZMB", "SSS", "FMB", "NAS", "EKB", "DKB", "NDS", "WNB", "SIS", "IML", "XML", "MLS", "WMH", "FOL", "FRC", "FRS", "FRF", "PAC", "PAS", "PAF", "OCC", "OCS", "OCF", "VRR", "NAR", "SSR", "PRR", "DKR", "ZOR", "FMR", "EKR", "ZMR", "AVR", "DKS", "SOS", "GLS")
linMeasTraitsHumanLog <- traitsHOWELLlog[,linMeasNames]
rownames(linMeasTraitsHumanLog) <- gsub(" ", "", rownames(linMeasTraitsHumanLog))
linMeasCovHumanLog <- covHumanLog[linMeasNames, linMeasNames]

#### write to file in a way compatible with revbayes ####

ratesTSV(linMeasCovHumanLog, "logHowellsRate_F.tsv")
meansNexus(linMeasTraitsHumanLog, "logHowellsTraits_F.nex")



##########################################
###### Get Neandertals in there Too ######
##########################################
# 
# d.orig <- read.csv("neandertal.csv")
# d <- d.orig
# d <- as.matrix(d[1:20, 3:29])
# covNean <- cov(d, use = "pairwise.complete.obs")
# colmeans <- sapply(1:length(d[1,]), function (x) mean(d[,x], na.rm = T))
# d.imp <- matrix(nrow = 20, ncol = 27, dimnames = list(rownames(d), colnames(d)))
# for (i in 1:20){ 
#   for (j in 1:27) { 
#     print(j)
#     if (is.na(d[i,j])) {d.imp[i,j] <- colmeans[j]} 
#     else {d.imp[i,j] <- d[i,j]}  
#   }
# }
# covNean <- cov(d.imp)
# 
# NeanHumanTraits <- intersect(rownames(covNean), rownames(covHuman))
# sharedNeanCov <- covNean[NeanHumanTraits, NeanHumanTraits]
# sharedHumanCov <- covHuman[NeanHumanTraits, NeanHumanTraits]

## write neandertal means in with humans and the appropriate rate matrix

# humNeanTraits <- linMeasTraitsHuman[,NeanHumanTraits]
# NEANDERTAL <- colmeans
# humNeanTraits <- rbind(humNeanTraits, NEANDERTAL)
# all(colnames(d) == colnames(humNeanTraits))
# rownames(humNeanTraits) <- gsub(" ", "", rownames(humNeanTraits))
# 
# sharedHumanCovr <- round(sharedHumanCov, digits = 4)
# isSymmetric(sharedHumanCovr)
# 
# #adding in the additional neandertal traits
# 
# neanTraits <- "
# + XX X 111.1
# + X X 93.4
# + X X X X X 98.9
# + X X X X X 112.5
# + X X X X X 119.4
# + XX X 116.9
# + X X X X X X X X 19.9
# + X X X X X 85.3
# + X X X X X 201.8
# + X X X X X X 197.0
# + X X 126.6
# + X X X X X X 148.7
# + X X X X X X X X 121.6
# + X X X X X X X X 134.0
# + X X X X 121.1
# + X X X X X 62.1
# + X X X X X X 36.4
# + X X X X X 44.1
# + X X X X X X 32.5
# + X X X X X X X 74.3
# + X X X X X X X 22.1
# + X X XXXX 111.5
# + X X XXXX 36.4
# + X X X X X X X X XXX 112.7
# + X X X X X X X 22.5
# + X X X X X 109.5
# + X X X X X X 28.1
# + X X X X X X X X 111.0
# + X X X X X X X X 57.3
# + X X X X X X X 108.5
# + X X X X X X X 18.6
# + X X X X X X X 54.7
# + X X 33.2
# + X X 40.7
# + X X X X X 118.4
# + X X X X X X X 107.3
# + X X X X 77.2"
# neanTraits <- gsub("X", "", neanTraits); neanTraits <- gsub("\n", "", neanTraits)
# neanTraits <- strsplit(neanTraits, split = " ")
# neanTraits <- neanTraits[[1]][which(neanTraits[[1]] != "")][-1]
# neanTraits[length(neanTraits)] <- paste0(neanTraits[length(neanTraits)], "+")
# neanTraits <- gsub('.{1}$', '', neanTraits)
# as.numeric(neanTraits)
# 
# 
# ntnames <- "BNL Basion-nasion length X X X X
# OCC Lambda-opisthion chord X X X X
# AVR Molar alveolus radius X X X
# SSR Subspinale radius X X X X X
# PRR Prosthion radius X X X X X
# BPL Basion-prosthion length X X X X
# FRS Nasion-bregma subtense X X X X X X X
# NPH Nasion-prosthion height X X X X X
# GOL Glabella-occipital length X X X X X X
# NOL Nasio-occipital length X X X X X X
# BBH Basion-bregma height X X X X
# XCB Maximum cranial breadth X X X X X X X
# XFB Maximum frontal breadth X X X X X X X X
# AUB Biauricular breadth X X X X X X
# ASB Biasterionic breadth X X X X X X
# NLH Nasal height X X X X X
# OBH Orbital height X X X X X X
# OBB Orbit breadth X X X X X X
# NLB Nasal breadth X X X X X
# MAB Palate breadth X X X X X X
# MDH Mastoid height X X X X X
# ZMB Bimaxillary breadth X X X X X X X
# SSS Subspinale radius X X X X X
# FMB Bifrontal breadth X X X X X X X X
# NAS Nasio-frontal subtense X X X X X X X X
# EKB Biorbital breadth X X X X X X X
# DKB Interorbital breadth X X X X X X
# FRC Nasion-bregma chord X X X X X X X
# FRF Nasion-subtense fraction X X X X X X X
# PAC Bregma-lambda chord X X X X X X
# PAS Bregma-lambda subtense X X X X X X
# PAF Bregma-subtense fraction X X X X X X
# OCS Lambda-opisthion subtense X X X X
# OCF Lambda-subtense fraction X X X X
# VRR Vertex radius X X X X X X
# NAR Nasion radius X X X X X X X
# ZMR"
# ntnames <- strsplit(gsub("X\n", "", ntnames), " ")
# ntnames <- ntnames[[1]][sapply(1:length(ntnames[[1]]), function (x) nchar(ntnames[[1]][x])) == 3]
# 
# 
# NeanHumanTraits <- intersect(ntnames, rownames(covHuman))
# sharedHumanCov <- covHuman[ntnames, ntnames]
# 
# humNeanTraits <- linMeasTraitsHuman[,ntnames]
# NEANDERTAL <- as.numeric(neanTraits)
# humNeanTraits <- rbind(humNeanTraits, NEANDERTAL)
# rownames(humNeanTraits) <- gsub(" ", "", rownames(humNeanTraits))
# humNeanTraits <- round(humNeanTraits, 2)
# sink("C:\\Users\\Nikolai\\Google Drive\\data\\humNeanMeans.nex", append=F)
# cat("#NEXUS\n\nBegin data;\n\tDimensions ntax=31 nchar=37;\n\tFormat datatype=Continuous missing=?;\n\tMatrix\n\n")
# for(i in 1:length(humNeanTraits[,1])) {cat(rownames(humNeanTraits)[i], humNeanTraits[i,], "\n")}
# cat(";\nEnd;")
# sink()
# 
# sharedHumanCovr <- round(sharedHumanCov, digits = 4)
# isSymmetric(sharedHumanCovr)
# sink("C:\\Users\\Nikolai\\Google Drive\\data\\humNeanRateMatrix.phy", append=F)
# cat("37\n")
# for(i in 1:length(sharedHumanCovr[,1])) {
#   cat(rownames(sharedHumanCovr)[i], paste0(as.vector(sharedHumanCovr[i, (1:(length(sharedHumanCovr[i,])-1))])))
#   cat(" ")
#   cat(sharedHumanCovr[i, (length(sharedHumanCovr[i,]))], "\n", sep = "")
# }
# sink()

## excluding female and indeterminate neandertals from the sample

d.orig <- read.csv("Neandertal_meas37.csv")
d <- d.orig
m <- d[d$SEX == "F",]
m <- m[,-(1:2)]
m <- log(m)
NEANDERTAL <- sapply(1:length(m[1,]), function (x) mean(m[,x], na.rm = T))
traitnames <- colnames(m)
sharedHumanCovLog <- covHumanLog[traitnames, traitnames]

humNeanTraitsLog <- linMeasTraitsHumanLog[,traitnames]
humNeanTraitsLog <- rbind(humNeanTraitsLog, NEANDERTAL)
rownames(humNeanTraitsLog) <- gsub(" ", "", rownames(humNeanTraitsLog))
humNeanTraitsLog <- round(humNeanTraitsLog, 3)

isSymmetric(sharedHumanCovLog)

#### write to file in a way compatible with revbayes ####

ratesTSV(sharedHumanCovLog, "logMaleNeandertalModHumRate_F.tsv")
meansNexus(humNeanTraitsLog, "logMaleNeandertalModHumTraits_F.nex")


# let's try doing both men and women in the same analysis, but as separate tips
# let's try doing both men and women in the same analysis, but as separate tips
# let's try doing both men and women in the same analysis, but as separate tips
# let's try doing both men and women in the same analysis, but as separate tips
# let's try doing both men and women in the same analysis, but as separate tips
# let's try doing both men and women in the same analysis, but as separate tips
# let's try doing both men and women in the same analysis, but as separate tips
# let's try doing both men and women in the same analysis, but as separate tips
# let's try doing both men and women in the same analysis, but as separate tips
# let's try doing both men and women in the same analysis, but as separate tips
# let's try doing both men and women in the same analysis, but as separate tips
# let's try doing both men and women in the same analysis, but as separate tips


d.orig <- read.csv("Howell.csv")
d <- d.orig
all <- d
linMeasNames <- c("GOL", "NOL", "BNL", "BBH", "XCB", "XFB", "STB", "ZYB", "AUB", "WCB", "ASB", "BPL", "NPH", "NLH", 
                  "OBH", "OBB", "JUB", "NLB", "MAB", "MDH", "MDB", "ZMB", "SSS", "FMB", "NAS", "EKB", "DKB", "NDS", 
                  "WNB", "SIS", "IML", "XML", "MLS", "WMH", "FOL", "FRC", "FRS", "FRF", "PAC", "PAS", "PAF", "OCC", 
                  "OCS", "OCF", "VRR", "NAR", "SSR", "PRR", "DKR", "ZOR", "FMR", "EKR", "ZMR", "AVR", "DKS", "SOS", "GLS")
all <- cbind(all[,c(1,2,3,4)], all[,linMeasNames])
all$Population <- paste0(all$Population, ".", all$Sex)
all[all == 0] <- NA
all
all <- all[complete.cases(all),]
all
pops <- sapply(1:length(unique(all$Population)), function (x) as.matrix((t(all[all$Population == unique(all$Population)[x],-(1:4)]))))

#get traits
traits <- list()
for(i in 1:length(linMeasNames)){ 
  print(i)
  trait <- list()
  trait[[1]] <- pops[[1]][i,]
  for (j in 2:56){
    trait[[j]] <- pops[[j]][i,]
  } 
  traits[[i]] <- trait
}

#make cov matrix
cov <- matrix(nrow = 57, ncol = 57)
for (i in 1:57){
  print(i)
  for(j in 1:57){
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
traitsHOWELL <- matrix(nrow=npop, ncol=57)
for(i in 1:npop){
  for(j in 1:57){
    traitsHOWELL[i,j] <- mean(pops[[i]][j,])
  }
}
rownames(traitsHOWELL) <- unique(all$Population)
colnames(traitsHOWELL) <- colnames(all)[-(1:4)]

rownames(cov) <- colnames(cov) <- colnames(all)[-(1:4)]
covHuman <- cov



# populations <- c("N JAPAN", "S JAPAN", "S MAORI", "N MAORI", "BUSHMAN", "AUSTRALI", "TASMANIA", "ARIKARA", "SANTA CR", "PERU", "BERG", "ZALAVAR", "AINU", "NORSE", "EASTER I", "MOKAPU")
#populations <- c("N JAPAN", "S JAPAN", "S MAORI", "N MAORI", "BUSHMAN", "AUSTRALI", "TASMANIA", "ARIKARA", "SANTA CR", "PERU")
#taken from roseman and weaver 2004, only BPC -> BPL (due to a typo? it's the basion-prosthion length )
linMeasNames <- c("GOL", "NOL", "BNL", "BBH", "XCB", "XFB", "STB", "ZYB", "AUB", "WCB", "ASB", "BPL", "NPH", "NLH", "OBH", "OBB", "JUB", "NLB", "MAB", "MDH", "MDB", "ZMB", "SSS", "FMB", "NAS", "EKB", "DKB", "NDS", "WNB", "SIS", "IML", "XML", "MLS", "WMH", "FOL", "FRC", "FRS", "FRF", "PAC", "PAS", "PAF", "OCC", "OCS", "OCF", "VRR", "NAR", "SSR", "PRR", "DKR", "ZOR", "FMR", "EKR", "ZMR", "AVR", "DKS", "SOS", "GLS")
linMeasTraitsHuman <- traitsHOWELL[,linMeasNames]
rownames(linMeasTraitsHuman) <- gsub(" ", "", rownames(linMeasTraitsHuman))

linMeasCovHuman <- covHuman[linMeasNames, linMeasNames]

#### write to file in a way compatible with revbayes ####

ratesTSV(linMeasCovHuman, "nonlogHowellsRate_MF.tsv")
meansNexus(linMeasTraitsHuman, "nonlogHowellsTraits_MF.nex")



#log transform

d.orig <- read.csv("Howell.csv")
d <- d.orig
all <- d
linMeasNames <- c("GOL", "NOL", "BNL", "BBH", "XCB", "XFB", "STB", "ZYB", "AUB", "WCB", "ASB", "BPL", "NPH", "NLH", 
                  "OBH", "OBB", "JUB", "NLB", "MAB", "MDH", "MDB", "ZMB", "SSS", "FMB", "NAS", "EKB", "DKB", "NDS", 
                  "WNB", "SIS", "IML", "XML", "MLS", "WMH", "FOL", "FRC", "FRS", "FRF", "PAC", "PAS", "PAF", "OCC", 
                  "OCS", "OCF", "VRR", "NAR", "SSR", "PRR", "DKR", "ZOR", "FMR", "EKR", "ZMR", "AVR", "DKS", "SOS", "GLS")
all <- cbind(all[,c(1,2,3,4)], all[,linMeasNames])
all
all[all == 0] <- NA
all <- all[complete.cases(all),]
all
all[,-c(1:4)] <- log(all[,-c(1:4)])
all$Population <- paste0(all$Population, ".", all$Sex)


pops <- sapply(1:length(unique(all$Population)), function (x) as.matrix((t(all[all$Population == unique(all$Population)[x],-(1:4)]))))


#get traits
traits <- list()
for(i in 1:57){ 
  print(i)
  trait <- list()
  trait[[1]] <- pops[[1]][i,]
  for (j in 2:56){
    trait[[j]] <- pops[[j]][i,]
  } 
  traits[[i]] <- trait
}

#make cov matrix
cov <- matrix(nrow = 57, ncol = 57)
for (i in 1:57){
  print(i)
  for(j in 1:57){
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
traitsHOWELLlog <- matrix(nrow=npop, ncol=57)
for(i in 1:npop){
  for(j in 1:57){
    traitsHOWELLlog[i,j] <- mean(pops[[i]][j,])
  }
}
rownames(traitsHOWELLlog) <- unique(all$Population)
colnames(traitsHOWELLlog) <- colnames(all)[-(1:4)]

rownames(cov) <- colnames(cov) <- colnames(all)[-(1:4)]
covHumanLog <- cov



# populations <- c("N JAPAN", "S JAPAN", "S MAORI", "N MAORI", "BUSHMAN", "AUSTRALI", "TASMANIA", "ARIKARA", "SANTA CR", "PERU", "BERG", "ZALAVAR", "AINU", "NORSE", "EASTER I", "MOKAPU")
#populations <- c("N JAPAN", "S JAPAN", "S MAORI", "N MAORI", "BUSHMAN", "AUSTRALI", "TASMANIA", "ARIKARA", "SANTA CR", "PERU")
#taken from roseman and weaver 2004, only BPC -> BPL (due to a typo? it's the basion-prosthion length )
linMeasNames <- c("GOL", "NOL", "BNL", "BBH", "XCB", "XFB", "STB", "ZYB", "AUB", "WCB", "ASB", "BPL", "NPH", "NLH", "OBH", "OBB", "JUB", "NLB", "MAB", "MDH", "MDB", "ZMB", "SSS", "FMB", "NAS", "EKB", "DKB", "NDS", "WNB", "SIS", "IML", "XML", "MLS", "WMH", "FOL", "FRC", "FRS", "FRF", "PAC", "PAS", "PAF", "OCC", "OCS", "OCF", "VRR", "NAR", "SSR", "PRR", "DKR", "ZOR", "FMR", "EKR", "ZMR", "AVR", "DKS", "SOS", "GLS")
linMeasTraitsHumanLog <- traitsHOWELLlog[,linMeasNames]
rownames(linMeasTraitsHumanLog) <- gsub(" ", "", rownames(linMeasTraitsHumanLog))
linMeasCovHumanLog <- covHumanLog[linMeasNames, linMeasNames]

#### write to file in a way compatible with revbayes ####

ratesTSV(linMeasCovHumanLog, "logHowellsRate_MF.tsv")
meansNexus(linMeasTraitsHumanLog, "logHowellsTraits_MF.nex")


# let's try doing both men and women in the same analysis, but as the same tips
# let's try doing both men and women in the same analysis, but as the same tips
# let's try doing both men and women in the same analysis, but as the same tips
# let's try doing both men and women in the same analysis, but as the same tips
# let's try doing both men and women in the same analysis, but as the same tips
# let's try doing both men and women in the same analysis, but as the same tips
# let's try doing both men and women in the same analysis, but as the same tips
# let's try doing both men and women in the same analysis, but as the same tips
# let's try doing both men and women in the same analysis, but as the same tips
# let's try doing both men and women in the same analysis, but as the same tips
# let's try doing both men and women in the same analysis, but as the same tips
# let's try doing both men and women in the same analysis, but as the same tips
# let's try doing both men and women in the same analysis, but as the same tips
# let's try doing both men and women in the same analysis, but as the same tips
# let's try doing both men and women in the same analysis, but as the same tips
# let's try doing both men and women in the same analysis, but as the same tips

d.orig <- read.csv("Howell.csv")
d <- d.orig
all <- d
linMeasNames <- c("GOL", "NOL", "BNL", "BBH", "XCB", "XFB", "STB", "ZYB", "AUB", "WCB", "ASB", "BPL", "NPH", "NLH", 
                  "OBH", "OBB", "JUB", "NLB", "MAB", "MDH", "MDB", "ZMB", "SSS", "FMB", "NAS", "EKB", "DKB", "NDS", 
                  "WNB", "SIS", "IML", "XML", "MLS", "WMH", "FOL", "FRC", "FRS", "FRF", "PAC", "PAS", "PAF", "OCC", 
                  "OCS", "OCF", "VRR", "NAR", "SSR", "PRR", "DKR", "ZOR", "FMR", "EKR", "ZMR", "AVR", "DKS", "SOS", "GLS")
all <- cbind(all[,c(1,2,3,4)], all[,linMeasNames])
all[all == 0] <- NA
all
all <- all[complete.cases(all),]
all
pops <- sapply(1:length(unique(all$Population)), function (x) as.matrix((t(all[all$Population == unique(all$Population)[x],-(1:4)]))))

#get traits
traits <- list()
for(i in 1:length(linMeasNames)){ 
  print(i)
  trait <- list()
  trait[[1]] <- pops[[1]][i,]
  for (j in 2:30){
    trait[[j]] <- pops[[j]][i,]
  } 
  traits[[i]] <- trait
}

#make cov matrix
cov <- matrix(nrow = 57, ncol = 57)
for (i in 1:57){
  print(i)
  for(j in 1:57){
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
traitsHOWELL <- matrix(nrow=npop, ncol=57)
for(i in 1:npop){
  for(j in 1:57){
    traitsHOWELL[i,j] <- mean(pops[[i]][j,])
  }
}
rownames(traitsHOWELL) <- unique(all$Population)
colnames(traitsHOWELL) <- colnames(all)[-(1:4)]

rownames(cov) <- colnames(cov) <- colnames(all)[-(1:4)]
covHuman <- cov



# populations <- c("N JAPAN", "S JAPAN", "S MAORI", "N MAORI", "BUSHMAN", "AUSTRALI", "TASMANIA", "ARIKARA", "SANTA CR", "PERU", "BERG", "ZALAVAR", "AINU", "NORSE", "EASTER I", "MOKAPU")
#populations <- c("N JAPAN", "S JAPAN", "S MAORI", "N MAORI", "BUSHMAN", "AUSTRALI", "TASMANIA", "ARIKARA", "SANTA CR", "PERU")
#taken from roseman and weaver 2004, only BPC -> BPL (due to a typo? it's the basion-prosthion length )
linMeasNames <- c("GOL", "NOL", "BNL", "BBH", "XCB", "XFB", "STB", "ZYB", "AUB", "WCB", "ASB", "BPL", "NPH", "NLH", "OBH", "OBB", "JUB", "NLB", "MAB", "MDH", "MDB", "ZMB", "SSS", "FMB", "NAS", "EKB", "DKB", "NDS", "WNB", "SIS", "IML", "XML", "MLS", "WMH", "FOL", "FRC", "FRS", "FRF", "PAC", "PAS", "PAF", "OCC", "OCS", "OCF", "VRR", "NAR", "SSR", "PRR", "DKR", "ZOR", "FMR", "EKR", "ZMR", "AVR", "DKS", "SOS", "GLS")
linMeasTraitsHuman <- traitsHOWELL[,linMeasNames]
rownames(linMeasTraitsHuman) <- gsub(" ", "", rownames(linMeasTraitsHuman))

linMeasCovHuman <- covHuman[linMeasNames, linMeasNames]

#### write to file in a way compatible with revbayes ####

ratesTSV(linMeasCovHuman, "nonlogHowellsRate_MwF.tsv")
meansNexus(linMeasTraitsHuman, "nonlogHowellsTraits_MwF.nex")



#log transform

d.orig <- read.csv("Howell.csv")
d <- d.orig
all <- d
linMeasNames <- c("GOL", "NOL", "BNL", "BBH", "XCB", "XFB", "STB", "ZYB", "AUB", "WCB", "ASB", "BPL", "NPH", "NLH", 
                  "OBH", "OBB", "JUB", "NLB", "MAB", "MDH", "MDB", "ZMB", "SSS", "FMB", "NAS", "EKB", "DKB", "NDS", 
                  "WNB", "SIS", "IML", "XML", "MLS", "WMH", "FOL", "FRC", "FRS", "FRF", "PAC", "PAS", "PAF", "OCC", 
                  "OCS", "OCF", "VRR", "NAR", "SSR", "PRR", "DKR", "ZOR", "FMR", "EKR", "ZMR", "AVR", "DKS", "SOS", "GLS")
all <- cbind(all[,c(1,2,3,4)], all[,linMeasNames])
all
all[all == 0] <- NA
all <- all[complete.cases(all),]
all
all[,-c(1:4)] <- log(all[,-c(1:4)])


pops <- sapply(1:length(unique(all$Population)), function (x) as.matrix((t(all[all$Population == unique(all$Population)[x],-(1:4)]))))


#get traits
traits <- list()
for(i in 1:57){ 
  print(i)
  trait <- list()
  trait[[1]] <- pops[[1]][i,]
  for (j in 2:30){
    trait[[j]] <- pops[[j]][i,]
  } 
  traits[[i]] <- trait
}

#make cov matrix
cov <- matrix(nrow = 57, ncol = 57)
for (i in 1:57){
  print(i)
  for(j in 1:57){
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
traitsHOWELLlog <- matrix(nrow=npop, ncol=57)
for(i in 1:npop){
  for(j in 1:57){
    traitsHOWELLlog[i,j] <- mean(pops[[i]][j,])
  }
}
rownames(traitsHOWELLlog) <- unique(all$Population)
colnames(traitsHOWELLlog) <- colnames(all)[-(1:4)]

rownames(cov) <- colnames(cov) <- colnames(all)[-(1:4)]
covHumanLog <- cov



# populations <- c("N JAPAN", "S JAPAN", "S MAORI", "N MAORI", "BUSHMAN", "AUSTRALI", "TASMANIA", "ARIKARA", "SANTA CR", "PERU", "BERG", "ZALAVAR", "AINU", "NORSE", "EASTER I", "MOKAPU")
#populations <- c("N JAPAN", "S JAPAN", "S MAORI", "N MAORI", "BUSHMAN", "AUSTRALI", "TASMANIA", "ARIKARA", "SANTA CR", "PERU")
#taken from roseman and weaver 2004, only BPC -> BPL (due to a typo? it's the basion-prosthion length )
linMeasNames <- c("GOL", "NOL", "BNL", "BBH", "XCB", "XFB", "STB", "ZYB", "AUB", "WCB", "ASB", "BPL", "NPH", "NLH", "OBH", "OBB", "JUB", "NLB", "MAB", "MDH", "MDB", "ZMB", "SSS", "FMB", "NAS", "EKB", "DKB", "NDS", "WNB", "SIS", "IML", "XML", "MLS", "WMH", "FOL", "FRC", "FRS", "FRF", "PAC", "PAS", "PAF", "OCC", "OCS", "OCF", "VRR", "NAR", "SSR", "PRR", "DKR", "ZOR", "FMR", "EKR", "ZMR", "AVR", "DKS", "SOS", "GLS")
linMeasTraitsHumanLog <- traitsHOWELLlog[,linMeasNames]
rownames(linMeasTraitsHumanLog) <- gsub(" ", "", rownames(linMeasTraitsHumanLog))
linMeasCovHumanLog <- covHumanLog[linMeasNames, linMeasNames]

#### write to file in a way compatible with revbayes ####

ratesTSV(linMeasCovHumanLog, "logHowellsRate_MwF.tsv")
meansNexus(linMeasTraitsHumanLog, "logHowellsTraits_MwF.nex")


