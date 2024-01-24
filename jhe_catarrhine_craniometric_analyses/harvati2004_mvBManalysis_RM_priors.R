library(geomorph)
library(Matrix)
library(phytools)
library(phangorn)

standardize <- F #make traits all have unit variance in RM
superStandardize <- T
noPCA <- T
PCA_covBG <- F
collapseHomo <- F
collapseSubsp <- F

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

data <- readland.nts("/Users/nikolai/data/allPD4.nts")
class(data)
str(data) #landmarks on the rows, dimensions (xyz) on the columns, species on the slices

#mean center all landmarks? landmark by landmark or vs the overall mean?

codes <- attributes(data)$dimnames[[3]]
codes_linnaen <- substr(codes, 1, 5); codes_linnaen <- sort(unique(codes_linnaen))
codes_linnaen <- cbind(codes_linnaen, sapply(1:length(codes_linnaen), function(x) sum(startsWith(codes, paste0(codes_linnaen[x], "M")))),
                       sapply(1:length(codes_linnaen), function(x) sum(startsWith(codes, paste0(codes_linnaen[x], "F")))))
colnames(codes_linnaen) <- c("code", "males", "females")
taxa <- read.csv("/Users/nikolai/data/Harvati2004_Taxa.csv", header = T)
sexMatch <- sapply(1:length(codes_linnaen[,1]), function(y) 
  which(sapply(1:length(taxa[,1]), function(x) 
    all.equal(as.numeric(codes_linnaen[y,2:3]), as.numeric(taxa[x,c("Males", "Females")])) == T))[1])
sexMatch <- unlist(replace(sexMatch, !sapply(sexMatch,length),NA))
codes_linnaen <- cbind(codes_linnaen, as.character(taxa[sexMatch,1]))
codes_linnaen[startsWith(codes_linnaen[,1], "PPPXX"),4] <- "Papio hamadryas papio" #process of elimination!
codes_linnaen[startsWith(codes_linnaen[,1], "HMN"),4] <- "Homo neanderthalensis"
codes_linnaen[startsWith(codes_linnaen[,1], "HMSAM"),4] <- "Homo sapiens (UP)"
codes_linnaen[startsWith(codes_linnaen[,1], "HMSUP"),4] <- "Homo sapiens (UP)"
codes_linnaen[startsWith(codes_linnaen[,1], "HMSAU"),4] <- "Homo sapiens (Australasian)" #right
codes_linnaen[startsWith(codes_linnaen[,1], "HMSTL"),4] <- "Homo sapiens (Australasian)" #right
codes_linnaen[startsWith(codes_linnaen[,1], "HMSEU"),4] <- "Homo sapiens (Eurasian)" #right
codes_linnaen[startsWith(codes_linnaen[,1], "HMSBG"),4] <- "Homo sapiens (Eurasian)" #right
codes_linnaen[startsWith(codes_linnaen[,1], "HMSDG"),4] <- "Homo sapiens (African)"
codes_linnaen[startsWith(codes_linnaen[,1], "HMSES"),4] <- "Homo sapiens (Greenland)"
codes_linnaen[startsWith(codes_linnaen[,1], "HMSSN"),4] <- "Homo sapiens (African)"
codes_linnaen[startsWith(codes_linnaen[,1], "HMS") & is.na(codes_linnaen[,4]),4] <- paste0(
  "Homo sapiens (Unknown ", 1:length(codes_linnaen[startsWith(codes_linnaen[,1], "HMS") & is.na(codes_linnaen[,4]),4]), ")")
if(collapseHomo){
  codes_linnaen[startsWith(codes_linnaen[,1], "HMS"), 4] <- "Homo sapiens"
}
subsp <- unique(codes_linnaen[,4])
if(collapseSubsp){
  codes_linnaen[,4] <- sapply(1:length(codes_linnaen[,4]), function(x) paste(strsplit(codes_linnaen[,4], " ")[[x]][1:2], collapse = " "))
}
species <- unique(codes_linnaen[,4])

#what is the 5th letter? subspecies? probably not, since TT is troglodytes troglodytes
# coded_genus <- substr(codes, 1, 2)
# genus_key <- cbind(c("HM", "MC", "MN", "PP", "PN", "GO"), c("Homo", "Macaca", "Mandrillus", "Papio", "Pan", "Gorilla"))
# sex <- substr(codes, 6, 6)
# linnaen <- as.character(unique(taxa[,1]))
# genus <- sapply(1:length(linnaen), function(x) strsplit(linnaen[x], split = " ")[[1]][1])
# species <- sapply(1:length(linnaen), function(x) strsplit(linnaen[x], split = " ")[[1]][2])
# subspecies <- sapply(1:length(linnaen), function(x) strsplit(linnaen[x], split = " ")[[1]][3])

#load in the molecular tree

#mean center all landmarks
for(i in 1:15){for(j in 1:3){data[i,j,] <- data[i,j,] - mean(data[i,j,])}}
t.data <- gpagen(data)$data
t.data$logCsize <- log(t.data$Csize)
t.data <- as.matrix(subset(x = t.data, select = colnames(t.data)[colnames(t.data) != "Csize"]))
covTotal <- cov(t.data) #total covariance matrix
if(PCA_covBG == T){
  d_PCA <- as.data.frame(t.data)
  ntraits_PCA <- length(d_PCA[1,])
  Population <- rep("species", nrow(d_PCA))
  for(i in 1:nrow(codes_linnaen)){Population[startsWith(rownames(d_PCA), codes_linnaen[,1][i])] <- codes_linnaen[,4][i]}
  d_PCA$Population <- Population
  pops_PCA <- sapply(1:length(unique(d_PCA$Population)), function (x) as.matrix((t(d_PCA[d_PCA$Population == unique(d_PCA$Population)[x],-(47)]))))
  
  #get traits
  traits_PCA <- list()
  for(i in 1:ntraits_PCA){ 
    print(i)
    trait <- list()
    trait[[1]] <- as.numeric(pops_PCA[[1]][i,])
    for (j in 2:27){
      trait[[j]] <- as.numeric(pops_PCA[[j]][i,])
    } 
    traits_PCA[[i]] <- trait
  }
  
  #make cov matrix
  covBG <- matrix(nrow = ntraits_PCA, ncol = ntraits_PCA)
  for (i in 1:ntraits_PCA){
    print(i)
    for(j in 1:ntraits_PCA){
      trait1 <- traits_PCA[[i]]
      trait2 <- traits_PCA[[j]]
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
      covBG[i,j] <- covariance
    }
  }
}


#do PCA on procrustes-transformed data using total covariance matrix
if(PCA_covBG != T){
  V <- eigen(covTotal)$vectors
  L <- eigen(covTotal)$values
} else{
  V <- eigen(covBG)$vectors
  L <- eigen(covBG)$values
}
VarExpl <- round(L / sum(L), 4) * 100
PC_Scores <- t(V) %*% t(t.data) 
if(normalize){
  PC_Scores <- PC_Scores / L^0.5 #set sd to 1
  PC_Scores <- PC_Scores - apply(PC_Scores, 1, mean) #mean center
}
PC_Scores <- PC_Scores[1:(nrow(PC_Scores)-7),]

if(!noPCA){
  d <- as.data.frame(t(PC_Scores))
} else {
  d <- as.data.frame(t.data)
}

ntraits <- length(d[1,])


Population <- rep("species", nrow(d))
for(i in 1:nrow(codes_linnaen)){Population[startsWith(rownames(d), codes_linnaen[,1][i])] <- codes_linnaen[,4][i]}
d$Population <- Population
if(!noPCA){
  pops <- sapply(1:length(unique(d$Population)), function (x) as.matrix((t(d[d$Population == unique(d$Population)[x],-(41)]))))
} else {
  pops <- sapply(1:length(unique(d$Population)), function (x) as.matrix((t(d[d$Population == unique(d$Population)[x],-(47)]))))
}
#get traits
traits <- list()
npop <- length(pops)
for(i in 1:ntraits){ 
  print(i)
  trait <- list()
  trait[[1]] <- as.numeric(pops[[1]][i,])
  for (j in 2:npop){
    trait[[j]] <- as.numeric(pops[[j]][i,])
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
if(noPCA){cov <- nearPD(cov)$mat}

#compute matrix of population means
traitsHARVATI <- matrix(nrow=npop, ncol=ntraits)
for(i in 1:npop){
  for(j in 1:ntraits){
    traitsHARVATI[i,j] <- mean(as.numeric(pops[[i]][j,]))
  }
}
rownames(traitsHARVATI) <- gsub(x = gsub(x = gsub(x = unique(d$Population), 
                                                  pattern =  " ", replacement =  "_"), 
                                         pattern = ")", replacement = ""),
                                pattern = "\\(", replacement = "")

colnames(traitsHARVATI) <- colnames(d)[-(ntraits+1)]

rownames(cov) <- colnames(cov) <- colnames(d)[-(ntraits+1)]
cov <- as.matrix(cov)


meansNexus <- function (characterMatrix, outputFilePath, traitsOnColumns = T){
  if (traitsOnColumns) {ntax = dim(characterMatrix)[1]; nchar = dim(characterMatrix)[2]}
  else {ntax = dim(characterMatrix)[2]; nchar = dim(characterMatrix)[1]}
  sink(outputFilePath, append=F)
  cat("#NEXUS\n\nBegin data;\n\tDimensions ntax=", ntax, "nchar=", nchar, ";\n\tFormat datatype=Continuous missing=?;\n\tMatrix\n\n")
  for(i in 1:length(characterMatrix[,1])) {cat(rownames(characterMatrix)[i], characterMatrix[i,], "\n")}
  cat(";\nEnd;")
  sink()
}

#write data for uninformative rate matrix prior
if(standardize){
  traitsHARVATI <- sapply(1:ncol(traitsHARVATI), function(trait) (traitsHARVATI[,trait] - mean(traitsHARVATI[,trait])) / sqrt(diag(cov)[trait]))
}
if(standardize){meansNexus(traitsHARVATI, outputFilePath = "~/data/Harvati_noPCA_homopops_standardized.nex")}

#write data for informative rate matrix prior
if(superStandardize){
  L <- t(chol(cov))
  Li <- solve(L)
  traitsHARVATI <- t(Li %*% t(traitsHARVATI))
  #mean center
  traitsHARVATI <- sapply(1:ncol(traitsHARVATI), function(trait) traitsHARVATI[,trait] - mean(traitsHARVATI[,trait]))
  meansNexus(traitsHARVATI, outputFilePath = "~/data/Harvati_noPCA_homopops_superStandardized.nex")
}


R <- rlkj(10); sds <- exp(rnorm(10, 1, 0.5)); R <- diag(sds) %*% R %*% diag(sds)
x <- mvtnorm::rmvnorm(1E5, sigma = R)
xc <- t(solve(t(chol(R))) %*% t(x))
cov(xc)
