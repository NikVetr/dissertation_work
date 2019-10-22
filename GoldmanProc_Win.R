d.orig <- read.csv("C:\\Users\\Nikolai\\Google Drive\\dataGD\\GoldmanBB.csv")
d <- d.orig

#find average of left and right sides, else use only side present; if bothg sides absent, mark NA
d$HML <- sapply(1:length(d$LHML), function(x) 
  if(!is.na(d$LHML[x]) & !is.na(d$RHML[x])) {(d$LHML[x] + d$RHML[x])/2}
  else if(!is.na(d$LHML[x])) {d$LHML[x]}
  else if(!is.na(d$RHML[x])) {d$RHML[x]}
  else (NA)
)
d$FHD <- sapply(1:length(d$LFHD), function(x) 
  if(!is.na(d$LFHD[x]) & !is.na(d$RFHD[x])) {(d$LFHD[x] + d$RFHD[x])/2}
  else if(!is.na(d$LFHD[x])) {d$LFHD[x]}
  else if(!is.na(d$RFHD[x])) {d$RFHD[x]}
  else (NA)
)
d$RML <- sapply(1:length(d$LRML), function(x) 
  if(!is.na(d$LRML[x]) & !is.na(d$RRML[x])) {(d$LRML[x] + d$RRML[x])/2}
  else if(!is.na(d$LRML[x])) {d$LRML[x]}
  else if(!is.na(d$RRML[x])) {d$RRML[x]}
  else (NA)
)
d$FML <- sapply(1:length(d$LFML), function(x) 
  if(!is.na(d$LFML[x]) & !is.na(d$RFML[x])) {(d$LFML[x] + d$RFML[x])/2}
  else if(!is.na(d$LFML[x])) {d$LFML[x]}
  else if(!is.na(d$RFML[x])) {d$RFML[x]}
  else (NA)
)
d$TML <- sapply(1:length(d$LTML), function(x) 
  if(!is.na(d$LTML[x]) & !is.na(d$RTML[x])) {(d$LTML[x] + d$RTML[x])/2}
  else if(!is.na(d$LTML[x])) {d$LTML[x]}
  else if(!is.na(d$RTML[x])) {d$RTML[x]}
  else (NA)
)

#get only traits in Roseman Data
d <- d[,c("Sex", "BIB", "HML", "FHD", "RML", "FML", "TML", "Brachial", "Crural", "NOTE")]

#discard individuals with missing data (impute later?)
d <- d[complete.cases(d),]

#use only men to avoid complicating effects of sexual dimosrphism (and agree with Roseman procedure)
men <- d[d$Sex == 0,]

#only use populations with more than a certain # of individuals, in this case arbitrarily set to 9
inclOverNine <- names(summary(men$NOTE)[summary(men$NOTE) > 9])[-35]
inclOverNine <- inclOverNine[-21]
data <- data.frame()
for(i in 1:length(inclOverNine)){
  data <- rbind(data, men[men$NOTE == inclOverNine[i],])
}

#use only traits of interest and population info
d <- data[,c("BIB", "HML", "FHD", "RML", "FML", "TML", "Brachial", "Crural", "NOTE")]

#compute the pooled within-group variance-covariance matrix
head(as.matrix(d))
unique(d$NOTE)[1]
sp <- sapply(1:length(unique(d$NOTE)), function (x) as.matrix((t(d[d$NOTE == unique(d$NOTE)[x],-(9)]))))
ntraits <- 8
nspecies <- 32

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
rownames(cov) <- colnames(cov) <- colnames(d)[-9]

#decompose covariance matrix to sds and correlations to reconstruct in rb
library(corpcor)
vars <- decompose.cov(cov)$v
sds <- vars^.5
corrs <- decompose.cov(cov)$r
corrsUT <- corrs[upper.tri(corrs)]
corrsUT

#get a vector representing the upper triangular elements of a matrix the way god intended
#i.e. in English, read left to right
upperTriLR <- function(matr){
  if(!is.matrix(matr)){print("ERROR, ERROR, NOT A MATRIX"); break()}
  UT <- vector()
  nrows <- dim(matr)[1]
  ncols <- dim(matr)[2]
  for (i in 1:(nrows-1)){
    for(j in (i+1):(ncols)){
      UT <- c(UT, matr[i,j])
    }
  }
  return(UT)
}

corrsUT <- upperTriLR(corrs)
cat(corrsUT)

#prepare data for import into rb
means <- read.csv("C:\\Users\\Nikolai\\Google Drive\\BodegaProject\\BBRosemanGoldmanData.csv", header = T)
names <- means$NAME
means <- means[,c("BIB", "HML", "FHD", "RML", "FML", "TML", "Brachial", "Crural", "LAT")]
rownames(means) <- names
means

meansNexus <- function (characterMatrix, outputFilePath, traitsOnColumns = T){
  if (traitsOnColumns) {ntax = dim(characterMatrix)[1]; nchar = dim(characterMatrix)[2]}
  else {ntax = dim(characterMatrix)[2]; nchar = dim(characterMatrix)[1]}
  sink(outputFilePath, append=F)
  cat("#NEXUS\n\nBegin data;\n\tDimensions ntax=", ntax, "nchar=", nchar, ";\n\tFormat datatype=Continuous missing=?;\n\tMatrix\n\n")
  for(i in 1:length(characterMatrix[,1])) {cat(rownames(characterMatrix)[i], characterMatrix[i,], "\n")}
  cat(";\nEnd;")
  sink()
}

#write data to nexus file
meansNexus(as.matrix(means), "C:\\Users\\nvetr\\goldMeans.nex")

### ERROR WITH IMPORT ###

#prepare dramatically reduced dataset for import into rb
means <- read.csv("/Users/nvetr/Google Drive/BodegaProject/BBRosemanGoldmanData.csv", header = T)
names <- means$NAME
means <- means[,c("BIB", "Crural", "LAT")]
rownames(means) <- names
means
meansNexus(as.matrix(means), "C:\\Users\\nvetr\\goldMeansRed2.nex")

reduCov <- cov[c("BIB", "Crural"), c("BIB", "Crural")]
library(corpcor)
rvars <- decompose.cov(reduCov)$v
rsds <- rvars^.5
rcorrs <- decompose.cov(reduCov)$r
rcorrsUT <- rcorrs[upper.tri(rcorrs)]
rcorrs


#plot bivariate data to show relationship when failing to account for phylogeny/population history
library(ggplot2)
ggplot(means, aes(x=LAT, y=BIB)) + geom_point(shape=1) + geom_smooth(method=lm, se=FALSE) + theme_bw() +
  labs(x="Latitude",y="Bi-illiac Breadth")

ggplot(means, aes(x=LAT, y=Crural)) + geom_point(shape=1) + geom_smooth(method=lm, se=FALSE) + theme_bw()+
  labs(x="Latitude",y="Crural Index")

  
