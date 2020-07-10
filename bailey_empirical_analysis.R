setwd("/Volumes/1TB/Bailey/")
d <- read.csv("bailey_data.csv", header = T, stringsAsFactors = F)
d_orig <-  d
d <- t(sapply(1:nrow(d), function(row_of_d) trimws(d[row_of_d,])))
colnames(d) <- colnames(d_orig)

#recode the data
d[d == ""] <- NA
d[d == " "] <- NA
d[d == "-"] <- NA
d[d == "*"] <- NA
d[d == "."] <- NA
d[d == "?"] <- NA
d[d == "--"] <- NA
d[d == ".."] <- NA
d[d == "worn"] <- NA
d[which(d == "Y", arr.ind = T)] <- "1+"
d[which(d == "X", arr.ind = T)] <- "1+"
d[which(d == "x", arr.ind = T)] <- "1+"
d[d == "retained di1 btwn I1s"] <- NA
d[which(d == "+", arr.ind = T)] <- "1+"
d[do.call(rbind, sapply(1:10, function(name) which(d == c("D", "D (asym)", "dis", "dis" , "DM", "m", "M", "MD", "med", "mes")[name], arr.ind = T)))] <- "1+"
d[which(d == "TRI", arr.ind = T)] <- 6
d[do.call(rbind, sapply(1:3, function(name) which(d == c("P", "P/R?", "R")[name], arr.ind = T)))] <- "1"
d[which(d == "A", arr.ind = T)] <- 0
d[d == "Lmissing"] <- NA

UPM3.BMxP_MR <- d[,"UPM3.BMxP"]
UPM3.BMxP_MR[which(nchar(UPM3.BMxP_MR) > 1)] <- sapply(1:length(UPM3.BMxP_MR[which(nchar(UPM3.BMxP_MR) > 1)]), function(indiv) substr(UPM3.BMxP_MR[which(nchar(UPM3.BMxP_MR) > 1)][indiv], 1, 1))
UPM3.BMxP_DR <- d[,"UPM3.BMxP"]
UPM3.BMxP_DR[which(nchar(UPM3.BMxP_DR) > 1)] <- sapply(1:length(UPM3.BMxP_DR[which(nchar(UPM3.BMxP_DR) > 1)]), function(indiv) substr(UPM3.BMxP_DR[which(nchar(UPM3.BMxP_DR) > 1)][indiv], 3, 3))
UPM3.BMxP_DR[which(UPM3.BMxP_DR == ",")] <- NA

UPM3.LMxP_MR <- d[,"UPM3.LMxP"]
UPM3.LMxP_MR[which(nchar(UPM3.LMxP_MR) > 1)] <- sapply(1:length(UPM3.LMxP_MR[which(nchar(UPM3.LMxP_MR) > 1)]), function(indiv) substr(UPM3.LMxP_MR[which(nchar(UPM3.LMxP_MR) > 1)][indiv], 1, 1))
UPM3.LMxP_DR <- d[,"UPM3.LMxP"]
UPM3.LMxP_DR[which(nchar(UPM3.LMxP_DR) > 1)] <- sapply(1:length(UPM3.LMxP_DR[which(nchar(UPM3.LMxP_DR) > 1)]), function(indiv) substr(UPM3.LMxP_DR[which(nchar(UPM3.LMxP_DR) > 1)][indiv], 3, 3))

d <- cbind(d, UPM3.BMxP_MR = UPM3.BMxP_MR)
d <- cbind(d, UPM3.BMxP_DR = UPM3.BMxP_DR)
d <- cbind(d, UPM3.LMxP_MR = UPM3.LMxP_MR)
d <- cbind(d, UPM3.LMxP_DR = UPM3.LMxP_DR)
d <- d[,-which(colnames(d) == "UPM3.LMxP")]
d <- d[,-which(colnames(d) == "UPM3.BMxP")]
d <- d[,-which(colnames(d) == "UPM4.ROT")]
d <- d[,-which(colnames(d) == "LP4.PAT")]



UPM4.LMxP_MR <- d[,"UPM4.LMxP"]
UPM4.LMxP_MR[which(nchar(UPM4.LMxP_MR) > 1)] <- sapply(1:length(UPM4.LMxP_MR[which(nchar(UPM4.LMxP_MR) > 1)]), function(indiv) substr(UPM4.LMxP_MR[which(nchar(UPM4.LMxP_MR) > 1)][indiv], 1, 1))
UPM4.LMxP_DR <- d[,"UPM4.LMxP"]
UPM4.LMxP_DR[which(nchar(UPM4.LMxP_DR) > 1)] <- sapply(1:length(UPM4.LMxP_DR[which(nchar(UPM4.LMxP_DR) > 1)]), function(indiv) substr(UPM4.LMxP_DR[which(nchar(UPM4.LMxP_DR) > 1)][indiv], 3, 3))
d <- cbind(d, UPM4.LMxP_MR = UPM4.LMxP_MR)
d <- cbind(d, UPM4.LMxP_DR = UPM4.LMxP_DR)
d <- d[,-which(colnames(d) == "UPM4.LMxP")]

UPM4.BMxP_MR <- d[,"UPM4.BMxP"]
UPM4.BMxP_MR[which(nchar(UPM4.BMxP_MR) > 1)] <- sapply(1:length(UPM4.BMxP_MR[which(nchar(UPM4.BMxP_MR) > 1)]), function(indiv) substr(UPM4.BMxP_MR[which(nchar(UPM4.BMxP_MR) > 1)][indiv], 1, 1))
UPM4.BMxP_DR <- d[,"UPM4.BMxP"]
UPM4.BMxP_DR[which(nchar(UPM4.BMxP_DR) > 1)] <- sapply(1:length(UPM4.BMxP_DR[which(nchar(UPM4.BMxP_DR) > 1)]), function(indiv) substr(UPM4.BMxP_DR[which(nchar(UPM4.BMxP_DR) > 1)][indiv], 3, 3))
d <- cbind(d, UPM4.BMxP_MR = UPM4.BMxP_MR)
d <- cbind(d, UPM4.BMxP_DR = UPM4.BMxP_DR)
d <- d[,-which(colnames(d) == "UPM4.BMxP")]

d[which(d == "1x2", arr.ind = T)] <- 1
d[which(d == "1X2", arr.ind = T)] <- 1
d[which(d == "2x1", arr.ind = T)] <- 1
d[d == "Rev"] <- NA
d[which(d[, "UM1.EE"] < 1, arr.ind = T), "UM1.EE"] <- 0
d[which(d[, "UM2.EE"] < 1, arr.ind = T), "UM2.EE"] <- 0
d[which(d[, "LM1.EE"] < 1, arr.ind = T), "LM1.EE"] <- 0
d[which(d[, "LM2.EE"] < 1, arr.ind = T), "LM2.EE"] <- 0
d[which(d[, "LM2.AF"] < 1, arr.ind = T), "LM2.AF"] <- 0
d[which(d == "<3", arr.ind = T)] <- "3-"
d[which(d == "2x2", arr.ind = T)] <- "2"
d[which(d == "2,2", arr.ind = T)] <- "2"
d[which(d == "Y,+", arr.ind = T)] <- "1+"
d[which(d == "", arr.ind = T)] <- 
d[which(d[, "UM2.EE"] < 1, arr.ind = T), "UM2.EE"] <- 0  

d[do.call(rbind, sapply(1:9, function(name) 
  which(d == c("agenesis" , "cariou", "could be agenesis need radiograph", "could be agenesis on right", 
  "left missing postmortem", "IN CRYPT", "in crypts", "post mort loss", "R missing")[name], arr.ind = T)))] <- NA
d[t(sapply(1:5, function(name) which(d == c("L missing", "missing check CT", "Lmissing, Rimpacted", 
                                            "missing M3s", "could be agenesis on right, left missing postmortem")[name], arr.ind = T)))] <- NA


d[do.call(rbind, sapply(1:3, function(name) which(d == c("P", "P/R", "R")[name], arr.ind = T)))] <- "1"
d[which(d == "AB", arr.ind = T)] <- 0
d[which(d == "R 45 M left", arr.ind = T)] <- NA
d[which(d == "2m", arr.ind = T)] <- 1
d[which(d == "tomes", arr.ind = T)] <- NA
d[which(d == "2Bif", arr.ind = T)] <- "1+"

d[which(d == "Z", arr.ind = T)] <- "1+"
d[which(d == "1M", arr.ind = T)] <- "1+"
d[which(d == "0-1", arr.ind = T)] <- "1-"
d[which(d == "2.2", arr.ind = T)] <- "2"
d[which(d == "pit", arr.ind = T)] <- "0"
d[which(d == "PIT", arr.ind = T)] <- "0"
d[which(d == ">2", arr.ind = T)] <- "2+"
d[which(d == ">1", arr.ind = T)] <- "1+"
d[which(d == "4  +", arr.ind = T)] <- "4+"
d[which(d == "1A", arr.ind = T)] <- "0"
d[which(d == "2*", arr.ind = T)] <- "2"

# d[which(d == "9", arr.ind = T)] 
d[which(d[, "LM1.EE"] == 9, arr.ind = T), "LM1.EE"] <- NA
d[which(d[, "LM2.4CUS"] == 0, arr.ind = T), "LM2.4CUS"] <- "5+"
d[which(d[, "LM3.4CUS"] == 0, arr.ind = T), "LM3.4CUS"] <- "5+"

d[do.call(rbind, sapply(1:3, function(name) which(d == c("agenesis", "impacted", "missing")[name], arr.ind = T)))] <- NA

d[d == "-"] <- NA
d[which(d == "-1", arr.ind = T)] <- "1"
d[which(d == "11", arr.ind = T)] <- "1"
d[which(d == "1*", arr.ind = T)] <- "1"
d[which(d == "2.0", arr.ind = T)] <- "2"
d[which(d == "4.0", arr.ind = T)] <- "4"
d[which(d == "3.0", arr.ind = T)] <- "3"
d[which(d == "0.0", arr.ind = T)] <- "0"
d[which(d == "NA", arr.ind = T)] <- NA
d[which(d == "3.5", arr.ind = T)] <- "3.4"
d[which(d == "2.5", arr.ind = T)] <- "2.3"
d[which(d == "3.1", arr.ind = T)] <- "3.4"
d[which(d == "1.0", arr.ind = T)] <- "1"

#check for remaining weird values
unique(as.vector(d[,-(1:2)]))
d[which(d == "11", arr.ind = T)] 
colnames(d)[which(d == "1.0", arr.ind = T)[,2]]

#ok, now let's code the ambiguity matrix

#looking at the possible options for each trait
poss <- sapply(3:139, function(trait) paste0(colnames(d)[trait], ": ", paste0(levels(d[,trait]), collapse = ", "), "\n"))
nf <- which(sapply(3:139, function(trait) class(d[,trait])) != "factor")
poss_nf <- sapply(nf, function(trait) paste0(colnames(d)[trait], ": ", paste0(unique(d[,trait]), collapse = ", "), "\n"))
poss[nf] <- poss_nf
cat(poss)

#missing data patterns
present <- as.matrix(!is.na(d[,-c(1:2)]))
present[present] <- 1

popcols <- RColorBrewer::brewer.pal(length(unique(d[,1])), "Dark2")
popcols <- popcols[as.factor(d[,1])]
row_i <- as.vector(matrix(1:nrow(present), nrow = nrow(present), ncol = ncol(present)))
col_j <- as.vector(matrix(1:ncol(present), nrow = nrow(present), ncol = ncol(present), byrow = T))
cols <- c("white", "black")
png("/Volumes/macOS/Users/nikolai/Pictures/bailey_missingdatavisualizaton.png", width = 800, height = 2500)
plot(1, 1, ylim = c(20,dim(present)[1]), xlim = c(1,dim(present)[2]), col = "white", bty = "n", xlab = "", ylab = "")
title(xlab = "trait index", ylab = "individual index", line = 2.25,  cex.lab = 2)
title(main = "missing data visualization\n colors indicate population", cex.main = 4, line = -3)
title(main = "black dots indicate trait presence, white dots trait absence", cex.main = 2, line = -5)
points(x = col_j, y = row_i, col = cols[as.vector(present)+1], cex = 0.5, pch = 15)
points(x = rep(-3, length(popcols)), y = 1:length(popcols), col = popcols, cex = 1.5, pch = 15, cex.main = 2)
dev.off()

par(mfrow = c(1,2))
plot(sort(apply(present, 1, mean), decreasing = T), type = "l", xlab = "individual index", ylab = "proportion traits present")
plot(sort(apply(present, 2, mean), decreasing = T), type = "l", xlab = "trait index", ylab = "proportion individuals present")

prop_indiv_filter <- 0
prop_traits_filter <- 0.2
present_filtered <- present[apply(present, 1, mean) > prop_indiv_filter, apply(present, 2, mean) > prop_traits_filter]
all(levels(d[,1])[as.integer(d[,1])] == d[,1])
npops <- length(unique(as.integer(d[,1])))
ntraits <- length(d[1,]) - 2
per_pop_ntraits_thresh <- 10
sum(sapply(1:ntraits, function(trait) all(sapply(1:npops, function(pop) sum(present[,trait][as.integer(d[,1]) == pop]))  > per_pop_ntraits_thresh)))

par(mfrow = c(1,1))
xloc = c(8,25,12,20,17,18,16,11)
for(i in 1:8){
  ngroups_threshold <- i
  y <- sapply(1:30, function(per_pop_ntraits_thresh) 
    sum(sapply(1:ntraits, function(trait) sum(sapply(1:npops, function(pop) sum(present[,trait][as.integer(d[,1]) == pop]))  > per_pop_ntraits_thresh) >= ngroups_threshold))
  )
  if(i == 1){
    plot(1:30, y, type = "l", xlab = "number of a traits observed in * populations threshold", ylab = "number of traits", xlim = c(1,30), ylim = c(0,140))
    text(x = xloc[i], y = y[xloc[i]]+3, labels = paste0("* = ", i))
  } else {
    lines(1:30, y)
    text(x = xloc[i], y = y[xloc[i]]+3, labels = paste0("* = ", i))
  }
}

# x <- (sapply(1:ntraits, function(trait) (sapply(1:npops, function(pop) sum(present[,trait][as.integer(d[,1]) == pop])))))
# apply(x > 8, 1, sum)

#filter out traits based in inflection point
ngroups_threshold <- 7
per_pop_ntraits_thresh <- 8
trait_inds_that_satisfy <- sapply(1:ntraits, function(trait) sum(sapply(1:npops, function(pop) sum(present[,trait][as.integer(d[,1]) == pop]))  > per_pop_ntraits_thresh) >= ngroups_threshold)
colnames(d[,-(1:2)])[!trait_inds_that_satisfy]

#filter out traits and prepare final dataset
traits_to_exclude_because_of_science <- 
