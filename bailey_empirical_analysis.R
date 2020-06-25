setwd("/Volumes/1TB/Bailey/")
d <- read.csv("bailey_data.csv", header = T, stringsAsFactors = T)
d_orig <-  d
levels(d$UI1.TD)
d[d == ""] <- NA
d[d == " "] <- NA
poss <- sapply(3:139, function(trait) paste0(colnames(d)[trait], ": ", paste0(levels(d[,trait]), collapse = ", "), "\n"))

nf <- which(sapply(3:139, function(trait) class(d[,trait])) != "factor")
poss_nf <- sapply(nf, function(trait) paste0(colnames(d)[trait], ": ", paste0(unique(d[,trait]), collapse = ", "), "\n"))
poss[nf] <- poss_nf

cat(poss)
