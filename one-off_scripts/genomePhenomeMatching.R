Howells <- read.csv("Howells_GeoCoordinates.csv")
Howells
#convert to +/- format
Lats <- sapply(1:length(Howells$Lat), function (i) as.numeric(substr(x = as.character(Howells$Lat[i]), 1, nchar(as.character(Howells$Lat[i]))-1)))
LatSigns <- sapply(1:length(Howells$Lat), function (i) substr(x = as.character(Howells$Lat[i]), start = nchar(as.character(Howells$Lat[i])), stop =  nchar(as.character(Howells$Lat[i]))))
LatSigns <- (LatSigns == "N") * 2 - 1
Howells$Lat <- Lats * LatSigns

Longs <- sapply(1:length(Howells$Long), function (i) as.numeric(substr(x = as.character(Howells$Long[i]), 1, nchar(as.character(Howells$Long[i]))-1)))
LongSigns <- sapply(1:length(Howells$Long), function (i) substr(x = as.character(Howells$Long[i]), start = nchar(as.character(Howells$Long[i])), stop =  nchar(as.character(Howells$Long[i]))))
LongSigns <- (LongSigns == "E") * 2 - 1
Howells$Long <- Longs * LongSigns


SGDP <- read.csv("SGDP_metadata.csv")
subSGDP <- SGDP[,c("Illumina_ID", "Population_ID", "Latitude", "Longitude")]

library(geosphere)
?distm (c(lon1, lat1), c(lon2, lat2), fun = distHaversine)

dists <- lapply(1:length(Howells$Pop), function (x) list(as.character(Howells$Pop[x])))
closest <- data.frame()
for(i in 1:length(dists)){
  print(i)
  HPop <- as.character(Howells$Pop[i])
  HLat <- Howells$Lat[i]
  HLong <- Howells$Long[i]
  dist <- vector()
  for (j in 1:length(subSGDP$Population_ID)) {
    SLat <- subSGDP$Latitude[j]
    SLong <- subSGDP$Longitude[j]
    dist[j] <- distm (c(HLong, HLat), c(SLong, SLat), fun = distHaversine)
  }
  ppDists <- data.frame(rep(HPop, length(dist)), subSGDP$Population_ID, dist, subSGDP$Illumina_ID)
  ppDists <- ppDists[order(ppDists$dist),]
  colnames(ppDists) <- c("HPop", "SPop", "Haversine.dist", "S.ID")
  dists[[i]][[2]] <- ppDists
  closest <- rbind(closest, ppDists[1,])

}
