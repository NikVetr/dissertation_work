library(RColorBrewer)
runs <- sapply(1:100, function(x) cumsum(rnorm(n = 1000, mean = 0, sd = 1)))
cols <- c(brewer.pal(9, "Set1"),brewer.pal(8, "Set2"),brewer.pal(12, "Set3"), brewer.pal(8, "Dark2"))
plot(x,y, xlim=c(min(runs),max(runs)), ylim=c(0,1000), xlab = i, type='n', yaxt = 'n', xaxt = "n", bty = 'n') 
for(i in 1:100){
  lines(runs[,i], 1:1000, type = "l", yaxt = 'n', xaxt = 'n', col = sample(cols, size = 1))
}
