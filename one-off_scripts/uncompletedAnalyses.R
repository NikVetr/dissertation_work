#this script looks in /output to find most recently completed jobs
#then creates a new job.text file by finding which of those jobs have not yet been completed

setwd("/Volumes/2TB/500/full/mvBM_sims_PB_noiseless_500//")
jobs <- as.character(read.table("job.txt")$V1)
jobs <- as.character(read.table("job_pb.txt")$V1)
jobs <- as.character(read.table("job1_100.txt")$V1)
jobs <- c(jobs, as.character(read.table("job101_200.txt")$V1))
jobs <- c(jobs, as.character(read.table("job201_300.txt")$V1))
jobs <- c(jobs, as.character(read.table("job301_400.txt")$V1))
jobs <- c(jobs, as.character(read.table("job401_500.txt")$V1))

noisyness <- strsplit(jobs[1], "_")[[1]]; noisyness <- noisyness[startsWith(noisyness, "nois")]

#by looking at output
outputDir <- "output"
done <- list.files(outputDir)
done <- done[sapply(1:length(done), function(x) endsWith(done[x], "log"))] #get just the log files
# done <- done[order(sapply(1:length(done), function(x) file.mtime(paste0(outputDir, "/",done[x]))), decreasing = T)[-(1:8)]] #remove most recent 8 scripts
done <- sapply(1:length(done), function(x) gsub("params", "script", done[x])) #remove params, replace with script
done <- paste0("scripts/", sapply(1:length(done), function(x) strsplit(done[x], split = ".log")), ".Rev") #replace .log with .Rev, add "scripts/" to front
done <- sapply(1:length(done), function(x) gsub("Strength", paste0("Strength_", noisyness), done[x])) #add noisyness

todo <- setdiff(jobs, done)
sink("job_leftover.txt", append = F)
for(i in 1:length(todo)){
  cat(todo[i], "\n", sep = "")
}
sink()

#by looking at ppss_dir
done <- list.files("ppss_dir/job_log/")
done <- sapply(1:length(done), function(x) gsub("scripts_", "scripts/", done[x])) 
done <- sapply(1:length(done), function(x) gsub("_Rev", ".Rev", done[x])) 
todo <- setdiff(jobs, done)
sink("job_leftover.txt", append = F)
for(i in 1:length(todo)){
  cat(todo[i], "\n", sep = "")
}
sink()
