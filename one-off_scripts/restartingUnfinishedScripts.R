# identify completed scripts and create new script with uncompleted scripts

all <- as.vector(as.matrix(read.table("/Users/nikolai/mvBM_sims_PB/job_pb4_5.txt")))
done <- list.files("/Users/nikolai/mvBM_sims_PB/output")
done <- strsplit(done, split = "_")
cleanDone <- vector()
for (i in 1:length(done)){
  if(any(done[[i]][9] == c("1.log", "2.log"))){
    cleanDone[length(cleanDone) + 1] <- gsub("params", "script", gsub("log", "Rev", paste0("scripts/", paste0(done[[i]][1:9], collapse = "_"))))
  }
}
notDone <- setdiff(all, cleanDone)
sink(file = "/Users/nikolai/mvBM_sims_PB/job_pb45_finish.txt", append = F)
for(i in 1:length(notDone)){
  cat(paste0(notDone[i], "\n"), sep = "")
}
sink()
