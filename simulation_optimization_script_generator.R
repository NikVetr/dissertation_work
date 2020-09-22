setwd("/Volumes/1TB/Bailey/")

lotsaDATA <- T
rep_i <- 1
optim_rep_i <- 1
nreps <- 1:500
if(lotsaDATA){
  nreps <- 501:1000
}
n_optim_reps <- 4
for(rep_i in nreps){
  for(optim_rep_i in 1:n_optim_reps){
    temp_script <- readLines("~/scripts/dissertation_work/bailey_simulations_step1_metascript.R", warn = F)
    if(lotsaDATA){temp_script[31] <- "setwd(\"/Volumes/1TB/Bailey_LOTSADATA/\")"}
    if(lotsaDATA){temp_script[13] <- "ignoreImputedData <- T"}
    if(lotsaDATA){temp_script[2045] <- "current_params <- recompileMeansCorrsThreshes(par, dat)"}
    temp_script[1] <- paste0("rep_i <- ", rep_i)
    temp_script[2] <- paste0("optim_rep_i <- ", optim_rep_i)
    temp_script[1599] <- paste0("init_cor <- diag(d_traits)")
    writeLines(temp_script, paste0("scripts_sim/script_replicate_", rep_i, "_optim_rep_", optim_rep_i, ".R"))
  }
}

for(rep_i in nreps){
  for(optim_rep_i in 1:n_optim_reps){
    sink(file = paste0("jobs_bailey_sims", ifelse(lotsaDATA, "_lotsaDATA", "") ,".txt"), append = !(rep_i == 1 & optim_rep_i == 1))
    cat(paste0("cd /Volumes/1TB/Bailey; /usr/local/bin/RScript scripts_sim/script_replicate_", rep_i, "_optim_rep_", optim_rep_i, ".R\n"))
    sink()
  }
}


cat("parallel --jobs 8 --sshloginfile /Volumes/macOS/Users/nikolai/instances_noDog --eta --results bailey_sims_screenout --files < jobs_bailey_sims.txt")
