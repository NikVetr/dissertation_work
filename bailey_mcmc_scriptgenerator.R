setwd("/Volumes/1TB/Bailey/")

lotsaDATA <- F
nreps <- 1:1000
if(lotsaDATA){
  nreps <- 501:1000
}
n_mcmc_reps <- 4
for(rep_i in nreps){
  for(mcmc_rep_i in 1:n_mcmc_reps){
    temp_script <- readLines("~/scripts/dissertation_work/bailey_mcmc_metascript.R", warn = F)
    temp_script[3] <- paste0("rep_i <- ", rep_i)
    temp_script[4] <- paste0("mcmc_rep_i <- ", mcmc_rep_i)
    writeLines(temp_script, paste0("scripts_sim/script_replicate_", rep_i, "_mcmc_rep_", mcmc_rep_i, ".R"))
  }
}

for(rep_i in nreps){
  for(mcmc_rep_i in 1:n_mcmc_reps){
    sink(file = paste0("jobs_bailey_mcmc_sims", ifelse(lotsaDATA, "_lotsaDATA", "") ,".txt"), append = !(rep_i == 1 & mcmc_rep_i == 1))
    cat(paste0("cd /Volumes/1TB/Bailey; /usr/local/bin/RScript scripts_sim/script_replicate_", rep_i, "_mcmc_rep_", mcmc_rep_i, ".R\n"))
    sink()
  }
}

cat("parallel --jobs 8 --sshloginfile /Volumes/macOS/Users/nikolai/instances_noDog --eta --results bailey_sims_mcmc_screenout --files < jobs_bailey_mcmc_sims.txt")

#condition on true means and correlations
nreps <- 1:500
if(lotsaDATA){
  nreps <- 501:1000
}
n_mcmc_reps <- 4
for(rep_i in nreps){
  for(mcmc_rep_i in 1:n_mcmc_reps){
    temp_script <- readLines("~/scripts/dissertation_work/bailey_mcmc_metascript_truemeans-corrs.R", warn = F)
    temp_script[3] <- paste0("rep_i <- ", rep_i)
    temp_script[4] <- paste0("mcmc_rep_i <- ", mcmc_rep_i)
    writeLines(temp_script, paste0("scripts_sim/truevals_script_replicate_", rep_i, "_mcmc_rep_", mcmc_rep_i, ".R"))
  }
}

for(rep_i in nreps){
  for(mcmc_rep_i in 1:n_mcmc_reps){
    sink(file = paste0("jobs_bailey_mcmc_sims_truevals", ifelse(lotsaDATA, "_lotsaDATA", "") ,".txt"), append = !(rep_i == 1 & mcmc_rep_i == 1))
    cat(paste0("cd /Volumes/1TB/Bailey; /usr/local/bin/RScript scripts_sim/truevals_script_replicate_", rep_i, "_mcmc_rep_", mcmc_rep_i, ".R\n"))
    sink()
  }
}

cat("parallel --jobs 8 --sshloginfile /Volumes/macOS/Users/nikolai/instances_noDog --eta --results bailey_sims_mcmc_screenout --files < jobs_bailey_mcmc_sims_truevals.txt")


#specify reduced set of just 100 reps
nreps <- c(1:100, 501:600, 1001:1100)
for(rep_i in nreps){
  for(mcmc_rep_i in 1:n_mcmc_reps){
    sink(file = paste0("jobs_bailey_mcmc_sims_shortset.txt"), append = !(rep_i == 1 & mcmc_rep_i == 1))
    if(rep_i < 1001){
      cat(paste0("cd /Volumes/1TB/Bailey; /usr/local/bin/RScript scripts_sim/script_replicate_", rep_i, "_mcmc_rep_", mcmc_rep_i, ".R\n"))
    } else {
      cat(paste0("cd /Volumes/1TB/Bailey; /usr/local/bin/RScript scripts_sim/truevals_script_replicate_", rep_i - 1000, "_mcmc_rep_", mcmc_rep_i, ".R\n"))
    }
    sink()
  }
}
cat("parallel --jobs 8 --sshloginfile /Volumes/macOS/Users/nikolai/instances_noDog --eta --results bailey_sims_mcmc_shortset_screenout --files < jobs_bailey_mcmc_sims_shortset.txt")
