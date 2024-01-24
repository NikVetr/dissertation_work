#Processing output of simulations with BONSAI

# devtools::build_vignettes()

library(bonsai)

# setwd("~/repos/bonsai/bonsai")

fileName <- "howellsLogUnifMCMC_"


# mvBM simulations
output_dir  <- paste0("~/repos/bonsai/howellsLogUnifMCMCvMC3/", fileName)
input_files <- c(paste0("~/output/", fileName, "1.log"),
                 paste0("~/output/", fileName, "2.log"),
                 paste0("~/output/", fileName, "1.trees"),
                 paste0("~/output/", fileName, "2.trees"))

input_files <- c(paste0("~/output/", fileName, "1.log"),
                 paste0("~/output/", fileName, "2.log"),
                 paste0("~/output/", fileName, "1.trees"),
                 paste0("~/output/", fileName, "2.trees"),
                 paste0("~/output/howellsLogUnifMC3_3c_Sims1.trees"),
                 paste0("~/output/howellsLogUnifMC3_3c_Sims2.trees"),
                 paste0("~/output/howellsLogUnifMC3_3c_Sims1.log"),
                 paste0("~/output/howellsLogUnifMC3_3c_Sims2.log")
                 )



name <- paste0("howellsEmpiricalAnalysis_", fileName)

project <- bonsai(name, output_dir, posterior_directories=input_files,
                  clade_posterior_threshold=0.1)

# project <- bonsai(name, output_dir, posterior_directories=input_files,
#                   ignore_list=c("Posterior","Likelihood","Prior","br_lens"),
#                   clade_posterior_threshold=0.1)

# compile the markdown
project$makeRMarkdown()


# omnibus summaries
project$omnibus_diagnostics_posterior
project$omnibus_diagnostics_posterior_summary
project$plotCompareTrees()
project$posterior_compare_trees[2,1][[1]][[1]]




