#Processing output of simulations with BONSAI

# devtools::build_vignettes()

library(bonsai)

# setwd("~/repos/bonsai/bonsai")


traitNumIter <- c(2,4,8,16,32,64)
traitNumIter <- traitNumIter[6:1]

degreeMisspecification <- c("perfect", "interregional", "chimpHuman")


# k iterates over degree of misspecification

counterVar <- 0

for (j in 1:6) {

for (k in 1:length(degreeMisspecification)) {
    
    for(i in 1:5){
      
      counterVar <- counterVar + 1
      print(paste(j, i, k, counterVar))
      
      ntraits <- traitNumIter[j]
      rateMatrixSpecification <- degreeMisspecification[k]
      replicateNum <- i
      
      fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_", "replicate_", replicateNum)
      
      
      # mvBM simulations
      output_dir  <- paste0("~/repos/bonsai/mvBM_Simulations/", fileName)
      input_files <- c(paste0("~/Google Drive/mvBM_sims/output/", fileName, "_params_run_1.log"),
                       paste0("~/Google Drive/mvBM_sims/output/", fileName, "_params_run_2.log"),
                       paste0("~/Google Drive/mvBM_sims/output/", fileName, "_trees_run_1.trees"),
                       paste0("~/Google Drive/mvBM_sims/output/", fileName, "_trees_run_2.trees"))
      name <- paste0("mvBM_sims_", fileName)
      
      project <- bonsai(name, output_dir, posterior_directories=input_files,
                        clade_posterior_threshold=0.1)
      
      # project <- bonsai(name, output_dir, posterior_directories=input_files,
      #                   ignore_list=c("Posterior","Likelihood","Prior","br_lens"),
      #                   clade_posterior_threshold=0.1)
      
      # compile the markdown
      project$makeRMarkdown()
      
      
    }
  }
}
# omnibus summaries
project$omnibus_diagnostics_posterior
project$omnibus_diagnostics_posterior_summary
project$plotCompareTrees()
project$posterior_compare_trees[2,1][[1]][[1]]




