setwd(dir = "/Volumes/1TB/500/full/uvBM_sims_PB_noiseless_500/") 
noisyness <- strsplit(getwd(), "_")[[1]][startsWith(strsplit(getwd(), "_")[[1]], "nois")]
traitNumIter <- c(2,4,8,16,32,64,128)
degreeCorrelation <- c("2", "4", "6", "8", "9")

for(i in 501:1500){
  for (k in 1:length(degreeCorrelation)) {
    for (j in 1:length(traitNumIter)) {
    
      ntraits <- traitNumIter[j]
      rateMatrixCorrelation <- degreeCorrelation[k]
      replicateNum <- i
      
      fileName <- paste0("simulations_", ntraits, "traits_perfectRateMatrix_", rateMatrixCorrelation, "correlationStrength_", noisyness, "_replicate_", replicateNum)

      #check if I named the file correctly
      if(!file.exists(paste0("output/", fileName, "_params_run_1.log"))){
        if(!file.exists(paste0("output/", fileName, "_trees_run_1.trees"))){
          oldFileName <- paste0("simulations_", ntraits, "traits_perfectRateMatrix_", rateMatrixCorrelation, "correlationStrength_replicate_", replicateNum)
          if(file.exists(paste0("output/", oldFileName, "_params_run_1.log"))){
            cat(paste0("changing name of file ", i, "\n"))
            file.rename(from = paste0("output/", oldFileName, "_params_run_1.log"), to = paste0("output/", fileName, "_params_run_1.log"))
            file.rename(from = paste0("output/", oldFileName, "_params_run_2.log"), to = paste0("output/", fileName, "_params_run_2.log"))
            file.rename(from = paste0("output/", oldFileName, "_trees_run_1.trees"), to = paste0("output/", fileName, "_trees_run_1.trees"))
            file.rename(from = paste0("output/", oldFileName, "_trees_run_2.trees"), to = paste0("output/", fileName, "_trees_run_2.trees"))
          }
        }
      }
    }
  }
}    

#replace awkwardly restarted runs
for(i in 501:1500){
  for (k in 1:length(degreeCorrelation)) {
    for (j in 1:length(traitNumIter)) {
      
      ntraits <- traitNumIter[j]
      rateMatrixCorrelation <- degreeCorrelation[k]
      replicateNum <- i
      
      fileName <- paste0("simulations_", ntraits, "traits_perfectRateMatrix_", rateMatrixCorrelation, "correlationStrength_", noisyness, "_replicate_", replicateNum)
      oldFileName <- paste0("simulations_", ntraits, "traits_perfectRateMatrix_", rateMatrixCorrelation, "correlationStrength_replicate_", replicateNum)

      if(file.exists(paste0("output/", oldFileName, "_params_run_1.log"))){
        cat(paste0("changing name of file ", i, "\n"))
        file.rename(from = paste0("output/", oldFileName, "_params_run_1.log"), to = paste0("output/", fileName, "_params_run_1.log"))
        file.rename(from = paste0("output/", oldFileName, "_params_run_2.log"), to = paste0("output/", fileName, "_params_run_2.log"))
        file.rename(from = paste0("output/", oldFileName, "_trees_run_1.trees"), to = paste0("output/", fileName, "_trees_run_1.trees"))
        file.rename(from = paste0("output/", oldFileName, "_trees_run_2.trees"), to = paste0("output/", fileName, "_trees_run_2.trees"))
      }
    }
  }
}    

