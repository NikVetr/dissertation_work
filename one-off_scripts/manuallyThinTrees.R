#script to go through MCMC_output and manually thin trees


traitNumIter <- c(2,4,8,16,32,64)

degreeMisspecification <- c("perfect", "interregional", "chimpHuman")


# k iterates over degree of misspecification
for (k in 1:length(degreeMisspecification)) {
  
  for (j in 1:6) {
    
    for(i in 1:5){
      
      cat(c(j, i, k, "... "))
      
      ntraits <- traitNumIter[j]
      rateMatrixSpecification <- degreeMisspecification[k]
      replicateNum <- i
      
      fileName <- paste0("simulations_", ntraits, "traits_", rateMatrixSpecification, "RateMatrix_", "replicate_", replicateNum)
      
      # continuous character diagnostics
      log1 <- read.table(file = paste0("A:\\mvBM_sims\\output\\", fileName, "_params_run_1.log"), header = T)
      #log1 <- log1[2501:12500,] #discard burn-in
      log1 <- log1[seq(0, 12500, by = 5),] #thin further
      write.table(x = log1, sep = "\t", row.names = F, append = FALSE, quote = FALSE, file = paste0("A:\\mvBM_sims\\output_thinned\\", fileName, "_params_run_1.log"))
      
      log2 <- read.table(file = paste0("A:\\mvBM_sims\\output\\", fileName, "_params_run_2.log"), header = T)
      #log2 <- log2[2501:12500,] #discard burn-in
      log2 <- log2[seq(0, 12500, by = 5),] #thin further
      write.table(x = log2, sep = "\t", row.names = F, append = FALSE, quote = FALSE, file = paste0("A:\\mvBM_sims\\output_thinned\\", fileName, "_params_run_2.log"))
      

      trees1 <- read.table(file = paste0("A:\\mvBM_sims\\output\\", fileName, "_trees_run_1.trees"))
      trees1_thinned <- trees1[seq(from = 1, to = length(trees1[,1]), by = 5),]
      trees1_thinned <- as.matrix(trees1_thinned)
      sink(paste0("A:\\mvBM_sims\\output_thinned\\", fileName, "_trees_run_1.trees"), append = F)
      for(i in 1:length(trees1_thinned[,1])){
        cat(trees1_thinned[i,], sep = "\t")
        cat("\n")
      }
      sink()
      
      trees2 <- read.table(file = paste0("A:\\mvBM_sims\\output\\", fileName, "_trees_run_2.trees"))
      trees2_thinned <- trees2[seq(from = 1, to = length(trees1[,1]), by = 5),]
      trees2_thinned <- as.matrix(trees2_thinned)
      sink(paste0("A:\\mvBM_sims\\output_thinned\\", fileName, "_trees_run_2.trees"), append = F)
      for(i in 1:length(trees2_thinned[,1])){
        cat(trees2_thinned[i,], sep = "\t")
        cat("\n")
      }
      sink()
      
    }
  }
}




trees1 <- read.table(file = "C:\\Users\\Nikolai\\Google Drive\\outputGD\\howellsLogUnifMC3_3c_run2.trees")

trees1_thinned <- trees1[seq(from = 1, to = length(trees1[,1]), by = 40),]

trees1_thinned <- as.matrix(trees1_thinned)

sink("C:\\Users\\Nikolai\\Google Drive\\outputGD\\howellsLogUnifMC3_3c_run2_thinned.trees", append = F)
for(i in 1:length(trees1_thinned[,1])){
  cat(trees1_thinned[i,], sep = "\t")
  cat("\n")
}
sink()

#################

trees1 <- read.table(file = "C:\\Users\\Nikolai\\Google Drive\\outputGD\\howellsLogUnifMC3_3c_run3.trees")

trees1_thinned <- trees1[seq(from = 1, to = length(trees1[,1]), by = 40),]

trees1_thinned <- as.matrix(trees1_thinned)

sink("C:\\Users\\Nikolai\\Google Drive\\outputGD\\howellsLogUnifMC3_3c_run3_thinned.trees", append = F)
for(i in 1:length(trees1_thinned[,1])){
  cat(trees1_thinned[i,], sep = "\t")
  cat("\n")
}
sink()

#################

trees1 <- read.table(file = "C:\\Users\\Nikolai\\Google Drive\\outputGD\\howellsDirichMC3_3c_run4.trees")

trees1_thinned <- trees1[seq(from = 1, to = length(trees1[,1]), by = 40),]

trees1_thinned <- as.matrix(trees1_thinned)

sink("C:\\Users\\Nikolai\\Google Drive\\outputGD\\howellsDirichMC3_3c_run4_thinned.trees", append = F)
for(i in 1:length(trees1_thinned[,1])){
  cat(trees1_thinned[i,], sep = "\t")
  cat("\n")
}
sink()

#################

trees1 <- read.table(file = "C:\\Users\\Nikolai\\Google Drive\\outputGD\\howellsDirichMC3_3c_run5.trees")

trees1_thinned <- trees1[seq(from = 1, to = length(trees1[,1]), by = 40),]

trees1_thinned <- as.matrix(trees1_thinned)

sink("C:\\Users\\Nikolai\\Google Drive\\outputGD\\howellsDirichMC3_3c_run5_thinned.trees", append = F)
for(i in 1:length(trees1_thinned[,1])){
  cat(trees1_thinned[i,], sep = "\t")
  cat("\n")
}
sink()

