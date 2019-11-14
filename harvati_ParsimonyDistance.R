setwd(dir = "/Volumes/1TB/Harvati/")

traitNumIter <- c(46)
nreps <- 500

# k iterates over degree of misspecification

for (j in 1:length(traitNumIter)) { #iterates over trait number
  
  for(i in 1:nreps){ #iterates over replicates
    
    replicateNum <- i
    ntraits <- traitNumIter[j]
    fileName <- paste0("simulations_", ntraits, "traits_", "replicate_", replicateNum)
    
    if(i == 1 & j == 1){
      file.remove("jobs_parsimony_distance.txt")
    }
    
    # write the  scripts
    tempScript <- readLines("metascript_parsimony_distance.R")
    tempScript[1] <- paste0("setwd(\"", getwd(), "/\") ")
    tempScript[2] <- paste0(tempScript[2], replicateNum)
    tempScript[3] <- paste0(tempScript[3], ntraits)
    writeLines(tempScript, paste0("scripts_parsimony_distance/", fileName, "_ParsimonyDistanceAnalysis.R"))
    
    sink(file = "jobs_parsimony_distance.txt", append = T)
    cat(paste0("cd ", getwd(), "; ",  "/usr/local/bin/RScript ", "scripts_parsimony_distance/", fileName, "_ParsimonyDistanceAnalysis.R\n"))
    sink()
  
    
  }
}

#terminal command; can use system() but screen output is messy
cat(paste0("cd ", getwd(), "; parallel --jobs 7 --sshloginfile /Volumes/macOS/Users/nikolai/instances_all --eta --results terminal_output_parsimony_distance --files < jobs_parsimony_distance.txt"))

#analyses that first do a PC transform on continuous character data
for (j in 1:length(traitNumIter)) { #iterates over trait number
  
  for(i in 1:nreps){ #iterates over replicates
    
    replicateNum <- i
    ntraits <- traitNumIter[j]
    fileName <- paste0("simulations_", ntraits, "traits_", "replicate_", replicateNum)
    
    if(i == 1 & j == 1){
      file.remove("jobs_parsimony_PCA.txt")
    }
    
    # write the  scripts
    tempScript <- readLines("metascript_parsimony_PCA.R")
    tempScript[1] <- paste0("setwd(\"", getwd(), "/\") ")
    tempScript[2] <- paste0(tempScript[2], replicateNum)
    tempScript[3] <- paste0(tempScript[3], ntraits)
    writeLines(tempScript, paste0("scripts_parsimony_distance/", fileName, "_ParsimonyDistanceAnalysis_PCA.R"))
    
    sink(file = "jobs_parsimony_PCA.txt", append = T)
    cat(paste0("cd ", getwd(), "; ",  "/usr/local/bin/RScript ", "scripts_parsimony_distance/", fileName, "_ParsimonyDistanceAnalysis_PCA.R\n"))
    sink()
    
    
  }
}

#terminal command; can use system() but screen output is messy
cat(paste0("cd ", getwd(), "; parallel --jobs 7 --sshloginfile /Volumes/macOS/Users/nikolai/instances_all --eta --results terminal_output_parsimony_distance_PCA --files < jobs_parsimony_PCA.txt"))
