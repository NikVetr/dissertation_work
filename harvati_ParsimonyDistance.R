setwd(dir = "/Volumes/1TB/Harvati/")

traitNumIter <- c(2^(0:4)*(46-1)+1)
nreps <- 100

# k iterates over degree of misspecification

for (j in 2:length(traitNumIter)) { #iterates over trait number
  
  for(i in 1:nreps){ #iterates over replicates
    
    replicateNum <- i
    ntraits <- traitNumIter[j]
    fileName <- paste0("simulations_", ntraits, "traits_", "replicate_", replicateNum)
    
    if(i == 1 & j == 1){
      file.remove("jobs_parsimony_distance_manytraits.txt")
    }
    
    # write the  scripts
    tempScript <- readLines("/Volumes/macOS/Users/nikolai/scripts/dissertation_work/metascript_parsimony_distance.R")
    tempScript[1] <- paste0("setwd(\"", getwd(), "/\") ")
    tempScript[2] <- paste0(tempScript[2], replicateNum)
    tempScript[3] <- paste0(tempScript[3], ntraits)
    writeLines(tempScript, paste0("scripts_parsimony_distance/", fileName, "_ParsimonyDistanceAnalysis.R"))
    
    sink(file = "jobs_parsimony_distance_manytraits.txt", append = T)
    cat(paste0("cd ", getwd(), "; ",  "/usr/local/bin/RScript ", "scripts_parsimony_distance/", fileName, "_ParsimonyDistanceAnalysis.R\n"))
    sink()
  
    
  }
}

#terminal command; can use system() but screen output is messy
cat(paste0("cd ", getwd(), "; parallel --jobs 7 --sshloginfile /Volumes/macOS/Users/nikolai/instances_all --eta --results terminal_output_parsimony_distance_manytraits --files < jobs_parsimony_distance_manytraits.txt"))

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


####################################
## Restarting Runs Failed by Crash #
####################################
traitNumIter <- c(2^(0:4)*(46-1)+1)
nreps <- 100

# k iterates over degree of misspecification

counter <- 0
for (j in 2:length(traitNumIter)) { #iterates over trait number
  
  for(i in 1:nreps){ #iterates over replicates
    
    counter <- counter + 1
    replicateNum <- i
    ntraits <- traitNumIter[j]
    fileName <- paste0("simulations_", ntraits, "traits_", "replicate_", replicateNum)
    
    if(counter == 1){
      file.remove("jobs_parsimony_distance_manytraits_restart.txt")
    }
    
    #check for completion
    notYetDone <- !file.exists(file = paste0("output_distance/", fileName, "_upgma.txt"))

      # write the  scripts
    if(notYetDone){
      print(paste0(counter, " ", fileName))
      tempScript <- readLines("/Volumes/macOS/Users/nikolai/scripts/dissertation_work/metascript_parsimony_distance.R")
      tempScript[1] <- paste0("setwd(\"", getwd(), "/\") ")
      tempScript[2] <- paste0(tempScript[2], replicateNum)
      tempScript[3] <- paste0(tempScript[3], ntraits)
      writeLines(tempScript, paste0("scripts_parsimony_distance/", fileName, "_ParsimonyDistanceAnalysis.R"))
    
      sink(file = "jobs_parsimony_distance_manytraits_restart.txt", append = T)
      cat(paste0("cd ", getwd(), "; ",  "/usr/local/bin/RScript ", "scripts_parsimony_distance/", fileName, "_ParsimonyDistanceAnalysis.R\n"))
      sink()
      
    }
    
  }
}

#terminal command; can use system() but screen output is messy
cat(paste0("cd ", getwd(), "; parallel --jobs 4 --sshloginfile /Volumes/macOS/Users/nikolai/instances_all --eta --results terminal_output_parsimony_distance_manytraits --files < jobs_parsimony_distance_manytraits_restart.txt"))
