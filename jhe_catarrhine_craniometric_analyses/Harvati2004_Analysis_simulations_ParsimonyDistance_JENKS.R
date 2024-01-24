setwd(dir = "/Volumes/1TB/Harvati/")

nreps <- 100
traitNumIter <- c((46-1)*2^(0:4)+1)
jenks_exp <- c(2,4,6)

counter <- 0
for(i in 1:nreps){
  
  cat(paste0("\n", i, ":"))
  
  for(j in 1:length(traitNumIter)){
    
    cat(paste0(" ", j, " "))
    
    ntraits <- traitNumIter[j]
    replicateNum <- i
    
    fileName <- paste0("simulations_", ntraits, "traits_", "replicate_", replicateNum)
    
    #get true tree just in case lol
    trueTree <- read.nexus(file = paste0("trees/", fileName, "_tree.nex"))
    tips <- trueTree$tip.label
    ntips <- length(tips)
    
    counter <- counter + 1
    if(counter == 1){
      file.remove("jobs_pars_jenks.txt")
    }
    
    for(transData in c("raw", "PCs")){
      
      for(jenks_exp_cats in 1:length(jenks_exp)){
        
        # write the  scripts
        tempScript <- readLines("/Volumes/macOS/Users/nikolai/scripts/dissertation_work/metascript_parsimony_distance_jenks.R", warn = F)
        tempScript[1] <- paste0("setwd(\"", getwd(), "/\") ")
        tempScript[2] <- paste0(tempScript[2], replicateNum)
        tempScript[3] <- paste0(tempScript[3], ntraits)
        tempScript[4] <- paste0(tempScript[4], jenks_exp[jenks_exp_cats])
        tempScript[5] <- paste0(tempScript[5], "\"", transData, "\"")
        writeLines(tempScript, paste0("scripts_parsimony_jenks/", fileName, "_", transData, "_jenks_", expCat, "_expcats_ParsimonyAnalysis.R"))
        
        sink(file = "jobs_pars_jenks.txt", append = T)
        cat(paste0("cd ", getwd(), "; ",  "/usr/local/bin/RScript ", "scripts_parsimony_jenks/", fileName, "_", transData, "_jenks_", expCat, "_expcats_ParsimonyAnalysis.R\n"))
        sink()
        
      }
    }
  }
}

cat(paste0("cd ", getwd(), "; parallel --jobs 8 --sshloginfile /Volumes/macOS/Users/nikolai/instances_justPigMonkey --eta --results terminal_output_parsimony_jenks --files < jobs_pars_jenks.txt"))
