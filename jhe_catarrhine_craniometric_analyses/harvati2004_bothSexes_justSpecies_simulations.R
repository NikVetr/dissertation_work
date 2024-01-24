setwd("/Volumes/1TB/Harvati")

library(ape)
library(coda)
library(rwty)
library(phangorn)
library(lattice)
library(latticeExtra)
library(MCMCpack)
library(mvMORPH)
library(phytools)
library(RevGadgets)

meansNexus <- function (characterMatrix, outputFilePath, traitsOnColumns = T){
  if (traitsOnColumns) {ntax = dim(characterMatrix)[1]; nchar = dim(characterMatrix)[2]}
  else {ntax = dim(characterMatrix)[2]; nchar = dim(characterMatrix)[1]}
  sink(outputFilePath, append=F)
  cat("#NEXUS\n\nBegin data;\n\tDimensions ntax=", ntax, "nchar=", nchar, ";\n\tFormat datatype=Continuous missing=?;\n\tMatrix\n\n")
  for(i in 1:length(characterMatrix[,1])) {cat(rownames(characterMatrix)[i], characterMatrix[i,], "\n")}
  cat(";\nEnd;")
  sink()
}

trees_m <- read.tree("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_males_noHomopops_justSpecies_uninf_mvBM_fixTreeInfRM_aggrProps_c1.trees")
trees_f <- read.tree("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_females_noHomopops_justSpecies_uninf_mvBM_fixTreeInfRM_aggrProps_c1.trees")

params_m_1 <- read.table("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_males_noHomopops_justSpecies_uninf_mvBM_fixTreeInfRM_aggrProps_c1.log", header = T)
params_m_2 <- read.table("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_males_noHomopops_justSpecies_uninf_mvBM_fixTreeInfRM_aggrProps_c2.log", header = T)
params_m_3 <- read.table("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_males_noHomopops_justSpecies_uninf_mvBM_mixTreeInfRM_c1.log", header = T)
params_m_4 <- read.table("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_males_noHomopops_justSpecies_uninf_mvBM_mixTreeInfRM_c2.log", header = T)
params_m_all <- rbind(params_m_1, params_m_2, params_m_3, params_m_4)
# sort(coda::effectiveSize(params_m_all))

params_f_1 <- read.table("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_females_noHomopops_justSpecies_uninf_mvBM_fixTreeInfRM_aggrProps_c1.log", header = T)
params_f_2 <- read.table("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_females_noHomopops_justSpecies_uninf_mvBM_fixTreeInfRM_aggrProps_c2.log", header = T)
params_f_3 <- read.table("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_females_noHomopops_justSpecies_uninf_mvBM_fixTreeInfRM_c1.log", header = T)
params_f_4 <- read.table("/Volumes/1TB/Harvati_Empirical/output/justSpecies/harvati_PCA_females_noHomopops_justSpecies_uninf_mvBM_fixTreeInfRM_c2.log", header = T)
params_f_all <- rbind(params_f_1, params_f_2, params_f_3, params_f_4)
# sort(coda::effectiveSize(params_f_all))

recomposeRateMatrix <- function(params, sds_keyword = "relative_sd", corrs_keyword = "correlation", ntraits = NA){
  #get the standard deviations
  sds <- params[grep(sds_keyword, colnames(params))]
  if(!all(diff(as.numeric(sapply(1:length(sds), function(x) strsplit(colnames(sds)[[x]], split = "\\.")[[1]][2]))) == 1)){
    stop("error something has gone terrible wrong")
  }
  sds <- diag(as.numeric(sds))
  
  if(is.na(ntraits)){
    ntraits <- nrow(sds)
  }
  
  #get the correlations
  corrmat <- diag(ntraits)
  upper_tri_corrs <- params[grep(corrs_keyword, colnames(params))]
  if(!all(diff(as.numeric(sapply(1:length(upper_tri_corrs), function(x) strsplit(colnames(upper_tri_corrs)[[x]], split = "\\.")[[1]][2]))) == 1)){
    stop("error something done goofed")
  }
  corrmat[lower.tri(corrmat)] <- as.numeric(upper_tri_corrs)
  corrmat <- corrmat + t(corrmat) - diag(ntraits)
  
  #recompose and return
  mvBM_R <- sds %*% corrmat %*% sds
  return(mvBM_R)
}

matchNodes = function(phy){
  
  # get some useful info
  num_tips = length(phy$tip.label)
  num_nodes = phy$Nnode
  tip_indexes = 1:num_tips
  node_indexes = num_tips + num_nodes:1
  
  node_map = data.frame(R=1:(num_tips + num_nodes), Rev=NA, visits=0)
  current_node = phy$Nnode + 2
  k = 1
  t = 1
  
  while(TRUE) {
    
    if ( current_node <= num_tips ) {
      node_map$Rev[node_map$R == current_node] = t
      current_node = phy$edge[phy$edge[,2] == current_node,1]
      t = t + 1
    } else {
      
      if ( node_map$visits[node_map$R == current_node] == 0 ) {
        node_map$Rev[node_map$R == current_node] = node_indexes[k]
        k = k + 1
      }
      node_map$visits[node_map$R == current_node] = node_map$visits[node_map$R == current_node] + 1
      
      if ( node_map$visits[node_map$R == current_node] == 1 ) {
        # go right
        current_node = phy$edge[phy$edge[,1] == current_node,2][2]
      } else if ( node_map$visits[node_map$R == current_node] == 2 ) {
        # go left
        current_node = phy$edge[phy$edge[,1] == current_node,2][1]
      } else if ( node_map$visits[node_map$R == current_node] == 3 ) {
        # go down
        if (current_node == num_tips + 1) {
          break
        } else {
          current_node = phy$edge[phy$edge[,2] == current_node,1]
        }
      }
    }
    
  }
  
  return(node_map[,1:2])
  
}
# trees_f[[1]]$edge.length[order(trees_f[[1]]$edge[,2])][order(matchNodes(trees_f[[1]])[,2])]

getBranchRates <- function(params, tree, keyword = "branch_rates_"){
  br <- params[grep(keyword, colnames(params))]
  if(!all(diff(as.numeric(sapply(1:length(br), function(x) strsplit(colnames(br)[[x]], split = "\\.")[[1]][2]))) == 1)){
    stop("error something has gone terrible wrong")
  }
  br <- as.numeric(br)
  
  #get table of tree edges
  tr_edge <- cbind(edge_ind = 1:length(tree$edge.length), tree$edge)
  tr_edge <- rbind(tr_edge, c(NA, NA, length(tree$tip.label)+1))
  tr_edge <- tr_edge[order(tr_edge[,3]),] #reorder by tip index
  edge_indices_of_desc_tip <- tr_edge[,1]
    
  #construct table of relations between indices
  matcher <- cbind(matchNodes(tree), R_edge_ind = edge_indices_of_desc_tip)
  #reorder to privelege rb
  matcher <- matcher[order(matcher$Rev),]
  #get br's in there
  matcher <- cbind(matcher, rb_br = c(br, NA))
  #reorder by phylo's edge lengths
  matcher <- matcher[order(matcher$R_edge_ind),]
  reordered_brs <- matcher$rb_br
  if(is.na(reordered_brs[length(reordered_brs)])){
    reordered_brs <- reordered_brs[-length(reordered_brs)]
  } else {
    stop("something's afoot with the indexing here")
  }
  return(reordered_brs)
}

n_rep <- 100
nrun <- 2
ntraits_m <- 20
ntraits_f <- 19
perc_var_m_all <- rep(0, n_rep)
perc_var_f_all <- rep(0, n_rep)
for(i in 1:n_rep){
  cat(paste0(i, " "))
  fileName <- paste0("simulations_99PCA_bothSexes_justSpecies_replicate_", i)
  
  m_sample_index <- sample(size = 1, 1:length(nrow(params_m_all)))
  f_sample_index <- sample(size = 1, 1:length(nrow(params_f_all)))
  
  params_m <- params_m_all[m_sample_index,]
  params_f <- params_f_all[f_sample_index,]
  
  tree <- trees_m[[1]]
  tree$edge.length <- tree$edge.length * getBranchRates(params_m, tree)
  
  rateMatrix_m <- recomposeRateMatrix(params_m)
  rateMatrix_f <- recomposeRateMatrix(params_f)
  rateMatrix_mf <- (rateMatrix_m + rateMatrix_f) / 2
  rateMatrix <- kronecker(matrix(c(1,0.9,0.9,1), 2, 2), rateMatrix_mf)
  
  traits <- mvSIM(tree = tree, nsim = 1, model = "BM1", param = list(ntraits = dim(rateMatrix)[1], sigma=rateMatrix, mu= rep(0, times = dim(rateMatrix)[1])))
  traits_m <- traits[,1:(ncol(traits)/2)]
  traits_f <- traits[,(ncol(traits)/2+1):ncol(traits)]
  
  es_m <- eigen(rateMatrix_m)
  es_f <- eigen(rateMatrix_f)
  
  perc_var_m_all[i] <- perc_var_m <- (cumsum(es_m$values) / sum(es_m$values))[ntraits_m]
  perc_var_f_all[i] <- perc_var_f <- (cumsum(es_f$values) / sum(es_f$values))[ntraits_f]
  
  PCs_m <- t((t(es_m$vectors) %*% t(traits_m)))
  PCs_m <-sapply(1:length(es_m$values), function(ev) PCs_m[,ev] / es_m$values[ev]^0.5)
  PCs_m <- PCs_m[,1:ntraits_m]
  
  PCs_f <- t((t(es_f$vectors) %*% t(traits_f)))
  PCs_f <-sapply(1:length(es_f$values), function(ev) PCs_f[,ev] / es_f$values[ev]^0.5)
  PCs_f <- PCs_f[,1:ntraits_f]
  
  traitsPath_m <- paste0("data/bothSexes_justSpecies/", fileName, "_m_traits.nex")
  traitsPath_f <- paste0("data/bothSexes_justSpecies/", fileName, "_f_traits.nex")
  
  meansNexus(PCs_f, traitsPath_f)
  meansNexus(PCs_m, traitsPath_m)
  
  if(i == 1){
    file.remove("jobs_justSpecies_bothSexes_sims.txt")
  }
  
  # write the master scripts
  for (l in 1:nrun) {
    
    tempScript <- readLines("harvati2004_infPrior_mvBM_metascript.txt")
    tempScript[1] <- paste0("setwd(\"", getwd(), "/\") ")
    tempScript[2] <- paste0(tempScript[2], l, ")")
    tempScript[3] <- paste0(tempScript[3], l)
    tempScript[4] <- paste0("replicateNum <- ", i)
    tempScript[8] <- gsub(x = tempScript[8], pattern = "placeholder_f.nex", replacement = traitsPath_f)
    tempScript[9] <- gsub(x = tempScript[9], pattern = "placeholder_m.nex", replacement = traitsPath_m)
    writeLines(tempScript, paste0("scripts/bothSexes_justSpecies/", fileName, "_script_run_", l, ".Rev"))
    
    sink(file = "jobs_justSpecies_bothSexes_sims.txt", append = T)
    cat(paste0("cd ", getwd(), "; ",  "/Users/nikolai/repos/revbayes-development/projects/cmake/rb scripts/bothSexes_justSpecies/", fileName, "_script_run_", l, ".Rev\n"))
    sink()
  }
  
}
    
#rewrite rb scripts to make br_lens visible

for(i in 1:n_rep){
  
  cat(paste0(i, " "))
  fileName <- paste0("simulations_99PCA_bothSexes_justSpecies_replicate_", i)
  traitsPath_m <- paste0("data/bothSexes_justSpecies/", fileName, "_m_traits.nex")
  traitsPath_f <- paste0("data/bothSexes_justSpecies/", fileName, "_f_traits.nex")
  
  if(i == 1){
    file.remove("jobs_justSpecies_bothSexes_sims_visibleBLs.txt")
  }
  
  # write the master scripts
  for (l in 1:nrun) {
    
    tempScript <- readLines("harvati2004_infPrior_mvBM_metascript_visibleBLs.txt")
    tempScript[1] <- paste0("setwd(\"", getwd(), "/\") ")
    tempScript[2] <- paste0(tempScript[2], l, ")")
    tempScript[3] <- paste0(tempScript[3], l)
    tempScript[4] <- paste0("replicateNum <- ", i)
    tempScript[8] <- gsub(x = tempScript[8], pattern = "placeholder_f.nex", replacement = traitsPath_f)
    tempScript[9] <- gsub(x = tempScript[9], pattern = "placeholder_m.nex", replacement = traitsPath_m)
    writeLines(tempScript, paste0("scripts/bothSexes_justSpecies/", fileName, "_script_visibleBLs_run_", l, ".Rev"))
    
    sink(file = "jobs_justSpecies_bothSexes_sims_visibleBLs.txt", append = T)
    cat(paste0("cd ", getwd(), "; ",  "/Users/nikolai/repos/revbayes-development/projects/cmake/rb scripts/bothSexes_justSpecies/", fileName, "_script_visibleBLs_run_", l, ".Rev\n"))
    sink()
  }
}

#terminal command; can use system() but screen output is messy
cat(paste0("cd ", getwd(), "; parallel --jobs 8 --sshloginfile /Volumes/macOS/Users/nikolai/instances_noDog --eta --results terminal_output_bSjS --files < jobs_justSpecies_bothSexes_sims.txt"))

#for visible BLs runs
cat(paste0("cd ", getwd(), "; parallel --jobs 8 --sshloginfile /Volumes/macOS/Users/nikolai/instances_noDog_noMonkey --eta --results terminal_output_bSj_visBLsS --files < jobs_justSpecies_bothSexes_sims_visibleBLs.txt"))


#MCMC diagnostics
metaScript <- readLines("diagTemp.R", warn = F)
dir.create("diagnosticScripts")
nrep <- 500
for(i in 1:nrep){
  tempScript <- metaScript
  tempScript[1] <- paste0("setwd(dir = \"", getwd(), "/\") ")
  tempScript[2] <- paste0(tempScript[2], ntraits)
  tempScript[3] <- paste0(tempScript[3], i)
  writeLines(tempScript, paste0("diagnosticScripts/", i, ".R"))
}

sink(file = "diagnosticJobs.txt", append = F)

for(i in 1:nrep){
  cat(paste0("cd ", getwd(), "; ",  "/usr/local/bin/RScript ", "diagnosticScripts/", i, ".R\n"))
}   

sink()

#terminal command; can use system() but screen output is messy
cat(paste0("cd ", getwd(), "; parallel --jobs 7 --sshloginfile /Volumes/macOS/Users/nikolai/instances_all --eta --results terminal_output_diags --files < diagnosticJobs.txt"))


# instantiate matrices to read diagnostic output back in
traitNumIter <- c(46)

diagnostics <- matrix(ncol = 5 + 376, nrow = nrep*14)
colnames(diagnostics) <- c("diagnostic", "chain", "ntraits", "replicateNum","Posterior","Likelihood","TL",
                           "collessIndex","pair1","pair2","pair3","pair4","pair5","pair6","pair7","pair8","pair9","pair10","pair11",
                           "pair12","pair13","pair14","pair15","pair16","pair17","pair18","pair19","pair20","pair21","pair22","pair23",
                           "pair24","pair25","pair26","pair27","pair28","pair29","pair30","pair31","pair32","pair33","pair34","pair35",
                           "pair36","pair37","pair38","pair39","pair40","pair41","pair42","pair43","pair44","pair45","pair46","pair47",
                           "pair48","pair49","pair50","pair51","pair52","pair53","pair54","pair55","pair56","pair57","pair58","pair59",
                           "pair60","pair61","pair62","pair63","pair64","pair65","pair66","pair67","pair68","pair69","pair70","pair71",
                           "pair72","pair73","pair74","pair75","pair76","pair77","pair78","pair79","pair80","pair81","pair82","pair83",
                           "pair84","pair85","pair86","pair87","pair88","pair89","pair90","pair91","pair92","pair93","pair94","pair95",
                           "pair96","pair97","pair98","pair99","pair100","pair101","pair102","pair103","pair104","pair105","pair106",
                           "pair107","pair108","pair109","pair110","pair111","pair112","pair113","pair114","pair115","pair116","pair117",
                           "pair118","pair119","pair120","pair121","pair122","pair123","pair124","pair125","pair126","pair127","pair128",
                           "pair129","pair130","pair131","pair132","pair133","pair134","pair135","pair136","pair137","pair138","pair139",
                           "pair140","pair141","pair142","pair143","pair144","pair145","pair146","pair147","pair148","pair149","pair150",
                           "pair151","pair152","pair153","pair154","pair155","pair156","pair157","pair158","pair159","pair160","pair161",
                           "pair162","pair163","pair164","pair165","pair166","pair167","pair168","pair169","pair170","pair171","pair172",
                           "pair173","pair174","pair175","pair176","pair177","pair178","pair179","pair180","pair181","pair182","pair183",
                           "pair184","pair185","pair186","pair187","pair188","pair189","pair190","pair191","pair192","pair193","pair194",
                           "pair195","pair196","pair197","pair198","pair199","pair200","pair201","pair202","pair203","pair204","pair205",
                           "pair206","pair207","pair208","pair209","pair210","pair211","pair212","pair213","pair214","pair215","pair216",
                           "pair217","pair218","pair219","pair220","pair221","pair222","pair223","pair224","pair225","pair226","pair227",
                           "pair228","pair229","pair230","pair231","pair232","pair233","pair234","pair235","pair236","pair237","pair238",
                           "pair239","pair240","pair241","pair242","pair243","pair244","pair245","pair246","pair247","pair248","pair249",
                           "pair250","pair251","pair252","pair253","pair254","pair255","pair256","pair257","pair258","pair259","pair260",
                           "pair261","pair262","pair263","pair264","pair265","pair266","pair267","pair268","pair269","pair270","pair271",
                           "pair272","pair273","pair274","pair275","pair276","pair277","pair278","pair279","pair280","pair281","pair282",
                           "pair283","pair284","pair285","pair286","pair287","pair288","pair289","pair290","pair291","pair292","pair293",
                           "pair294","pair295","pair296","pair297","pair298","pair299","pair300","pair301","pair302","pair303","pair304",
                           "pair305","pair306","pair307","pair308","pair309","pair310","pair311","pair312","pair313","pair314","pair315",
                           "pair316","pair317","pair318","pair319","pair320","pair321","pair322","pair323","pair324","pair325","pair326",
                           "pair327","pair328","pair329","pair330","pair331","pair332","pair333","pair334","pair335","pair336","pair337",
                           "pair338","pair339","pair340","pair341","pair342","pair343","pair344","pair345","pair346","pair347","pair348",
                           "pair349","pair350","pair351","rfDist1",
                           "kfDist1","clade1","clade2","clade3","clade4","clade5","clade6","clade7","clade8","clade9","clade10","clade11",
                           "clade12","clade13","clade14","clade15","clade16","clade17","clade18","clade19","clade20")

comptrees <- matrix(ncol = 4, nrow = nrep)
colnames(comptrees) <- c("ntraits", "replicateNum", "correlation", "p.value")

counter <- 0
for(i in 1:nrep){
  print(i)    
  
  counter <- counter + 1
  
  ntraits <- 46
  replicateNum <- i
  
  fileName <- paste0("simulations_", ntraits, "traits_", "replicate_", replicateNum)
  
  load(file = paste0("diagnostics/geweke/", fileName, "_c1.txt"))
  load(file = paste0("diagnostics/geweke/", fileName, "_c2.txt"))
  load(file = paste0("diagnostics/geweke/", fileName, "_cb.txt"))
  
  load(file = paste0("diagnostics/ess/", fileName, "_c1.txt"))
  load(file = paste0("diagnostics/ess/", fileName, "_c2.txt"))
  load(file = paste0("diagnostics/ess/", fileName, "_cb.txt"))
  
  load(file = paste0("diagnostics/heiwel/", fileName, "_c1.txt"))
  load(file = paste0("diagnostics/heiwel/", fileName, "_c2.txt"))
  load(file = paste0("diagnostics/heiwel/", fileName, "_cb.txt"))
  
  if(is.vector(heiwel)){heiwel <- rbind(heiwel, heiwel)}; if(is.vector(heiwel1)){heiwel1 <- rbind(heiwel1, heiwel1)}; if(is.vector(heiwel2)){heiwel2 <- rbind(heiwel2, heiwel2)}
  
  load(file = paste0("diagnostics/psrf/", fileName, "_cb.txt"))    
  if(is.vector(psrf)){psrf <- rbind(psrf, psrf); print("changing...")}
  
  diagnostics[(((counter-1)*14+1):(counter*14)),] <- rbind(
    c(deparse(substitute(geweke1)),1,ntraits,replicateNum, geweke1[1:377]),
    c(deparse(substitute(ess1)),1,ntraits,replicateNum, ess1[1:377]),
    c(paste0(deparse(substitute(heiwel1)), "_CVM.p.val"),1,ntraits,replicateNum, heiwel1[1,1:377]),
    c(paste0(deparse(substitute(heiwel1)), "_halfwidth.passed"),1,ntraits,replicateNum, heiwel1[2,1:377]),
    
    c(deparse(substitute(geweke2)),2,ntraits,replicateNum, geweke2[1:377]),
    c(deparse(substitute(ess2)),2,ntraits,replicateNum, ess2[1:377]),
    c(paste0(deparse(substitute(heiwel2)), "_CVM.p.val"),2,ntraits,replicateNum, heiwel2[1,1:377]),
    c(paste0(deparse(substitute(heiwel2)), "_halfwidth.passed"),2,ntraits,replicateNum, heiwel2[2,1:377]),
    
    c(deparse(substitute(psrf_point)),"both",ntraits,replicateNum, psrf[1,1:377]),
    c(deparse(substitute(psrf_upper)),"both",ntraits,replicateNum, psrf[2,1:377]),
    c(deparse(substitute(geweke)),"both",ntraits,replicateNum, geweke[1:377]),
    c(deparse(substitute(ess)),"both",ntraits,replicateNum, ess[1:377]),
    c(paste0(deparse(substitute(heiwel)), "_CVM.p.val"),"both",ntraits,replicateNum, heiwel[1,1:377]),
    c(paste0(deparse(substitute(heiwel)), "_halfwidth.passed"),"both",ntraits,replicateNum, heiwel[2,1:377])
  )
  
  comptrees[counter,] <- unlist(read.table(file = paste0("diagnostics/compareTree/comptrees_for_", fileName, ".txt")))
}

diagnostics <- as.data.frame(diagnostics)
diagnostics[,6:381] <- sapply(6:381, function(x) as.numeric(as.character(diagnostics[,x])))
save(diagnostics, file = "diagnostics.txt")
save(comptrees, file = "CompareTreesCorrsPs")


#loading diagnostics back in

print(paste0("# of original failures: ", length(read.delim("job_MCMC_Fail.txt", header = F)[[1]])/2))
load(file = "diagnostics.txt")
load(file = "CompareTreesCorrsPs")

#look at ESS and R2
essBothChains <- diagnostics[diagnostics$diagnostic == "ess",]
thresholdESS <- 1000
thresholdR2 <- 0.95
failures <- c(which(as.numeric(as.character(essBothChains$Posterior)) < thresholdESS), 
              which(as.numeric(as.character(essBothChains$Likelihood)) < thresholdESS),
              which(as.numeric(as.character(essBothChains$TL)) < thresholdESS),
              which(as.numeric(as.character(essBothChains$rfDist1)) < thresholdESS),
              which(as.numeric(as.character(essBothChains$kfDist1)) < thresholdESS),
              which(as.numeric(as.character(essBothChains$kfDist1)) < thresholdESS),
              unique(unlist(sapply(1:length(which(startsWith(colnames(essBothChains), "pair"))), 
                                   function(x) which(as.numeric(as.character(essBothChains[,which(startsWith(colnames(essBothChains), "pair"))[x]])) < thresholdESS)))),
              which(comptrees[,3]^2 < thresholdR2))
failures <- sort(unique(failures))
failures <- essBothChains[failures,3:5]
failures <- failures[2*(1:(length(failures[,1])/2)),]
summary(failures$rateMatrixCorrelation)
