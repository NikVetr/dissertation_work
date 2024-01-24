setwd(dir = "/Volumes/2TB/mvBM_sims_PB_extra/output_pb_noisy_mvn_128")
list.files()
file.rename(from="simulations_128traits_perfectRateMatrix_2correlationStrength_replicate_1_params_run_1.log",to="simulations_128traits_perfectRateMatrix_2correlationStrength_replicate_1_params_run_1_test.log")

setwd(dir = "/Volumes/2TB/mvBM_sims_PB_extra/output_pb_noisy_mvn_128")
fileNames <- list.files()
sapply(1:length(fileNames), function(x) file.rename(from = fileNames[x], to = sub(pattern="correlationStrength",replacement="correlationStrength_mvBM_noisy", fileNames[x])))
sapply(1:length(fileNames), function(x) file.rename(from = fileNames[x], to = sub(pattern="correlationStrength_mvBM_noisy_mvBM_noisy",replacement="correlationStrength_noisy", fileNames[x])))

setwd(dir = "/Volumes/2TB/uvBM_sims_PB_noiseless/trees")
fileNames <- list.files()
sapply(1:length(fileNames), function(x) file.rename(from = fileNames[x], to = sub(pattern="_noiseless_uvBM",replacement="", fileNames[x])))
sapply(1:length(fileNames), function(x) file.rename(from = fileNames[x], to = sub(pattern="correlationStrength",replacement="correlationStrength_noiseless", fileNames[x])))
fileNames <- list.files()
sapply(1:length(fileNames), function(x) file.rename(from = fileNames[x], to = sub(pattern="3correlationStrength_noisy",replacement="3correlationStrength", fileNames[x])))
fileNames <- list.files()
sapply(1:length(fileNames), function(x) file.rename(from = fileNames[x], to = sub(pattern="3correlationStrength",replacement="3correlationStrength_noiseless", fileNames[x])))

setwd(dir = "/Volumes/2TB/mvBM_sims_PB_noiseless/output")
fileNames <- list.files()
sapply(1:length(fileNames), function(x) file.rename(from = fileNames[x], to = sub(pattern="correlationStrength",replacement="correlationStrength_noiseless", fileNames[x])))


setwd(dir = "/Volumes/2TB/mvBM_sims_PB_extra/output_noiseless_pb_uvn_128")
fileNames <- list.files()
sapply(1:length(fileNames), function(x) file.rename(from = fileNames[x], to = sub(pattern="correlationStrength",replacement="correlationStrength_noiseless", fileNames[x])))
