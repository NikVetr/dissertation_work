setwd("/Volumes/1TB/Bailey/")

data_files <- list.files("data")
data <- lapply(1:length(data_files), function(x) "placeholder")
for(i in 1:length(data)){
  load(paste0("data/", data_files[i]))
  data[[i]] <- current_params
}

data <- data[grep(data_files, pattern = "NJ")]

mean_of_means <- data[[1]]$means
for(i in 2:length(data)){
  mean_of_means <- mean_of_means + data[[i]]$means
}
mean_of_means <- mean_of_means / length(data)


means <- sapply(1:length(data), function(run) as.vector(data[[run]]$means))
pairs(means, col = rgb(0,0,0,0.1))

cors <- sapply(1:length(data), function(run) as.vector(data[[run]]$cor[upper.tri(data[[run]]$cor)]))
pairs(cors, col = rgb(0,0,0,0.3))

threshes <- sapply(1:length(data), function(run) as.vector(data[[run]]$threshold_mat)[as.vector(data[[run]]$threshold_mat > 0 & as.vector(data[[run]]$threshold_mat) != Inf)])
pairs(threshes, col = rgb(0,0,0,0.3))


mean_of_cors <- data[[1]]$cor
for(i in 2:length(data)){
  mean_of_cors <- mean_of_cors + data[[i]]$cor
}
mean_of_cors <- mean_of_cors / length(data)

eigensystem <- eigen(mean_of_cors)

plot(cumsum(eigensystem$values), type = "l"); abline(v = which.max(cumsum(eigensystem$values)))
