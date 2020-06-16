
# #instantiate the appropriate dimension ordered freek; substitute nStates for 195
# for (i in 1:nStates) {
#   Q_temp195[i] = rep(0,nStates)
#   if (i != 1 & i != nStates) {Q_temp[i][i-1] = 1.0; Q_temp[i][i+1] = 1.0} 
#   if (i == 1) {Q_temp[i][i+1] = 1.0}
#   if (i == nStates) {Q_temp[i][i-1] = 1.0}
# }
# 
# Q[dats] = fnFreeK(Q_temp195, rescaled = TRUE, v("scalingAndSquaring","scalingAndSquaringPade","scalingAndSquaringTaylor","uniformization","eigen")[2])
# 
# }

path <- "orderedFreekScripts/orderedFreekNUM.Rev"
temp <- readLines(path)
for(i in 1:500){
  print(i)
  newScript <- gsub(x = temp, pattern = "NUM", replacement = i)
  writeLines(newScript, con = gsub(x = path, pattern = "NUM", replacement = i))
}