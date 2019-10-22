tuningParam <- runif(1,0,1)
currentMatrix <- rlkjcorr(1,10,1)
lkjDraw <- rlkjcorr(1,10,1)
proposedMatrix <- currentMatrix*tuningParam + lkjDraw*(1-tuningParam)
all(eigen(proposedMatrix)$values > 0) #check for PSD

#proposal ratio?
#what matrixdo we need to add to proposedMatrix with given weights to get currentMatrix?
proposedLKJReverseDirection <- (currentMatrix - tuningParam*proposedMatrix) / (1-tuningParam)
eigen(proposedLKJReverseDirection)$values
mean(abs(currentMatrix - (proposedMatrix*tuningParam + proposedLKJReverseDirection*(1-tuningParam))))

t(chol(currentMatrix)) %*% chol(currentMatrix)
