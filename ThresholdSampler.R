#library(rethinking)

  nHave <- 0
  nTot <- 100
  pseq <- seq(0, 1, by = .00001)
  prior <- dbeta(x = pseq, shape1 = 1.01, shape2 = 1.01)
  posterior <- dbinom(nHave, nTot, pseq) * prior; posterior <- posterior / sum(posterior)
  plot(pseq, posterior, type = "l")
  freqDistr <- sample(x = pseq, size = 1e6, replace = T, prob = posterior)
  liabs <- qnorm(freqDistr, mean = 0, sd = 1)
  hist(liabs, freq = F, breaks = 30)
  simDens <- rnorm(n = 1e6, mean = mean(liabs), sd =  sd(liabs))
  dens(simDens, add = T, col = 2)

#given some frequency of presence/absence, returns a Bayesian estimate of the mean liability from weak beta prior
liabPost <- function(nHave, nTot, shape1 = 1.01, shape2 = 1.01) {
  pseq <- seq(0, 1, by = .00001)
  prior <- dbeta(x = pseq, shape1 = shape1, shape2 = shape2)
  posterior <- dbinom(nHave, nTot, pseq) * prior; posterior <- posterior / sum(posterior)
  freqDistr <- sample(x = pseq, size = 1e6, replace = T, prob = posterior)
  liabs <- qnorm(freqDistr, mean = 0, sd = 1)
  return(c(mean(liabs), sd(liabs)))
}

liabPost(4, 15)
