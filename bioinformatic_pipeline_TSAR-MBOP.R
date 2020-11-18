library(DiagrammeR)
library(DiagrammeRsvg)  # for conversion to svg
library(magrittr)
library(rsvg)  # for saving svg
graph_spec <- "digraph {graph [layout = dot, rankdir = TB, compound=true]


newrank='true'

# define the global styles of the nodes. We can override these in box if we wish
node [shape = rectangle, style = filled, fillcolor = White]

data1 [label = '@@1', shape = folder, fillcolor = Beige]

Filter [label =  '@@2']
FilteredData [label = '@@3', shape = folder, fillcolor = Beige]
optim [label = 'Iterative Optimization']
optimReg [label= 'Optimize Regularization Parameters']
optimMean [label = 'Optimize Population Means']
optimThreshold [label = 'Optimize Trait Thresholds']
optimCorrelation [label = 'Optimize Bivariate Correlations']
correctCor [label = 'Find Nearest PSD Correlation Matrix \n Regularize by Weighing with Identity']
  
estimates [label = 'Optimization Output', shape = folder]

combEst [label = 'Combine Estimates', shape = folder]

infPhylo [label = 'Infer Phylogeny, Conditional on\nEstimated Means + Correlations']
propTree [label = 'Make Proposal to \nTree Topology (15%)'] 
propRates [label = 'Make Proposal to \nTrait-Specific Rates (60%)'] 
propTL [label = 'Make Proposal to \nTree Length (7%)'] 
propBL [label = 'Make Proposal to \nBranch-Length Proportions (15%)'] 
propRateReg [label = 'Make Proposal to \nTrait-Rate Regularization (4%)']
mcmcOutput [label = 'MCMC Output', shape = folder]
mcmcDiag [label = 'MCMC Diagnostics']
combChains [label = 'Combine Chains', shape = folder]
processResults [label = 'Process and Visualize Results']

stateProbs [label = 'Compute State Conditionals Probabilities']
missingProbs [label = 'Compute State Missingness Probabilities']
combine [label = 'Combine with Bayes Theorem']
sample[label = 'Sample States']

gelmanRubin [label = 'Gelman-Rubin Convergence Diagnostic < 1.01']
ESS [label = 'Effective Sample Size > 1,000']
compareTreesR2 [label = 'Bipartition Probabilities R&#xb2; > 0.99']


subgraph clusterStep1 {
  label = <<B>Step 1</B>>
  fontsize = 20
  graph[style=dashed]
  
  subgraph clusterOptim {
  fontsize = 15
  label = <<I>Optimization</I> <br/> (Loop 1:12)>
  graph[style=solid]
  edge [dir='back', minlen = 1.25]
  rank = same optimReg optimMean optimThreshold optimCorrelation
  optimReg -> optimMean -> optimThreshold -> optimCorrelation
  optimCorrelation -> optimThreshold -> optimMean -> optimReg
  }
  
  subgraph clusterImpute {
  fontsize = 15
  label = <<I>Stochastic Imputation</I> <br/> (Loop 5:12)>
  graph[style=solid]
  edge [dir='forward']
  rank = same stateProbs missingProbs combine sample
  stateProbs -> missingProbs -> combine -> sample
  }

  estimates
  #evaluateConcordance
  combEst
}


subgraph clusterStep2 {
  label = <<B>Step 2</B>>
  fontsize = 20
  graph[style=dashed]

  subgraph clusterMCMC {
  
  label = <<I>Markov Chain Monte Carlo</I> <br/> (10,000,000 Iterations)>
  fontsize = 15
  edge [dir='forward']
  graph[style=solid]
  rank = same propTree propRates propTL propBL propRateReg
  propRates -> propTL -> propBL -> propRateReg    
  propRateReg -> propBL -> propTL -> propRates
  propTree -> propRates -> propTree
  
  }

  subgraph clusterMCMCdiag {
    
  label = <<I>MCMC Diagnostics</I>>
  fontsize = 15
  edge [dir='forward']
  graph[style=solid]
  rank = same gelmanRubin ESS compareTreesR2
  
  
  }

mcmcOutput
mcmcDiag 
combChains

}

# edge definitions with the node IDs
{data1} -> Filter -> FilteredData -> optim
optim -> optimMean            [lhead = clusterOptim, arrowhead = none, label = '']
optimCorrelation -> stateProbs [lhead = clusterImpute, ltail = clusterOptim, minlen = 2]
sample -> optimReg [lhead = clusterOptim, ltail = clusterImpute, minlen = 2]
optim -> estimates 
optim -> estimates
optim -> estimates
optim -> estimates
# estimates -> evaluateConcordance
# estimates -> evaluateConcordance
# estimates -> evaluateConcordance
# estimates -> evaluateConcordance
# evaluateConcordance -> combEst
# evaluateConcordance -> combEst
# evaluateConcordance -> combEst
# evaluateConcordance -> combEst
estimates -> combEst
estimates -> combEst
estimates -> combEst
estimates -> combEst
combEst -> correctCor
correctCor -> infPhylo
infPhylo -> propTree            [lhead = clusterMCMC, arrowhead = none, label = '']
infPhylo -> mcmcOutput
infPhylo -> mcmcOutput
infPhylo -> mcmcOutput
infPhylo -> mcmcOutput
mcmcOutput -> mcmcDiag
mcmcOutput -> mcmcDiag
mcmcOutput -> mcmcDiag
mcmcOutput -> mcmcDiag
mcmcDiag -> combChains
mcmcDiag -> combChains
mcmcDiag -> combChains
mcmcDiag -> combChains
mcmcDiag -> gelmanRubin [lhead = clusterMCMCdiag, arrowhead = none, minlen = 2]
combChains -> processResults
propTL -> ESS [lhead = clusterMCMCdiag, ltail = clusterMCMC]
}



[1]: paste0('Raw Dataset \\n (722 individuals, 137 traits)')    
[2]: paste0('Data \\n Filtration')
[3]: paste0('Filtered Dataset \\n (684 individuals, 118 traits)')
"

grViz(graph_spec)
grViz(graph_spec) %>% export_svg %>% charToRaw %>% rsvg_pdf("~/dissertation/figures/tsarmbop.pdf")

