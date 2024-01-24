countObsStates <- function(states, maxStates = 7){
  sttab <- (table(states))
  obs_states <- rep(0, maxStates)
  obs_states[as.numeric(attr(sttab, which = "dimnames")$states) + 1] <- as.integer(sttab)
  return(obs_states)
}

countUnobsStates <- function(trait_inds, state_probs, maxStates = 7){
  conv_buff_probs <- t(sapply(1:length(state_probs), function(foo) softmax(c(state_probs[[foo]], rep(-Inf, maxStates - length(state_probs[[foo]]))))))
  unique_trait_inds <- unique(trait_inds)
  t(sapply(unique_trait_inds, function(trait) apply(conv_buff_probs[trait_inds == trait,], 2, sum)))
}

rbeta_flat_update_sample_matrix <- function(matrix_a, matrix_b){
  samp_probs <- sapply(1:nrow(matrix_a), function(rowi) 
    sapply(1:ncol(matrix_a), function(colj) 
      rbeta(1, shape1 = 1 + matrix_a[rowi, colj], shape2 = 1 + matrix_b[rowi, colj])
    )
  )
  return(t(samp_probs))
}

find_missing_probs <- function(n_obs_trait_state_per_tip, n_imputed_trait_state_per_tip, per_pop = F, per_trait = T, raw = F){
  if(raw){
    if(per_pop){
      return(lapply(1:length(n_obs_trait_state_per_tip), function(tip) n_imputed_trait_state_per_tip[[tip]] / (n_imputed_trait_state_per_tip[[tip]] + n_obs_trait_state_per_tip[[tip]])))
    } else {
      n_obs_trait_state_per_tip_total <- n_obs_trait_state_per_tip[[1]]
      for(i in 2:length(n_obs_trait_state_per_tip)){
        n_obs_trait_state_per_tip_total <- n_obs_trait_state_per_tip_total + n_obs_trait_state_per_tip[[i]]
      }
      n_imputed_trait_state_per_tip_total <- n_imputed_trait_state_per_tip[[1]]
      for(i in 2:length(n_imputed_trait_state_per_tip)){
        n_imputed_trait_state_per_tip_total <- n_imputed_trait_state_per_tip_total + n_imputed_trait_state_per_tip[[i]]
      }
      if(!per_trait){
        return(sapply(1:7, function(x) sum(n_imputed_trait_state_per_tip_total[,x]) / sum((n_imputed_trait_state_per_tip_total[,x] + n_obs_trait_state_per_tip_total[,x]))))
      } else {
        return(n_imputed_trait_state_per_tip_total / (n_imputed_trait_state_per_tip_total + n_obs_trait_state_per_tip_total)) 
      }
    }
  } else {
    #update beta distribution with binomial count
    if(per_pop){
      return(lapply(1:length(n_obs_trait_state_per_tip), function(tip) 
        rbeta_flat_update_sample_matrix(n_imputed_trait_state_per_tip[[tip]], n_obs_trait_state_per_tip[[tip]])))
    } else {
      n_obs_trait_state_per_tip_total <- n_obs_trait_state_per_tip[[1]]
      for(i in 2:length(n_obs_trait_state_per_tip)){
        n_obs_trait_state_per_tip_total <- n_obs_trait_state_per_tip_total + n_obs_trait_state_per_tip[[i]]
      }
      n_imputed_trait_state_per_tip_total <- n_imputed_trait_state_per_tip[[1]]
      for(i in 2:length(n_imputed_trait_state_per_tip)){
        n_imputed_trait_state_per_tip_total <- n_imputed_trait_state_per_tip_total + n_imputed_trait_state_per_tip[[i]]
      }
      if(!per_trait){
        tr_imp <- sapply(1:ncol(n_imputed_trait_state_per_tip_total), function(x) sum(n_imputed_trait_state_per_tip_total[,x]))
        tr_obs <- sapply(1:ncol(n_obs_trait_state_per_tip_total), function(x) sum(n_obs_trait_state_per_tip_total[,x]))
        return(rbeta_flat_update_sample_matrix(t(as.matrix(tr_imp)), t(as.matrix(tr_obs))))
      } else {
        return(rbeta_flat_update_sample_matrix(n_imputed_trait_state_per_tip_total, n_obs_trait_state_per_tip_total)) 
      }
    }
  }
}

maxStates <- 7

n_obs_trait_state_per_tip <- lapply(1:length(obs), function(tip)
  t(sapply(1:ncol(obs[[tip]]), function(trait) countObsStates(traits_indiv_discr[[tip]][obs[[tip]][,trait],trait], maxStates = maxStates))
  )
)

n_unobs_true_trait_state_per_tip <- lapply(1:length(obs), function(tip)
  t(sapply(1:ncol(obs[[tip]]), function(trait) countObsStates(traits_indiv_discr_true[[tip]][!obs[[tip]][,trait],trait], maxStates = maxStates))
  )
)

n_imputed_trait_state_per_tip <- lapply(1:length(imputed_probs), function(tip)
  countUnobsStates(imputed_probs[[tip]][[1]][,2], imputed_probs[[tip]][[2]], maxStates = maxStates)
)

n_imputed_sampled_trait_state_per_tip <- lapply(1:length(obs), function(tip)
  t(sapply(1:ncol(obs[[tip]]), function(trait) countObsStates(traits_indiv_discr[[tip]][!obs[[tip]][,trait],trait], maxStates = maxStates))
  )
)

# n_imputed_trait_state_per_tip <- n_unobs_true_trait_state_per_tip

find_missing_probs(n_obs_trait_state_per_tip, n_imputed_sampled_trait_state_per_tip, per_pop = F, raw = F)
find_missing_probs(n_obs_trait_state_per_tip, n_imputed_trait_state_per_tip, per_pop = F, raw = F)
find_missing_probs(n_obs_trait_state_per_tip, n_unobs_true_trait_state_per_tip, per_pop = F, raw = F)
plot(find_missing_probs(n_obs_trait_state_per_tip, n_imputed_trait_state_per_tip, per_pop = F, raw = T, per_trait = F), sort(unique(probs_miss), decreasing = T)); abline(0,1)

