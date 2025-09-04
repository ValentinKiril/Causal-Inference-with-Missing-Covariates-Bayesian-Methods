find_patterns = function(Q) {
  M = rbind(replicate(100*2^Q, sample(0:1, Q, replace = TRUE)))
  patterns = apply(M, 2, function(x) {return(paste0(x, collapse = ""))})
  return(unique(patterns))
}

sim_data_s1 = function(N=500, P=2, p_zero=0, Q=3, trt=NULL){

  ## Simulation parameters will keep these as constant truths throughout simulation runs

  #coef_X = c(1, sample(c(-1,1), P, replace=T )) 
  coef_X = c(1, -1, 1)
  bsl_miss = 0 ## baseline odds of missingness for each variable
  coef_M = c(bsl_miss, sample(c(-1,1), P, replace=T ))

  ## assuming propensity score coefficients same across patterns..may want to modify.
  Potential_bsl_trt = c(-0.5, 0, 0.5)
  Potential_Coef_1 = c(-0.5)
  Potential_Coef_2 = c(0.5)
  coef_A = matrix(c(0), ncol = 1+P+Q, nrow = length(find_patterns(Q)))
  row.names(coef_A) = find_patterns(Q)
  for (t in find_patterns(Q)[1:4]) {
    coef_A[t, ] = c(sample(Potential_bsl_trt, 1, replace=TRUE),
                  1, -1, ## coefficients on V; needs to be updated if changing P and Q
                  sample(Potential_Coef_1, Q, replace=TRUE)) ## coefficients on Q
  }
  for (t in find_patterns(Q)[4:8]) {
    coef_A[t, ] = c(sample(Potential_bsl_trt, 1, replace=TRUE),
                  1, -1, ## coefficients on V; needs to be updated if changing P and Q
                  sample(Potential_Coef_2, Q, replace=TRUE)) ## coefficients on Q
  }

  bsl_outcome = 1
  coef_Y = c(bsl_outcome, sample(c(-1,1), P+Q, replace=T ) )

  #------------------------------------------------------------------------------#
  #--------------------------Simulation Covariates     --------------------------#
  #------------------------------------------------------------------------------#
  ## simulate from joint covariate distribution f(v, x) = f(x | v) f(v)
  
  ## assume V's are independence and X's are independent conditional on V
  ##  will consider scenarios with correlation structure.
  V = cbind(replicate(P, rnorm(N)))
  X = cbind(replicate(Q, rnorm(N, mean = cbind(1,V) %*% coef_X  )))
  
  #------------------------------------------------------------------------------#
  #--------------------------Simulate Missingness   -----------------------------#
  #------------------------------------------------------------------------------#
  ## simulate missingness mechanism (assume not dependent on X itself)
  ##    for Aim 4, we will relax this.
  ## also assume missingness for each var are conditionally independent
  ##    will consider scenarios with correlation structure.
  
  M = cbind(replicate(Q, expr = {
    rbinom(n = N, size = 1, prob = plogis( cbind(1,V) %*% coef_M )) 
  } ))

  #------------------------------------------------------------------------------#
  #--------------------------Simulate Treatment   -------------------------------#
  #------------------------------------------------------------------------------#

  ## assuming treatment coefficient on X_j is same across patterns 
  ##    consider scenarios where this is not true.
  pattern_ids = apply(M, 1, function(row) paste0(row, collapse = '') )
  table(pattern_ids)
  
  unique_patterns = unique(pattern_ids)
  n_obs_patterns = length(table(pattern_ids))
  
  A = numeric(N)
  
  if(is.null(trt)){
    
    for(pattern in unique_patterns){
      
      in_pattern = (pattern_ids==pattern)
      
      N_pattern = sum(in_pattern)
      
      pattern_idx = as.logical(as.integer(unlist(strsplit(pattern,split = ''))))
      
      V_pattern = V[in_pattern, , drop=F ]
      X_pattern = X[in_pattern, pattern_idx , drop=F ]
      
      coef_A_pattern = coef_A[row.names(coef_A) == pattern, c(rep(TRUE, P+1), pattern_idx)]
      
      eta_pattern = cbind(1, V_pattern, X_pattern ) %*% coef_A_pattern
      
      A[in_pattern] = rbinom(N_pattern, 1, plogis( eta_pattern  ))
      
    }
  }else{
    A[1:N] = trt
  }
  
  #------------------------------------------------------------------------------#
  #--------------------------Simulate Outcome   ---------------------------------#
  #------------------------------------------------------------------------------#
  ## linear/additive covariate effects
  ##    may consider other scenarios with complex functional forms
  
  Y = rnorm(N, cbind(1, V, X) %*% coef_Y + 5*A )
  
  X_obs = X
  X_obs[M==0] = NA
  
  V = cbind(V, replicate(p_zero, rnorm(N)))
  
  data = data.frame(V=V, X=X_obs, M=M, A, Y)
  for (i in 1:P) {
    data[,i] = unlist(data[,i])
  }
  
  data$pattern = apply(data[, grepl("^M", colnames(data)),drop=F], 1, paste0, collapse='')
  
  return(data)
}
