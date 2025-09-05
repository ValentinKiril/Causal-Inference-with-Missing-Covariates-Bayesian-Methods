
bayes_psom_horse = function(data){
  # data = data[indices, ] ## for bootstrap if indices=1:N, then nothing happens.
  
  ## create missing data pattern indicator
  data$pattern = apply(data[, grepl("^M", colnames(data)),drop=F], 1, paste0, collapse='')
  
  ## create unique list of observed pattern/count how many there are
  pattern_list = unique(data$pattern)
  n_patterns = length(pattern_list)
  pattern_id = 1:n_patterns
  names(pattern_id)= pattern_list
  
  data$pattern_id = pattern_id[data$pattern]
  
  Q = sum(grepl("^M", colnames(data)))
  P = sum(grepl("^V", colnames(data)))
  
  xmat = data[ , grepl("^X", colnames(data))]
  xmat[is.na(xmat)] = 999
  
  xmat_split = split( xmat, data$pattern_id )
  xmat_split = lapply(xmat_split, function(x) x[ , apply(x,2, function(x) sum(x==999) )==0, drop=F ] )
  vmat_split = split( data[,grepl("^V", colnames(data))], data$pattern_id )
  y_split = split( data$Y, data$pattern_id)
  a_split = split( data$A, data$pattern_id)

  names(xmat_split) = paste0('xmat', 1:n_patterns)
  names(vmat_split) = paste0('vmat', 1:n_patterns)
  names(y_split) = paste0('y', 1:n_patterns)
  names(a_split) = paste0('a', 1:n_patterns)

  yvec = c()
  avec = c()
  for (i in 1:n_patterns) {
    yvec = append(yvec, y_split[[i]])
    avec = append(avec, a_split[[i]])
  }

  # Add indicator variable (lots of 1s) to design matrix of V
  for(i in 1:n_patterns) {
    vmat_split[[i]] = cbind(data.frame(Int = rep(1, nrow(vmat_split[[i]]))), 
                      vmat_split[[i]])
  }

  nvec = sapply(1:n_patterns, function(x) sum(data$pattern_id==x) )
  names(nvec) = paste0('N', 1:n_patterns)
  
  # Each row of M are indicators for inclusion of coefficient in pattern
  M = matrix(c(0), ncol = Q, nrow = max(data$pattern_id))
  for (i in 1:max(data$pattern_id)) {
    M[i,] = as.numeric(((data[data$pattern_id == i,])[,grepl("^M", colnames(data)),drop=F])[1,])
  }

  patmat = lapply(1:n_patterns, function(x){
    tt = data[data$pattern_id==x,grepl("^M", colnames(data)) ,drop=F] 
    res =  c(1:Q)[colSums((tt==1)) >0]
    ares = as.array(res, dim=length(res))
    return( ares  )
  })
  
  names(patmat) = paste0('qidx', 1:n_patterns)
  
  Qvec = lapply(patmat, function(x) length(x))
  names(Qvec) = paste0('Q',1:n_patterns)
  
  # Vector of indicators for selecting beta coefficients
  Indicator = c()
  for (i in 1:n_patterns) {
    Indicator = append(Indicator, rep(1, P+1))
    Indicator = append(Indicator, M[i, ])
  }
  
  # Total number of betas: for Q, for P, and intercepts
  Num_Beta = sum(M) + (P+1)*n_patterns

  # Selector for Mu parameter in STAN
  Mu_ind = c()
  for (i in 1:Num_Beta) {
    count = 0
    j=0
    while(count < i) {
      j = j+1
      if (Indicator[j]==1) {
        count <- count + 1
      }
    }
    Mu_ind[i] = j %% (1+P+Q)
  }
  Mu_ind[Mu_ind==0] = 1+P+Q

  # Selector for Tau parameter in STAN
  Tau_ind = c()
  for (i in 1:Num_Beta) {
    count = 0
    j=0
    while(count < i) {
      j = j+1
      if (Indicator[j]==1) {
        count <- count + 1
      }
    }
    Tau_ind[i] = floor(j / (1+P+Q)) + 1
  }
  Tau_ind[Tau_ind>n_patterns] = n_patterns

  Select_l = c(1)
  Select_u = c(nvec[[1]])
  for (j in 1:(n_patterns-1)) {
    Select_l[j+1] = sum(nvec[1:j]) + 1
    Select_u[j+1] = sum(nvec[1:(j+1)])
  }

  Delta_ind = c()
  for (i in 1:n_patterns) {
    Delta_ind = append(Delta_ind, rep(i, nvec[i]))
  }
  
  # Create diagonal block matrix of data (the covariates)
  Block_List <- list()
  for(i in 1:n_patterns) {
    Block_List = list.append(Block_List,
                    cbind(as.matrix(vmat_split[[i]]), as.matrix(xmat_split[[i]])))
  }
  Diag_Data = as.matrix(Matrix::bdiag(Block_List))

  stan_data1 =  list( Diag_Data=list(Diag_Data), c(N=nrow(data), P=P, Q=Q, J=n_patterns,
                        Num_Beta=Num_Beta, Mu_Ind=list(Mu_ind),
                        Tau_Ind=list(Tau_ind), Delta_Ind=list(Delta_ind),
                        y=list(yvec), Nvec=list(nvec), a=list(avec),
                        Select_l = list(Select_l), Select_u = list(Select_u)),
                        n = list(nvec), obs_per_pat = list( 15*(rowSums(M) + P+1)))
  stan_data = c(stan_data1, Prior=1)
  stan_data <- as.list(unlist(stan_data, recursive = FALSE)) ## flatten
  stan_data <- list.append(stan_data, a=avec, y=yvec)
      
  modres = sampling(stan_model_psom_horse, data = stan_data, chains = 1,
     	              warmup = 1000, iter = 2000, thin = 1,
                    pars = c('trt_effect'), verbose = F)


  output_results = summary(modres)$summary[ 'trt_effect' , 'mean' ]
  
  post_draws = rstan::extract(modres, pars = 'trt_effect')$trt_effect
  
  post_mean = mean(post_draws)
  post_lwr = quantile(post_draws, probs = .025)
  post_upr = quantile(post_draws, probs = .975)
  
  output_results = c(post_mean, post_lwr, post_upr)
  names(output_results) = c('mean', 'lwr', 'upr')
  
  return(output_results)
}

bayes_psps_horse = function(data, indices=1:N){

  data = data[indices, ]
  
  data$pattern = apply(data[, grepl("^M", colnames(data)),drop=F], 1, paste0, collapse='')
  
  ## create unique list of observed pattern/count how many there are
  pattern_list = unique(data$pattern)
  n_patterns = length(pattern_list)
  pattern_id = 1:n_patterns
  names(pattern_id)= pattern_list
  
  data$pattern_id = pattern_id[data$pattern]
  Q = sum(grepl("^M", colnames(data)))
  P = sum(grepl("^V", colnames(data)))

  xmat = data[ , grepl("^X", colnames(data))]
  xmat[is.na(xmat)] = 999
  
  xmat_split = split( xmat, data$pattern_id )
  xmat_split = lapply(xmat_split, function(x) x[ , apply(x,2, function(x) sum(x==999) )==0, drop=F ] )
  vmat_split = split( data[,grepl("^V", colnames(data))], data$pattern_id )
  y_split = split( data$Y, data$pattern_id)
  a_split = split( data$A, data$pattern_id)

  names(xmat_split) = paste0('xmat', 1:n_patterns)
  names(vmat_split) = paste0('vmat', 1:n_patterns)
  names(y_split) = paste0('y', 1:n_patterns)
  names(a_split) = paste0('a', 1:n_patterns)

  yvec = c()
  avec = c()
  for (i in 1:n_patterns) {
    yvec = append(yvec, y_split[[i]])
    avec = append(avec, a_split[[i]])
  }

  # Add indicator variable (lots of 1s) to design matrix of V
  for(i in 1:n_patterns) {
    vmat_split[[i]] = cbind(data.frame(Int = rep(1, nrow(vmat_split[[i]]))), 
                      vmat_split[[i]])
  }
  
  nvec = sapply(1:n_patterns, function(x) sum(data$pattern_id==x) )
  names(nvec) = paste0('N', 1:n_patterns)
  
  # Each row of M are indicators for inclusion of coefficient in pattern
  M = matrix(c(0), ncol = Q, nrow = max(data$pattern_id))
  for (i in 1:max(data$pattern_id)) {
    M[i,] = as.numeric(((data[data$pattern_id == i,])[,grepl("^M", colnames(data)),drop=F])[1,])
  }

  # Total number of betas: for Q, for P, and intercepts
  Num_Beta = sum(M) + (P+1)*n_patterns
  
  patmat = lapply(1:n_patterns, function(x){
    tt = data[data$pattern_id==x,grepl("^M", colnames(data)) ,drop=F] 
    res =  c(1:Q)[colSums((tt==1)) >0]
    ares = as.array(res, dim=length(res))
    return( ares  )
  })
  
  names(patmat) = paste0('qidx', 1:n_patterns)
  
  Qvec = lapply(patmat, function(x) length(x))
  names(Qvec) = paste0('Q',1:n_patterns)

  # Vector of indicators for selecting beta coefficients
  Indicator = c()
  for (i in 1:n_patterns) {
    Indicator = append(Indicator, rep(1, P+1))
    Indicator = append(Indicator, M[i, ])
  }
 
  # Selector for Mu parameter in STAN
  Mu_ind = c()
  for (i in 1:Num_Beta) {
    count = 0
    j=0
    while(count < i) {
      j = j+1
      if (Indicator[j]==1) {
        count <- count + 1
      }
    }
    Mu_ind[i] = j %% (1+P+Q)
  }
  Mu_ind[Mu_ind==0] = 1+P+Q

  # Selector for Tau parameter in STAN
  Tau_ind = c()
  for (i in 1:Num_Beta) {
    count = 0
    j=0
    while(count < i) {
      j = j+1
      if (Indicator[j]==1) {
        count <- count + 1
      }
    }
    Tau_ind[i] = floor(j / (1+P+Q)) + 1
  }
  Tau_ind[Tau_ind>n_patterns] = n_patterns

  # Create diagonal block matrix of data (the covariates)
  Block_List <- list()
  for(i in 1:n_patterns) {
    Block_List = list.append(Block_List,
                    cbind(as.matrix(vmat_split[[i]]), as.matrix(xmat_split[[i]])))
  }
  Diag_Data = as.matrix(Matrix::bdiag(Block_List))

  stan_data1 =  list( Diag_Data=list(Diag_Data), c(N=nrow(data), P=P, Q=Q, J=n_patterns,
                         Num_Beta=Num_Beta, a=list(avec), Mu_Ind=list(Mu_ind), 
                         Tau_Ind=list(Tau_ind)))
  stan_data = c(stan_data1, Prior=1)
  stan_data <- as.list(unlist(stan_data, recursive = FALSE)) ## flatten

  optim_res = rstan::optimizing(stan_model_psps_horse, data = stan_data,
                                as_vector = F)

  calc_weight = function(pscore, a)  (a*mean(data$A)/pscore) + ((1-a)*(1-mean(data$A))/(1-pscore))
  
  w1 = calc_weight(optim_res$par['pscore']$pscore, avec)
  
  n = nrow(data)
  
  delta_hat = lm( yvec ~ avec, weights = w1)$coefficients[2]
  
  return(delta_hat)
}

bayes_psom_lasso = function(data){
  # data = data[indices, ] ## for bootstrap if indices=1:N, then nothing happens.
  
  ## create missing data pattern indicator
  data$pattern = apply(data[, grepl("^M", colnames(data)),drop=F], 1, paste0, collapse='')
  
  ## create unique list of observed pattern/count how many there are
  pattern_list = unique(data$pattern)
  n_patterns = length(pattern_list)
  pattern_id = 1:n_patterns
  names(pattern_id)= pattern_list
  
  data$pattern_id = pattern_id[data$pattern]
  
  Q = sum(grepl("^M", colnames(data)))
  P = sum(grepl("^V", colnames(data)))
  
  xmat = data[ , grepl("^X", colnames(data))]
  xmat[is.na(xmat)] = 999
  
  xmat_split = split( xmat, data$pattern_id )
  xmat_split = lapply(xmat_split, function(x) x[ , apply(x,2, function(x) sum(x==999) )==0, drop=F ] )
  vmat_split = split( data[,grepl("^V", colnames(data))], data$pattern_id )
  y_split = split( data$Y, data$pattern_id)
  a_split = split( data$A, data$pattern_id)

  names(xmat_split) = paste0('xmat', 1:n_patterns)
  names(vmat_split) = paste0('vmat', 1:n_patterns)
  names(y_split) = paste0('y', 1:n_patterns)
  names(a_split) = paste0('a', 1:n_patterns)

  yvec = c()
  avec = c()
  for (i in 1:n_patterns) {
    yvec = append(yvec, y_split[[i]])
    avec = append(avec, a_split[[i]])
  }

  # Add indicator variable (lots of 1s) to design matrix of V
  for(i in 1:n_patterns) {
    vmat_split[[i]] = cbind(data.frame(Int = rep(1, nrow(vmat_split[[i]]))), 
                      vmat_split[[i]])
  }

  nvec = sapply(1:n_patterns, function(x) sum(data$pattern_id==x) )
  names(nvec) = paste0('N', 1:n_patterns)
  
  # Each row of M are indicators for inclusion of coefficient in pattern
  M = matrix(c(0), ncol = Q, nrow = max(data$pattern_id))
  for (i in 1:max(data$pattern_id)) {
    M[i,] = as.numeric(((data[data$pattern_id == i,])[,grepl("^M", colnames(data)),drop=F])[1,])
  }

  patmat = lapply(1:n_patterns, function(x){
    tt = data[data$pattern_id==x,grepl("^M", colnames(data)) ,drop=F] 
    res =  c(1:Q)[colSums((tt==1)) >0]
    ares = as.array(res, dim=length(res))
    return( ares  )
  })
  
  names(patmat) = paste0('qidx', 1:n_patterns)
  
  Qvec = lapply(patmat, function(x) length(x))
  names(Qvec) = paste0('Q',1:n_patterns)
  
  # Vector of indicators for selecting beta coefficients
  Indicator = c()
  for (i in 1:n_patterns) {
    Indicator = append(Indicator, rep(1, P+1))
    Indicator = append(Indicator, M[i, ])
  }
  
  # Total number of betas: for Q, for P, and intercepts
  Num_Beta = sum(M) + (P+1)*n_patterns

  # Selector for Mu parameter in STAN
  Mu_ind = c()
  for (i in 1:Num_Beta) {
    count = 0
    j=0
    while(count < i) {
      j = j+1
      if (Indicator[j]==1) {
        count <- count + 1
      }
    }
    Mu_ind[i] = j %% (1+P+Q)
  }
  Mu_ind[Mu_ind==0] = 1+P+Q

  # Selector for Tau parameter in STAN
  Tau_ind = c()
  for (i in 1:Num_Beta) {
    count = 0
    j=0
    while(count < i) {
      j = j+1
      if (Indicator[j]==1) {
        count <- count + 1
      }
    }
    Tau_ind[i] = floor(j / (1+P+Q)) + 1
  }
  Tau_ind[Tau_ind>n_patterns] = n_patterns

  Select_l = c(1)
  Select_u = c(nvec[[1]])
  for (j in 1:(n_patterns-1)) {
    Select_l[j+1] = sum(nvec[1:j]) + 1
    Select_u[j+1] = sum(nvec[1:(j+1)])
  }

  Delta_ind = c()
  for (i in 1:n_patterns) {
    Delta_ind = append(Delta_ind, rep(i, nvec[i]))
  }
  
  # Create diagonal block matrix of data (the covariates)
  Block_List <- list()
  for(i in 1:n_patterns) {
    Block_List = list.append(Block_List,
                    cbind(as.matrix(vmat_split[[i]]), as.matrix(xmat_split[[i]])))
  }
  Diag_Data = as.matrix(Matrix::bdiag(Block_List))

  stan_data1 =  list( Diag_Data=list(Diag_Data), c(N=nrow(data), P=P, Q=Q, J=n_patterns,
                        Num_Beta=Num_Beta, Mu_Ind=list(Mu_ind),
                        Tau_Ind=list(Tau_ind), Delta_Ind=list(Delta_ind),
                        y=list(yvec), Nvec=list(nvec), a=list(avec),
                        Select_l = list(Select_l), Select_u = list(Select_u)),
                        n = list(nvec), obs_per_pat = list( 15*(rowSums(M) + P+1)))
  stan_data = c(stan_data1, Prior=2)
  stan_data <- as.list(unlist(stan_data, recursive = FALSE)) ## flatten

  modres = sampling(stan_model_psom_lasso, data = stan_data, chains = 1,
     	              warmup = 1000, iter = 2000, thin = 1,
                    pars = c('trt_effect'), verbose = F)


  output_results = summary(modres)$summary[ 'trt_effect' , 'mean' ]
  
  post_draws = rstan::extract(modres, pars = 'trt_effect')$trt_effect
  
  post_mean = mean(post_draws)
  post_lwr = quantile(post_draws, probs = .025)
  post_upr = quantile(post_draws, probs = .975)
  
  output_results = c(post_mean, post_lwr, post_upr)
  names(output_results) = c('mean', 'lwr', 'upr')
  
  return(output_results)
}

bayes_psps_lasso = function(data, indices=1:N){

  data = data[indices, ]

  data$pattern = apply(data[, grepl("^M", colnames(data)),drop=F], 1, paste0, collapse='')
  
  ## create unique list of observed pattern/count how many there are
  pattern_list = unique(data$pattern)
  n_patterns = length(pattern_list)
  pattern_id = 1:n_patterns
  names(pattern_id)= pattern_list
  
  data$pattern_id = pattern_id[data$pattern]

  Q = sum(grepl("^M", colnames(data)))
  P = sum(grepl("^V", colnames(data)))

  xmat = data[ , grepl("^X", colnames(data))]
  xmat[is.na(xmat)] = 999
  
  xmat_split = split( xmat, data$pattern_id )
  xmat_split = lapply(xmat_split, function(x) x[ , apply(x,2, function(x) sum(x==999) )==0, drop=F ] )
  vmat_split = split( data[,grepl("^V", colnames(data))], data$pattern_id )
  y_split = split( data$Y, data$pattern_id)
  a_split = split( data$A, data$pattern_id)

  names(xmat_split) = paste0('xmat', 1:n_patterns)
  names(vmat_split) = paste0('vmat', 1:n_patterns)
  names(y_split) = paste0('y', 1:n_patterns)
  names(a_split) = paste0('a', 1:n_patterns)

  yvec = c()
  avec = c()
  for (i in 1:n_patterns) {
    yvec = append(yvec, y_split[[i]])
    avec = append(avec, a_split[[i]])
  }

  # Add indicator variable (lots of 1s) to design matrix of V
  for(i in 1:n_patterns) {
    vmat_split[[i]] = cbind(data.frame(Int = rep(1, nrow(vmat_split[[i]]))), 
                      vmat_split[[i]])
  }
  
  nvec = sapply(1:n_patterns, function(x) sum(data$pattern_id==x) )
  names(nvec) = paste0('N', 1:n_patterns)
  
  # Each row of M are indicators for inclusion of coefficient in pattern
  M = matrix(c(0), ncol = Q, nrow = max(data$pattern_id))
  for (i in 1:max(data$pattern_id)) {
    M[i,] = as.numeric(((data[data$pattern_id == i,])[,grepl("^M", colnames(data)),drop=F])[1,])
  }

  # Total number of betas: for Q, for P, and intercepts
  Num_Beta = sum(M) + (P+1)*n_patterns
  
  patmat = lapply(1:n_patterns, function(x){
    tt = data[data$pattern_id==x,grepl("^M", colnames(data)) ,drop=F] 
    res =  c(1:Q)[colSums((tt==1)) >0]
    ares = as.array(res, dim=length(res))
    return( ares  )
  })
  
  names(patmat) = paste0('qidx', 1:n_patterns)
  
  Qvec = lapply(patmat, function(x) length(x))
  names(Qvec) = paste0('Q',1:n_patterns)

  # Vector of indicators for selecting beta coefficients
  Indicator = c()
  for (i in 1:n_patterns) {
    Indicator = append(Indicator, rep(1, P+1))
    Indicator = append(Indicator, M[i, ])
  }
 
  # Selector for Mu parameter in STAN
  Mu_ind = c()
  for (i in 1:Num_Beta) {
    count = 0
    j=0
    while(count < i) {
      j = j+1
      if (Indicator[j]==1) {
        count <- count + 1
      }
    }
    Mu_ind[i] = j %% (1+P+Q)
  }
  Mu_ind[Mu_ind==0] = 1+P+Q

  # Selector for Tau parameter in STAN
  Tau_ind = c()
  for (i in 1:Num_Beta) {
    count = 0
    j=0
    while(count < i) {
      j = j+1
      if (Indicator[j]==1) {
        count <- count + 1
      }
    }
    Tau_ind[i] = floor(j / (1+P+Q)) + 1
  }
  Tau_ind[Tau_ind>n_patterns] = n_patterns

  # Create diagonal block matrix of data (the covariates)
  Block_List <- list()
  for(i in 1:n_patterns) {
    Block_List = list.append(Block_List,
                    cbind(as.matrix(vmat_split[[i]]), as.matrix(xmat_split[[i]])))
  }
  Diag_Data = as.matrix(Matrix::bdiag(Block_List))

  stan_data1 =  list( Diag_Data=list(Diag_Data), c(N=nrow(data), P=P, Q=Q, J=n_patterns,
                         Num_Beta=Num_Beta, a=list(avec), Mu_Ind=list(Mu_ind), 
                         Tau_Ind=list(Tau_ind)))
  stan_data = c(stan_data1, Prior=2)
  stan_data <- as.list(unlist(stan_data, recursive = FALSE)) ## flatten

  optim_res = rstan::optimizing(stan_model_psps_lasso, data = stan_data, as_vector = F)

  calc_weight = function(pscore, a)  (a*mean(data$A)/pscore) + ((1-a)*(1-mean(data$A))/(1-pscore))
  
  w1 = calc_weight(optim_res$par['pscore']$pscore, avec)
  
  n = nrow(data)
  
  delta_hat = lm( yvec ~ avec, weights = w1)$coefficients[2]
  
  return(delta_hat)
}

bayes_psom_ridge = function(data){
  # data = data[indices, ] ## for bootstrap if indices=1:N, then nothing happens.
  
  ## create missing data pattern indicator
  data$pattern = apply(data[, grepl("^M", colnames(data)),drop=F], 1, paste0, collapse='')
  
  ## create unique list of observed pattern/count how many there are
  pattern_list = unique(data$pattern)
  n_patterns = length(pattern_list)
  pattern_id = 1:n_patterns
  names(pattern_id)= pattern_list
  
  data$pattern_id = pattern_id[data$pattern]
  
  Q = sum(grepl("^M", colnames(data)))
  P = sum(grepl("^V", colnames(data)))
  
  xmat = data[ , grepl("^X", colnames(data))]
  xmat[is.na(xmat)] = 999
  
  xmat_split = split( xmat, data$pattern_id )
  xmat_split = lapply(xmat_split, function(x) x[ , apply(x,2, function(x) sum(x==999) )==0, drop=F ] )
  vmat_split = split( data[,grepl("^V", colnames(data))], data$pattern_id )
  y_split = split( data$Y, data$pattern_id)
  a_split = split( data$A, data$pattern_id)

  names(xmat_split) = paste0('xmat', 1:n_patterns)
  names(vmat_split) = paste0('vmat', 1:n_patterns)
  names(y_split) = paste0('y', 1:n_patterns)
  names(a_split) = paste0('a', 1:n_patterns)

  yvec = c()
  avec = c()
  for (i in 1:n_patterns) {
    yvec = append(yvec, y_split[[i]])
    avec = append(avec, a_split[[i]])
  }

  # Add indicator variable (lots of 1s) to design matrix of V
  for(i in 1:n_patterns) {
    vmat_split[[i]] = cbind(data.frame(Int = rep(1, nrow(vmat_split[[i]]))), 
                      vmat_split[[i]])
  }

  nvec = sapply(1:n_patterns, function(x) sum(data$pattern_id==x) )
  names(nvec) = paste0('N', 1:n_patterns)
  
  # Each row of M are indicators for inclusion of coefficient in pattern
  M = matrix(c(0), ncol = Q, nrow = max(data$pattern_id))
  for (i in 1:max(data$pattern_id)) {
    M[i,] = as.numeric(((data[data$pattern_id == i,])[,grepl("^M", colnames(data)),drop=F])[1,])
  }

  patmat = lapply(1:n_patterns, function(x){
    tt = data[data$pattern_id==x,grepl("^M", colnames(data)) ,drop=F] 
    res =  c(1:Q)[colSums((tt==1)) >0]
    ares = as.array(res, dim=length(res))
    return( ares  )
  })
  
  names(patmat) = paste0('qidx', 1:n_patterns)
  
  Qvec = lapply(patmat, function(x) length(x))
  names(Qvec) = paste0('Q',1:n_patterns)
  
  # Vector of indicators for selecting beta coefficients
  Indicator = c()
  for (i in 1:n_patterns) {
    Indicator = append(Indicator, rep(1, P+1))
    Indicator = append(Indicator, M[i, ])
  }
  
  # Total number of betas: for Q, for P, and intercepts
  Num_Beta = sum(M) + (P+1)*n_patterns

  # Selector for Mu parameter in STAN
  Mu_ind = c()
  for (i in 1:Num_Beta) {
    count = 0
    j=0
    while(count < i) {
      j = j+1
      if (Indicator[j]==1) {
        count <- count + 1
      }
    }
    Mu_ind[i] = j %% (1+P+Q)
  }
  Mu_ind[Mu_ind==0] = 1+P+Q

  # Selector for Tau parameter in STAN
  Tau_ind = c()
  for (i in 1:Num_Beta) {
    count = 0
    j=0
    while(count < i) {
      j = j+1
      if (Indicator[j]==1) {
        count <- count + 1
      }
    }
    Tau_ind[i] = floor(j / (1+P+Q)) + 1
  }
  Tau_ind[Tau_ind>n_patterns] = n_patterns

  Select_l = c(1)
  Select_u = c(nvec[[1]])
  for (j in 1:(n_patterns-1)) {
    Select_l[j+1] = sum(nvec[1:j]) + 1
    Select_u[j+1] = sum(nvec[1:(j+1)])
  }

  Delta_ind = c()
  for (i in 1:n_patterns) {
    Delta_ind = append(Delta_ind, rep(i, nvec[i]))
  }
  
  # Create diagonal block matrix of data (the covariates)
  Block_List <- list()
  for(i in 1:n_patterns) {
    Block_List = list.append(Block_List,
                    cbind(as.matrix(vmat_split[[i]]), as.matrix(xmat_split[[i]])))
  }
  Diag_Data = as.matrix(Matrix::bdiag(Block_List))

  stan_data1 =  list( Diag_Data=list(Diag_Data), c(N=nrow(data), P=P, Q=Q, J=n_patterns,
                        Num_Beta=Num_Beta, Mu_Ind=list(Mu_ind),
                        Tau_Ind=list(Tau_ind), Delta_Ind=list(Delta_ind),
                        y=list(yvec), Nvec=list(nvec), a=list(avec),
                        Select_l = list(Select_l), Select_u = list(Select_u)),
                        n = list(nvec), obs_per_pat = list( 15*(rowSums(M) + P+1)))
  stan_data = c(stan_data1, Prior=0)
  stan_data <- as.list(unlist(stan_data, recursive = FALSE)) ## flatten

  modres = sampling(stan_model_psom_ridge, data = stan_data, chains = 1,
     	              warmup = 1000, iter = 2000, thin = 1,
                    pars = c('trt_effect'), verbose = F)

  output_results = summary(modres)$summary[ 'trt_effect' , 'mean' ]
  
  post_draws = rstan::extract(modres, pars = 'trt_effect')$trt_effect
  
  post_mean = mean(post_draws)
  post_lwr = quantile(post_draws, probs = .025)
  post_upr = quantile(post_draws, probs = .975)
  
  output_results = c(post_mean, post_lwr, post_upr)
  names(output_results) = c('mean', 'lwr', 'upr')
  
  return(output_results)
}

bayes_psps_ridge = function(data, indices=1:N){

  data = data[indices, ]
  
  data$pattern = apply(data[, grepl("^M", colnames(data)),drop=F], 1, paste0, collapse='')
  
  ## create unique list of observed pattern/count how many there are
  pattern_list = unique(data$pattern)
  n_patterns = length(pattern_list)
  pattern_id = 1:n_patterns
  names(pattern_id)= pattern_list
  
  data$pattern_id = pattern_id[data$pattern]

  Q = sum(grepl("^M", colnames(data)))
  P = sum(grepl("^V", colnames(data)))

  xmat = data[ , grepl("^X", colnames(data))]
  xmat[is.na(xmat)] = 999
  
  xmat_split = split( xmat, data$pattern_id )
  xmat_split = lapply(xmat_split, function(x) x[ , apply(x,2, function(x) sum(x==999) )==0, drop=F ] )
  vmat_split = split( data[,grepl("^V", colnames(data))], data$pattern_id )
  y_split = split( data$Y, data$pattern_id)
  a_split = split( data$A, data$pattern_id)

  names(xmat_split) = paste0('xmat', 1:n_patterns)
  names(vmat_split) = paste0('vmat', 1:n_patterns)
  names(y_split) = paste0('y', 1:n_patterns)
  names(a_split) = paste0('a', 1:n_patterns)

  yvec = c()
  avec = c()
  for (i in 1:n_patterns) {
    yvec = append(yvec, y_split[[i]])
    avec = append(avec, a_split[[i]])
  }

  # Add indicator variable (lots of 1s) to design matrix of V
  for(i in 1:n_patterns) {
    vmat_split[[i]] = cbind(data.frame(Int = rep(1, nrow(vmat_split[[i]]))), 
                      vmat_split[[i]])
  }
  
  nvec = sapply(1:n_patterns, function(x) sum(data$pattern_id==x) )
  names(nvec) = paste0('N', 1:n_patterns)
  
  # Each row of M are indicators for inclusion of coefficient in pattern
  M = matrix(c(0), ncol = Q, nrow = max(data$pattern_id))
  for (i in 1:max(data$pattern_id)) {
    M[i,] = as.numeric(((data[data$pattern_id == i,])[,grepl("^M", colnames(data)),drop=F])[1,])
  }

  # Total number of betas: for Q, for P, and intercepts
  Num_Beta = sum(M) + (P+1)*n_patterns
  
  patmat = lapply(1:n_patterns, function(x){
    tt = data[data$pattern_id==x,grepl("^M", colnames(data)) ,drop=F] 
    res =  c(1:Q)[colSums((tt==1)) >0]
    ares = as.array(res, dim=length(res))
    return( ares  )
  })
  
  names(patmat) = paste0('qidx', 1:n_patterns)
  
  Qvec = lapply(patmat, function(x) length(x))
  names(Qvec) = paste0('Q',1:n_patterns)

  # Vector of indicators for selecting beta coefficients
  Indicator = c()
  for (i in 1:n_patterns) {
    Indicator = append(Indicator, rep(1, P+1))
    Indicator = append(Indicator, M[i, ])
  }
 
  # Selector for Mu parameter in STAN
  Mu_ind = c()
  for (i in 1:Num_Beta) {
    count = 0
    j=0
    while(count < i) {
      j = j+1
      if (Indicator[j]==1) {
        count <- count + 1
      }
    }
    Mu_ind[i] = j %% (1+P+Q)
  }
  Mu_ind[Mu_ind==0] = 1+P+Q

  # Selector for Tau parameter in STAN
  Tau_ind = c()
  for (i in 1:Num_Beta) {
    count = 0
    j=0
    while(count < i) {
      j = j+1
      if (Indicator[j]==1) {
        count <- count + 1
      }
    }
    Tau_ind[i] = floor(j / (1+P+Q)) + 1
  }
  Tau_ind[Tau_ind>n_patterns] = n_patterns

  # Create diagonal block matrix of data (the covariates)
  Block_List <- list()
  for(i in 1:n_patterns) {
    Block_List = list.append(Block_List,
                    cbind(as.matrix(vmat_split[[i]]), as.matrix(xmat_split[[i]])))
  }
  Diag_Data = as.matrix(Matrix::bdiag(Block_List))

  stan_data1 =  list( Diag_Data=list(Diag_Data), c(N=nrow(data), P=P, Q=Q, J=n_patterns,
                         Num_Beta=Num_Beta, a=list(avec), Mu_Ind=list(Mu_ind), 
                         Tau_Ind=list(Tau_ind)))
  stan_data = c(stan_data1, Prior=0)
  stan_data <- as.list(unlist(stan_data, recursive = FALSE)) ## flatten

  optim_res = rstan::optimizing(stan_model_psps_ridge, data = stan_data, as_vector = F)

  calc_weight = function(pscore, a)  (a*mean(data$A)/pscore) + ((1-a)*(1-mean(data$A))/(1-pscore))
  
  w1 = calc_weight(optim_res$par['pscore']$pscore, avec)
  
  n = nrow(data)
  
  delta_hat = lm( yvec ~ avec, weights = w1)$coefficients[2]
  
  return(delta_hat)
}
