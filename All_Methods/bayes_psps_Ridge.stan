data {
  int<lower=0> P;
  int<lower=0> Q;
  int<lower=0> N;
  int<lower=0> J;

  int<lower=0> Prior;
  
  int<lower=1> Num_Beta;
  
  array[Num_Beta] int Mu_Ind;
  array[Num_Beta] int Tau_Ind;
  matrix[N, Num_Beta] Diag_Data;

  array[N] int<lower=0, upper=1> a;

}

parameters {
  
  vector[Q+P+1] mu; 
  vector[Num_Beta] eps;
  array[J] real tau;

}


transformed parameters {
  
  vector[Num_Beta] gamma;
  vector<lower=0, upper=1>[N] eta;

  for (q in 1:Num_Beta) {
    gamma[q] = mu[Mu_Ind[q]] + tau[Tau_Ind[q]] * eps[q];
  }

  eta = inv_logit( Diag_Data * gamma);

}

model {
 
  mu ~ normal(0, 10);
  eps ~ normal(0, 1);
  tau ~ normal(0, 1);

  a ~ bernoulli(eta);
  
}



generated quantities{
  
  vector[N] pscore;
  
  pscore = eta;
  
}
