data {
  int<lower=0> P;
  int<lower=0> Q;
  int<lower=0> N;
  int<lower=0> J;

  int<lower=0> Prior;

  int<lower=1> Num_Beta;
  
  array[Num_Beta] int Mu_Ind;
  array[Num_Beta] int Tau_Ind;
  array[N] int Delta_Ind;
  matrix[N, Num_Beta] Diag_Data;
  array[J] int Nvec;
  array[J] int Select_l;
  array[J] int Select_u;

  vector[N] y;
  vector[N] a;
  //vector[J] obs_per_pat;
  //vector[J] n;

}

parameters {
  
  array[J] real<lower=0> phi; // outcome sd
  array[J] real delta; // treatment effect
  
  vector[Q+P+1] mu;
  vector[Num_Beta] eps;
  array[J] real tau;
  array[Num_Beta] real sigma;

}


transformed parameters {
  
  vector[Num_Beta] gamma;

  for (q in 1:Num_Beta) {
    gamma[q] = mu[Mu_Ind[q]] + tau[Tau_Ind[q]] * sigma[q] * eps[q];
  }

}

model {
  
  phi ~ normal(0,1);
  delta ~ normal(0,10);
  eps ~ normal(0, 1);
  
  mu ~ normal(0,10);
  tau ~ cauchy(0,1);
  sigma ~ cauchy(0,1);

  for (j in 1:J) {
    y[Select_l[j]:Select_u[j]] ~ normal(delta[j]*a[Select_l[j]:Select_u[j]] + Diag_Data[Select_l[j]:Select_u[j], ] * gamma, phi[j]);
  }
  
}



generated quantities{
  
  real Nr = N;
  
  real trt_effect = 0;

  for (j in 1:J) {
    trt_effect = trt_effect + Nvec[j] * delta[j]/Nr;
  }
  
}
