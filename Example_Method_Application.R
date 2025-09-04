library(rlist)
library(dplyr)
library(tidyr)
library(rstan)

dir = "C:\\Users\\amval\\.vscode-R\\ForGitHub"
setwd(dir)

source("data_generate.R")
source("All_Methods/methods.R")

#------------------------------------------------------------------------------#
#--------------------------------- Generate Data ------------------------------#
#------------------------------------------------------------------------------#
set.seed(123)

data = sim_data_s1()

stan_model_psps_horse = rstan::stan_model("All_Methods/bayes_psps_Horse.stan")
stan_model_psom_horse = rstan::stan_model("All_Methods/bayes_psom_Horse.stan")
stan_model_psps_lasso = rstan::stan_model("All_Methods/bayes_psps_LASSO.stan")
stan_model_psom_lasso = rstan::stan_model("All_Methods/bayes_psom_LASSO.stan")
stan_model_psps_ridge = rstan::stan_model("All_Methods/bayes_psps_Ridge.stan")
stan_model_psom_ridge = rstan::stan_model("All_Methods/bayes_psom_Ridge.stan")

#------------------------------------------------------------------------------#
#-------------------------- 2. Run Methods on Data   --------------------------#
#------------------------------------------------------------------------------#

## for each method store 1) point estimate, 2) int. lower, & 3) int. upper
bayes_horse_psom_res = matrix(nrow=1, ncol=3) ## pattern-specific horse outcome model
bayes_horse_psps_res = matrix(nrow=1, ncol=3) ## pattern-specific horse pscore model
bayes_lasso_psom_res = matrix(nrow=1, ncol=3) ## pattern-specific lasso outcome model
bayes_lasso_psps_res = matrix(nrow=1, ncol=3) ## pattern-specific lasso pscore model
bayes_ridge_psom_res = matrix(nrow=1, ncol=3) ## pattern-specific ridge outcome model
bayes_ridge_psps_res = matrix(nrow=1, ncol=3) ## pattern-specific ridge pscore model


### --------  run Bayesian pattern-specific (Horse) p-score method -------###

bayes_horse_psps_output = bayes_psps_horse(data)
bayes_horse_psps_res[1, 1] = bayes_horse_psps_output['mean']
bayes_horse_psps_res[1, 2] = bayes_horse_psps_output['lwr']
bayes_horse_psps_res[1, 3] = bayes_horse_psps_output['upr']

### ---------  run Bayesian pattern-specific outcome model method ---------###
  
bayes_horse_psom_output = bayes_psom_horse(data)
bayes_horse_psom_res[1, 1] = bayes_horse_psom_output['mean']
bayes_horse_psom_res[1, 2] = bayes_horse_psom_output['lwr']
bayes_horse_psom_res[1, 3] = bayes_horse_psom_output['upr']

### --------  run Bayesian pattern-specific (LASSO) p-score method -------###

bayes_lasso_psps_output = bayes_psps_lasso(data)
bayes_lasso_psps_res[1, 1] = bayes_lasso_psps_output['mean']
bayes_lasso_psps_res[1, 2] = bayes_lasso_psps_output['lwr']
bayes_lasso_psps_res[1, 3] = bayes_lasso_psps_output['upr']
  
### ---------  run Bayesian pattern-specific outcome model method ---------###

bayes_lasso_psom_output = bayes_psom_lasso(data)
bayes_lasso_psom_res[1, 1] = bayes_lasso_psom_output['mean']
bayes_lasso_psom_res[1, 2] = bayes_lasso_psom_output['lwr']
bayes_lasso_psom_res[1, 3] = bayes_lasso_psom_output['upr']

### --------  run Bayesian pattern-specific (Ridge) p-score method -------###

bayes_ridge_psps_output = bayes_psps_ridge(data)
bayes_ridge_psps_res[1, 1] = bayes_ridge_psps_output['mean']
bayes_ridge_psps_res[1, 2] = bayes_ridge_psps_output['lwr']
bayes_ridge_psps_res[1, 3] = bayes_ridge_psps_output['upr']

### ---------  run Bayesian pattern-specific outcome model method ---------###

bayes_ridge_psom_output = bayes_psom_ridge(data)
bayes_ridge_psom_res[1, 1] = bayes_ridge_psom_output['mean']
bayes_ridge_psom_res[1, 2] = bayes_ridge_psom_output['lwr']
bayes_ridge_psom_res[1, 3] = bayes_ridge_psom_output['upr']


## ----------------------------------------------------------------------- ##
## --------------------------- Package Results --------------------------- ##
## ----------------------------------------------------------------------- ##

res = rbind( data.frame(method='Horseshoe_Outcome', bayes_horse_psom_res),
            data.frame(method='Horseshoe_IPW', bayes_horse_psps_res),
            data.frame(method='LASSO_Outcome', bayes_lasso_psom_res),
            data.frame(method='LASSO_IPW', bayes_lasso_psps_res),
            data.frame(method='Ridge_Outcome', bayes_ridge_psom_res),
            data.frame(method='Ridge_IPW', bayes_ridge_psps_res))

names(res) = c('method','point','lwr', 'upr')

# Round to 3 digits
for (i in 2:4) {
  for (j in 1:6) {
    res[j, i] = round(res[j, i], digits = 3)
  }
}

# The true treatment effect is 5.000
res