

data {
  int<lower=0> N_obs;
  real obs[N_obs];
}

parameters {
  real<lower=0> sigma_pro;
  real<lower=0> sigma_obs;
  real tru_est[N_obs];
}

model {
  sigma_pro ~ cauchy(0, 5);
  sigma_obs ~ cauchy(0, 5);
  
  // Probability of initial true state
  tru_est[1] ~ normal(0, sigma_pro);
  
  // Probability of true state
  for (n in 2:N_obs)
  tru_est[n] ~ normal(tru_est[n-1], sigma_pro);
  
  // Likeilhood of obs
  for (n in 1:N_obs)
  obs[n] ~ normal(tru_est[n], sigma_obs);
}


