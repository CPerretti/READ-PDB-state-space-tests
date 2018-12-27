

data {
  //int<lower=0> N_obs;
  int<lower=0> N_mis;
  real<lower=0> rland_obs[N_obs];
  //real<lower=0> cland4obs[N_obs];
  real<lower=0> cland4mis[N_mis];
}

parameters {
  real<lower=0> sigma;
  real<lower=0> rland_mis[N_mis];
  real a;
  real b;
}

model {
  //sigma ~ cauchy(0, 5);
  //a ~ normal(0, 1e4);
  //b ~ normal(0, 1e4);
  
  // Likelihood of observed
  // for (n in 1:N_obs)
  // rland_obs[n] ~ normal(cland4obs[n], sigma);
  
  // Likeilhood of missing
  for (n in 1:N_mis)
  rland_mis[n] ~ normal(cland4mis[n], sigma);
}


