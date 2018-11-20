#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(observations_i);
  
  // Parameters
  PARAMETER(log_sigma_pro);
  PARAMETER(log_sigma_obs);
  PARAMETER_VECTOR(estimates_i);
  
  Type sigma_pro = exp(log_sigma_pro);
  Type sigma_obs = exp(log_sigma_obs);
  
  // Objective function
  Type jnll = 0;
  
  
  // Probability of initial observation
    jnll -= dnorm(observations_i(0),
                  estimates_i(0),
                  sigma_obs,
                  true);
  
  
  // Process model likelihood
  int n_i = estimates_i.size();
  for( int i = 1; i < n_i; i++) {
    jnll -= dnorm(estimates_i(i),
                  estimates_i(i-1),
                  sigma_pro, true);
  }
  
  // Observation model likelihood 
  for( int i = 1; i < n_i; i++) {
      jnll -= dnorm(observations_i(i),
                    estimates_i(i),
                    sigma_obs, true);
  }
  
  // Reporting
  ADREPORT(sigma_pro);
  ADREPORT(sigma_obs);
  ADREPORT(estimates_i);
  
  return jnll;
}
