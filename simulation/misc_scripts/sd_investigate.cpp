#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(y);
  
  // Parameters
  PARAMETER(log_sigma);
  
  
  Type sigma = exp(log_sigma);
  
  // Number of observations
  int n_y = y.size();

  // Likelihood for intial obsevation
  Type jnll = 0;
  // jnll -= dnorm(y(0), 
  //               Type(0), sigma, true);
  
  
  // Likliehood for all other observations
  for (int i = 1; i < n_y; i++) {
    jnll -= dnorm(y(i) - y(i-1),
                  Type(0), sigma, true);
  }
  
  // Reporting
  ADREPORT(sigma);
  
  return jnll;
}
