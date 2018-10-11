# See if we can replicate the simulate feature for a fitted model
n_ages <- 
for i in 2:n_t {
  N[1, i] <- w[1, i-1] * p[, i-1] %*% N[, i-1] * exp(rnorm(1))  
}
