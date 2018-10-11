# See if we can replicate the simulate feature for a fitted model 
# (start simple, not running yet)

# Use herring data and model fit to set up simulation


# Run simulation
n_a <- 8 # number of age-classes
n_t <- 20 # length of time series

for i in 2:n_t {
  N[1, i] <- r_func(w[, i-1] * p[, i-1] %*% N[, i-1]) * exp(rnorm(1))
  N[-c(1, n_a), i] <- N[-c(n_a-1, n_a), i-1] *
                      exp(f[-c(n_a-1, n_a), i-1]) * 
                      exp(m[-c(n_a-1, n_a), i-1]) * exp(rnorm(1))
  N[n_a, i] <- 
}
