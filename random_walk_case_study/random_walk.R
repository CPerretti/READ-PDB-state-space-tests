# Why the process error problem?
# Can we replicate it in the simplest state-space model?
# Try a random walk.
library(TMB)
library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(321) # for reproducibility

# Random walk
nT <- 100
true_states_i <- vector(mode = "numeric", length = nT)
true_states_i[1] <- rnorm(n = 1) # intial condition
for (i in 2:nT) {
  true_states_i[i] <- true_states_i[i-1] + rnorm(n = 1) # step forward with process error
}

observations_i <- true_states_i + rnorm(n = nT) # add observation error

# plot it 
plot(true_states_i, xlab = "Year", ylab = "Value", type = "l")
points(observations_i)

# Plot many realizations of the original model to get a sense of 
# the variability (no obs error)
nReps <- 1000
true_states_i <- matrix(data = NA, nrow = nT, ncol = nReps)
true_states_i[1,] <- rnorm(n = nReps) # intial condition
for (i in 2:100) {
  true_states_i[i,] <- true_states_i[i-1,] + rnorm(n = nReps) # step forward with process error
}

df2plot <-
  true_states_i %>%
  as.data.frame() %>%
  dplyr::mutate(year = 1:nT) %>%
  tidyr::gather(variable, value, -year)



# Variance of a random walk scales as var * t where t is the length of the time series.
# (Because the innovations are independent so their variance simply sums over time.)
# If var = 1 then at year 50 the var of possible states = 50.
# If var = 1 corresponds to the log scale, then on the natural scale at year
# 50 var = exp(50). I.e. it blows up.

# Assume this is log-scale
ggplot(df2plot, aes(x = year, y = value, color = variable)) +
  geom_line() +
  theme(legend.position = "none") +
  ggtitle("Log-scale things seems reasonable")

# Then this is natural scale
ggplot(df2plot, aes(x = year, y = exp(value), color = variable)) +
  geom_line() +
  theme(legend.position = "none") +
  ggtitle("Natural scale blow up")