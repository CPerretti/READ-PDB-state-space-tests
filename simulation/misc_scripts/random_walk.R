# Why the process error problem?
# Try a random walk.
library(TMB)
library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(321) # for reproducibility

# Plot many realizations of a random walk to get a sense of
# the variability (no obs error)
nT <- 100
nReps <- 1000
true_states_i <- matrix(data = NA, nrow = nT, ncol = nReps)
true_states_i[1,] <- rnorm(n = nReps) # intial condition
for (i in 2:100) {
  true_states_i[i,] <- true_states_i[i-1,] + rnorm(n = nReps) # step forward with process error
}


# Variance of a random walk scales as var * t where t is the length of the time series.
# (Because the innovations are independent so their variance simply sums over time.)
# If var = 1 then at year 100 the var of possible states = 1 * 100.
# If var = 1 corresponds to the log scale, then on the natural scale at year
# 100 var = exp(100). I.e. it blows up.

# Let's take a look
df2plot <-
  true_states_i %>%
  as.data.frame() %>%
  dplyr::mutate(year = 1:nT) %>%
  tidyr::gather(variable, value, -year)

# Assume this is log-scale
ggplot(df2plot, aes(x = year, y = value, color = variable)) +
  geom_line() +
  theme(legend.position = "none") +
  ggtitle("Log-scale things seem reasonable")

# Then this is natural scale
ggplot(df2plot, aes(x = year, y = exp(value), color = variable)) +
  geom_line() +
  theme(legend.position = "none") +
  ggtitle("But this is what it implies for the natural scale")

