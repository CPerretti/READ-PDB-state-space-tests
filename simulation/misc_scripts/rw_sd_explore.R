
# load required libraries
library(ggplot2)
library(TMB)
library(dplyr)

set.seed(123) # for reproducibility

n_t <- 1e6 # length of time series

# Create rw
e <- vector(mode = "numeric", length = n_t) # vector for errors
sigma <- 1 # sd of white noise component
e[1] <- rnorm(1, sd = sigma) # error in first year
phi <- 1 # autocorrelation parameter
for (i in 2:n_t) e[i] <- phi * e[i-1] + rnorm(1, sd = sigma) # create errors

Year <- 1:length(e)
# Plot it
ggplot(data.frame(Year = Year, e = e),
       aes(x = Year)) +
  geom_line(aes(y = e)) +
  theme_bw() +
  theme(axis.title = element_text(size = 14))

print(empirical_sd <- sd(e))
print(analytical_sd <- sigma / (1 - phi^2)^0.5)
      
      
      