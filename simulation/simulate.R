# See if we can replicate the simulate feature for a fitted model 
# (start simple)

# Use atl herring model fit to set up simulation
load("../atlherring_example/output/fit.Rdata")


# Run simulation
nA <- ncol(fit$data$propF) # number of age-classes
nT <- fit$data$noYears # length of time series

N <- matrix(data = NA, # N container
            nrow = nA, 
            ncol = nT, 
            dimnames = list(c(1:nA), NULL)) 

N[, 1] <-  exp(fit$pl$logN[,1]) # initial condition from fitted model

# Set F (need to replicate some elements to match ModelConf)
f <- exp(fit$pl$logF[(fit$conf$keyLogFsta[1,] + 1),])
# Set M
m <- t(fit$data$natMor)
# Set N process sd
sdN <- exp(fit$pl$logSdLogN[(fit$conf$keyVarLogN + 1)])

# simulate
for (i in 2:nT) {
  N[1, i] <- N[1, i-1] * exp(rnorm(1))
  N[-c(1, nA), i] <- N[-c(nA-1, nA), i-1] *
                      exp(-f[-c(nA-1, nA), i-1]) * 
                      exp(-m[-c(nA-1, nA), i-1]) * exp(rnorm(1, sd = ))
  N[nA, i] <- N[nA-1, i-1] * 
               exp(-f[nA-1, i-1]) *  
               exp(-m[nA-1, i-1]) +
               N[nA, i-1] * 
               exp(-f[nA, i-1]) *  
               exp(-m[nA, i-1]) * exp(rnorm(1))
  
}



# Plot abundance-at-age
library(dplyr)
df2plot <- # setup data for plot
  N %>%
  t() %>%
  as.data.frame() %>%
  dplyr::mutate(year = fit$data$years) %>%
  tidyr::gather(age, N, -year)

library(ggplot2) # plot it
ggplot(data = df2plot, aes(x = year, y = N)) +
  geom_point() +
  geom_line() +
  facet_wrap(~age, scales = "free")
