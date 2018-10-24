# testing simulation of data

# Required packages
library(dplyr) 
library(ggplot2)
library(tidyr)
library(stockassessment)

# Use atl herring fit to set up simulation
load("../atlherring_example/output/fit.Rdata")

# use default configuration and parameter intitial guesses for this comparison
conf <- defcon(fit$data)
par <- defpar(fit$data, conf)

# fit data twice, once with sim.condRE = TRUE (default) and second time = FALSE
fitT <- sam.fit(fit$data, conf, par, sim.condRE = TRUE)
fitF <- sam.fit(fit$data, conf, par, sim.condRE = FALSE)

# check get same parameter estimates
fitT$pl$logN - fitF$pl$logN

# simulate data from both cases
simdatT <- stockassessment:::simulate.sam(fitT, nsim = 10, full.data = FALSE)
simdatF <- stockassessment:::simulate.sam(fitF, nsim = 10, full.data = FALSE)

# get the N at age from both cases
df_sim <- data.frame() 
for (i in 1:length(simdatT)) {
  df_sim <-
    rbind(df_sim, data.frame(replicate = as.factor(i),
                             year = fit$data$years,
                             NT = t(exp(simdatT[[i]]$logN)),
                             NF = t(exp(simdatF[[i]]$logN))))
}
df2plotSim <-
  df_sim %>%
  cbind(data.frame(fit = fit$pl$logN %>% exp %>% t)) %>% # original fit
  tidyr::gather(variable, value, -year, -replicate) %>%
  tidyr::separate(variable, into = c("variable", "age"))

# plot simulate.sam replicates
ggplot(df2plotSim %>% dplyr::filter(variable %in% c("NT")),
       aes(x = year, y = value, color = replicate)) +
  geom_line() + 
  facet_wrap(~age, scales = "free") +
  ylab("Abundance (1000's)") +
  ggtitle("sim.condRE = TRUE") # all 10 replicates identical
ggplot(df2plotSim %>% dplyr::filter(variable %in% c("NF")),
       aes(x = year, y = value, color = replicate)) +
  geom_line() + 
  facet_wrap(~age, scales = "free") +
  ylab("Abundance (1000's)") +
  ggtitle("sim.condRE = FALSE") # all 10 replicates wildly different - some blow up
ggplot(df2plotSim %>% dplyr::filter(replicate == 1, variable %in% c("fit", "NT")),
       aes(x = year, y = value, color = variable)) +
  geom_line() + 
  facet_wrap(~age, scales = "free") +
  ylab("Abundance (1000's)") +
  ggtitle("fit vs NT") # fit and NT are slightly different - why? had expected these to be identical

