# Replicate SAM model fit using the equations written here

# Required packages
library(plyr) # always load before dplyr
library(dplyr)
library(ggplot2)
library(tidyr)
library(stockassessment)

# Load user-defined functions
source("1-functions.R")

## Run simulation #########################################
# Use atl herring fit to set up simulation
load("../atlherring_example/output/fitHer.Rdata")

# How many simulation replicates to do
nRep <- 100

# Generate simulation replicates
simOut <- list()
for (i in 1:nRep) {
  simOut[[i]] <- sim(fit = fitHer)  
}



## Plo true vs observed vs *herring* fit #################

## (1) N-at-age (1000s)
# plotN(N = simOut$N, 
#       fit = fitHer)

## (2) Catch (mt)
# plotC(Cobs_mt = simOut$Cobs_mt, 
#       Ctru_mt = simOut$Ctru_mt, 
#       fit = fitHer)

## (3) Survey (1000s)
# plotS(Sobs_N = simOut$Sobs_N, 
#       Stru_N = simOut$Stru_N, 
#       fit = fitHer)

## Plot some simulations from simulate.sam()
#plotSimSAM(fitHer, nsim = 10, seed = NULL)


## Fit sam to a simulation ################################
fitSim <- list()
for (i in 1:nRep) {
  # Prep simulation data for read.ices()
  prepSimData(Sobs_N = simOut[[i]]$Sobs_N, 
              fit = fitHer, 
              Cobs_N = simOut[[i]]$Cobs_N) 
  
  # Read in data, set initial params and configuration
  setupOut <- setupModel()
  
  # Fit sam
  fitSim[[i]] <- sam.fit(setupOut$dat, setupOut$conf, 
                    setupOut$par, sim.condRE = FALSE)
}



## Plot true vs observed vs fit to observed ################
# (1) N-at-age (1000s)
# plotN(N = simOut$N, 
#       fit = fitSim)

# (2) Catch (mt)
# plotC(Cobs_mt = simOut$Cobs_mt, 
#       Ctru_mt = simOut$Ctru_mt, 
#       fit = fitSim)

# (3) Survey (1000s)
# plotS(Sobs_N = simOut$Sobs_N, 
#       Stru_N = simOut$Stru_N, 
#       fit = fitSim)


## Plot fit vs true parameter values #######################
fitSimAll <- data.frame()
for (i in 1:nRep) {
  fitSimAll <-
    rbind(fitSimAll,
        data.frame(variable = paste("sdLogN", 1:length(fitHer$pl[[4]]), sep = "."),
                     #paste(names(fitHer$pl)[[4]], 1:length(fitHer$pl[[4]]), sep = "."),
                   tru = 0.33 * exp(fitHer$pl[[4]]),
                   est = exp(fitSim[[i]]$pl[[4]]),
                   sd  = exp(fitSim[[i]]$plsd[[4]]),
                   replicate = i),
        data.frame(variable = paste(names(fitHer$pl)[[5]], 1:length(fitHer$pl[[5]]), sep = "."),
                   tru = fitHer$pl[[5]],
                   est = fitSim[[i]]$pl[[5]],
                   sd  = fitSim[[i]]$plsd[[5]],
                   replicate = i))
}

df2plot <-
  fitSimAll %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(tru = unique(tru),
                   est_mean = mean(est),
                   est_se = sd(est)/sqrt(nRep))

ggplot(df2plot, aes(y = variable)) +
  geom_point( aes(x = tru), color = "red", size  = 3) +
  geom_point( aes(x = est_mean), color = "blue", size = 3) +
  geom_errorbarh(aes(xmin = est_mean - 1.96*est_se,
                     xmax = est_mean + 1.96*est_se), height = 0.1, color = "blue") +
  theme_bw() +
  ylab("") +
  xlab("Value") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

# Fit sam to a simulate.sam replicate
# Need to resimulate with full.data = TRUE to get output that sam.fit() can use
# set.seed(123)
#simOut2fit <- stockassessment:::simulate.sam(fitHer, nsim = 10, full.data = TRUE)
#sam.fit(data = simOut2fit[[10]], conf = fitHer$conf, par = defpar(fitHer$data, fitHer$conf))

