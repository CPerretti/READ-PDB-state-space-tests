## Perform SAM simulation tests

# Debt:
# (x) Clean plot functions so it just takes simOut[[i]]
# (x) Fix second plot iteration
# ( ) Duplicate F's in output to match the config key
# (x) Moke the plots below into functions
# To do:
# ( ) Try parallelizing

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
nRep <- 10#0

# Generate simulation replicates
simOut <- list()
for (i in 1:nRep) {
  simOut[[i]] <- sim(fit = fitHer)  
}



## Plot an example true vs observed vs *herring* fit ########

## (1) N-at-age (1000s)
plotN(simOut = simOut[[1]],
      fit = fitHer)

## (2) F-at-age
plotF(simOut = simOut[[1]],
      fit = fitHer)

## (3) Catch (mt)
plotC(simOut = simOut[[1]],
      fit = fitHer)

## (4) Survey (1000s)
plotS(simOut = simOut[[1]],
      fit = fitHer)

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



## Plot example true vs observed vs fit to observed ########
## (1) N-at-age (1000s)
plotN(simOut = simOut[[1]],
      fit = fitSim[[1]])

## (2) F-at-age
plotF(simOut = simOut[[1]],
      fit = fitSim[[1]])

## (3) Catch (mt)
plotC(simOut = simOut[[1]],
      fit = fitSim[[1]])

## (4) Survey (1000s)
plotS(simOut = simOut[[1]],
      fit = fitSim[[1]])


## Plot fit vs true parameter values #######################

# Plot parmeters true vs fit
plotPars(fitSim, simOut)

# Plot error for N and F
err_logNF <- calcTsError(fitSim, simOut)
plotTsError(err_logNF)
  
         
        
# Fit sam to a simulate.sam replicate
# Need to resimulate with full.data = TRUE to get output that sam.fit() can use
# set.seed(123)
#simOut2fit <- stockassessment:::simulate.sam(fitHer, nsim = 10, full.data = TRUE)
#sam.fit(data = simOut2fit[[10]], conf = fitHer$conf, par = defpar(fitHer$data, fitHer$conf))

