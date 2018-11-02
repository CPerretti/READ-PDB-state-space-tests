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

# Generate a simulation replicate
simOut <- sim(fit = fitHer)


## Plo true vs observed vs *herring* fit #################

## (1) N-at-age (1000s)
plotN(N = simOut$N, 
      fit = fitHer)

## (2) Catch (mt)
plotC(Cobs_mt = simOut$Cobs_mt, 
      Ctru_mt = simOut$Ctru_mt, 
      fit = fitHer)

## (3) Survey (1000s)
plotS(Sobs_N = simOut$Sobs_N, 
      Stru_N = simOut$Stru_N, 
      fit = fitHer)

## Plot some simulations from simulate.sam()
#plotSimSAM(fitHer, nsim = 10, seed = NULL)


## Fit sam to a simulation ################################

# Prep simulation data for read.ices()
prepSimData(Sobs_N = simOut$Sobs_N, 
            fit = fitHer, 
            Cobs_N = simOut$Cobs_N) 

# Read in data, set initial params and configuration
setupOut <- setupModel()

# Fit sam
fitSim <- sam.fit(setupOut$dat, setupOut$conf, 
                  setupOut$par, sim.condRE = FALSE)


## Plot true vs observed vs fit to observed ################
# (1) N-at-age (1000s)
plotN(N = simOut$N, 
      fit = fitSim)

# (2) Catch (mt)
plotC(Cobs_mt = simOut$Cobs_mt, 
      Ctru_mt = simOut$Ctru_mt, 
      fit = fitSim)

# (3) Survey (1000s)
plotS(Sobs_N = simOut$Sobs_N, 
      Stru_N = simOut$Stru_N, 
      fit = fitSim)




# Fit sam to a simulate.sam replicate
# Need to resimulate with full.data = TRUE to get output that sam.fit() can use
# set.seed(123)
#simOut2fit <- stockassessment:::simulate.sam(fitHer, nsim = 10, full.data = TRUE)
#sam.fit(data = simOut2fit[[10]], conf = fitHer$conf, par = defpar(fitHer$data, fitHer$conf))

