# Compare summary statistics of simulations to the
# real data

# Required packages
library(plyr) # always load before dplyr
library(dplyr)
library(ggplot2)
library(tidyr)
library(stockassessment)
library(parallel)

# Load user-defined functions
source("../1-functions.R")

## Run simulation #########################################
# Use atl herring fit to set up simulation
# example_dir <- "atlherring_example"
# load(paste0("../", example_dir, "/output/fitHer.Rdata"))
# fitReal <- fitHer

# Use North Sea cod fit to set up simulation
example_dir <- "nscod_example"
load(paste0("../../", example_dir, "/fitNScod.Rdata"))
fitReal <- fitNScod

#set.seed(321) # for reproducibility

# How many simulation replicates to do

nRep <- 100

# Generate simulation replicates
simOut <- list()
for (i in 1:nRep) {
  simOut[[i]] <- sim(fit = fitReal)  
}


## Plot an example true vs observed vs *real data* fit ########

# ## (1) N-at-age (1000s)
# plotN(simOut = simOut[[1]],
#       fit = fitReal)
# 
# ## (2) F-at-age
# plotF(simOut = simOut[[1]],
#       fit = fitReal) 
# 
# ## (3) Catch (mt)
# plotC(simOut = simOut[[1]],
#       fit = fitReal)
# 
# ## (4) Survey (1000s)
# plotS(simOut = simOut[[1]],
#       fit = fitReal)


## Plot quantiles vs real data for all plots.
plotN_quants(simOut = simOut, fit = fitReal)

