## Perform SAM simulation tests
 
# Install github version of package to use default package
#devtools::install_github("fishfollower/SAM/stockassessment")
# Install local version of package with changes
#devtools::install_local("../../SAM/stockassessment/", force = TRUE)


# Required packages
library(plyr) # always load before dplyr
library(dplyr)
library(ggplot2)
library(tidyr)
library(stockassessment2) # with changes
library(stockassessment)  # default
library(parallel)

# Load user-defined functions
source("1-functions.R")


## Run simulation #########################################
# Use atl herring fit to set up simulation
# example_dir <- "atlherring_example"
# load(paste0("../", example_dir, "/output/fitHer.Rdata"))
# fitReal <- fitHer

# Use North Sea cod fit to set up simulation
example_dir <- "nscod_example"
load(paste0("../", example_dir, "/fitNScod.Rdata"))
fitReal <- fitNScod

fitReal$conf$constRecBreaks <- numeric(0) # Needed for new SAM

#set.seed(321) # for reproducibility
scenarios <- c("random", "fixed", "none")
nRep <- 10 # Number of simulation replicates
noScaledYears <- 10

# Study output container
container <- expand.grid(replicate = 1:nRep, scenario = scenarios, stringsAsFactors = F)
container$simOut          <- vector("list", length = nrow(container))
container$setupMod_random <- vector("list", length = nrow(container))
container$setupMod_fixed  <- vector("list", length = nrow(container))
container$setupMod_none   <- vector("list", length = nrow(container))


# Setup configurations for each model
keyLogScale <- fitReal$conf$keyLogFsta
keyLogScale[keyLogScale > -1] <- 0:(length(keyLogScale[keyLogScale > -1])-1)
nAs <- sum(keyLogScale[1,] > -1)

confLogScale_random <- 
  list(logScaleType = "random",
       keyLogScale = keyLogScale,
       keyVarLogScale = rep(0, nAs),
       noScaledYears = noScaledYears,
       keyScaledYears = (max(fitReal$data$years) - noScaledYears + 1) : 
         max(fitReal$data$years))

confLogScale_fixed <-
  list(logScaleType = "fixed",
       keyLogScale = keyLogScale,
       noScaledYears = noScaledYears,
       keyScaledYears = (max(fitReal$data$years) - noScaledYears + 1) : 
         max(fitReal$data$years),
       keyParScaledYA =  matrix(data = rep(0, noScaledYears * ncol(fitReal$data$propF)),
                                nrow = noScaledYears)) # One parameter all years

confLogScale_none <- list(logScaleType = "none")


# Generate simulation replicates
for (i in 1:nrow(container)) {
  container$simOut[[i]] <- sim(fit = fitReal,
                               keyLogScale = keyLogScale,
                               noScaledYears = noScaledYears,
                               container_i = container[i,])
}

## Plot an example true vs observed vs *real data* fit ########

# ## (1) N-at-age (1000s)
# plotN(simOut = container$simOut[[1]],
#       fit = fitReal)
# 
# ## (2) F-at-age
# plotF(simOut = container$simOut[[1]],
#       fit = fitReal)
# 
# ## (3) Catch (mt)
# plotC(simOut = container$simOut[[1]],
#       fit = fitReal)
# 
# ## (4) Survey (1000s)
# plotS(simOut = container$simOut[[1]],
#       fit = fitReal)

## Plot some simulations from simulate.sam()
#plotSimSAM(fitReal, nsim = 10, seed = NULL)


## Fit sam to a simulation ################################

for (i in 1:nrow(container)) {
  # Prep simulation data for read.ices()
  prepSimData(Sobs_N = container$simOut[[i]]$Sobs_N, 
              fit = fitReal, 
              Cobs_N = container$simOut[[i]]$Cobs_N) 
  
  # Read in data, set initial params and configuration
  container$setupMod_random[[i]] <- setupModel(conf = fitReal$conf, 
                                               example_dir = example_dir, 
                                               noScaledYears = noScaledYears,
                                               confLogScale = confLogScale_random)
  
  container$setupMod_fixed[[i]] <- setupModel(conf = fitReal$conf,
                                              example_dir = example_dir,
                                              noScaledYears = noScaledYears,
                                              confLogScale = confLogScale_fixed)
  
  container$setupMod_none[[i]] <- setupModel(conf = fitReal$conf,
                                             example_dir = example_dir,
                                             noScaledYears = noScaledYears,
                                             confLogScale = confLogScale_none)
  
  # setupOut[[i]]$par$logFpar <- fitReal$pl$logFpar # For debugging
  
}

# Fit model to replicates in parallel
# fitSimTest1 <- sam.fit(setupOut_fixed[[1]]$dat,
#                       setupOut_fixed[[1]]$conf,
#                       setupOut_fixed[[1]]$par)

# fitSimTest2 <- sam.fit_cp(setupOut_random[[1]]$dat,
#                           setupOut_random[[1]]$conf,
#                           setupOut_random[[1]]$par)


cl <- makeCluster(detectCores() - 1) #setup nodes for parallel <<<CONTIUE HERE<<<<<<
#load stockassessment and functions to each node
clusterEvalQ(cl, {library(stockassessment); source("1-functions.R")}) 

# Estimate with misreporting as random effect
fitSim_random <- parLapply(cl, 
                           setupOut_random, 
                           function(x){try(sam.fit_cp(x$dat, x$conf, x$par))})

# Estimate with misreporting as fixed effect
fitSim_fixed <- parLapply(cl, 
                          setupOut_fixed, 
                          function(x){try(sam.fit(x$dat, x$conf, x$par))})

# Assume no misreporting
fitSim_none <- parLapply(cl, 
                         setupOut_none, 
                         function(x){try(sam.fit(x$dat, x$conf, x$par))})

stopCluster(cl) #shut down nodes


## Error handling #####
# Exclude TMB fails
fitSimAccept_random <- fitSim_random[!(sapply(fitSim_random, class) == "try-error")]
simOutAccept_random <- simOut[!(sapply(fitSim_random, class) == "try-error")]
fitSimAccept_fixed <- fitSim_fixed[!(sapply(fitSim_fixed, class) == "try-error")]
simOutAccept_fixed <- simOut[!(sapply(fitSim_fixed, class) == "try-error")]
fitSimAccept_none <- fitSim_none[!(sapply(fitSim_none, class) == "try-error")]
simOutAccept_none <- simOut[!(sapply(fitSim_none, class) == "try-error")]
# Exclude non-convergences
x <- unlist(sapply(fitSimAccept_random, function (x) x[[6]][3]))
fitSimAccept_random <- fitSimAccept_random[x != 1]
simOutAccept_random <- simOutAccept_random[x != 1]  
nRepAccept_random <- length(fitSimAccept_random)
x <- unlist(sapply(fitSimAccept_fixed, function (x) x[[6]][3]))
fitSimAccept_fixed <- fitSimAccept_fixed[x != 1]
simOutAccept_fixed <- simOutAccept_fixed[x != 1]
nRepAccept_fixed <- length(fitSimAccept_fixed)
x <- unlist(sapply(fitSimAccept_none, function (x) x[[6]][3]))
fitSimAccept_none <- fitSimAccept_none[x != 1]
simOutAccept_none <- simOutAccept_none[x != 1]
nRepAccept_none <- length(fitSimAccept_none)

# Save output
# suffix <- paste0(Sys.time(), ".Rdata")
# save(list = "fitSim", file = paste0("./output/fitSim", suffix))
# save(list = "fitSimAccept", file = paste0("./output/fitSimAccept", suffix))
# save(list = "simOut", file = paste0("./output/simOut", suffix))
# save(list = "simOutAccept", file = paste0("./output/simOutAccept", suffix))

## Plot example true vs observed vs fit to observed ########
## (1) N-at-age (1000s)
plotN(simOut = simOutAccept_random[[1]],
      fit = fitSimAccept_random[[1]])

## (2) F-at-age
plotF(simOut = simOutAccept_random[[1]],
      fit = fitSimAccept_random[[1]])

## (3) Catch (mt)
plotC(simOut = simOutAccept_random[[1]],
      fit = fitSimAccept_random[[1]])

## (4) Survey (1000s)
plotS(simOut = simOutAccept_random[[1]],
      fit = fitSimAccept_random[[1]])


## Plot fit vs true parameter values #######################

# Plot error of time series estimates
errRe_random    <- calcReTsError(fitSimAccept_random, 
                                 simOutAccept_random,
                                 confLogScale_random)
errCSSB_random  <- calcCSSBError(fitSimAccept_random, simOutAccept_random)
err_random      <- rbind(errRe_random, errCSSB_random)

errRe_fixed <- calcReTsError(fitSimAccept_fixed, 
                             simOutAccept_fixed, 
                             confLogScale_fixed)
errCSSB_fixed  <- calcCSSBError(fitSimAccept_fixed, simOutAccept_fixed)
err_fixed <- rbind(errRe_fixed, errCSSB_fixed)

errRe_none <- calcReTsError(fitSimAccept_none, 
                            simOutAccept_none,
                            confLogScale_none)
errCSSB_none <- calcCSSBError(fitSimAccept_none, simOutAccept_none)
err_none <- rbind(errRe_none, errCSSB_none)

plotTsError(err_random, noYears = fitSimAccept_random[[1]]$data$noYears)
plotTsError(err_fixed, noYears = fitSimAccept_fixed[[1]]$data$noYears)
plotTsError(err_none, noYears = fitSimAccept_none[[1]]$data$noYears)

# plotTsMeanError(err_random, nRepAccept_random)
# plotTsMeanError(err_fixed, nRepAccept_fixed)
# plotTsMeanError(err_none, nRepAccept_none)

# Plot parameters true vs fit
plotPars(fitSimAccept_random, simOutAccept_random)
plotPars(fitSimAccept_fixed, simOutAccept_fixed)
plotPars(fitSimAccept_none, simOutAccept_none)


