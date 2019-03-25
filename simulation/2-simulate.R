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
scenarios <- c("random", "fixed", "none") # Right now "fixed" is uniform random <<<
nRep <- 4#10 # Number of simulation replicates
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
# fitSimTest1 <- sam.fit(container$setupMod_fixed[[1]]$dat,
#                        container$setupMod_fixed[[1]]$conf,
#                        container$setupMod_fixed[[1]]$par)

# fitSimTest2 <- sam.fit_cp(container$setupMod_random[[1]]$dat,
#                           container$setupMod_random[[1]]$conf,
#                           container$setupMod_random[[1]]$par)


cl <- makeCluster(detectCores() - 1) #setup nodes for parallel
#load stockassessment and functions to each node
clusterEvalQ(cl, {library(stockassessment); source("1-functions.R")}) 

# Estimate with misreporting as random effect
container$fitSim_random <- parLapply(cl, 
                                     container$setupMod_random, 
                                     function(x){try(sam.fit_cp(x$dat, x$conf, x$par))})

# Estimate with misreporting as fixed effect
container$fitSim_fixed <- parLapply(cl, 
                                    container$setupMod_fixed, 
                                    function(x){try(sam.fit(x$dat, x$conf, x$par))})

# Assume no misreporting
container$fitSim_none <- parLapply(cl, 
                                   container$setupMod_none, 
                                   function(x){try(sam.fit(x$dat, x$conf, x$par))})

stopCluster(cl) #shut down nodes


## Error handling #####
# Only include replicates where all three models fit successfully
containerAccept <- # Exclude TMB fails
  container[sapply(container$fitSim_random, class) != "try-error" &
            sapply(container$fitSim_fixed, class)  != "try-error" &
            sapply(container$fitSim_none, class)   != "try-error", ]

containerAccept <- # Also exclude non-covergences
  containerAccept[unlist(sapply(containerAccept$fitSim_random, 
                                function (x) x[[6]][3])) != 1 &
                  unlist(sapply(containerAccept$fitSim_fixed, 
                                function (x) x[[6]][3])) != 1 &
                  unlist(sapply(containerAccept$fitSim_none, 
                                function (x) x[[6]][3])) != 1, ]

# Save output
# suffix <- paste0(Sys.time(), ".Rdata")
# save(list = "container", file = paste0("./output/container", suffix))
# save(list = "containerAccept", file = paste0("./output/containerAccept", suffix))

## Plot example true vs observed vs fit to observed ########
## (1) N-at-age (1000s)
plotN(simOut = containerAccept$simOut[[1]],
      fit = containerAccept$fitSim_random[[1]])

## (2) F-at-age
plotF(simOut = containerAccept$simOut[[1]],
      fit = containerAccept$fitSim_random[[1]])

## (3) Catch (mt)
plotC(simOut = containerAccept$simOut[[1]],
      fit = containerAccept$fitSim_random[[1]])

## (4) Survey (1000s)
plotS(simOut = containerAccept$simOut[[1]],
      fit = containerAccept$fitSim_random[[1]])


## Plot fit vs true parameter values #######################

# Calculate fit error
errRe_random    <- calcReTsError(containerAccept$fitSim_random[[1]], 
                                 containerAccept$simOut[[1]],
                                 confLogScale_random)

# <<< CONTIUE HERE by combining both types of errors and all models into a single error matrix
containerAccept$errRe_random <- vector("list", length = nrow(containerAccept))
for (i in 1:length(containerAccept)) {
  containerAccept[i, "errRe_random"] <- calcReTsError(containerAccept$fitSim_random[[i]], 
                                                      containerAccept$simOut[[i]],
                                                      confLogScale_random)  
}


errCSSB_random  <- calcCSSBError(containerAccept$fitSim_random, 
                                 containerAccept$simOut)
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

# Plot time series error
plotTsError(err_random, noYears = containerAccept$fitSim_random[[1]]$data$noYears)
plotTsError(err_fixed, noYears = fitSimAccept_fixed[[1]]$data$noYears)
plotTsError(err_none, noYears = fitSimAccept_none[[1]]$data$noYears)

# plotTsMeanError(err_random, nRepAccept_random)
# plotTsMeanError(err_fixed, nRepAccept_fixed)
# plotTsMeanError(err_none, nRepAccept_none)

# Plot parameters true vs fit
plotPars(fitSimAccept_random, simOutAccept_random)
plotPars(fitSimAccept_fixed, simOutAccept_fixed)
plotPars(fitSimAccept_none, simOutAccept_none)


