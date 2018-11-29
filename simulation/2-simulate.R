## Perform SAM simulation tests

# Debt:
# - Duplicate F's in output to match the config key
# - In simulate code, move exp locations to places that make sense
# To do:



# Required packages
library(plyr) # always load before dplyr
library(dplyr)
library(ggplot2)
library(tidyr)
library(stockassessment)
library(parallel)

# Load user-defined functions
source("1-functions.R")

## Run simulation #########################################
# Use atl herring fit to set up simulation
load("../atlherring_example/output/fitHer.Rdata")

# How many simulation replicates to do
set.seed(321) # for reproducibility
nRep <- 100

# Generate simulation replicates
simOut <- list()
for (i in 1:nRep) {
  simOut[[i]] <- sim(fit = fitHer)  
}



## Plot an example true vs observed vs *herring* fit ########

## (1) N-at-age (1000s)
plotN(simOut = simOut[[1]],
      fit = fitHer)

plot(1:length(simOut[[1]]$logN[8,]), simOut[[1]]$logN[8,], type = "l")

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
setupOut <- list()
for (i in 1:nRep) {
  # Prep simulation data for read.ices()
  prepSimData(logSobs_N = simOut[[i]]$logSobs_N, 
              fit = fitHer, 
              logCobs_N = simOut[[i]]$logCobs_N) 
  
  # Read in data, set initial params and configuration
  setupOut[[i]] <- setupModel()
}


# test
# fitSim <- list()
# fitSim[[1]] <- sam.fit(setupOut[[1]]$dat, 
#                        setupOut[[1]]$conf, setupOut[[1]]$par)

# Fit model to replicates in parallel
cl <- makeCluster(detectCores() - 1) #setup nodes for parallel
clusterEvalQ(cl, {library(stockassessment)}) #load stockassessment to each node
fitSim <- parLapply(cl, setupOut, 
                    function(x){try(sam.fit(x$dat, x$conf, x$par))})
stopCluster(cl) #shut down nodes

## Error handling #####
# Exclude TMB fails
fitSimAccept <- fitSim[!(sapply(fitSim, class) == "try-error")]
simOutAccept <- simOut[!(sapply(fitSim, class) == "try-error")]
# Exclude non-convergences
x <- unlist(sapply(fitSimAccept, function (x) x[[6]][3]))
fitSimAccept <- fitSimAccept[x != 1]
simOutAccept <- simOutAccept[x != 1]  
nRepAccept <- length(fitSimAccept)

## Plot example true vs observed vs fit to observed ########
## (1) N-at-age (1000s)
plotN(simOut = simOutAccept[[1]],
      fit = fitSimAccept[[1]])

## (2) F-at-age
plotF(simOut = simOutAccept[[1]],
      fit = fitSimAccept[[1]])

## (3) Catch (mt)
plotC(simOut = simOutAccept[[1]],
      fit = fitSimAccept[[1]])

## (4) Survey (1000s)
plotS(simOut = simOutAccept[[1]],
      fit = fitSimAccept[[1]])


## Plot fit vs true parameter values #######################

# Plot parameters true vs fit
plotPars(fitSimAccept, simOutAccept)

# Plot error for N and F
err_logNF <- calcTsError(fitSimAccept, simOutAccept)
plotTsError(err_logNF)

err_logNF_mean <-
  err_logNF %>%
  dplyr::group_by(variable, age) %>%
  dplyr::summarise(error_mean = mean(error),
                   error_se   = sd(error) / sqrt(nRepAccept))
df2plot <-
  err_logNF %>%
  dplyr::left_join(err_logNF_mean)

ggplot(df2plot, aes(x = error)) +
  geom_vline(aes(xintercept = 0), color = "black") +
  geom_vline(aes(xintercept = error_mean), color = "blue") +
  geom_vline(aes(xintercept = error_mean - 1.96 * error_se), 
             color = "blue", linetype = "dashed") +
  geom_vline(aes(xintercept = error_mean + 1.96 * error_se), 
             color = "blue", linetype = "dashed") +
  facet_wrap(variable~age, scales = "free")
  

# Save output
suffix <- paste0(Sys.time(), ".Rdata")
save(list = "fitSim", file = paste0("./output/fitSim", suffix))
save(list = "fitSimAccept", file = paste0("./output/fitSimAccept", suffix))
save(list = "simOut", file = paste0("./output/simOut", suffix))
save(list = "simOutAccept", file = paste0("./output/simOutAccept", suffix))
        
# Fit sam to a simulate.sam replicate
# Need to resimulate with full.data = TRUE to get output that sam.fit() can use
# set.seed(123)
#simOut2fit <- stockassessment:::simulate.sam(fitHer, nsim = 10, full.data = TRUE)
#sam.fit(data = simOut2fit[[10]], conf = fitHer$conf, par = defpar(fitHer$data, fitHer$conf))

