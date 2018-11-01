# Replicate SAM model fit using the equations written here

# Required packages
library(plyr) # always load before dplyr
library(dplyr)
library(ggplot2)
library(tidyr)
library(stockassessment)

# Load user-defined functions
source("1-functions.R")

## Set up simulation parameters ###########################

# Use atl herring fit to set up simulation
load("../atlherring_example/output/fit.Rdata")

# Set F (need to replicate some elements to match ModelConf)
f <- exp(fit$pl$logF[(fit$conf$keyLogFsta[1,] + 1),])
# Set M
m <- t(fit$data$natMor)

# Calcuate total mortality
z <- f + m

# Set up matrix to record N-at-age
nA <- ncol(fit$data$propF) # number of age-classes
nT <- fit$data$noYears # length of time series

logN <- matrix(data = NA, # N container
               nrow = nA, 
               ncol = nT, 
               dimnames = list(paste0("simulated.", c(1:nA)), fit$data$years)) 

logN[, 1] <- fit$pl$logN[,1] # initial condition from fitted model

# Calculate the process errors that were estimated in the fit so we can 
# exactly replicate the fit
errPro_exact <- matrix(data = NA,
                       nrow = nA, 
                       ncol = nT-1)

errPro_exact[1, ] <- fit$pl$logN[1, 2:nT] - fit$pl$logN[1, 1:(nT-1)]
errPro_exact[-c(1, nA), ] <-  fit$pl$logN[-c(1, nA), 2:nT] -
                              (fit$pl$logN[-c(nA-1, nA), 1:(nT-1)] -
                              z[-c(nA-1, nA), 1:(nT-1)])
errPro_exact[nA, ] <- fit$pl$logN[nA, 2:nT] -
                      log(exp(fit$pl$logN[nA-1, 1:(nT-1)]) *
                      exp(-z[nA-1, 1:(nT-1)]) +
                      exp(fit$pl$logN[nA, 1:(nT-1)]) *
                      exp(-z[nA, 1:(nT-1)]))
# Calculate a new relization of process errors with sd's from the fit
# Set N process sd
errPro <- matrix(data = NA,
                 nrow = nA, 
                 ncol = nT-1)
sdLogN <- exp(fit$pl$logSdLogN[(fit$conf$keyVarLogN + 1)])
for (i in 1:(nT-1)) { # Create process error (N-at-age)
  errPro[, i] <-  rnorm(n = nA, sd = sdLogN)
}

#errPro <- errPro_exact # Use if you want the fit N

## Simulate process model #################################
# N fit
# Simulate N-at-age
for (i in 2:nT) {
  logN[1, i] <- logN[1, i-1] + errPro[1, i-1]
  logN[-c(1, nA), i] <- logN[-c(nA-1, nA), i-1] - 
                        z[-c(nA-1, nA), i-1] + 
                        errPro[-c(1, nA), i-1]
  logN[nA, i] <- log(exp(logN[nA-1, i-1]) * exp(-z[nA-1, i-1]) +
              exp(logN[nA, i-1]) * 
              exp(-z[nA, i-1])) + errPro[nA, i-1]
}

N <- exp(logN)

## Simulate observation model #############################
# Observation error (for catch and surveys)
errObs <- array(data = NA, # error container (3-d: age x year x fleet)
                dim = c(nA, nT, fit$data$noFleets),
                dimnames = list(paste0("error.", c(1:nA)), 
                                fit$data$years, 
                                attr(fit$data,"fleetNames")))
# Need to replicate some sd's to match config file
index <- as.vector(t(fit$conf$keyVarObs + 1))
index[index == 0] <- NA
sdLogObs <- exp(fit$pl$logSdLogObs[index])

# Make observation error (can only do uncorrelated error right now)
  for (j in 1:nT) { # all surveys in a year
    errObs[, j, ] <- rnorm(n = length(sdLogObs), sd = sdLogObs)
  }

# Check that sds are correctly assigned to each age x survey combo
# sdInput <- c()
# for (i in 1:fit$data$noFleets) {
#   sdInput <- c(sdInput, apply(errObs[, , i], 1, sd))  
# }
# 
# cbind(sdInput, sdLogObs)


# Simulate Catch

# Calculate catch on N-scale (1000s)
logCtru_N <- log(f / z * (1 - exp(-z)) * N)
logCobs_N <- logCtru_N + errObs[, , "Residual catch"]

Ctru_N <- exp(logCtru_N)
Cobs_N <- exp(logCobs_N)

# Convert to MT (1000s * kg = mt)
Ctru_mt <- Ctru_N * t(fit$data$catchMeanWeight)
Cobs_mt <- Cobs_N * t(fit$data$catchMeanWeight)


# Simulate Survey (SIMULATIONS HAVE DATA IN ALL YEARS RIGHT NOW)
logS_N <- array(data = NA, # survey container (3-d: age x year x survey)
              dim = c(nA, nT, fit$data$noFleets),
              dimnames = list(paste0("simulated.", c(1:nA)), 
                              fit$data$years, 
                              attr(fit$data,"fleetNames")))

logSq <- matrix(data = NA, # survey q-at-age matrix
                nrow = nrow(fit$conf$keyLogFpar), 
                ncol = ncol(fit$conf$keyLogFpar))
logSq[which(fit$conf$keyLogFpar != -1)] <- # fill with herring fit values
  fit$pl$logFpar[fit$conf$keyLogFpar[fit$conf$keyLogFpar != -1] + 1]
Sq <- exp(logSq)

surveyIndex <- # some fleets are fishermen not surveys
  (1:fit$data$noFleets)[fit$data$fleetTypes == 2] 
for (i in surveyIndex) {
logS_N[, , i] <- log(Sq[i,] * exp(-z * fit$data$sampleTimes[i]) * N[,]) +
                errObs[, , i]
}

S_N <- exp(logS_N)


## Make plots for comparison of simulated vs original fit ##################
# (1) N-at-age
plotN(N, fit)

# (2) Catch (mt)
plotC(Cobs_mt, Ctru_mt, fit)

# (3) Survey #<<<< CHANGE THIS PLOT SO IT HAS TRUE AND OBSERVED AND CHECK UNITS>>>>>>
plotS(S_N, fit)

## Plot some simulations from simulate.sam()
#plotSimSAM(fit, nsim = 10, seed = NULL)


## Fit sam to a simulation ################################

# Prep simulation data for read.ices()
prepSimData(S_N, fit, Cobs_N) 

# Read in data
cn <- read.ices("./sim_data/catch.dat") # catch abundace-at-age
cw <- read.ices("../atlherring_example/data/Herrcw.dat") # catch mean weight-at-age
dw <- cw # discards mean weight-at-age (using catch mean weight-at-age for now)
lw <- cw # landings mean weight-at-age (using catch mean weight-at-age for now)
pf <- read.ices("../atlherring_example/data/Herrpf.dat") # proportion of f before spawning
lf <- pf; lf[,] <- 1 # fraction of catch that is landed (set to 1 for now)
mo <- read.ices("../atlherring_example/data/Herrmo.dat") # maturity-at-age ogive
nm <- read.ices("../atlherring_example/data/Herrnm.dat") # natural mortality-at-age
pm <- read.ices("../atlherring_example/data/Herrpm.dat") # proportion of m before spawning
sw <- read.ices("../atlherring_example/data/Herrsw.dat") # stock weight-at-age (kg)
surveys <- read.ices("./sim_data/surveys.dat") #surveys
#surveys <- read.ices("../atlherring_example/data/Herrsurvey_BigSep_NoAcoust.dat") #surveys

# setup the data as needed for SAM
dat <- setup.sam.data(surveys = surveys,
                      residual.fleet = cn,
                      prop.mature = mo,
                      stock.mean.weight = sw,
                      dis.mean.weight = dw,
                      land.mean.weight = lw,
                      land.frac = lf,
                      prop.f = pf,
                      prop.m = pm,
                      natural.mortality = nm,
                      catch.mean.weight = cw)

# Load model configuration file
conf <- loadConf(dat = dat, file = "../atlherring_example/ModelConf_original.txt")

par <- defpar(dat, conf) # some default starting values

fitSim <- sam.fit(dat, conf, par, sim.condRE = FALSE) # fit the model

## Plot simulation vs fit to simulation
# (1) N-at-age
plotN(N, fitSim)

# (2) Catch
plotC(Cobs_mt, Ctru_mt, fitSim)
catchplot(fitSim)

# (3) Survey
plotS(S_N, fitSim)





#fit2sim <- sam.fit(data = fit$data, conf = fit$conf, 
#                   par = defpar(fit$data, fit$conf))

# Try to fit sam to one of the simulate.sam replicates
# Need to resimulate with full.data = TRUE to get output that sam.fit() can use
# set.seed(123)
# simOut2fit <- stockassessment:::simulate.sam(fit, nsim = 10, full.data = TRUE)
# sam.fit(data = simOut2fit[[1]], conf = fit$conf, par = defpar(fit$data, fit$conf))

