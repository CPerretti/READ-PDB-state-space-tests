# Replicate SAM model fit using the equations written here

# Required packages
library(plyr) # always load before dplyr
library(dplyr)
library(ggplot2)
library(tidyr)
library(stockassessment)

# Use atl herring fit to set up simulation
load("../atlherring_example/output/fit.Rdata")


# Set up simulation parameters

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

errPro <- errPro_exact # Use if you want the fit N

## Process model ##################################################
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

## Observation model ##############################################
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


# Catch fit
logC <- matrix(data = NA, # Catch container
               nrow = nA, 
               ncol = nT, 
               dimnames = list(paste0("simulated.", c(1:nA)), 
                               fit$data$years))

 
# Calculate catch index
logC[,] <- log(f / z * (1 - exp(-z)) * N[,] * 
           t(fit$data$catchMeanWeight)) + errObs[, , "Residual catch"]

C <- exp(logC)

# Survey fits (SIMULATIONS HAVE DATA IN ALL YEARS RIGHT NOW i.e. SummerNMFS has
# data from 1965-2017)
logS <- array(data = NA, # survey container (3-d: age x year x survey)
              dim = c(nA, nT, fit$data$noFleets),
              dimnames = list(paste0("simulated.", c(1:nA)), 
                              fit$data$years, 
                              attr(fit$data,"fleetNames")))

logSq <- matrix(data = NA, # survey q-at-age matrix
                nrow = nrow(fit$conf$keyLogFpar), 
                ncol = ncol(fit$conf$keyLogFpar))
logSq[which(fit$conf$keyLogFpar != -1)] <- # fill with fit values
  fit$pl$logFpar[fit$conf$keyLogFpar[fit$conf$keyLogFpar != -1] + 1]
Sq <- exp(logSq)

surveyIndex <- # some fleets are fishermen not surveys
  (1:fit$data$noFleets)[fit$data$fleetTypes == 2] 
for (i in surveyIndex) {
logS[, , i] <- log(Sq[i,] * exp(-z * fit$data$sampleTimes[i]) * N[,]) +
                errObs[, , i]
}

S <- exp(logS)


## Make plots for comparison of simulated vs original fit ##################
## (1) N-at-age
# Setup N-at-age to plot
df2plotN <- 
  N %>%
  t() %>%
  as.data.frame() %>%
  cbind(data.frame(fit = fit$pl$logN %>% exp %>% t)) %>%
  dplyr::mutate(year = fit$data$years) %>%
  tidyr::gather(variable, N, -year) %>%
  tidyr::separate(variable, c("source", "age"))

# Plot N-at-age (all ages should match exactly when using errPro_exact)
ggplot(data = df2plotN,
       aes(x = year, y = N, color = source)) +
  geom_line() +
  facet_wrap(~age, scales = "free") +
  ylab("Abundance (1000's)")

## (2) Catch
# Setup Catch to plot
df_fitCatch <- 
  data.frame(variable = names(fit$sdrep$value),
             value = fit$sdrep$value) %>%
  dplyr::filter(variable == "logCatch") %>%
  dplyr::rename(logCatch = value) %>%
  dplyr::mutate(fit.total = exp(logCatch),
                year = fit$data$years) %>%
  dplyr::select(year, fit.total)
  
df2plotC <- 
  C %>%
  t() %>%
  as.data.frame() %>%
  dplyr::mutate(simulated.total = rowSums(.),
                year = fit$data$years) %>%
  dplyr::left_join(df_fitCatch) %>%
  tidyr::gather(variable, Catch, -year) %>%
  tidyr::separate(variable, c("source", "age"))

# Plot Catch (should exactly match in total subplot when using errPro_exact)
p <- 
  ggplot(data = df2plotC,
       aes(x = year, y = Catch, color = source)) +
  geom_line() +
  facet_wrap(~age, scales = "free") +
  ylab("Catch (mt)")

p

## (3) Survey
# Set up Survey data to plot
df_simS <- # convert 3-d array to long data.frame
  as.data.frame(S) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::mutate(year_fleet = row.names(.)) %>%
  tidyr::separate(col = year_fleet, 
                  into = c("year", "fleetNames"), 
                  extra = "merge") %>%
  tidyr::gather(variable, value, -year, -fleetNames) %>%
  tidyr::separate(variable, c("source", "age")) %>%
  dplyr::mutate(year = as.numeric(year),
                age = as.numeric(age),
                observed = NA)


# Organize SAM survey fits
fleets <- unique(fit$data$aux[,"fleet"])
idx <- fit$data$aux[,"fleet"]%in%fleets 
df_fitS <-
  data.frame(year = fit$data$aux[idx,"year"], 
             observed = exp(fit$data$logobs[idx]), 
             value = exp(fit$obj$report(c(fit$sdrep$par.fixed,
                                        fit$sdrep$par.random))$predObs[idx]), 
             age = fit$data$aux[idx,"age"], 
             fleetNames = attr(fit$data,"fleetNames")[fit$data$aux[idx,"fleet"]],
             source = "fit")

df2plotS <- # combine simulated Survey and fit Survey
  bind_rows(df_simS, df_fitS)

# Plot Survey (should match exactly when using errPro_exact)
p <-
  ggplot(data = df2plotS %>% dplyr::filter(age>1, fleetNames != "Residual catch"),
       aes(x = year)) +
  geom_point(aes(y = observed)) +
  geom_line(aes(y = value, color = source)) +
  facet_wrap(age~fleetNames, scales = "free_y", nrow = 7)

p


## Compare my simulation to simulate.sam() ####
# SAM simulate feature
set.seed(123) # for reproduciblilty and to match with full.data = TRUE (below)
simOut <- stockassessment:::simulate.sam(fit, nsim = 10, full.data = FALSE)

df_simOut <- data.frame() 
for (i in 1:length(simOut)) {
  df_simOut <-
    rbind(df_simOut, data.frame(replicate = as.factor(i),
                              year = fit$data$years,
                              N = t(exp(simOut[[i]]$logN)),
                              F = t(exp(simOut[[i]]$logF))))
}
df2plotSimOut <-
  df_simOut %>%
  cbind(data.frame(fit = fit$pl$logN %>% exp %>% t)) %>% # original fit
  tidyr::gather(variable, value, -year, -replicate) %>%
  tidyr::separate(variable, into = c("variable", "age"))

# plot simulate.sam replicates
p <-
  ggplot(df2plotSimOut %>% dplyr::filter(variable %in% c("N")),
       aes(x = year, y = value, group = replicate)) +
  geom_line() + 
  facet_wrap(~age, scales = "free") +
  ylab("Abundance (1000's)")
#p


# Fit sam to a simulation (UNDER CONSTRUCTION!)
# Convert survey data to SAM input format
surveys <-
  lapply(X = plyr::alply(S[, , surveyIndex], 3, .dims = TRUE), t) %>%
  lapply(function(x) {colnames(x) <- sub(colnames(x), # Change colnames
                                        pattern = "simulated.", 
                                        replacement = ""); x}) %>%
  lapply(function(x) x[,-which(colSums(x, na.rm = TRUE) == 0)]) # Remove ages not in survey

for(i in 1:length(surveys)) { # set survey times to match real data
  attr(surveys[[i]], "time", fit$data$sampleTimes[surveyIndex[i]] + c(-0.25, 0.25))
}

# Export the simulated surveys to a text file, then import it
# using read.ices().
# First set file header
write.table(rbind("US Atlantic Herring Survey Data", 100 + length(surveys2)), 
            file = "surveys.dat", sep = " ", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
# Then append surveys
for (i in 1:length(surveys2)) { # loop over surveys
  min_year <- min(as.numeric(rownames(surveys2[[i]])))
  max_year <- max(as.numeric(rownames(surveys2[[i]])))
  header_survey <-
    rbind(names(surveys2)[i], 
          paste(min_year, max_year),
          paste(1, 1, fit$data$sampleTimes[surveyIndex[i]] - 0.25, 
                      fit$data$sampleTimes[surveyIndex[i]] + 0.25),
          paste(min(colnames(surveys2[[i]])), max(colnames(surveys2[[1]]))))
  a_survey <- cbind(1, as.data.frame(surveys2[[i]]))
  
  write.table(header_survey, file = "surveys.dat", sep = " ", 
              row.names = FALSE, col.names = FALSE, quote = FALSE,
              append = TRUE)
  
  write.table(a_survey, file = "surveys.dat", sep = " ", 
              row.names = FALSE, col.names = FALSE, quote = FALSE,
              append = TRUE)
}

surveys <- read.ices("surveys.dat")
#fit2sim <- sam.fit(data = fit$data, conf = fit$conf, 
#                   par = defpar(fit$data, fit$conf))

# Try to fit sam to one of the simulate.sam replicates
# Need to resimulate with full.data = TRUE to get output that sam.fit() can use
# set.seed(123)
# simOut2fit <- stockassessment:::simulate.sam(fit, nsim = 10, full.data = TRUE)
# sam.fit(data = simOut2fit[[1]], conf = fit$conf, par = defpar(fit$data, fit$conf))

