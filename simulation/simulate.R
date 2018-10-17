# See if we can replicate the simulate feature for a fitted model 
# (start simple)

# Use atl herring model fit to set up simulation
load("../atlherring_example/output/fit.Rdata")


# Set up simulation parameters

# Set F (need to replicate some elements to match ModelConf)
f <- exp(fit$pl$logF[(fit$conf$keyLogFsta[1,] + 1),])
# Set M
m <- t(fit$data$natMor)
# Set N process sd
sdN <- exp(fit$pl$logSdLogN[(fit$conf$keyVarLogN + 1)])

# Set up matrix to record N-at-age
nA <- ncol(fit$data$propF) # number of age-classes
nT <- fit$data$noYears # length of time series

N <- matrix(data = NA, # N container
            nrow = nA, 
            ncol = nT, 
            dimnames = list(paste0("simulated.", c(1:nA)), fit$data$years)) 

N[, 1] <-  exp(fit$pl$logN[,1]) # initial condition from fitted model

# Calculate the process errors that were estimated in the fit so we can 
# exactly replicate the fit
errPro <- matrix(data = NA,
               nrow = nA, 
               ncol = nT-1)
errPro[1, ] <- exp(fit$pl$logN[1, 2:nT]) / exp(fit$pl$logN[1, 1:(nT-1)])
errPro[-c(1, nA), ] <-  exp(fit$pl$logN[-c(1, nA), 2:nT]) / 
                       (exp(fit$pl$logN[-c(nA-1, nA), 1:(nT-1)]) *
                        exp(-f[-c(nA-1, nA), 1:(nT-1)]) * 
                        exp(-m[-c(nA-1, nA), 1:(nT-1)]))
errPro[nA, ] <- exp(fit$pl$logN[nA, 2:nT]) / 
                (exp(fit$pl$logN[nA-1, 1:(nT-1)]) *
                   exp(-f[nA-1, 1:(nT-1)]) *  
                   exp(-m[nA-1, 1:(nT-1)]) +
                   exp(fit$pl$logN[nA, 1:(nT-1)]) *
                   exp(-f[nA, 1:(nT-1)]) *  
                   exp(-m[nA, 1:(nT-1)]))

#### Process model ####
# Simulate using process errors from fit (should match exactly)
for (i in 2:nT) {
  N[1, i] <- N[1, i-1] * errPro[1, i-1] #* exp(rnorm(1, sd = sdN[1]))
  N[-c(1, nA), i] <- N[-c(nA-1, nA), i-1] *
                     exp(-f[-c(nA-1, nA), i-1]) * 
                     exp(-m[-c(nA-1, nA), i-1]) * 
                     errPro[-c(1, nA), i-1]
                     #exp(rnorm(nA-2, sd = sdN[-c(1, nA)]))
  N[nA, i] <- (N[nA-1, i-1] *
              exp(-f[nA-1, i-1]) *  
              exp(-m[nA-1, i-1]) +
              N[nA, i-1] * 
              exp(-f[nA, i-1]) *  
              exp(-m[nA, i-1])) *
              errPro[nA, i-1]
              #exp(rnorm(1, sd = sdN[nA]))
}


# Plot process N-at-age
library(dplyr) # setup data for plot
df2plotN <- 
  N %>%
  t() %>%
  as.data.frame() %>%
  cbind(data.frame(fit = fit$pl$logN %>% exp %>% t)) %>%
  dplyr::mutate(year = fit$data$years) %>%
  tidyr::gather(variable, N, -year) %>%
  tidyr::separate(variable, c("source", "age"))

library(ggplot2) # plot it (all ages should match exactly)
ggplot(data = df2plotN,
       aes(x = year, y = N, color = source)) +
  geom_line() +
  facet_wrap(~age, scales = "free") +
  ylab("Abundance (1000's)")


#### Observation model ####
# Calculate catch index

C <- matrix(data = NA, # Catch container
            nrow = nA, 
            ncol = nT, 
            dimnames = list(paste0("simulated.", c(1:nA)), fit$data$years))


  

C[,] <- f/(f + m) * (1 - exp(-(f + m))) * N[,] * # exp(errObs[]) *
        t(fit$data$catchMeanWeight)

# Calculate survey indices (under construction!)

# Survey container (3-d: age x year x survey)
S <- array(data = NA, 
           dim = c(nA, nT, fit$data$noFleets),
           dimnames = list(paste0("simulated.", c(1:nA)), 
                           fit$data$years, 
                           attr(fit$data,"fleetNames")))

logSq <- matrix(data = NA, # Survey q-at-age matrix
                nrow = nrow(fit$conf$keyLogFpar), 
                ncol = ncol(fit$conf$keyLogFpar))
logSq[which(fit$conf$keyLogFpar != -1)] <- # fill with fit values
  fit$pl$logFpar[fit$conf$keyLogFpar[fit$conf$keyLogFpar != -1] + 1]
Sq <- exp(logSq)

surveyIndex <- # some fleets are fishermen not surveys
  (1:fit$data$noFleets)[fit$data$fleetTypes == 2] 

for (i in surveyIndex) {
S[, , i] <- Sq[i,] * 
            exp(-(f + m) * fit$data$sampleTimes[i]) * N[,] #* exp(errObs[]) 
}


# Extract survey fits and observations (under construction!)
fleets <- unique(fit$data$aux[,"fleet"])
idx <- fit$data$aux[,"fleet"] %in%fleets
  
fit$obj$report(c(fit$sdrep$par.fixed,fit$sdrep$par.random))$predObs[idx]


x <- fit$obj$report()
x2 <- fit$obj$report(c(fit$sdrep$par.fixed,fit$sdrep$par.random))


#### Make plots for comparison of simulated vs fit ####
# Set up Catch data to plot
fitCatch <- 
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
  dplyr::left_join(fitCatch) %>%
  tidyr::gather(variable, Catch, -year) %>%
  tidyr::separate(variable, c("source", "age"))

# plot Catch data (should exactly match in total subplot)
ggplot(data = df2plotC,
       aes(x = year, y = Catch, color = source)) +
  geom_line() +
  facet_wrap(~age, scales = "free") +
  ylab("Catch (mt)")



