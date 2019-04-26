## Perform SAM simulation tests
 
#() Catch advice metric
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
library(RColorBrewer)

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
scenarios <- c("uniform random", 
               "random walk",
               "fixed",
               "no misreporting")
nRep <- 100 # Number of simulation replicates
noScaledYears <- 10

# Output container
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
  list(logScaleType = "random walk",
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
       keyParScaledYA =  matrix(data = rep(c(0), 
                                           each = noScaledYears * ncol(fitReal$data$propF)),
                                nrow = noScaledYears)) # one param all years

confLogScale_none <- list(logScaleType = "no misreporting")


# Generate simulation replicates
for (i in 1:nrow(container)) {
  container$simOut[[i]] <- sim(fit = fitReal,
                               keyLogScale = keyLogScale,
                               noScaledYears = noScaledYears,
                               container_i = container[i,])
}

## Plot an example true vs observed vs *real data* fit ########

# # ## (1) N-at-age (1000s)
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
  
  
  #container$setupMod_fixed[[i]]$par$logFpar <- fitReal$pl$logFpar # For debugging
  
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
                                              #map = list("logScale" = factor(c(NA, 0)))))})

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
# non_converged_random <- which(unlist(sapply(container$fitSim_random,
#                             function (x) x[[6]][3])) == 1)
# non_converged_fixed <- which(unlist(sapply(container$fitSim_fixed,
#                   function (x) x[[6]][3])) == 1)
# non_converged_none <- which(unlist(sapply(container$fitSim_none,
#                   function (x) x[[6]][3])) == 1)

# Also exclude runs with unrealistic N estimates 
# (fit N 1000 times true N). This occurs because occasionaly high N,
# low F, low Q, results in a decent fit, but this would be rejected
# by an analyst given the large scale difference.
i2include <- vector()
for (i in 1:nrow(containerAccept)) {
  if (mean(exp(containerAccept$fitSim_random[[i]]$pl$logN) / 
       exp(containerAccept$simOut[[i]]$trueParams$pl$logN)) < 1000 |
      mean(exp(containerAccept$fitSim_fixed[[i]]$pl$logN) / 
           exp(containerAccept$simOut[[i]]$trueParams$pl$logN)) < 1000 |
      mean(exp(containerAccept$fitSim_none[[i]]$pl$logN) / 
           exp(containerAccept$simOut[[i]]$trueParams$pl$logN)) < 1000) {
    i2include <- c(i2include, i)
  }
}
containerAccept <- containerAccept[i2include,]
 

## Perform retro runs
# For accepted runs only. This is already parallelized over peels ("year") so
# it needs to be looped over simulations.
containerAccept$retro_random <- vector("list", length = nrow(containerAccept))
containerAccept$retro_fixed  <- vector("list", length = nrow(containerAccept))
containerAccept$retro_none   <- vector("list", length = nrow(containerAccept))
for (i in 1:nrow(containerAccept)) {
  containerAccept$retro_random[[i]] <- tryCatch(retro_cp(containerAccept$fitSim_random[[i]], year = 5),
                                                error = function(e) "error", 
                                                warning=function(w) "non-converge")
  if (is.character(containerAccept$retro_random[[i]])) next
  containerAccept$retro_fixed[[i]]  <- tryCatch(retro(containerAccept$fitSim_fixed[[i]], year = 5), 
                                                error = function(e) "error", 
                                                warning=function(w) "non-converge")
  if (is.character(containerAccept$retro_fixed[[i]])) next
  containerAccept$retro_none[[i]]   <- tryCatch(retro(containerAccept$fitSim_none[[i]], year = 5), 
                                                error = function(e) "error", 
                                                warning=function(w) "non-converge")
  if (is.character(containerAccept$retro_none[[i]])) next
}


## Plot example true vs observed vs fit to observed ########
## (1) N-at-age (1000s)
  # plotN(simOut = containerAccept$simOut[[1]],
  #       fit = containerAccept$fitSim_none[[1]])

## (2) F-at-age
# plotF(simOut = containerAccept$simOut[[1]],
#       fit = containerAccept$fitSim_none[[1]])

## (3) Catch (mt)
# plotC(simOut = containerAccept$simOut[[1]],
#       fit = containerAccept$fitSim_random[[1]])

## (4) Survey (1000s)
# plotS(simOut = containerAccept$simOut[[1]],
#       fit = containerAccept$fitSim_random[[1]])


## Plot error statistics #######################

# Calculate Mohn's rho
# Exclude replicates where any peel didn't converge
# Dimensions are number of models * number of rows in containerAccept
df_mohn_random <- data.frame(`R(age 1)` = numeric(length = nrow(containerAccept)), 
                             `SSB` = numeric(length = nrow(containerAccept)), 
                             `Fbar(4-6)` = numeric(length = nrow(containerAccept)))
df_mohn_fixed <- data.frame(`R(age 1)` = numeric(length = nrow(containerAccept)), 
                             `SSB` = numeric(length = nrow(containerAccept)), 
                             `Fbar(4-6)` = numeric(length = nrow(containerAccept)))
df_mohn_none <- data.frame(`R(age 1)` = numeric(length = nrow(containerAccept)), 
                           `SSB` = numeric(length = nrow(containerAccept)), 
                           `Fbar(4-6)` = numeric(length = nrow(containerAccept)))
df_mohn_random$model <- "random walk"
df_mohn_fixed$model <- "fixed"
df_mohn_none$model <- "no misreporting"
df_mohn_random$scenario <- containerAccept$scenario
df_mohn_fixed$scenario  <- containerAccept$scenario
df_mohn_none$scenario   <- containerAccept$scenario
df_mohn_random$replicate <- containerAccept$replicate
df_mohn_fixed$replicate <- containerAccept$replicate
df_mohn_none$replicate <- containerAccept$replicate

for(i in 1:nrow(containerAccept)) {
  if (is.character(containerAccept$retro_random[[i]]) ||
      is.character(containerAccept$retro_fixed[[i]]) ||
      is.character(containerAccept$retro_none[[i]])){
    df_mohn_random[i,1:3] <- NA
    df_mohn_fixed[i,1:3]  <- NA
    df_mohn_none[i,1:3]   <- NA
  } else {
    df_mohn_random[i,1:3] <- stockassessment::mohn(containerAccept$retro_random[[i]])
    df_mohn_fixed[i,1:3] <- stockassessment::mohn(containerAccept$retro_fixed[[i]])
    df_mohn_none[i,1:3] <- stockassessment::mohn(containerAccept$retro_none[[i]])
  }
}

df_mohn <- 
  rbind(df_mohn_random, df_mohn_fixed, df_mohn_none) %>%
  tidyr::gather(variable, mohn_rho, -scenario, -replicate, -model) %>%
  dplyr::mutate(model  = factor(model, levels = c("no misreporting", 
                                                  "fixed",
                                                  "random walk")),
                scenario  = factor(scenario, levels = c("no misreporting", 
                                                        "fixed",
                                                        "random walk",
                                                        "uniform random")),
                abs_mohn = abs(mohn_rho)) %>%
  dplyr::group_by(model, scenario, variable) %>%
  dplyr::summarise(abs_mohn_mean = mean(abs_mohn, na.rm = T),
                   abs_mohn_se   = sd(abs_mohn, na.rm = T)/sum(!is.na(abs_mohn)))

# Plot Mohn's rho results
#plotMohn(df_mohn)
colors2use <- RColorBrewer::brewer.pal(3, "Dark2")
ggplot(df_mohn,
       aes(x = model, y = abs_mohn_mean, color = model)) +
  geom_point() +
  geom_errorbar(aes(ymin = abs_mohn_mean - 1.96 * abs_mohn_se,
                    ymax = abs_mohn_mean + 1.96 * abs_mohn_se),
                width = 0.2) +
  geom_hline(aes(yintercept = 0), size = 0.2) +
  facet_grid(scenario~variable, scales = "free_y") +
  theme_bw() +
  guides(color=guide_legend(title="Estimation model")) +
  scale_color_manual(values = colors2use) +
  ylab("Mean absolute mohn's rho") +
  xlab("Estimation model")
  
# Calculate F40% for fitted models and true models
# F40% is the F that results in the SSB/R being 40% of the
# unfished SSB/R.
#ypr_cp(containerAccept$fitSim_none[[1]])

fit <- containerAccept$fitSim_none[[1]]
simOut <- containerAccept$simOut[[1]]
idxno <- which(fit$data$years==max(fit$data$years))
aveYears <- 5
ave.sw<-colMeans(fit$data$stockMeanWeight[(idxno-aveYears+1):idxno,,drop=FALSE])
ave.cw<-colMeans(fit$data$catchMeanWeight[(idxno-aveYears+1):(idxno-1),,drop=FALSE])
ave.pm<-colMeans(fit$data$propMat[(idxno-aveYears+1):idxno,,drop=FALSE])
ave.nm<-colMeans(fit$data$natMor[(idxno-aveYears+1):idxno,,drop=FALSE])
ave.sel_est <- colMeans(faytable(fit)[(idxno-aveYears+1):idxno,,drop=FALSE]/
                         apply(faytable(fit)[(idxno-aveYears+1):idxno,,drop=FALSE],1,max))
ave.sel_tru <- colMeans(t(exp(simOut$trueParams$pl$logF))[(idxno-aveYears+1):idxno,,drop=FALSE]/
                         apply(t(exp(simOut$trueParams$pl$logF))[(idxno-aveYears+1):idxno,,drop=FALSE],1,max) )

naa_est <- vector("numeric", length(fitReal$conf$minAge:fitReal$conf$maxAge))
naa_est[1] <- 1
naa_tru <- vector("numeric", length(fitReal$conf$minAge:fitReal$conf$maxAge))
naa_tru[1] <- 1
full_f <- seq(0, 10, 0.01)
maxFage_est <- which(ave.sel == 1)
maxFage_tru <- which(ave.sel == 1)
ssbr_est <- vector("numeric", length = length(full_f))
ssbr_tru <- vector("numeric", length = length(full_f))
for (h in 1:length(full_f)) {
  for (i in 2:fit$conf$maxAge)  {
    naa_est[i] <- naa_est[i-1] * exp(-(full_f[h] * ave.sel_est[i-1] + ave.nm[i-1]))
    naa_tru[i] <- naa_tru[i-1] * exp(-(full_f[h] * ave.sel_tru[i-1] + ave.nm[i-1]))
    if (i == fit$conf$maxAge & fit$conf$maxAgePlusGroup) { # if plus group
      naa_est[i] <- naa_est[i] * 1/(1-exp(-(full_f[h] * ave.sel[maxFage] + ave.nm[maxFage_est])))
      naa_tru[i] <- naa_tru[i] * 1/(1-exp(-(full_f[h] * ave.sel[maxFage] + ave.nm[maxFage_tru])))
    }
  }
  ssbr_est[h] <- sum(naa_est * ave.sw * ave.pm)
  ssbr_tru[h] <- sum(naa_tru * ave.sw * ave.pm)
}

plot(full_f, abs(ssbr_tru/ssbr_tru[1] - .4), type = "l", col ="black")
lines(full_f, abs(ssbr_est/ssbr_est[1] - .4), type = "l", col ="blue")
f40ind_tru <- which.min(abs(ssbr_tru/ssbr_tru[1] - .4))
f40ind_est <- which.min(abs(ssbr_est/ssbr_est[1] - .4))
f40_tru <- full_f[f40ind_tru] * ave.sel_tru
f40_est <- full_f[f40ind_est] * ave.sel_est
z40_tru <- f40_tru + ave.nm
z40_est <- f40_est + ave.nm

catch40_mt_tru <- exp(simOut$trueParams$pl$logN[,idxno]) * 
                    (1 - exp(-z40_tru)) * f40_tru/z40_tru * ave.cw
catch40_mt_est <- exp(fit$pl$logN[,idxno]) * 
                    (1 - exp(-z40_est)) * f40_est/z40_est * ave.cw #<< MAKE INTO A FUNCTION
  
# Calculate fit error
containerAccept$err_random <- vector("list", length = nrow(containerAccept))
containerAccept$err_fixed  <- vector("list", length = nrow(containerAccept))
containerAccept$err_none   <- vector("list", length = nrow(containerAccept))
for (i in 1:nrow(containerAccept)) {
  errRe_random <- calcReTsError(containerAccept$fitSim_random[[i]], 
                                containerAccept$simOut[[i]],
                                confLogScale_random)
  errRe_fixed <- calcReTsError(containerAccept$fitSim_fixed[[i]],
                                containerAccept$simOut[[i]],
                                confLogScale_fixed)
  errRe_none <- calcReTsError(containerAccept$fitSim_none[[i]], 
                                containerAccept$simOut[[i]],
                                confLogScale_none)
  
  errCSSB_random <- calcCSSBError(containerAccept$fitSim_random[[i]], 
                                  containerAccept$simOut[[i]])
  errCSSB_fixed <- calcCSSBError(containerAccept$fitSim_fixed[[i]], 
                                  containerAccept$simOut[[i]])
  errCSSB_none <- calcCSSBError(containerAccept$fitSim_none[[i]], 
                                  containerAccept$simOut[[i]])
  
  containerAccept$err_random[[i]] <- rbind(errRe_random, errCSSB_random)
  containerAccept$err_fixed[[i]]  <- rbind(errRe_fixed, errCSSB_fixed)
  containerAccept$err_none[[i]]   <- rbind(errRe_none, errCSSB_none)
}


# Save output
suffix <- paste0(Sys.time(), ".Rdata")
save(list = "container", file = paste0("./output/container", suffix))
save(list = "containerAccept", file = paste0("./output/containerAccept", suffix))

# Plot time series error
plotTsError(containerAccept) #<< CONTINUE TO FIGURE OUT WHY THE LAST PLOT IS TRUE

# Plot parameters true vs fit
plotPars(containerAccept,
         models2plot = c("random walk", "fixed", "no misreporting"))



