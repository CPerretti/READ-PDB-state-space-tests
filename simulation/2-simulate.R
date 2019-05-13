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

# ## (1) N-at-age (1000s)
plotN(simOut = container$simOut[[2]],
      fit = fitReal)

## (2) F-at-age
plotF(simOut = container$simOut[[2]],
      fit = fitReal)

## (3) Catch (mt)
plotC(simOut = container$simOut[[2]],
      fit = fitReal)

## (4) Survey (1000s)
plotS(simOut = container$simOut[[2]],
      fit = fitReal)

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
       exp(containerAccept$simOut[[i]]$trueParams$pl$logN)) < 1000 &
      mean(exp(containerAccept$fitSim_fixed[[i]]$pl$logN) / 
           exp(containerAccept$simOut[[i]]$trueParams$pl$logN)) < 1000 &
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
# (1) N-at-age (1000s)
plotN(simOut = containerAccept$simOut[[2]],
      fit = containerAccept$fitSim_none[[2]])

# (2) F-at-age
plotF(simOut = containerAccept$simOut[[1]],
      fit = containerAccept$fitSim_none[[1]])

# (3) Catch (mt)
plotC(simOut = containerAccept$simOut[[2]],
      fit = containerAccept$fitSim_random[[2]])

# (4) Survey (1000s)
plotS(simOut = containerAccept$simOut[[1]],
      fit = containerAccept$fitSim_random[[1]])


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
                scenario = paste(scenario, "scenario"),
                scenario  = factor(scenario, levels = c("no misreporting scenario", 
                                                        "fixed scenario",
                                                        "random walk scenario",
                                                        "uniform random scenario")),
                abs_mohn = abs(mohn_rho)
                ) #%>%
  #dplyr::group_by(model, scenario, variable) %>%
  #dplyr::summarise(abs_mohn_mean = mean(abs_mohn, na.rm = T),
  #                 abs_mohn_se   = sd(abs_mohn, na.rm = T)/sqrt(sum(!is.na(abs_mohn))))


# Plot Mohn's rho results
colors2use <- RColorBrewer::brewer.pal(3, "Dark2")
ggplot(df_mohn %>% 
         dplyr::mutate(abs_mohn_plusgroup = ifelse(abs_mohn >= 1, 1, abs_mohn)),
       aes(x = abs_mohn_plusgroup, color = model, fill = model)) +
  geom_histogram(alpha=.5, position="identity") +
  geom_rug(data = df_mohn %>% 
             dplyr::group_by(model, scenario) %>%
             dplyr::summarise(median_abs_mohn = median(abs_mohn, na.rm = T)),
           aes(x = median_abs_mohn, color = model),
           size = 3) +
  facet_grid(scenario~variable, scales = "free_y") +
  theme_bw() +
  guides(color=guide_legend(title="Estimation model")) +
  guides(fill=guide_legend(title="Estimation model")) +
  scale_color_manual(values = colors2use) +
  scale_fill_manual(values = colors2use) +
  scale_x_continuous(breaks = seq(0,1,0.25),
                     labels = c(0, 0.25, 0.5, 0.75, "1+")) +
  ylab("Frequency") +
  xlab("Mohn's rho") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 12))
  
# Calculate F40% for fitted models and true models
# F40% is the F that results in the SSB/R being 40% of the
# unfished SSB/R.

df_errCatchAdvice <- data.frame()
for (i in 1:nrow(containerAccept)) {
  errCatchAdvice_random <- 
    calcCatchAdviceError(containerAccept$fitSim_random[[i]], 
                         containerAccept$simOut[[i]],
                         confLogScale_random)
  errCatchAdvice_fixed <- 
    calcCatchAdviceError(containerAccept$fitSim_fixed[[i]], 
                         containerAccept$simOut[[i]],
                         confLogScale_fixed)
  errCatchAdvice_none <- 
    calcCatchAdviceError(containerAccept$fitSim_none[[i]], 
                         containerAccept$simOut[[i]],
                         confLogScale_none)
  df_errCatchAdvice <-
    rbind(df_errCatchAdvice,
          rbind(errCatchAdvice_random, 
                errCatchAdvice_fixed, 
                errCatchAdvice_none) %>%
              dplyr::mutate(scenario = containerAccept$scenario[[i]],
                            replicate = containerAccept$replicate[[i]]))
}
df_errCatchAdvice$model <- factor(df_errCatchAdvice$model, 
                                  levels = c("no misreporting", 
                                             "fixed",
                                             "random walk"))
df_errCatchAdvice$scenario <- factor(paste(df_errCatchAdvice$scenario, "scenario"), 
                                     levels = c("no misreporting scenario", 
                                                "fixed scenario",
                                                "random walk scenario",
                                                "uniform random scenario"))

# Plot error in catch advice
colors2use <- RColorBrewer::brewer.pal(3, "Dark2")
ggplot(df_errCatchAdvice %>%
        dplyr::group_by(model, variable, scenario) %>%
        dplyr::summarise(mape    = mean(abs_error_pc),
                         mape_se = sd(abs_error_pc)/sqrt(length(abs_error_pc))),
       aes(x = model, y = mape, color = model)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mape - 1.96 * mape_se,
                    ymax = mape + 1.96 * mape_se),
                width = 0.2) +
  geom_hline(aes(yintercept = 0), size = 0.2) +
  facet_wrap(~scenario) +
  theme_bw() +
  guides(color=guide_legend(title="Estimation model")) +
  scale_color_manual(values = colors2use) +
  ylab("Catch advice error (MAPE)") +
  xlab("Estimation model") +
  #ggtitle("Catch advice error (recommended catch @ F40%)") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 12))

# Catch advice error table
df_errCatchAdvice %>% 
  dplyr::group_by(model, scenario) %>%
  dplyr::summarise(median_error = median(error, na.rm = T),
                   se_error = sd(error, na.rm = T) / sqrt(sum(!is.na(error))))

# histogram of catch error
ggplot(df_errCatchAdvice %>% 
         dplyr::mutate(error_plusgroup = ifelse(error >= 10000, 10000, error),
                       error_plusgroup = ifelse(error <= -10000, -10000, error_plusgroup)),
       aes(x = error_plusgroup, color = model, fill = model)) +
  geom_histogram(alpha=.5, position="identity") +
  geom_rug(data = df_errCatchAdvice %>% 
             dplyr::group_by(model, scenario) %>%
             dplyr::summarise(median_error = median(error, na.rm = T)),
           aes(x = median_error, color = model), size = 3) +
  facet_wrap(~scenario, scales = "free_y") +
  theme_bw() +
  guides(color=guide_legend(title="Estimation model")) +
  guides(fill=guide_legend(title="Estimation model")) +
  scale_color_manual(values = colors2use) +
  scale_fill_manual(values = colors2use) +
  scale_x_continuous(breaks = c(-10000, -5000, 0, 5000, 10000),
                     labels = c("-10000+", -5000, 0, 5000, "10000+")) +    
  ylab("Frequency") +
  xlab("Catch advice error (fit - true) (metric tons)") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 12))

# histogram of absolute catch error
ggplot(df_errCatchAdvice %>% 
         dplyr::mutate(abs_error_plusgroup = ifelse(abs_error >= 100000, 100000, abs_error)),
       aes(x = abs_error_plusgroup, color = model, fill = model)) +
  geom_histogram(alpha=.5, position="identity") +
  geom_rug(data = df_errCatchAdvice %>% 
             dplyr::group_by(model, scenario) %>%
             dplyr::summarise(median_abs_error = median(abs_error, na.rm = T)),
           aes(x = median_abs_error, color = model), size = 3) +
  facet_wrap(~scenario, scales = "free_y") +
  theme_bw() +
  guides(color=guide_legend(title="Estimation model")) +
  guides(fill=guide_legend(title="Estimation model")) +
  scale_color_manual(values = colors2use) +
  scale_fill_manual(values = colors2use) +
  scale_x_continuous(breaks = c(0, 25000, 50000, 75000, 100000),
                     labels = c(0, 25000, 50000, 75000, "100000+")) +    
  ylab("Frequency") +
  xlab("Catch advice absolute error, abs(fit - true), (metric tons)") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 12))
  
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
plotTsError(containerAccept)

# Plot parameters true vs fit
plotPars(containerAccept,
         models2plot = c("random walk", "fixed", "no misreporting"))



