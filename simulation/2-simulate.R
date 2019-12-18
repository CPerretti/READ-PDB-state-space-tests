## Perform SAM simulation tests
 

# Install github version of package to use default package
#devtools::install_github("fishfollower/SAM/stockassessment")
# Install local version of package with changes
#devtools::install_local("../../SAM/stockassessment/", force = TRUE)

set.seed(321) # for reproducibility

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

scenarios <- c("uniform random", 
               "random walk",
               "fixed",
               "no misreporting")
nRep <- 300 # Number of simulation replicates
noScaledYears <- 10


# Output container
sim_label <- expand.grid(replicate = 1:nRep, scenario = scenarios, stringsAsFactors = F)
simOut          <- vector("list", length = nrow(sim_label))
setupMod_random <- vector("list", length = nrow(sim_label))
setupMod_fixed  <- vector("list", length = nrow(sim_label))
setupMod_none   <- vector("list", length = nrow(sim_label))


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

scaled_years <- confLogScale_fixed$keyScaledYears

# Generate simulation replicates
for (i in 1:nrow(sim_label)) {
  simOut[[i]] <- sim(fit = fitReal,
                     keyLogScale = keyLogScale,
                     noScaledYears = noScaledYears,
                     scenario = sim_label$scenario[i])
}

## Plot an example true vs observed vs fit ########

# ## (1) N-at-age (1000s)
plotN(simOut = simOut[[1]],
      fit = fitReal)

## (2) F-at-age
plotF(simOut = simOut[[1]],
      fit = fitReal)

## (3) Catch (mt)
plotC(simOut = simOut[[1]],
      fit = fitReal)

## (4) Survey (1000s)
plotS(simOut = simOut[[1]],
      fit = fitReal)

## Plot variables in single plot
ind_nomis <- which(sim_label$scenario == "no misreporting")
simOut2plot <- simOut[[ind_nomis[3]]] 
plotAll(simOut = simOut2plot)

## Plot some simulations from simulate.sam()
#plotSimSAM(fitReal, nsim = 10, seed = NULL)


## Fit sam to a simulation ################################

for (i in 1:nrow(sim_label)) {
  # Prep simulation data for read.ices()
  prepSimData(Sobs_N = simOut[[i]]$Sobs_N, 
              fit = fitReal, 
              Cobs_N = simOut[[i]]$Cobs_N) 
  
  # Read in data, set initial params and configuration
  setupMod_random[[i]] <- setupModel(conf = fitReal$conf, 
                                     example_dir = example_dir, 
                                     noScaledYears = noScaledYears,
                                     confLogScale = confLogScale_random)
  
  setupMod_fixed[[i]] <- setupModel(conf = fitReal$conf,
                                    example_dir = example_dir,
                                    noScaledYears = noScaledYears,
                                    confLogScale = confLogScale_fixed)
  
  setupMod_none[[i]] <- setupModel(conf = fitReal$conf,
                                   example_dir = example_dir,
                                   noScaledYears = noScaledYears,
                                   confLogScale = confLogScale_none)
  
  
  #setupMod_fixed[[i]]$par$logFpar <- fitReal$pl$logFpar # For debugging
  
}


cl <- makeCluster(detectCores() - 1) #setup nodes for parallel
#load stockassessment and functions to each node
clusterEvalQ(cl, {library(stockassessment); source("1-functions.R")}) 

# Estimate with misreporting as random effect
fitSim_random <- parLapply(cl, setupMod_random, 
                           function(x){try(sam.fit_cp(x$dat, x$conf, x$par))})

# Estimate with misreporting as fixed effect
fitSim_fixed <- parLapply(cl, setupMod_fixed, 
                          function(x){try(sam.fit(x$dat, x$conf, x$par))})
                              #map = list("logScale" = factor(c(NA, 0)))))})

# Assume no misreporting
fitSim_none <- parLapply(cl, setupMod_none, 
                          function(x){try(sam.fit(x$dat, x$conf, x$par))})

stopCluster(cl) #shut down nodes



## Error handling #####
# Exclude runs with unrealistic N estimates 
# (fit N 1000 times true N). This occurs because occasionaly high N,
# low F, low Q, results in a decent fit, but this would be rejected
# by an analyst given the large scale difference.
# Only include replicates where all three models fit successfully
ind2keep0 <- which(sapply(fitSim_random, class) != "try-error" &
                     sapply(fitSim_fixed, class)  != "try-error" &
                     sapply(fitSim_none, class)   != "try-error")
                     # also exclude non-convergences
ind2keep1 <- ind2keep0[unlist(sapply(fitSim_random[ind2keep0],
                                   function (x) x[[6]][3])) != 1 &
                       unlist(sapply(fitSim_fixed[ind2keep0],
                                   function (x) x[[6]][3])) != 1 &
                       unlist(sapply(fitSim_none[ind2keep0],
                                   function (x) x[[6]][3])) != 1]

ind2keep2 <- vector()
for (i in ind2keep1) {
  if (mean(exp(fitSim_random[[i]]$pl$logN) / 
           exp(simOut[[i]]$trueParams$pl$logN)) < 1000 &
      mean(exp(fitSim_fixed[[i]]$pl$logN) / 
           exp(simOut[[i]]$trueParams$pl$logN)) < 1000 &
      mean(exp(fitSim_none[[i]]$pl$logN) / 
           exp(simOut[[i]]$trueParams$pl$logN)) < 1000) {
    ind2keep2 <- c(ind2keep2, i)
  }
}

sim_labelAccept <- sim_label[ind2keep2,]
simOutAccept          <- simOut[ind2keep2]
setupMod_randomAccept <- setupMod_random[ind2keep2]
setupMod_fixedAccept  <- setupMod_fixed[ind2keep2]
setupMod_noneAccept   <- setupMod_none[ind2keep2]
fitSim_randomAccept   <- fitSim_random[ind2keep2]
fitSim_fixedAccept    <- fitSim_fixed[ind2keep2]
fitSim_noneAccept     <- fitSim_none[ind2keep2]

# Save fits and setups
suffix <- paste0(Sys.Date(), ".Rdata")
save(list = c("sim_labelAccept", 
              "simOutAccept", 
              "setupMod_randomAccept", 
              "setupMod_fixedAccept", 
              "setupMod_noneAccept",
              "fitSim_randomAccept", 
              "fitSim_fixedAccept", 
              "fitSim_noneAccept",
              "scaled_years",
              "confLogScale_random",
              "confLogScale_fixed",
              "confLogScale_none"), 
     file = paste0("./output/setupAndFits", suffix))

#load("output/setupAndFits2019-12-05.Rdata")


## Perform retro runs
# For accepted runs only. This is already parallelized over peels ("year") so
# it needs to be looped over simulations. Store only Mohn's rho to limit
# memory usage.
df_mohn_random <- data.frame(`R(age 1)` = numeric(length = nrow(sim_labelAccept)), 
                             `SSB` = numeric(length = nrow(sim_labelAccept)), 
                             `Fbar(4-6)` = numeric(length = nrow(sim_labelAccept)))
df_mohn_fixed <- data.frame(`R(age 1)` = numeric(length = nrow(sim_labelAccept)), 
                            `SSB` = numeric(length = nrow(sim_labelAccept)), 
                            `Fbar(4-6)` = numeric(length = nrow(sim_labelAccept)))
df_mohn_none <- data.frame(`R(age 1)` = numeric(length = nrow(sim_labelAccept)), 
                           `SSB` = numeric(length = nrow(sim_labelAccept)), 
                           `Fbar(4-6)` = numeric(length = nrow(sim_labelAccept)))
df_mohn_random$model <- "random walk"
df_mohn_fixed$model <- "fixed"
df_mohn_none$model <- "no misreporting"
df_mohn_random$scenario <- sim_labelAccept$scenario
df_mohn_fixed$scenario  <- sim_labelAccept$scenario
df_mohn_none$scenario   <- sim_labelAccept$scenario
df_mohn_random$replicate <- sim_labelAccept$replicate
df_mohn_fixed$replicate <- sim_labelAccept$replicate
df_mohn_none$replicate <- sim_labelAccept$replicate

for(i in 1:nrow(sim_labelAccept)) {
  print(paste0("retro run ", i, " of ", nrow(sim_labelAccept)))
    retro_random <- tryCatch(retro_cp(fitSim_randomAccept[[i]], year = 7),
                                  error = function(e) "error", 
                                  warning=function(w) "non-converge")
    if (is.character(retro_random)) {
      retro_fixed <- "skip" 
    } else {
      retro_fixed  <- tryCatch(retro(fitSim_fixedAccept[[i]], year = 7), 
                                  error = function(e) "error", 
                                  warning=function(w) "non-converge")
    }
    if (is.character(retro_fixed)) {
      retro_none <- "skip"
    } else {
      retro_none   <- tryCatch(retro(fitSim_noneAccept[[i]], year = 7), 
                               error = function(e) "error", 
                               warning=function(w) "non-converge")  
    }
    

  if (is.character(retro_random) ||
      is.character(retro_fixed) ||
      is.character(retro_none)){
    df_mohn_random[i,1:3] <- NA
    df_mohn_fixed[i,1:3]  <- NA
    df_mohn_none[i,1:3]   <- NA
  } else {
    df_mohn_random[i,1:3] <- stockassessment::mohn(retro_random)
    df_mohn_fixed[i,1:3] <- stockassessment::mohn(retro_fixed)
    df_mohn_none[i,1:3] <- stockassessment::mohn(retro_none)
  }
}

df_mohn <- 
  rbind(df_mohn_random, df_mohn_fixed, df_mohn_none) %>%
  tidyr::gather(variable, mohn_rho, -scenario, -replicate, -model) %>%
  dplyr::filter(!is.nan(mohn_rho), !is.infinite(mohn_rho), !is.na(mohn_rho)) %>%
  dplyr::mutate(model  = factor(model, levels = c("no misreporting", 
                                                  "fixed",
                                                  "random walk")),
                scenario = paste(scenario, "scenario"),
                scenario  = factor(scenario, levels = c("no misreporting scenario", 
                                                        "fixed scenario",
                                                        "random walk scenario",
                                                        "uniform random scenario")),
                abs_mohn = abs(mohn_rho),
                variable = ifelse(variable == "Fbar.4.6.", "F", variable),
                variable = ifelse(variable == "R.age.1.", "Recruitment", variable))

# Save mohn's rho calculations
suffix <- paste0(Sys.Date(), ".Rdata")
save(list = c("df_mohn"), 
     file = paste0("./output/df_mohn", suffix))




# Calculate fit error
err <- data.frame()
for (i in 1:nrow(sim_labelAccept)) {
  errRe_random <- calcReTsError(fitSim_randomAccept[[i]], 
                                simOutAccept[[i]],
                                confLogScale_random)
  errRe_fixed <- calcReTsError(fitSim_fixedAccept[[i]],
                               simOutAccept[[i]],
                               confLogScale_fixed)
  errRe_none <- calcReTsError(fitSim_noneAccept[[i]], 
                              simOutAccept[[i]],
                              confLogScale_none)
  
  errCSSB_random <- calcCSSBError(fitSim_randomAccept[[i]], 
                                  simOutAccept[[i]])
  errCSSB_fixed <- calcCSSBError(fitSim_fixedAccept[[i]], 
                                 simOutAccept[[i]])
  errCSSB_none <- calcCSSBError(fitSim_noneAccept[[i]], 
                                simOutAccept[[i]])
  
  err_random <- rbind(errRe_random, errCSSB_random)
  err_fixed  <- rbind(errRe_fixed, errCSSB_fixed)
  err_none   <- rbind(errRe_none, errCSSB_none)
  
  err <-
    rbind(err, 
          {rbind(data.frame(err_random, model = "random walk"),
                 data.frame(err_fixed,  model = "fixed"),
                 data.frame(err_none,   model = "no misreporting")) %>%
              dplyr::mutate(replicate = sim_labelAccept$replicate[i],
                            scenario  = as.factor(paste(sim_labelAccept$scenario[i], "scenario")),
                            scenario  = factor(scenario, levels = c("no misreporting scenario",
                                                                    "fixed scenario",
                                                                    "random walk scenario",
                                                                    "uniform random scenario")),
                            model     = factor(model, levels = c("no misreporting", 
                                                                 "fixed",
                                                                 "random walk")))})
}

suffix <- paste0(Sys.Date(), ".Rdata")
save(list = "err", 
     file = paste0("./output/err", suffix))
#load("./output/err")

## Plot example true vs observed vs fit to observed ########
# (1) N-at-age (1000s)
plotN(simOut = simOutAccept[[1]],
      fit = fitSim_noneAccept[[1]])

# (2) F-at-age
plotF(simOut = simOutAccept[[1]],
      fit = fitSim_noneAccept[[1]])

# (3) Catch (mt)
plotC(simOut = simOutAccept[[1]],
      fit = fitSim_randomAccept[[1]])

# (4) Survey (1000s)
plotS(simOut = simOutAccept[[1]],
      fit = fitSim_randomAccept[[1]])


## Plot error statistics #######################

# Plot Mohn's rho results
colors2use <- RColorBrewer::brewer.pal(3, "Dark2")
ggplot(df_mohn %>%
         dplyr::rename(`Estimation model` = model) %>%
         dplyr::group_by(`Estimation model`, scenario, variable) %>%
         dplyr::summarise(abs_mohn_mean = mean(abs_mohn),
                          abs_mohn_se   = sd(abs_mohn)/sqrt(length(abs_mohn)))) +
  aes(x = `Estimation model`, y = abs_mohn_mean, color = `Estimation model`) +
  geom_point() +
  geom_errorbar(aes(ymin = abs_mohn_mean - 1.96 * abs_mohn_se, 
                    ymax = abs_mohn_mean + 1.96 * abs_mohn_se), 
                width = 0.3) +
  facet_grid(scenario ~ variable, scales = "free") +
  theme_bw() +
  scale_color_manual(values = colors2use) +
  scale_fill_manual(values = colors2use) +
  xlab("Estimation model") +
  ylab("Mean absolute Mohn's rho")
  
ggsave("./figures/mean_mohn.png", width = 8.5, height = 7)

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
for (i in 1:nrow(sim_labelAccept)) {
  errCatchAdvice_random <- 
    calcCatchAdviceError(fitSim_randomAccept[[i]], 
                         simOutAccept[[i]],
                         confLogScale_random)
  errCatchAdvice_fixed <- 
    calcCatchAdviceError(fitSim_fixedAccept[[i]], 
                         simOutAccept[[i]],
                         confLogScale_fixed)
  errCatchAdvice_none <- 
    calcCatchAdviceError(fitSim_noneAccept[[i]], 
                         simOutAccept[[i]],
                         confLogScale_none)
  df_errCatchAdvice <-
    rbind(df_errCatchAdvice,
          rbind(errCatchAdvice_random, 
                errCatchAdvice_fixed, 
                errCatchAdvice_none) %>%
              dplyr::mutate(scenario = sim_labelAccept$scenario[[i]],
                            replicate = sim_labelAccept$replicate[[i]]))
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
# mape dot plot
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
  xlab("Estimation model")

ggsave("./figures/catch_advice_err_mean.png", width = 7, height = 6)

# Catch advice error table
catch_error_table <-
  df_errCatchAdvice %>% 
  dplyr::group_by(scenario, model) %>%
  dplyr::summarise(mean_error_pc = median(error_pc, na.rm = T) %>% round(1),
                   se_error = sd(error, na.rm = T) / sqrt(sum(!is.na(error)))) %>%
  dplyr::select(-se_error) %>%
  spread(model, mean_error_pc)

# histogram of catch error
ggplot(df_errCatchAdvice, #%>% 
         #dplyr::mutate(error_plusgroup = ifelse(error >= 10000, 10000, error),
         #               error_plusgroup = ifelse(error <= -10000, -10000, error_plusgroup)),
       aes(x = error_pc_f40, color = model, fill = model)) +
  geom_histogram(alpha=.5, position="identity") +
  geom_rug(data = df_errCatchAdvice %>% 
             dplyr::group_by(model, scenario) %>%
             dplyr::summarise(median_error_pc = median(error_pc_f40, na.rm = T)),
           aes(x = median_error_pc, color = model), size = 3) +
  facet_wrap(~scenario, scales = "free_y") +
  theme_bw() +
  guides(color=guide_legend(title="Estimation model")) +
  guides(fill=guide_legend(title="Estimation model")) +
  scale_color_manual(values = colors2use) +
  scale_fill_manual(values = colors2use) +
  # scale_x_continuous(breaks = c(-10000, -5000, 0, 5000, 10000),
  #                    labels = c("-10000+", -5000, 0, 5000, "10000+")) +    
  ylab("Frequency") +
  xlab("Catch advice percent error (fit - true)") +
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


# Plot time series error
plotTsError(err, scaled_years = scaled_years)

# Plot parameters true vs fit
plotPars(fitSim_randomAccept, fitSim_fixedAccept, fitSim_noneAccept,
         simOutAccept, sim_labelAccept,
         models2plot = c("random walk", "fixed", "no misreporting"))

# Plot simulations
dat2plot0 <- data.frame(scenario = sim_labelAccept$scenario[1],
                       replicate = sim_labelAccept$replicate[1],
                       year = colnames(simOutAccept[[1]]$N),
                       Cobs_total = colSums(simOutAccept[[1]]$Cobs_mt),
                       N_total = colSums(simOutAccept[[1]]$N),
                       IBTS_Q1_total = colSums(simOutAccept[[1]]$logStru_N[1:5,,2]),
                       IBTS_Q3_total = colSums(simOutAccept[[1]]$logStru_N[1:4,,3]))

for (i in 1:nrow(sim_labelAccept)) {
  dat2plot0 <-
    rbind(dat2plot0,
      data.frame(scenario = sim_labelAccept$scenario[i],
                 replicate = sim_labelAccept$replicate[i],
                 year = colnames(simOutAccept[[i]]$logN),
                 Cobs_total = colSums(simOutAccept[[i]]$Cobs_mt),
                 N_total = colSums(simOutAccept[[i]]$N),
                 IBTS_Q1_total = colSums(simOutAccept[[i]]$logStru_N[,,2], na.rm = T),
                 IBTS_Q3_total = colSums(simOutAccept[[i]]$logStru_N[,,3], na.rm = T)))
}

dat2plot1 <-
  dat2plot0 %>%
  gather(variable, value, -scenario, -replicate, -year) %>%
  group_by(scenario, variable, year) %>%
  summarise(ci_lo = quantile(value, 0.25, na.rm = T),
            ci_mi = quantile(value, 0.50, na.rm = T),
            ci_hi = quantile(value, 0.75, na.rm = T)) %>%
  mutate(year = as.numeric(year),
         ci_mi = ifelse(ci_mi == 0, NA, ci_mi))

ggplot(dat2plot1, aes(x = year)) +
  geom_line(aes(y=ci_mi)) +
  #geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.2) +
  facet_grid(variable~scenario, scales = "free_y") +
  theme_bw()

ggsave("./figures/median_sims.png", height = 7, width = 5)
           
           
# Plot mean error in catch and ssb vs year
err2plot_CSSB <-
  err %>% 
  #dplyr::filter(model %in% c("fixed", "random walk")) %>%
  dplyr::filter(variable %in% c("catch",#"catch_observed",
                                "ssb", "F", "N")) %>%
  dplyr::select(model, scenario, year, variable, age, replicate, error, error_pc, abs_error_pc) %>%
  dplyr::group_by(model, scenario, year, variable) %>%
  dplyr::summarise(error_pc_mean = mean(error_pc, na.rm = T),
                   mape  = mean(abs_error_pc, na.rm = T),
                   nObs = length(abs_error_pc),
                   error_pc_hi   = error_pc_mean + 1.96 * sd(error_pc, na.rm = T)/sqrt(nObs),
                   error_pc_lo  = error_pc_mean - 1.96 * sd(error_pc, na.rm = T)/sqrt(nObs),
                   mape_hi   = mape + 1.96 * sd(abs_error_pc, na.rm = T)/sqrt(nObs),
                   mape_lo   = mape - 1.96 * sd(abs_error_pc, na.rm = T)/sqrt(nObs)) %>%
  dplyr::rename(`Estimation model` = model) %>%
  dplyr::mutate(variable = ifelse(variable == "catch", "Catch", variable),
                variable = ifelse(variable == "ssb", "SSB", variable))
p <-
  ggplot(err2plot_CSSB %>% 
          dplyr::filter(`Estimation model` %in% c("fixed", "random walk"))
         ) +
         aes(x = year, color = `Estimation model`, fill = `Estimation model`) +
  geom_line(aes(y = mape)) +
  geom_ribbon(aes(ymin = mape_lo, ymax = mape_hi), color = NA, alpha = 0.3) +
  geom_hline(yintercept = 0) +
  facet_grid(scenario~variable, scales = "free_y") +
  theme_bw() +
  xlab("Year") +
  ylab("Estimation error (MAPE)") +
  scale_color_manual(values = colors2use[2:3]) +
  scale_fill_manual(values = colors2use[2:3]) +
  theme(axis.title   = element_text(size = 14),
        plot.title   = element_text(size = 16),
        strip.text   = element_text(size = 9)) +
  scale_x_continuous(breaks=c(1970, 1990, 2010))

print(p)

#ggsave(plot = p, "./figures/est_error_all.png", height = 8, width = 7)
ggsave(plot = p, "./figures/est_error_fixedrw.png", height = 8, width = 7)

p <-
  ggplot(err2plot_CSSB) +
  aes(x = year, color = `Estimation model`, fill = `Estimation model`) +
  geom_line(aes(y = error_pc_mean)) +
  geom_ribbon(aes(ymin = error_pc_lo, ymax = error_pc_hi),
              color = NA, alpha = 0.3) +
  geom_hline(yintercept = 0) +
  facet_grid(scenario~variable, scales = "free_y") +
  theme_bw() +
  xlab("Year") +
  ylab("Estimation error (MPE)") +
  scale_color_manual(values = colors2use[1:3]) +
  scale_fill_manual(values = colors2use[1:3]) +
  theme(axis.title   = element_text(size = 14),
        plot.title   = element_text(size = 16),
        strip.text   = element_text(size = 9)) +
  scale_x_continuous(breaks=c(1970, 1990, 2010))

print(p)

ggsave(plot = p, "./figures/mpe_all.png", height = 8, width = 7)


# Plot average error for each scenario
err2plot_Scale_overall <-
  err %>%  
  dplyr::filter(variable %in% c("Scale")) %>%
  dplyr::select(model, scenario, year, variable, age, replicate, 
                error_pc, abs_error_pc, abs_error) %>%
  dplyr::group_by(model, scenario, year, variable) %>%
  dplyr::summarise(error_pc_mean = mean(error_pc, na.rm = T),
                   mape  = mean(abs_error_pc, na.rm = T),
                   mae   = mean(abs_error),
                   nObs  = length(abs_error_pc),
                   error_pc_hi = error_pc_mean + 1.96 * sd(error_pc, na.rm = T)/sqrt(nObs),
                   error_pc_lo = error_pc_mean - 1.96 * sd(error_pc, na.rm = T)/sqrt(nObs),
                   mape_hi   = mape + 1.96 * sd(abs_error_pc, na.rm = T)/sqrt(nObs),
                   mape_lo   = mape - 1.96 * sd(abs_error_pc, na.rm = T)/sqrt(nObs),
                   mae_hi   = mae + 1.96 * sd(mae, na.rm = T)/sqrt(nObs),
                   mae_lo   = mae - 1.96 * sd(mae, na.rm = T)/sqrt(nObs)) %>%
  dplyr::rename(`Estimation model` = model)

p <-
  ggplot(err2plot_Scale_overall) +
  aes(x = year, color = `Estimation model`, 
      fill = `Estimation model`) +
  geom_line(aes(y = mape)) +
  geom_ribbon(aes(ymin = mape_lo, ymax = mape_hi), color = NA, alpha = 0.3) +
  geom_hline(yintercept = 0) +
  facet_grid(~scenario) +
  theme_bw() +
  xlab("Year") +
  ylab("Estimaton error of scale parameter (MAPE)") +
  scale_color_manual(values = colors2use[2:3]) +
  scale_fill_manual(values = colors2use[2:3]) +
  scale_x_continuous(breaks=c(2006, 2009, 2012))

print(p)

ggsave(plot = p, "./figures/scale_error.png", height = 5, width = 8)

# Calculate 95% interval of scale parameter in random walk scenario
err %>%
  filter(scenario == "random walk scenario",
         variable == "Scale") %>%
  summarize(mean_scale = mean(tru))
  

# Calculate average estimation error of fixed and random walk model
# over all variables, ages, and years with misreporting
err %>% 
  dplyr::filter(model %in% c("fixed", "random walk")) %>%
  dplyr::filter(scenario %in% c("no misreporting scenario",
                                "uniform random scenario")) %>%
  dplyr::filter(variable %in% c("catch", "ssb", "F", "N"),
                year %in% confLogScale_random$keyScaledYears) %>%
  dplyr::select(model, scenario, year, variable, age, 
                replicate, error, error_pc, abs_error_pc) %>%
  dplyr::group_by(model, 
                  variable) %>%
  dplyr::summarise(mape  = mean(abs_error_pc, na.rm = T))


# Plot example fit vs true
err2plot_examples <-
  err %>%
  mutate(fit_975 = exp(log(fit) + 1.96 * sdLog),
         fit_025 = exp(log(fit) - 1.96 * sdLog),
         replicate = paste("replicate", replicate),
         age = paste("age-", age)) %>%
  filter(year %in% scaled_years)
p <-
  ggplot(err2plot_examples %>%
           filter(variable == c("F")) %>%
           filter(replicate %in% unique(replicate)[1]), 
         aes(x = year)) +
  geom_line(aes(y = fit, color = model)) +
  geom_line(aes(y = tru, color = "true")) +
  geom_ribbon(aes(ymin = fit_975, ymax = fit_025, fill = model), 
              alpha = 0.3) +
  facet_grid(scenario~age, scales = "free_y") +
  theme_bw() +
  xlab("Year") +
  ylab("F estimate") +
  #scale_y_continuous(breaks=seq(1,11,2)) +
  scale_x_continuous(breaks = seq(2005, 2014, 4)) +
  scale_color_manual(values = c(colors2use[1:3], "black")) +
  scale_fill_manual(values = c(colors2use[c(2,1,3)]), guide = "none") +
  theme(axis.title   = element_text(size = 14),
        plot.title   = element_text(size = 16),
        strip.text   = element_text(size = 12),
        legend.title = element_blank())
print(p)

# Calculate the average width of the confidence interval for
# each variable in each estimation model
ci_width <-
  err %>%
  mutate(fit_975 = exp(log(fit) + 1.96 * sdLog),
         fit_025 = exp(log(fit) - 1.96 * sdLog),
         cv = fit/exp(sdLog),
         replicate = paste("replicate", replicate),
         age = paste("age-", age)) %>%
  filter(year %in% scaled_years,
         variable != "catch_observed") %>%
  group_by(scenario, model, variable) %>%
  summarise(#median_ci_width = (median(fit_975 - fit_025, na.rm = T) %>%
            #  round(2)),
            median_cv = median(cv, na.rm = T) %>% round(2)) %>%
  spread(variable, median_cv) %>%
  arrange(scenario)

# Calculate whether a retro-adjustment would occur
# given the 90% CI rule

# (1) attach to df_mohn the terminal year value and 90%CI for SSB and F and R
# (2) calculate rho-adjusted SSB, F and R
# (3) Check if it falls outside the 90% CI

# Pull out terminal year values and calculate 90% CI's
terminal_fixed <- 
  as.data.frame(matrix(nrow = 3 * nrow(sim_labelAccept), ncol = 7)) %>%
  rename(variable = V1, log_value = V2, log_se = V3, model = V4, 
         scenario = V5, replicate = V6,  year = V7)
terminal_none <- terminal_fixed
terminal_random <- terminal_fixed

c <- 1
for (i in 1:nrow(sim_labelAccept)) {
  terminal_fixed[c:(c+2),] <-
    data.frame(variable = names(fitSim_fixedAccept[[i]]$sdrep$value),
               log_value = fitSim_fixedAccept[[i]]$sdrep$value,
               log_se = fitSim_fixedAccept[[i]]$sdrep$sd,
               model = "fixed",
               scenario = sim_labelAccept$scenario[[i]],
               replicate = sim_labelAccept$replicate[[i]], stringsAsFactors = F) %>%
    filter(variable %in% c("logssb", "logfbar", "logR")) %>%
    mutate(year = rep(fitSim_fixedAccept[[i]]$data$years, 3)) %>%
    filter(year == max(year))
  terminal_none[c:(c+2),] <-
    data.frame(variable = names(fitSim_noneAccept[[i]]$sdrep$value),
               log_value = fitSim_noneAccept[[i]]$sdrep$value,
               log_se = fitSim_noneAccept[[i]]$sdrep$sd,
               model = "no misreporting",
               scenario = sim_labelAccept$scenario[[i]],
               replicate = sim_labelAccept$replicate[[i]], stringsAsFactors = F) %>%
    filter(variable %in% c("logssb", "logfbar", "logR")) %>%
    mutate(year = rep(fitSim_noneAccept[[i]]$data$years, 3)) %>%
    filter(year == max(year))
  terminal_random[c:(c+2),] <-
    data.frame(variable = names(fitSim_fixedAccept[[i]]$sdrep$value),
               log_value = fitSim_fixedAccept[[i]]$sdrep$value,
               log_se = fitSim_fixedAccept[[i]]$sdrep$sd,
               model = "random walk",
               scenario = sim_labelAccept$scenario[[i]],
               replicate = sim_labelAccept$replicate[[i]], stringsAsFactors = F) %>%
    filter(variable %in% c("logssb", "logfbar", "logR")) %>%
    mutate(year = rep(fitSim_randomAccept[[i]]$data$years, 3)) %>%
    filter(year == max(year))
  c <- c+3
}
  terminal_all <-
    rbind(terminal_none, terminal_fixed, terminal_random) %>%
    mutate(value = exp(log_value),
           value_90hi = exp(log_value + 1.645 * log_se),
           value_90lo = exp(log_value - 1.645 * log_se),
           variable = ifelse(variable == "logssb", "SSB", variable),
           variable = ifelse(variable == "logfbar", "F", variable),
           variable = ifelse(variable == "logR", "Recruitment", variable),
           scenario = paste(scenario, "scenario")) %>%
    select(model, scenario, replicate, variable, value, value_90hi, value_90lo)



  df_retro_adjust <-
    df_mohn %>%
    left_join(terminal_all) %>%
    mutate(value_rho_adjusted = value /(1 + mohn_rho),
           adjustment_needed = ifelse((value_rho_adjusted > value_90hi | 
                                         value_rho_adjusted < value_90lo),
                                      "yes",
                                      "no"))
  df_retro_summary <-
    df_retro_adjust %>%
    group_by(scenario, model) %>%
    summarise(percent_adjusted = 
                round(100*sum(adjustment_needed == "yes")/n(), 1)) %>%
    spread(model, percent_adjusted) %>%
    select(scenario, `no misreporting`, fixed, `random walk`)
  
  
  
  
  
  

