# Replicate SAM model fit using the equations written here

# To do:
# (1) Clean plot functions so it just takes simOut[[i]]
# (2) Duplicate F's in output to match the config key

# Required packages
library(plyr) # always load before dplyr
library(dplyr)
library(ggplot2)
library(tidyr)
library(stockassessment)

# Load user-defined functions
source("1-functions.R")

## Run simulation #########################################
# Use atl herring fit to set up simulation
load("../atlherring_example/output/fitHer.Rdata")

# How many simulation replicates to do
nRep <- 100

# Generate simulation replicates
simOut <- list()
for (i in 1:nRep) {
  simOut[[i]] <- sim(fit = fitHer)  
}



## Plot an example true vs observed vs *herring* fit ########

## (1) N-at-age (1000s)
plotN(N = simOut[[1]]$N,
      fit = fitHer)

## (2) F-at-age
plotF(simOut = simOut[[1]],
      fit = fitHer)

## (3) Catch (mt)
plotC(Cobs_mt = simOut[[1]]$Cobs_mt,
      Ctru_mt = simOut[[1]]$Ctru_mt,
      fit = fitHer)

## (4) Survey (1000s)
plotS(Sobs_N = simOut[[1]]$Sobs_N,
      Stru_N = simOut[[1]]$Stru_N,
      fit = fitHer)

## Plot some simulations from simulate.sam()
#plotSimSAM(fitHer, nsim = 10, seed = NULL)


## Fit sam to a simulation ################################
fitSim <- list()
for (i in 1:nRep) {
  # Prep simulation data for read.ices()
  prepSimData(Sobs_N = simOut[[i]]$Sobs_N, 
              fit = fitHer, 
              Cobs_N = simOut[[i]]$Cobs_N) 
  
  # Read in data, set initial params and configuration
  setupOut <- setupModel()
  
  # Fit sam
  fitSim[[i]] <- sam.fit(setupOut$dat, setupOut$conf, 
                    setupOut$par, sim.condRE = FALSE)
}



## Plot true vs observed vs fit to observed ################
# (1) N-at-age (1000s)
# plotN(N = simOut$N, 
#       fit = fitSim)

# (2) F-at-age
# plotF(simOut = simOut[[1]],
#       fit = fitHer)

# (2) Catch (mt)
# plotC(Cobs_mt = simOut$Cobs_mt, 
#       Ctru_mt = simOut$Ctru_mt, 
#       fit = fitSim)

# (3) Survey (1000s)
# plotS(Sobs_N = simOut$Sobs_N, 
#       Stru_N = simOut$Stru_N, 
#       fit = fitSim)


## Plot fit vs true parameter values #######################

# Plot fixed effects
parsFixed <- which(names(fitSim[[1]]$pl) %in% names(fitSim[[1]]$obj$par))
df_parsOut <- data.frame()
for (h in parsFixed) {
  for (i in 1:nRep) {
    df_parsOut <-
      rbind(df_parsOut,
            data.frame(variable = paste(names(fitSim[[1]]$pl)[h], 
                                        1:length(fitSim[[1]]$pl[[h]]), sep = "."),
                       tru = simOut[[i]]$trueParams$pl[[h]],
                       est = fitSim[[i]]$pl[[h]],
                       sd  = fitSim[[i]]$plsd[[h]],
                       replicate = i))
  }  
}


df2plot <-
  df_parsOut %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(tru = unique(tru),
                   est_mean = mean(est),
                   est_se = sd(est)/sqrt(nRep))

ggplot(df2plot, aes(y = variable)) +
  geom_point( aes(x = tru), color = "red", size  = 3) +
  geom_point( aes(x = est_mean), color = "blue", size = 3) +
  geom_errorbarh(aes(xmin = est_mean - 1.96*est_se,
                     xmax = est_mean + 1.96*est_se), color = "blue") +
  theme_bw() +
  ylab("") +
  xlab("Value") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

# Plot error for N and F (CHECK TO MAKE SURE F IS SPECIFIED TO AGE CORRECTLY<< (it's not yet))
indNF <- which(names(fitSim[[1]]$pl) %in% c("logN", "logF"))
err_logNF <- data.frame()
for (h in indNF) {
  for (i in 1:nRep) {
    rownames(fitSim[[i]]$pl[[h]]) <- paste0("fit.", 1:nrow(fitSim[[i]]$pl[[h]]))
    Sd <- exp(fitSim[[i]]$plsd[[h]])
    rownames(Sd) <- paste0("Sd.", 1:nrow(fitSim[[i]]$pl[[h]]))
    err_logNF <-
      rbind(err_logNF,
            fitSim[[i]]$pl[[h]] %>%
              t() %>%
              exp %>%
              as.data.frame() %>%
              cbind(Sd %>%
                    t() %>%
                    as.data.frame()) %>%
              cbind(data.frame(simOut[[i]]$trueParams$pl[[h]] %>% t %>% exp)) %>%
              dplyr::mutate(year = as.numeric(fitSim[[i]]$data$years)) %>%
              tidyr::gather(variable, N, -year) %>%
              dplyr::mutate(variable = gsub(x = variable, pattern = "X", replacement = "tru.")) %>%
              tidyr::separate(variable, c("source", "age")) %>%
              tidyr::spread(source, N) %>%
              dplyr::mutate(age = as.numeric(age),
                            error = (fit - tru),
                            decile = ceiling(10 * pnorm(q    = log(tru), 
                                                        mean = log(fit), 
                                                        sd   = log(Sd))),
                            replicate = i,
                            variable = ifelse(names(fitSim[[i]]$pl[h]) == "logN", "N",
                                              if(names(fitSim[[i]]$pl[h]) == "logF") "F")))
  }
}

# Calculate median error
err_logNFannual <-
  err_logNF %>%
  dplyr::group_by(variable, age, year) %>%
  dplyr::summarise(median_error = median(error))

err_logNFoverall <-
  err_logNF %>%
  dplyr::group_by(variable, age) %>%
  dplyr::summarise(median_error = median(error))

# Plot the median error over all years
ggplot(err_logNFoverall %>% 
         dplyr::filter(variable == "F"),
       aes(x = age, y = median_error)) +
  geom_line() +
  theme_bw() +
  ylab("Median error (1000's)") +
  xlab("Age") +
  ggtitle("N estimation error") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)) #+
  #facet_wrap(~variable, scales = "free_y")

# Plot median error in each year for each age
ggplot(err_logNFannual %>% 
         dplyr::filter(variable == "F") %>%
         dplyr::mutate(as.numeric(year)),
       aes(x = year, y = median_error)) +
  geom_line() +
  theme_bw() +
  ylab("Median error (1000s)") +
  xlab("Year") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)) +
  facet_wrap(~age) +
  ggtitle("Estimation error vs year")

# Plot a few example age-1 fit vs tru time series
ggplot(err_logNF %>%
         dplyr::select(-Sd, -decile) %>%
         dplyr::filter(variable == "F",
                       age == 1,
                       replicate %in% 1:10) %>%
         tidyr::gather(source, value, 
                       -year, -age, -replicate, -variable),
       aes(x = year, y = value, color = source)) +
  geom_line() +
  facet_wrap(~replicate)

# Plot decile coverage for each age and variable
ggplot(err_logNF,
         aes(x = decile)) +
  geom_histogram(binwidth=1, colour="white") +
  geom_hline(yintercept = ncol(fitSim[[1]]$pl$logN) * nRep / 10, color = "dark grey") +
  theme_bw() +
  facet_grid(variable~age)
         
         

# Fit sam to a simulate.sam replicate
# Need to resimulate with full.data = TRUE to get output that sam.fit() can use
# set.seed(123)
#simOut2fit <- stockassessment:::simulate.sam(fitHer, nsim = 10, full.data = TRUE)
#sam.fit(data = simOut2fit[[10]], conf = fitHer$conf, par = defpar(fitHer$data, fitHer$conf))

