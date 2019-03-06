## Perform SAM simulation tests
 
# Install github version of package to use default package
#devtools::install_github("fishfollower/SAM/stockassessment")
# Install local version of package to make changes
#devtools::install_local("../../SAM/stockassessment/", force = TRUE)

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
# example_dir <- "atlherring_example"
# load(paste0("../", example_dir, "/output/fitHer.Rdata"))
# fitReal <- fitHer

# Use North Sea cod fit to set up simulation
example_dir <- "nscod_example"
load(paste0("../", example_dir, "/fitNScod.Rdata"))
fitReal <- fitNScod

set.seed(321) # for reproducibility

# How many simulation replicates to do

nRep <- 10
noScaledYears <- 10
logScale <- #c(log(2))#, log(4))#
    log(runif(n = noScaledYears * ncol(fitReal$data$propF), min = 1, max = 3)) # sequence of misreporting

# Generate simulation replicates
simOut <- list()
for (i in 1:nRep) {
  simOut[[i]] <- sim(fit = fitReal, noScaledYears = noScaledYears, logScale = logScale)  
}

## Plot an example true vs observed vs *real data* fit ########

# ## (1) N-at-age (1000s)
# plotN(simOut = simOut[[1]],
#       fit = fitReal)
# 
# ## (2) F-at-age
# plotF(simOut = simOut[[1]],
#       fit = fitReal) 
# 
# ## (3) Catch (mt)
# plotC(simOut = simOut[[1]],
#       fit = fitReal)
# 
# ## (4) Survey (1000s)
# plotS(simOut = simOut[[1]],
#       fit = fitReal)

## Plot some simulations from simulate.sam()
#plotSimSAM(fitReal, nsim = 10, seed = NULL)


## Fit sam to a simulation ################################
setupOut <- list()
for (i in 1:nRep) {
  # Prep simulation data for read.ices()
  prepSimData(Sobs_N = simOut[[i]]$Sobs_N, 
              fit = fitReal, 
              Cobs_N = simOut[[i]]$Cobs_N) 
  
  # Read in data, set initial params and configuration
  setupOut[[i]] <- setupModel(conf = fitReal$conf, 
                              example_dir = example_dir, 
                              noScaledYears = noScaledYears)
  #setupOut[[i]]$par$logScale <- logScale #Just to set map on true values.
}

# Fit model to replicates in parallel <<<< MAKE SURE LOGSCALE IS BEING APPLIED CORRECTLY (DOES ALL AGES MAKE SENSE??)
fitSimTest <- sam.fit_cp(setupOut[[1]]$dat, setupOut[[1]]$conf, setupOut[[1]]$par,
                         map = list("logScale" = factor(cbind(matrix(data = NA,
                                                                     nrow = nrow(setupOut[[1]]$par$logN),
                                                                     ncol = ncol(setupOut[[1]]$par$logN) - noScaledYears),
                                                              matrix(data = 1:(nrow(setupOut[[1]]$par$logN)*noScaledYears),
                                                                     nrow = nrow(setupOut[[1]]$par$logN),
                                                                     ncol = noScaledYears)))))

fitSim_random <- list()
fitSim_noMis <- list()

for (i in 1:length(setupOut)) {
  fitSim_random[[i]] <- try(sam.fit_cp(setupOut[[i]]$dat, setupOut[[i]]$conf, setupOut[[i]]$par,
                                   map = list("logScale" = factor(cbind(matrix(data = NA,
                                                                               nrow = nrow(setupOut[[1]]$par$logN),
                                                                               ncol = ncol(setupOut[[1]]$par$logN) - noScaledYears),
                                                                        matrix(data = 1:(nrow(setupOut[[1]]$par$logN)*noScaledYears),
                                                                               nrow = nrow(setupOut[[1]]$par$logN),
                                                                               ncol = noScaledYears))))))
  fitSim_noMis[[i]] <- try(sam.fit_cp(setupOut[[i]]$dat, setupOut[[i]]$conf, setupOut[[i]]$par,
                                    map = list("logScale" = factor(matrix(data = NA,
                                                                          nrow = nrow(setupOut[[1]]$par$logN),
                                                                          ncol = ncol(setupOut[[1]]$par$logN))))))
}
  
  
# cl <- makeCluster(detectCores() - 1) #setup nodes for parallel
# clusterEvalQ(cl, {library(stockassessment); source("1-functions.R")}) #load stockassessment to each node
# # Estimate with misreporting as random effect
# fitSim_random <- parLapply(cl, setupOut,
#                     function(x){try(sam.fit_cp(x$dat, x$conf, x$par,
#                                                map = list("logScale" = factor(cbind(matrix(data = NA, 
#                                                                                            nrow = nrow(x$par$logF), 
#                                                                                            ncol = ncol(x$par$logF) - x$conf$noScaledYears),
#                                                                                     matrix(data = 1:(nrow(x$par$logF)*x$conf$noScaledYears),
#                                                                                            nrow = nrow(x$par$logF), 
#                                                                                            ncol = x$conf$noScaledYears))))))})
# 
# # Estimate with misreporting as fixed effect
# # start.time <- Sys.time()
# # fitSim_fixed <- parLapply(cl, setupOut,
# #                     function(x){try(sam.fit_cp(x$dat, x$conf, x$par))})
# # end.time <- Sys.time()
# #time.taken_fixed <- end.time - start.time
# # Assume no misreporting
# fitSim_noMis<- parLapply(cl, setupOut,
#                          function(x){try(sam.fit_cp(x$dat, x$conf, x$par,
#                                             map = list("logScale" = factor(rep(NA, length(x$par$logScale))))))})
# stopCluster(cl) #shut down nodes


## Error handling #####
# Exclude TMB fails
fitSimAccept_random <- fitSim_random[!(sapply(fitSim_random, class) == "try-error")]
simOutAccept_random <- simOut[!(sapply(fitSim_random, class) == "try-error")]
#fitSimAccept_fixed <- fitSim_fixed[!(sapply(fitSim_fixed, class) == "try-error")]
#simOutAccept_fixed <- simOut[!(sapply(fitSim_fixed, class) == "try-error")]
fitSimAccept_noMis <- fitSim_noMis[!(sapply(fitSim_noMis, class) == "try-error")]
simOutAccept_noMis <- simOut[!(sapply(fitSim_noMis, class) == "try-error")]
# Exclude non-convergences
x <- unlist(sapply(fitSimAccept_random, function (x) x[[6]][3]))
fitSimAccept_random <- fitSimAccept_random[x != 1]
simOutAccept_random <- simOutAccept_random[x != 1]  
nRepAccept_random <- length(fitSimAccept_random)
# x <- unlist(sapply(fitSimAccept_fixed, function (x) x[[6]][3]))
# fitSimAccept_fixed <- fitSimAccept_fixed[x != 1]
# simOutAccept_fixed <- simOutAccept_fixed[x != 1]
# nRepAccept_fixed <- length(fitSimAccept_fixed)
x <- unlist(sapply(fitSimAccept_noMis, function (x) x[[6]][3]))
fitSimAccept_noMis <- fitSimAccept_noMis[x != 1]
simOutAccept_noMis <- simOutAccept_noMis[x != 1]  
nRepAccept_noMis <- length(fitSimAccept_noMis)

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
errNF_random <- calcNFTsError(fitSimAccept_random, simOutAccept_random)
errCSSB_random  <- calcCSSBError(fitSimAccept_random, simOutAccept_random)
err_random <- rbind(errNF_random, errCSSB_random)

# errNF_fixed <- calcNFTsError(fitSimAccept_fixed, simOutAccept_fixed)
# errCSSB_fixed  <- calcCSSBError(fitSimAccept_fixed, simOutAccept_fixed)
# err_fixed <- rbind(errNF_fixed, errCSSB_fixed)

errNF_noMis <- calcNFTsError(fitSimAccept_noMis, simOutAccept_noMis)
errCSSB_noMis  <- calcCSSBError(fitSimAccept_noMis, simOutAccept_noMis)
err_noMis <- rbind(errNF_noMis, errCSSB_noMis)

plotTsError(err_random, noYears = fitSimAccept_random[[1]]$data$noYears)
#plotTsError(err_fixed, noYears = fitSimAccept_fixed[[1]]$data$noYears)
plotTsError(err_noMis, noYears = fitSimAccept_noMis[[1]]$data$noYears)

# plotTsMeanError(err_random, nRepAccept_random)
# plotTsMeanError(err_fixed, nRepAccept_fixed)
# plotTsMeanError(err_noMis, nRepAccept_noMis)

# Plot parameters true vs fit
plotPars(fitSimAccept_random, simOutAccept_random)
#plotPars(fitSimAccept_fixed, simOutAccept_fixed)
plotPars(fitSimAccept_noMis, simOutAccept_noMis)

# Plot observed catch vs true catch vs estimated catch in N (not MT) because
# the model fit is in N, and we want to check to see how well it estiamtes the
# bias which occurs there.
# plotC_N(simOut = simOutAccept[[1]],
#       fit = fitSimAccept[[1]])

## Replicate using simulate.sam() ######################
# fitRealAdjusted <- fitReal
# fitRealAdjusted$pl$logSdLogFsta <- # Adjust (reduce) process error
#   #(c(0.1, rep(0.33, length(fitReal$pl$logSdLogFsta)-1)) * exp(fitReal$pl$logSdLogFsta)) %>%
#   (c(1, rep(1, length(fitReal$pl$logSdLogFsta)-1)) * exp(fitReal$pl$logSdLogFsta)) %>%
#   log
# 
# fitRealAdjusted$pl$logSdLogN[(fitRealAdjusted$conf$keyVarLogN + 1)] <-
#   fitRealAdjusted$pl$logSdLogN[(fitRealAdjusted$conf$keyVarLogN + 1)] %>%
#   exp %>% "*"(1) %>% log
# 
# fitRealAdjusted$pl$logSdLogObs <-
#  fitRealAdjusted$pl$logSdLogObs %>% exp %>% "*"(1) %>% log
# 
# set.seed(12)
# simOutSAM <- stockassessment:::simulate.sam(fitRealAdjusted, nsim = nRep, full.data = TRUE)
# set.seed(12) # Need to do this in order to get N & F and the data needed to run sam.fit to match
# simOutSAM4error <- stockassessment:::simulate.sam(fitRealAdjusted, 
#                                                   nsim = nRep, 
#                                                   full.data = FALSE)
# 
# # make sure the observations are the same
# #for (i in 1:nRep) print(all(simOutSAM[[i]]$logobs == simOutSAM4error[[i]]$logobs))
# 
# cl <- makeCluster(detectCores() - 1) #setup nodes for parallel
# clusterExport(cl, c("fitReal"))
# clusterEvalQ(cl, {library(stockassessment)}) #load stockassessment to each node
# fitSimSAM <- parLapply(cl, simOutSAM,
#                        function(x){try(sam.fit(data = x, conf = fitReal$conf, 
#                                                par = defpar(dat = x, 
#                                                             conf = fitReal$conf)))})
# stopCluster(cl) #shut down nodes
# 
# ## Error handling #####
# for(i in 1:nRep) simOutSAM[[i]]$trueParams <- fitRealAdjusted
# # Exclude TMB fails
# fitSimSAMAccept <- fitSimSAM[!(sapply(fitSimSAM, class) == "try-error")]
# simOutSAMAccept <- simOutSAM[!(sapply(fitSimSAM, class) == "try-error")]
# simOutSAM4errorAccept <- simOutSAM4error[!(sapply(fitSimSAM, class) == "try-error")]
# # Exclude non-convergences
# x <- unlist(sapply(fitSimSAMAccept, function (x) x[[6]][3]))
# fitSimSAMAccept <- fitSimSAMAccept[x != 1]
# simOutSAMAccept <- simOutSAMAccept[x != 1]
# simOutSAM4errorAccept <- simOutSAM4errorAccept[x != 1]
# nRepSAMAccept <- length(fitSimSAMAccept)
# 
# # Plot parameter error
# plotPars(fitSimSAMAccept, simOutSAMAccept)
# 
# plotPars(fitSimAccept, simOutAccept) #comapre to mine
# 
# # Calculate random effect error
# errNFSAM <- calcNFTsErrorSAM(fitSimSAMAccept, simOutSAM4errorAccept)
# plotTsError(errNFSAM)
# 
# plotTsError(err) # compare to mine
# # Plot one
# plotN(simOut = simOutSAM4errorAccept[[3]],
#       fit = fitSimSAMAccept[[3]])
# 
# plotN(simOut = simOutAccept[[3]], #compare to mine
#       fit = fitSimAccept[[3]])
#
# # plotF(simOut = simOutSAM4errorAccept[[1]],
# #      fit = fitSimSAMAccept[[1]])



# Some plots to inspect simulaiton output from the two sources
# mean(c(sapply(simOutSAM4error, function(x) x$logN[3,])))
# mean(c(sapply(simOut, function(x) x$trueParams$pl$logN[3,])))
# 
# mean(c(sapply(simOutSAM, function(x) x$logobs[x$aux[, "fleet"] == 1])))
# mean(c(sapply(simOut, function(x) x$logCobs_N)))

# mean(c(sapply(fitSimAccept, function(x) x$data$logobs[x$data$aux[, "fleet"] == 1])))
# mean(c(sapply(simOutAccept, function(x) x$logCobs_N)))
# 
# c(sapply(fitSimAccept, function(x) x$data$logobs[x$data$aux[, "fleet"] == 2])) 
# c(sapply(fitSimSAMAccept, function(x) x$data$logobs[x$data$aux[, "fleet"] == 2]))
# 
# mean(c(sapply(simOutSAM, function(x) x$logobs[x$aux[, "fleet"] == 2])))
# mean(c(sapply(simOut, function(x) x$logSobs_N[,,2][!is.na(x$logSobs_N[,,2])])))
# 
# mean(c(sapply(simOutSAM, function(x) x$logobs[x$aux[, "fleet"] == 3])))
# mean(c(sapply(simOut, function(x) x$logSobs_N[,,3][!is.na(x$logSobs_N[,,3])])))
# 
# mean(c(sapply(simOutSAM, function(x) x$logobs[x$aux[, "fleet"] == 4])))
# mean(c(sapply(simOut, function(x) x$logSobs_N[,,4][!is.na(x$logSobs_N[,,4])])))

#simOut2plot <- sapply(simOut, function(x) x$logSobs_N[6,,3][!is.na(x$logSobs_N[6,,3])])
#simOut2plot <- sapply(simOut, function(x) x$logCobs_N[6,])
#simOut2plot <- sapply(simOut, function(x) x$trueParams$pl$logN[3,])
# simOut2plot <- sapply(simOut, function(x) x$trueParams$pl$logF[1,])

#simOut2plotSAM <- sapply(simOutSAM, function(x) x$logobs[x$aux[, "fleet"] == 3 & x$aux[, "age"] == 6])
#simOut2plotSAM <- sapply(simOutSAM, function(x) x$logobs[x$aux[, "fleet"] == 1 & x$aux[, "age"] == 6])
#simOut2plotSAM <- sapply(simOutSAM4error, function(x) x$logN[3,])
# simOut2plotSAM <- sapply(simOutSAM4error, function(x) x$logF[1,])
# 
# 
# df2plot <-
#   data.frame(simOut2plot, source = "mine", index = 1:nrow(simOut2plot)) %>%
#   tidyr::gather(variable, value, -source, -index) %>%
#   rbind(data.frame(simOut2plotSAM, source = "SAM", index = 1:nrow(simOut2plot)) %>%
#           tidyr::gather(variable, value, -source, -index)) %>%
#   dplyr::mutate(variable = as.numeric(substr(variable, 2, 10))) %>%
#   dplyr::rename(replicate = variable) %>%
#   dplyr::group_by(source, index) %>%
#   dplyr::summarise(value_mean = mean(value),
#                    value_hi   = quantile(value, .95),
#                    value_low  = quantile(value, .05))
# 
# ggplot(df2plot,# %>% dplyr::filter(index > 10),
#        aes(x = index)) +
#   geom_line(aes(y = value_mean)) +
#   facet_wrap(~source)

#
# mean(sapply(fitSimAccept, function(x) cor(x$pl$logF[2,], x$pl$logN[2,])))
# mean(sapply(fitSimSAMAccept, function(x) cor(x$pl$logF[2,], x$pl$logN[2,])))
#
#
# plot(x<-simOutSAM4error[[2]]$logobs)
# plot(y<-fitSim[[2]]$data$logobs)
# 
# hist(sapply(simOutSAM4error, function(x) mean(x$logF %>% c())))
# hist(sapply(simOut, function(x) mean(x$trueParams$pl$logF[1:2,] %>% c())))
# 
# sapply(simOutSAM4error, function(x) apply(x$logF, 1, FUN = var))
# sapply(simOut, function(x) apply(x$trueParams$pl$logF[1:2,], 1, FUN = var))
# 
# 
# hist(sapply(simOutSAM4error, function(x) mean(x$logN %>% c())))
# hist(sapply(simOut, function(x) mean(x$trueParams$pl$logN %>% c())))
# 
# hist(simOut[[1]]$logCobs_N)
# hist(simOutSAM4error[[1]]$logobs[simOutSAM[[1]]$aux[,"fleet"] == 1])
# 
# hist(simOut[[2]]$logSobs_N)
# hist(simOutSAM4error[[2]]$logobs[simOutSAM[[1]]$aux[,"fleet"] != 1])

# Check if MLE is higher when using true logScale params
# par2use <- fitSimAccept[[1]]$opt$par
# par2use[which(names(par2use) == "logScale")] <- simOutAccept[[1]]$trueParams$pl$logScale
# par2use[which(names(par2use) == "logFpar")] <- simOutAccept[[1]]$trueParams$pl$logFpar
# par2use[which(names(par2use) == "logSdLogFsta")] <- simOutAccept[[1]]$trueParams$pl$logSdLogFsta
# par2use[which(names(par2use) == "logSdLogN")] <- simOutAccept[[1]]$trueParams$pl$logSdLogN
# par2use[which(names(par2use) == "logSdLogObs")] <- simOutAccept[[1]]$trueParams$pl$logSdLogObs
# fitSimAccept[[1]]$obj$fn(fitSimAccept[[1]]$opt$par)
# fitSimAccept[[1]]$obj$fn(par2use)


