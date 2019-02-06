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

#set.seed(321) # for reproducibility

# How many simulation replicates to do

nRep <- 5#0

# Generate simulation replicates
simOut <- list()
for (i in 1:nRep) {
  simOut[[i]] <- sim(fit = fitReal)  
}

## Plot an example true vs observed vs *real data* fit ########

## (1) N-at-age (1000s)
plotN(simOut = simOut[[1]],
      fit = fitReal)

## (2) F-at-age
plotF(simOut = simOut[[1]],
      fit = fitReal) 

## (3) Catch (mt)
plotC(simOut = simOut[[1]],
      fit = fitReal)

## (4) Survey (1000s)
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
  setupOut[[i]] <- setupModel(conf = fitReal$conf, example_dir = example_dir)
  
  #setupOut[[i]]$par$logFpar <- simOut[[1]]$trueParams$pl$logFpar #Just to set map on true values.
}

# Fit model to replicates in parallel
#sam.fit(setupOut[[1]]$dat, setupOut[[1]]$conf, setupOut[[1]]$par)
cl <- makeCluster(detectCores() - 1) #setup nodes for parallel
clusterEvalQ(cl, {library(stockassessment)}) #load stockassessment to each node
fitSim <- parLapply(cl, setupOut,
                    function(x){try(sam.fit(x$dat, x$conf, x$par#,
                                            #map = list("logFpar" = factor(rep(NA, length(x$par$logFpar))))
                                            ))})
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

# Save output
suffix <- paste0(Sys.time(), ".Rdata")
save(list = "fitSim", file = paste0("./output/fitSim", suffix))
save(list = "fitSimAccept", file = paste0("./output/fitSimAccept", suffix))
save(list = "simOut", file = paste0("./output/simOut", suffix))
save(list = "simOutAccept", file = paste0("./output/simOutAccept", suffix))

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

# Plot error of time series estimates
errNF <- calcNFTsError(fitSimAccept, simOutAccept)
errC  <- calcCatchError(fitSimAccept, simOutAccept)
err <- rbind(errNF, errC)

plotTsError(err)
plotTsMeanError(err)

# Plot parameters true vs fit
plotPars(fitSimAccept, simOutAccept)

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
