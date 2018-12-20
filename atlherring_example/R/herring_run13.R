library(stockassessment)
library(dplyr)

#cn <- read.ices("../data/Herrcn.dat") # catch abundace-at-age
cn <- read.ices("../../simulation/sim_data/catch.dat") # catch abundace-at-age
cw <- read.ices("../data/Herrcw.dat") # catch mean weight-at-age
dw <- cw # discards mean weight-at-age (using catch mean weight-at-age for now)
lw <- cw # landings mean weight-at-age (using catch mean weight-at-age for now)
pf <- read.ices("../data/Herrpf.dat") # proportion of f before spawning
lf <- pf; lf[,] <- 1 # fraction of catch that is landed (set to 1 for now)
mo <- read.ices("../data/Herrmo.dat") # maturity-at-age ogive
nm <- read.ices("../data/Herrnm.dat") # natural mortality-at-age
pm <- read.ices("../data/Herrpm.dat") # proportion of m before spawning
sw <- read.ices("../data/Herrsw.dat") # stock weight-at-age (kg)

cn <- read.ices("../data/Herrcn.dat") # catch abundace-at-age
#cn <- read.ices("../../simulation/sim_data/catch.dat") # catch abundace-at-age
#surveys <- read.ices("../data/Herrsurvey_BigSep_NoAcoust.dat") #surveys
surveys <- read.ices("../../simulation/sim_data/surveys.dat") #surveys

# setup the data as needed for SAM
dat_atl <- setup.sam.data(surveys = surveys,
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
conf <- loadConf(dat = dat_atl, file = "../ModelConf_simple.txt")

par <- defpar(dat_atl, conf) # some default starting values

fitHer <- sam.fit(dat_atl, conf, par, sim.condRE = FALSE) # fit the model

save(list = "fitHer", file = "../output/fitHerSimple.Rdata")

modelTable <- modeltable(fitHer) # AIC and # of params


# Make plots
ssbplot(fitHer)
fbarplot(fitHer)
recplot(fitHer)
catchplot(fitHer)


# res <- residuals(fitHer)
# corplot(res)
# stockassessment::srplot(fitHer)
# 
# pdf(file = "parameter_plot.pdf")
# parplot(fitHer)
# dev.off()
# 
# 
# #### Simulate from the fitted model and re-fit ####
# simdat <- simulate(fitHer, seed = 1, nsim = 1)[[1]]
# 
# simfit <- sam.fit(simdat, conf, par)
# 
# ssbplot(fitHer)
# ssbplot(simfit, add = TRUE, col = "blue")
# 
# fbarplot(fitHer, partial = FALSE)
# fbarplot(simfit, partial = FALSE, add = TRUE, col = "blue")
# 
# recplot(fitHer)
# recplot(simfit, add = TRUE, col = "blue")
# 
# catchplot(fitHer)
# catchplot(simfit, add = TRUE, col = "blue")
# 
# 
# # Run lots of simulations and fits
# simlist <- simulate(fitHer, seed=1, nsim=10)
# library(parallel)
# no_cores <- detectCores() - 1 #how many cores can we use
# if( no_cores>2 ) no_cores <- 2 # Cran check does not allow us to use more than two 
# cl <- makeCluster(no_cores) #set up some number of nodes
# 
# clusterExport(cl, c("conf", "par")) #send these objects to each node
# clusterEvalQ(cl, {library(stockassessment)}) #load the package to each node
# simfitslist <- parLapply(cl, simlist, function(x){sam.fit(x, conf, par)}) #do sam.fit to element
# stopCluster(cl) #shut it down
# 
# ssbplot(fitHer)#the original data
# trash <- lapply(simfitslist, function(x){ssbplot(x, ci=FALSE, add=TRUE)})
# 
# fbarplot(fitHer, partial = FALSE)
# trash <- lapply(simfitslist, function(x){fbarplot(x, ci=FALSE, partial = FALSE, add=TRUE)})
# 
# recplot(fitHer)
# trash <- lapply(simfitslist, function(x){recplot(x, ci=FALSE, add=TRUE)})
# 
# catchplot(fitHer)
# trash <- lapply(simfitslist, function(x){catchplot(x, ci=FALSE, add=TRUE)})
# 
# # Extract parameter estimates from the true model
# # remove estiamtes of missing values since they are not 
# # in the simulated dataset and we want the two to match up
# df_true <- 
#   summary(fitHer$sdrep) %>%
#   as.data.frame() %>%
#   dplyr::mutate(variable = row.names(.)) %>%
#   dplyr::filter(!grepl('missing', variable))
# 
# # Extract parameter estimates from the first replicate of the
# # simulated data
# df_replicate1 <-
#   summary(simfitslist[[6]]$sdrep) %>%
#   as.data.frame() %>%
#   dplyr::mutate(variable = row.names(.))
# 
# # Make sure the variables match
# all(df_true$variable == df_replicate1$variable)
# 
# 
# # Calculate proportion of true parameters that fall within the 
# # 95% interval of the estimates
# mean(df_true$Estimate < df_replicate1$Estimate + 1.96*df_replicate1$`Std. Error` &
#      df_true$Estimate > df_replicate1$Estimate - 1.96*df_replicate1$`Std. Error`)
# 
# 

