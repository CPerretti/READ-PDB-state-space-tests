library(stockassessment)
library(dplyr)


# CL:
# I suggest not using setwd approach. 
# Instead assume the files are downlaoded from github with that directory structure.
# In RStudio, under the Session tab, select Set Working Directory, then To Source File Location.
# Use relative references from the source file location. 
# For example, if all files in same directory, can just read and write the files.
# If R code in example/R directory and data in example/data, then read.ices("../data/Herrcn.dat").
# This approach will allow all of us to download from github and run the programs without having to change setwd.
# CP:
# Done

cn <- read.ices("../data/Herrcn.dat") # catch abundace-at-age
cw <- read.ices("../data/Herrcw.dat") # catch mean weight-at-age
dw <- cw # discards mean weight-at-age (using catch mean weight-at-age for now)
lw <- cw # landings mean weight-at-age (using catch mean weight-at-age for now)
pf <- read.ices("../data/Herrpf.dat") # proportion of f before spawning
lf <- matrix(NA, nrow = nrow(pf), ncol = ncol(pf))
lf[,] <- 1 # fraction of catch that is landed (set to 1 for now)
mo <- read.ices("../data/Herrmo.dat") # maturity-at-age ogive
nm <- read.ices("../data/Herrnm.dat") # natural mortality-at-age
pm <- read.ices("../data/Herrpm.dat") # proportion of m before spawning
sw <- read.ices("../data/Herrsw.dat") # stock weight-at-age (kg)
surveys <- read.ices("../data/Herrsurvey_BigSep_NoAcoust.dat") #surveys

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
conf <- loadConf(dat = dat_atl, file = "../ModelConf.txt")


par <- defpar(dat_atl, conf) # some default starting values

fit <- sam.fit(dat_atl, conf, par) # fit the model

save(list = "fit", file = "../output/fit.Rdata")

modelTable <- modeltable(fit) # AIC and # of params


# Make plots
ssbplot(fit)
fbarplot(fit)
recplot(fit)
catchplot(fit)


res <- residuals(fit)
corplot(res)
stockassessment::srplot(fit)

pdf(file = "parameter_plot.pdf")
parplot(fit)
dev.off()


#### Simulate from the fitted model and re-fit ####
simdat <- simulate(fit, seed = 1, nsim = 1)[[1]]

simfit <- sam.fit(simdat, conf, par)

ssbplot(fit)
ssbplot(simfit, add = TRUE, col = "blue")

fbarplot(fit, partial = FALSE)
fbarplot(simfit, partial = FALSE, add = TRUE, col = "blue")

recplot(fit)
recplot(simfit, add = TRUE, col = "blue")

catchplot(fit)
catchplot(simfit, add = TRUE, col = "blue")


# Run lots of simulations and fits
simlist <- simulate(fit, seed=1, nsim=10)
library(parallel)
no_cores <- detectCores() - 1 #how many cores can we use
if( no_cores>2 ) no_cores <- 2 # Cran check does not allow us to use more than two 
cl <- makeCluster(no_cores) #set up some number of nodes

clusterExport(cl, c("conf", "par")) #send these objects to each node
clusterEvalQ(cl, {library(stockassessment)}) #load the package to each node
simfitslist <- parLapply(cl, simlist, function(x){sam.fit(x, conf, par)}) #do sam.fit to element
stopCluster(cl) #shut it down

ssbplot(fit)#the original data
trash <- lapply(simfitslist, function(x){ssbplot(x, ci=FALSE, add=TRUE)})

fbarplot(fit, partial = FALSE)
trash <- lapply(simfitslist, function(x){fbarplot(x, ci=FALSE, partial = FALSE, add=TRUE)})

recplot(fit)
trash <- lapply(simfitslist, function(x){recplot(x, ci=FALSE, add=TRUE)})

catchplot(fit)
trash <- lapply(simfitslist, function(x){catchplot(x, ci=FALSE, add=TRUE)})

# Extract parameter estimates from the true model
# remove estiamtes of missing values since they are not 
# in the simulated dataset and we want the two to match up
df_true <- 
  summary(fit$sdrep) %>%
  as.data.frame() %>%
  dplyr::mutate(variable = row.names(.)) %>%
  dplyr::filter(!grepl('missing', variable))

# Extract parameter estimates from the first replicate of the
# simulated data
df_replicate1 <-
  summary(simfitslist[[6]]$sdrep) %>%
  as.data.frame() %>%
  dplyr::mutate(variable = row.names(.))

# Make sure the variables match
all(df_true$variable == df_replicate1$variable)


# Calculate proportion of true parameters that fall within the 
# 95% interval of the estimates
mean(df_true$Estimate < df_replicate1$Estimate + 1.96*df_replicate1$`Std. Error` &
     df_true$Estimate > df_replicate1$Estimate - 1.96*df_replicate1$`Std. Error`)



