library(stockassessment)

setwd("/Users/charlesperretti/Projects/state-space-tests/atlherring_example/")


cn <- read.ices("Herrcn.dat") # catch abundace-at-age
cw <- read.ices("Herrcw.dat") # catch mean weight-at-age
dw <- cw # discards mean weight-at-age (using catch mean weight-at-age for now)
lw <- cw # landings mean weight-at-age (using catch mean weight-at-age for now)
pf <- read.ices("Herrpf.dat") # proportion of f before spawning
lf <- pf; lf[,] <- 1 # MAYBE proportion of landings before spawning??? (set to 1 for now)
mo <- read.ices("Herrmo.dat") # maturity-at-age ogive
nm <- read.ices("Herrnm.dat") # natural mortality-at-age
pm <- read.ices("Herrpm.dat") # proportion of m before spawning
sw <- read.ices("Herrsw.dat") # stock weight-at-age (kg)
surveys <- read.ices("Herrsurvey_BigSep_NoAcoust.dat") #surveys

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
conf <- loadConf(dat = dat_atl, file = "ModelConf.txt")

par <- defpar(dat_atl, conf) # some default starting values

fit <- sam.fit(dat_atl, conf, par) # fit the model

modelTable <- modeltable(fit) # AIC and # of params

ssbplot(fit)
recplot(fit)
tsbplot(fit)
catchplot(fit)

parplot(fit)

obscov(fit)

corplot(fit)

#res <- residuals(fit)
#plot(res)

