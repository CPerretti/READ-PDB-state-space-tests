library(stockassessment)

datdirect <- "/Users/charlesperretti/Projects/state-space-tests/atlherring_example/"

cn <- read.ices(paste(datdirect,"Herrcn.dat",sep="/")) # catch abundace-at-age
cw <- read.ices(paste(datdirect,"Herrcw.dat",sep="/")) # catch mean weight-at-age
dw <- cw # discards mean weight-at-age (using catch mean weight-at-age for now)
lw <- cw # landings mean weight-at-age (using catch mean weight-at-age for now)
pf <- read.ices(paste(datdirect,"Herrpf.dat",sep="/")) # proportion of f before spawning
lf <- pf; lf[,] <- 1 # MAYBE proportion of landings before spawning??? (set to 1 for now)
mo <- read.ices(paste(datdirect,"Herrmo.dat",sep="/")) # maturity-at-age ogive
nm <- read.ices(paste(datdirect,"Herrnm.dat",sep="/")) # natural mortality-at-age
pm <- read.ices(paste(datdirect,"Herrpm.dat",sep="/")) # proportion of m before spawning
sw <- read.ices(paste(datdirect,"Herrsw.dat",sep="/")) # stock weight-at-age (kg)
surveys <- read.ices(paste(datdirect,"Herrsurvey_BigSep_NoAcoust.dat",sep="/")) #surveys

# setup the data as needed for SAM
dat_atl <- setup.sam.data(surveys = surveys,
                          residual.fleet = cn,
                          prop.mature = mo,
                          stock.mean.weight = sw,
                          dis.mean.weight=dw,
                          land.mean.weight=lw,
                          land.frac=lf,
                          prop.f = pf,
                          prop.m = pm,
                          natural.mortality = nm,
                          catch.mean.weight = cw)

conf <- defcon(dat_atl) # default configuration

par<-defpar(dat_atl, conf) # some default starting values

fit<-sam.fit(dat_atl, conf, par) # fit the model
