## WARNING THIS CODE CURRENTLY CRASHES R
library(stockassessment)

datdirect <- "/Users/charlesperretti/Projects/state-space-tests/atlherring_example/"

cn <- read.ices(paste(datdirect,"Herrcn.dat",sep="/"))
cw <- read.ices(paste(datdirect,"Herrcw.dat",sep="/"))
mo <- read.ices(paste(datdirect,"Herrmo.dat",sep="/"))
nm <- read.ices(paste(datdirect,"Herrnm.dat",sep="/"))
pf <- read.ices(paste(datdirect,"Herrpf.dat",sep="/"))
pm <- read.ices(paste(datdirect,"Herrpm.dat",sep="/"))
sw <- read.ices(paste(datdirect,"Herrsw.dat",sep="/"))
surveys <- read.ices(paste(datdirect,"Herrsurvey_BigSep_NoAcoust.dat",sep="/"))

# setup the data as needed for SAM
dat_atl <- setup.sam.data(surveys = surveys,
                          residual.fleet = cn,
                          prop.mature = mo,
                          stock.mean.weight = sw,
                          prop.f = pf,
                          prop.m = pm,
                          natural.mortality = nm,
                          catch.mean.weight = cw)

conf <- defcon(dat_atl) # default configuration

par<-defpar(dat_atl, conf) # some default starting values

fit<-sam.fit(dat_atl, conf, par) # fit the model
