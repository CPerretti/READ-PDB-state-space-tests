
# filestoget <- c("cn.dat", "cw.dat", "dw.dat", "lf.dat", "lw.dat", 
#                 "mo.dat", "nm.dat", "pf.dat", "pm.dat", "sw.dat", 
#                 "survey.dat")
# url <- "https://raw.githubusercontent.com/fishfollower/SAM/master/stockassessment/tests/nscod/"
# d <- lapply(filestoget, function(f)download.file(paste(url,f,sep=""), f))

library(stockassessment)
# note: For my code to run I deleted 2015 data from the survey files so everything ends in 2014.
cn <- read.ices("./data/cn.dat") # catch-at-age (thousands)
cw <- read.ices("./data/cw.dat") # catch weight-at-age (kg)
dw <- read.ices("./data/dw.dat") #
lf <- read.ices("./data/lf.dat") #
lw <- read.ices("./data/lw.dat")
mo <- read.ices("./data/mo.dat")
nm <- read.ices("./data/nm.dat")
pf <- read.ices("./data/pf.dat")
pm <- read.ices("./data/pm.dat")
sw <- read.ices("./data/sw.dat")
surveys <- read.ices("./data/survey.dat")

dat <- setup.sam.data(surveys=surveys,
                      residual.fleet=cn, 
                      prop.mature=mo, 
                      stock.mean.weight=sw, 
                      catch.mean.weight=cw, 
                      dis.mean.weight=dw,
                      land.mean.weight=lw,
                      prop.f=pf, 
                      prop.m=pm, 
                      natural.mortality=nm, 
                      land.frac=lf)



conf <- defcon(dat)

# Turn off ar1 correlation of F so it matches my simulation
conf$corFlag <- 0


par <- defpar(dat,conf)

fitNScod <- sam.fit(dat,conf,par, sim.condRE = FALSE) 

#save("fitNScod", file = "./fitNScod.Rdata")


ssbplot(fitNScod)

fbarplot(fitNScod)

recplot(fitNScod)

catchplot(fitNScod)

# res <- residuals(fitNScod)
# plot(res)
# 
# resp <- procres(fitNScod)
# plot(resp)
# 
# retro <- retro(fitNScod, year=10)
# plot(retro)
