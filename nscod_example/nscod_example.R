
# filestoget <- c("cn.dat", "cw.dat", "dw.dat", "lf.dat", "lw.dat", 
#                 "mo.dat", "nm.dat", "pf.dat", "pm.dat", "sw.dat", 
#                 "survey.dat")
# url <- "https://raw.githubusercontent.com/fishfollower/SAM/master/stockassessment/tests/nscod/"
# d <- lapply(filestoget, function(f)download.file(paste(url,f,sep=""), f))

library(stockassessment)
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

# Fix for my simulation script because commerical landings
# are shorter than survey landings
# mo <- mo[rownames(mo) < 2015,]
#nm <- nm[rownames(nm) < 2015,]
#pf <- pf[rownames(pf) < 2015,]
#pm <- pm[rownames(pm) < 2015,]
#sw <- sw[rownames(sw) < 2015,]
#surveys[[1]] <- surveys[[1]][rownames(surveys[[1]]) < 2015,]

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

# Make a couple of changes to conf so it matches simulation
conf$corFlag <- 0


par <- defpar(dat,conf)

#par$logFpar <- rep(0,9)

fitNScod <- sam.fit(dat,conf,par) 

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
