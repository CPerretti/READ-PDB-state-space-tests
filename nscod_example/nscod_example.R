
# filestoget <- c("cn.dat", "cw.dat", "dw.dat", "lf.dat", "lw.dat", 
#                 "mo.dat", "nm.dat", "pf.dat", "pm.dat", "sw.dat", 
#                 "survey.dat")
# url <- "https://raw.githubusercontent.com/fishfollower/SAM/master/stockassessment/tests/nscod/"
# d <- lapply(filestoget, function(f)download.file(paste(url,f,sep=""), f))

library(stockassessment)
cn <- read.ices("cn.dat") # catch-at-age (thousands)
cw <- read.ices("cw.dat") # catch weight-at-age (kg)
dw <- read.ices("dw.dat") #
lf <- read.ices("lf.dat") #
lw <- read.ices("lw.dat")
mo <- read.ices("mo.dat")
nm <- read.ices("nm.dat")
pf <- read.ices("pf.dat")
pm <- read.ices("pm.dat")
sw <- read.ices("sw.dat")
surveys <- read.ices("survey.dat")

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

# conf$fbarRange <- c(2,6)
# conf$corFlag <- 1
# conf$keyLogFpar <- 
#   matrix(c(-1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
#            -1,    0,    1,    2,    3,    4,    5,    6,   -1,
#            -1,    7,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
#            8,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1), 
#          nrow=4, byrow=TRUE)


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
