# This script should break R in sam.fit() because dw, lf and lw are not supplied.

filestoget <- c("cn.dat", "cw.dat", "dw.dat", "lf.dat", "lw.dat", 
                "mo.dat", "nm.dat", "pf.dat", "pm.dat", "sw.dat", 
                "survey.dat")
url <- "https://raw.githubusercontent.com/fishfollower/SAM/master/stockassessment/tests/nsher/"
d <- lapply(filestoget, function(f)download.file(paste(url,f,sep=""), f))

library(stockassessment)
cn <- read.ices("cn.dat") # catch-at-age (thousands)
cw <- read.ices("cw.dat") # catch weight-at-age (kg)
#dw <- read.ices("dw.dat") #
#lf <- read.ices("lf.dat") #
#lw <- read.ices("lw.dat")
mo <- read.ices("mo.dat")
nm <- read.ices("nm.dat")
pf <- read.ices("pf.dat")
pm <- read.ices("pm.dat")
sw <- read.ices("sw.dat")
surveys <- read.ices("survey.dat")

dat <- setup.sam.data(surveys=surveys,
                      residual.fleet=cn, 
                      prop.mature=mo, 
                      stock.mean.weight=sw, 
                      catch.mean.weight=cw, 
                      #dis.mean.weight=dw,
                      #land.mean.weight=lw,
                      prop.f=pf, 
                      prop.m=pm, 
                      natural.mortality=nm)#, 
                      #land.frac=lf)

conf <- defcon(dat)

par <- defpar(dat,conf)

fit <- sam.fit(dat,conf,par) 

