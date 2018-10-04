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

#Creates graphic axis names
#getnames<-function(confa,agesa,looptoa,looptob,namesa) 
#confa is the section of configuration file wanted
#agesa is vector of ages c(1,2,...)
#looptoa is number of unique states, couplings, or whatever in the given conf
#looptob is number of rows in the conf file that are relevant (e.g., first row of keyLogFpar is for fishery, and so looptob isn't nrows, it's number of surveys)
#names are exactly that, the names for the graphic (the fxn ties these names to proper age couplings)
source("functions.R")
ages <- seq(1, ncol(cn), 1)
qplotnames <- getnames(confa = conf$keyLogFpar[2:nrow(conf$keyLogFpar),],
                       agesa = ages,
                       looptoa = (max(conf$keyLogFpar)+1),
                       looptob = (length(names(surveys))),
                       namesa = names(surveys))
Fpronames <- getnames(confa = t(as.matrix(conf$keyVarF[1,])),
                      agesa = ages,
                      looptoa = (max(conf$keyVarF)+1),
                      looptob = 1,
                      namesa = "SD LogF process")
Npronames <- getnames(confa=t(as.matrix(conf$keyVarLogN)),
                      agesa=ages,
                      looptoa=(max(conf$keyVarLogN)+1),
                      looptob=1,
                      namesa="SD LogN process")
Obsnames <- getnames(confa=conf$keyVarObs, 
                     agesa=ages, 
                     looptoa=(max(conf$keyVarObs)+1),
                     looptob=nrow(conf$keyVarObs), 
                     namesa=c("Catch",names(surveys)))
Obsnames <- paste("SD", Obsnames)
varplotnames <- c(Fpronames,Npronames,Obsnames)
Fstanames <- getnames(confa=t(as.matrix(conf$keyLogFsta[1,])),
                      agesa=ages,
                      looptoa=(max(conf$keyLogFsta)+1),
                      looptob=1,namesa="F")

if(conf$noScaledYears>0){
  catchmultnames<-getnames(confa=as.matrix(conf$keyParScaledYA),
                           agesa=ages,
                           looptoa=(max(conf$keyParScaledYA)+1),
                           looptob=(max(conf$keyParScaledYA)+1),
                           namesa=rep("Scalar",max(conf$keyParScaledYA)+1))
  #to assign years to catch scalar names for plots
  catchmults<-data.frame(t(conf$keyScaledYears),conf$keyParScaledYA)
  names(catchmults)<-c("Year",ages)
  for(m in 1:(max(conf$keyParScaledYA)+1)) {
    myearnum<-c()
    for(r in 1:nrow(conf$keyParScaledYA)){
      myears<-sapply(data.frame(conf$keyParScaledYA)[r,],function(x) any(x==(m-1)))
      if(myears[m]) {
        myearnum<-c(myearnum,catchmults[r,1])
      }
    }
    if(m==1) {
      yearlabels<-paste(min(myearnum),max(myearnum),sep="to")
    } else {
      yearlabels<-c(yearlabels,paste(min(myearnum),max(myearnum),sep="to"))
    }
  }
  catchmultnames<-paste(catchmultnames,yearlabels,sep=" ")
  if(length(catchmultnames)==2*(max(conf$keyParScaledYA)+1)) {
    catchmultnames<-catchmultnames[1:(max(conf$keyParScaledYA)+1)]
  }
} else {
  catchmultnames<-NA
}

#make the plots and save to pdf in Run sub-directory
run <- "Run13_2017_changefbar"
plotfxn(afit = fit,
        datdirect = getwd(),
        run = run,
        confa = conf,
        qplotnamesa = qplotnames,
        varplotnamesa = varplotnames,
        fstanamesa = Fstanames,
        catchmultnamesa = catchmultnames,
        retroyrs = 5) #as.vector(seq(2011,2014)))

#does likelihood profile and makes some plots; puts results in run directory
#may not be working correctly
likeliprof = FALSE
if(likeliprof){ likeliprofile(fit = fit, datdirect = getwd(), run = run) }

forec = FALSE
if(forec){
  fit.forecast<-forecast(fit, nosim=1000, 
                         overwriteSelYears = seq(2012,2017),
                         fscale = rep(NA,50),
                         catchval = rep(NA,50),
                         fval = rep(0.39,50),
                         label = "F40Forecast",
                         deterministic = F,
                         ave.years = seq(2012,2017),
                         rec.years = (min(fit$data$years):max(fit$data$years)))
  windows(15,10)
  plot(fit.forecast)
  #mean(fit.forecast[[40]]$ssb)
}

if(FALSE){
  select<-matrix(nrow=nrow(faytable(fit)),ncol=ncol(faytable(fit)))
  for(s in 1:nrow(faytable(fit))){
    select[s,]<-faytable(fit)[s,]/max(faytable(fit)[s,])
  }
  
}
