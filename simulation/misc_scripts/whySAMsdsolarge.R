# trying to determine what is causing process errors to be so large for herring

library(dplyr)
library(ggplot2)
library(tidyr)
library(stockassessment)

# Use atl herring fit to set up simulation
load("../atlherring_example/output/fitHer.Rdata")

# get process errors
f <- exp(fitHer$pl$logF[(fitHer$conf$keyLogFsta[1,] + 1),])
m <- t(fitHer$data$natMor)
z <- f + m
sdLogN <- exp(fitHer$pl$logSdLogN[(fitHer$conf$keyVarLogN + 1)])
nA <- ncol(fitHer$data$propF) # number of age-classes
nT <- fitHer$data$noYears # length of time series
errPro <- matrix(data = NA,
                 nrow = nA, 
                 ncol = nT-1)
errPro[1, ] <- exp(fitHer$pl$logN[1, 2:nT]) / exp(fitHer$pl$logN[1, 1:(nT-1)])
errPro[-c(1, nA), ] <-  exp(fitHer$pl$logN[-c(1, nA), 2:nT]) / 
  (exp(fitHer$pl$logN[-c(nA-1, nA), 1:(nT-1)]) *
     exp(-z[-c(nA-1, nA), 1:(nT-1)]))
errPro[nA, ] <- exp(fitHer$pl$logN[nA, 2:nT]) / 
  (exp(fitHer$pl$logN[nA-1, 1:(nT-1)]) *
     exp(-z[nA-1, 1:(nT-1)]) +
     exp(fitHer$pl$logN[nA, 1:(nT-1)]) *
     exp(-z[nA, 1:(nT-1)]))

errProLogN <- log(errPro)
#errProLogN

# compute standard deviations of process errors at age and compare to SAM estimates
derivedsdLogN <-apply(errProLogN, 1, sd)
sddf <- data.frame(age = 1:nA,
                   derivedsdLogN = derivedsdLogN,
                   SAMsdLogN = sdLogN)
sddf
# why are SAM estimates larger than actual variability among ages?
