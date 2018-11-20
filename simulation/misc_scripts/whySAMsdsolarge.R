# trying to determine what is causing process errors to be so large for herring

library(dplyr)
library(ggplot2)
library(tidyr)
library(stockassessment)

# Use atl herring fit to set up simulation
load("../../atlherring_example/output/fitHer.Rdata")

# get process errors
f <- exp(fitHer$pl$logF[(fitHer$conf$keyLogFsta[1,] + 1),])
m <- t(fitHer$data$natMor)
z <- f + m
sdLogN <- exp(fitHer$pl$logSdLogN[(fitHer$conf$keyVarLogN + 1)])
nA <- ncol(fitHer$data$propF) # number of age-classes
nT <- fitHer$data$noYears # length of time series
errProLogN <- matrix(data = NA,
                     nrow = nA, 
                     ncol = nT-1)
errProLogN[1, ] <- fitHer$pl$logN[1, 2:nT] - fitHer$pl$logN[1, 1:(nT-1)]
errProLogN[-c(1, nA), ] <-  fitHer$pl$logN[-c(1, nA), 2:nT] -
  fitHer$pl$logN[-c(nA-1, nA), 1:(nT-1)] + z[-c(nA-1, nA), 1:(nT-1)]
errProLogN[nA, ] <- fitHer$pl$logN[nA, 2:nT] -  
  log(exp(fitHer$pl$logN[nA-1, 1:(nT-1)]) *
     exp(-z[nA-1, 1:(nT-1)]) +
     exp(fitHer$pl$logN[nA, 1:(nT-1)]) *
     exp(-z[nA, 1:(nT-1)]))


sd(errProLogN[2:8,]) # sd of log process errors across ages 2-8 - still doesn't match SAM sd

# compute standard deviations of process errors at age and compare to SAM estimates
derivedsdLogN <-apply(errProLogN, 1, sd)
sddf <- data.frame(age = 1:nA,
                   derivedsdLogN = derivedsdLogN,
                   SAMsdLogN = sdLogN)
sddf
# why are SAM estimates larger than actual variability among ages?

# distribution of process errors
pe <- errProLogN %>%
  data.frame(.) %>%
  mutate(age = rownames(.)) %>%
  gather(key = year, value = procerr, 1:52)

ggplot(filter(pe, age != 1), aes(x=procerr)) +
  geom_density() +
  theme_bw()

ggplot(filter(pe, age != 1), aes(x=age, y=procerr)) +
  geom_boxplot() +
  geom_jitter(width=0.2, color="blue", alpha=0.2) +
  theme_bw()
