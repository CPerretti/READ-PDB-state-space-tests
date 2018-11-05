# sim3.r
# testing to see if can fake out stockassessment to simulate data with different process error

library(dplyr) 
library(ggplot2)
library(stockassessment)

# Use atl herring fit to set up simulation
load("../../atlherring_example/output/fit.Rdata")

# use herring configuration and default parameter intitial guesses 
conf <- fit$conf
par <- defpar(fit$data, conf)

# fit data with sim.condRE = FALSE so that simulations will use process error
fitF <- sam.fit(fit$data, conf, par, sim.condRE = FALSE)

# simulate data with estimated process error to show huge variability in N and F at age 5
simdatF <- stockassessment:::simulate.sam(fitF, nsim = 10, full.data = FALSE)

# make dataframe of base case
df_sim <- data.frame() 
for (i in 1:length(simdatF)) {
  df_sim <-
    rbind(df_sim, data.frame(replicate = as.factor(i),
                             case = "Base",
                             year = fit$data$years,
                             N5 = exp(simdatF[[i]]$logN[5,]),
                             F5 = exp(simdatF[[i]]$logF[5,])))
}

# modify fitF to reduce process error to extremely small values
modfitF <- fitF
modfitF$pl$logSdLogN <- rep(-100, 2)
modfitF$pl$logSdLogFsta <- rep(-100, 4)

# simulate data with reduced process error
simdatF <- stockassessment:::simulate.sam(modfitF, nsim = 10, full.data = FALSE)

# add to dataframe 
for (i in 1:length(simdatF)) {
  df_sim <-
    rbind(df_sim, data.frame(replicate = as.factor(i),
                             case = "Modified",
                             year = fit$data$years,
                             N5 = exp(simdatF[[i]]$logN[5,]),
                             F5 = exp(simdatF[[i]]$logF[5,])))
}

# compare Base and Modified
ggplot(df_sim, aes(x=year, y=N5, color=case)) +
  geom_line() +
  facet_wrap(replicate~., scales="free_y") +
  ggtitle("N age 5") +
  theme_bw()

ggplot(df_sim, aes(x=year, y=F5, color=case)) +
  geom_line() +
  facet_wrap(replicate~., scales="free_y") +
  ggtitle("F age 5") +
  theme_bw()

# all 10 replicated of Modified equilibrate at same value 
# reducing process error to -100 has removed random walk, but went too far

# modify fitF to reduce process error proportionally
mod2fitF <- fitF
mod2fitF$pl$logSdLogN <- fitF$pl$logSdLogN / 100
mod2fitF$pl$logSdLogFsta <- fitF$pl$logSdLogFsta / 100

# simulate data with reduced process error
simdatF <- stockassessment:::simulate.sam(mod2fitF, nsim = 10, full.data = FALSE)

# add to dataframe 
for (i in 1:length(simdatF)) {
  df_sim <-
    rbind(df_sim, data.frame(replicate = as.factor(i),
                             case = "Mod2",
                             year = fit$data$years,
                             N5 = exp(simdatF[[i]]$logN[5,]),
                             F5 = exp(simdatF[[i]]$logF[5,])))
}

# compare Base and Modified
ggplot(df_sim, aes(x=year, y=N5, color=case)) +
  geom_line() +
  facet_wrap(replicate~., scales="free_y") +
  ggtitle("N age 5") +
  theme_bw()

ggplot(df_sim, aes(x=year, y=F5, color=case)) +
  geom_line() +
  facet_wrap(replicate~., scales="free_y") +
  ggtitle("F age 5") +
  theme_bw()

# nope still getting too much process error in Mod2 results
