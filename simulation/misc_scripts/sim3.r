# sim3.r
# testing to see if can fake out stockassessment to simulate data with different process error

library(dplyr) 
library(ggplot2)
library(stockassessment)

# Use atl herring fit to set up simulation
load("../atlherring_example/output/fitHer.Rdata")

# use herring configuration and default parameter intitial guesses 
conf <- fitHer$conf
par <- defpar(fitHer$data, conf)

# fit data with sim.condRE = FALSE so that simulations will use process error
fitF <- sam.fit(fitHer$data, conf, par, sim.condRE = FALSE)

# simulate data with estimated process error to show huge variability in N and F at age 5
simdatF <- stockassessment:::simulate.sam(fitF, nsim = 10, full.data = FALSE)

# check for diffs in N at age in first year
df_N <- data.frame() 
for (i in 1:length(simdatF)) {
  df_N <-
    rbind(df_N, data.frame(replicate = as.factor(i),
                           case = "Base",
                           age = 1:dim(simdatF[[i]]$logN)[1],  
                           Nyr1 = exp(simdatF[[i]]$logN[, 1]),
                           Nyr30 = exp(simdatF[[i]]$logN[, 30])))
}
df_tab <- df_N %>%
  group_by(age) %>%
  summarize(MeanNyr1 = mean(Nyr1), sdNyr1 = sd(Nyr1),
            MeanNyr30 = mean(Nyr30), sdNyr30 = sd(Nyr30))
df_tab # df_N has same values for all replicated of Nyr1 (thus sdNyr1 = 0), not true for Nyr30

# make dataframe of base case
df_sim <- data.frame() 
for (i in 1:length(simdatF)) {
  df_sim <-
    rbind(df_sim, data.frame(replicate = as.factor(i),
                             case = "Base",
                             year = fitHer$data$years,
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
                             year = fitHer$data$years,
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

# CP: logSd is negative so dividing by 100 increases the sd
# CP: mod2fitF$pl$logSdLogN <- fitF$pl$logSdLogN / 100 
# CP: mod2fitF$pl$logSdLogFsta <- fitF$pl$logSdLogFsta / 100
# CP: instead of dividing subtract some amount
mod2fitF$pl$logSdLogN <- fitF$pl$logSdLogN - 3 
mod2fitF$pl$logSdLogFsta <- fitF$pl$logSdLogFsta - 3

# simulate data with reduced process error
simdatF <- stockassessment:::simulate.sam(mod2fitF, nsim = 10, full.data = FALSE)

# add to dataframe 
for (i in 1:length(simdatF)) {
  df_sim <-
    rbind(df_sim, data.frame(replicate = as.factor(i),
                             case = "Mod2",
                             year = fitHer$data$years,
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

# nope still getting too much process error in Mod2 results (cp: works with new code)
