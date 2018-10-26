# Let's try to figure out whether these runs are realistic. If we establish
# a threshold for when a run is unrealistic then we can plot the
# proportion of replicates where a certain threshold is met:
library(dplyr)
library(tidyr)
library(ggplot2)
library(stockassessment)
# Use atl herring fit to set up simulation
load("../atlherring_example/output/fit.Rdata")

nA <- ncol(fit$data$propF) # number of age-classes

## Compare my simulation to simulate.sam() ####
# SAM simulate feature
set.seed(123) # for reproduciblilty and to match with full.data = TRUE (below)
simOut <- stockassessment:::simulate.sam(fit, nsim = 1000, full.data = FALSE)

df_simOut <- data.frame() 
for (i in 1:length(simOut)) {
  df_simOut <-
    rbind(df_simOut, data.frame(replicate = as.factor(i),
                                year = fit$data$years,
                                N = t(exp(simOut[[i]]$logN)),
                                F = t(exp(simOut[[i]]$logF))))
}
df2plotSimOut <-
  df_simOut %>%
  cbind(data.frame(fit = fit$pl$logN %>% exp %>% t)) %>% # original fit
  tidyr::gather(variable, value, -year, -replicate) %>%
  tidyr::separate(variable, into = c("variable", "age"))

# plot simulate.sam replicates
p <-
  ggplot(df2plotSimOut %>% dplyr::filter(variable %in% c("N")),
         aes(x = year, y = value, group = replicate)) +
  geom_line() + 
  facet_wrap(~age, scales = "free") +
  ylab("Abundance (1000's)")


# Calculate mean/median of each simulation and compare to the original
# fit.
fitN <- exp(fit$pl$logN)
colnames(fitN) <- fit$data$years
rownames(fitN) <- paste(1:nA)
df_meanNFit <-
  fitN %>%
  t() %>%
  as.data.frame() %>%
  dplyr::mutate(year = rownames(.),
                variable = "N") %>%
  tidyr::gather(age, value, -year, -variable) %>%
  dplyr::group_by(age, variable) %>%
  dplyr::summarise(meanNFit = mean(value))

df2plot <-
  df2plotSimOut %>%
  dplyr::filter(variable == "N") %>%
  dplyr::group_by(age, variable, replicate) %>%
  dplyr::summarise(meanValue = mean(value)) %>%
  dplyr::left_join(df_meanNFit) %>%
  dplyr::group_by(age, variable) %>%
  dplyr::summarise(`>1x` = mean(meanValue > 1*meanNFit), # proportion of runs above threshold
                   `>10x` = mean(meanValue > 10*meanNFit),
                   `>100x` = mean(meanValue > 100*meanNFit),
                   `>1000x` = mean(meanValue > 1000*meanNFit)) %>%
  tidyr::gather(threshold, value, -age, -variable) %>%
  dplyr::mutate(threshold = as.factor(threshold), # reverse levels for plot
                threshold = factor(threshold, levels = rev(levels(threshold)))) %>%
  dplyr::mutate(metric = "Mean")

ggplot(df2plot, aes(x = threshold, y = value)) +
  geom_point() +
  facet_wrap(~age) +
  xlab("Increase in mean of simulation compared to mean of fit") +
  ylab("Proportion of simulations") +
  ggtitle("How biased are the simulations compared to the fit?") +
  ylim(0,1)

# repeating exercise using median instead of mean
df_medianNFit <-
  fitN %>%
  t() %>%
  as.data.frame() %>%
  dplyr::mutate(year = rownames(.),
                variable = "N") %>%
  tidyr::gather(age, value, -year, -variable) %>%
  dplyr::group_by(age, variable) %>%
  dplyr::summarise(medianNFit = median(value))

df2plotmed <-
  df2plotSimOut %>%
  dplyr::filter(variable == "N") %>%
  dplyr::group_by(age, variable, replicate) %>%
  dplyr::summarise(medianValue = median(value)) %>%
  dplyr::left_join(df_medianNFit) %>%
  dplyr::group_by(age, variable) %>%
  dplyr::summarise(`>1x` = mean(medianValue > 1*medianNFit), # proportion of runs above threshold
                   `>10x` = mean(medianValue > 10*medianNFit),
                   `>100x` = mean(medianValue > 100*medianNFit),
                   `>1000x` = mean(medianValue > 1000*medianNFit)) %>%
  tidyr::gather(threshold, value, -age, -variable) %>%
  dplyr::mutate(threshold = as.factor(threshold), # reverse levels for plot
                threshold = factor(threshold, levels = rev(levels(threshold)))) %>%
  dplyr::mutate(metric = "Median")

df2plotcomb <- rbind(df2plot, df2plotmed)

ggplot(df2plotcomb, aes(x = threshold, y = value, color = metric)) +
  geom_point() +
  facet_wrap(~age) +
  xlab("Increase in metric of simulation compared to metric of fit") +
  ylab("Proportion of simulations") +
  ggtitle("How biased are the simulations compared to the fit?") +
  ylim(0,1)