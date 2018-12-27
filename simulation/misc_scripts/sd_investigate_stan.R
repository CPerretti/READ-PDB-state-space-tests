
# Code to illustrate how the variance of the state estimates
# changes with observation error level.
# Using a Bayesian approach.
library(dplyr)
library(ggplot2)
library(rstan)

rstan_options(auto.write = TRUE)

# Load catch data
set.seed(321)
N_obs <- 30
true_states_i <- vector(mode = "numeric", length = N_obs)
true_states_i[1] <- rnorm(n = 1) # intial condition
sd_pro <- 1
for (i in 2:N_obs) {
  true_states_i[i] <- true_states_i[i-1] + rnorm(n = 1, sd = sd_pro) # step forward with process error
}

sd_obs <- .1
observations_i <- true_states_i + rnorm(n = length(true_states_i), sd = sd_obs)

df_landings <-
  data.frame(rland = rep(NA, N_obs),
             cland = observations_i,
             Year  = 1:N_obs)
  # readxl::read_excel("../../../summer_flounder/Fluke_Catch_History.xlsx") %>%
  # dplyr::filter(Year %in% c(1970:2016)) %>%
  # dplyr::rename(rland = `Recreational Landings`,
  #               cland = `Commercial Landings`,
  #               rdisc = `Recreational Discards`,
  #               cdisc = `Commercial Discards`) %>%
  # dplyr::mutate(ccatch = `cland` + `cdisc`,
  #               rcatch = `rland` + `rdisc`,
  #               Year = as.integer(Year)) %>%
  # dplyr::select(`cland`, `rland`, Year)



# Plot landings
df2plot <- 
  df_landings %>%
  tidyr::gather(variable, value, -Year)
ggplot(data = df2plot, aes(x = Year, y = value, color = variable)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  ylab("Landings (MT)")

# Print stan file
writeLines(readLines("sd_investigate_stan.stan"))


#y <- df_landings$`Recreational Landings`
#sigma <- c(15, 10, 16, 11,  9, 11, 10, 18)

#stan_rdump(c("J", "y", "sigma"), file="eight_schools.data.R")
#data <- read_rdump("eight_schools.data.R")

rland_obs <- df_landings$rland[!is.na(df_landings$rland)]
N_obs <- length(rland_obs)
N_mis <- sum(is.na(df_landings$rland))
cland <- df_landings$cland

# Fit basic model with just mean
# fit <- stan(file = "impute_recland.stan", 
#             data = list(rland_obs = rland_obs,
#                         N_obs = N_obs,
#                         N_mis = N_mis), 
#             iter = 1e4)
# 
# print(fit)

# Fit model that uses commercial landings to predict recreational
fit <- stan(file = "sd_investigate_stan.stan", 
            data = list(rland_obs = rland_obs,
                        N_obs = N_obs,
                        N_mis = N_mis,
                        cland4obs = cland[which(!is.na(df_landings$rland))],
                        cland4mis = cland[which(is.na(df_landings$rland))]),
            iter = 1e4)

print(fit)

#shinystan::launch_shinystan(fit)

df_fit <- 
  summary(fit, pars = c("rland_mis"), probs = c(0.05, 0.95))$summary %>%
  as.data.frame()

df_landings$`rland_imputed` <- c(df_fit$mean, rep(NA, N_obs))
df_landings$`rland_imputed_low` <- c(df_fit$`5%`, rep(NA, N_obs)) 
df_landings$`rland_imputed_high` <- c(df_fit$`95%`, rep(NA, N_obs)) 

df_landings$rland_ratio <- c(df_landings$cland[1:N_mis] * (2/3), rep(NA, N_obs))

# Plot landings with imputed values and CI
df2plot <- 
  df_landings %>%
  tidyr::gather(variable, value, -Year)
ggplot(data = df2plot %>% 
         dplyr::filter(!(variable %in% c("rland_imputed_low", 
                                         "rland_imputed_high")))) +
  #geom_point(aes(x = Year, y = value, color = variable)) +
  geom_line(aes(x = Year, y = value, color = variable)) +
  geom_ribbon(data = df2plot %>% 
                dplyr::filter(variable %in% c("rland_imputed_low", 
                                              "rland_imputed_high")) %>%
                tidyr::spread(variable, value),
              aes(x = Year, ymin = rland_imputed_low, ymax = rland_imputed_high),
              alpha = 0.5,
              fill = "light blue") +
  theme_bw() +
  ylab("Landings (MT)")




