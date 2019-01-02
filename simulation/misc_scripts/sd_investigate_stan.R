
# Code to illustrate how the variance of the state estimates
# changes with observation error level.
# Using a Bayesian approach.
library(dplyr)
library(ggplot2)
library(rstan)

rstan_options(auto.write = TRUE)

# Load catch data
set.seed(321)
N_obs <- 100
sigma_pro <- 1
sigma_obs <- 5
true_states_i <- vector(mode = "numeric", length = N_obs)
true_states_i[1] <- rnorm(n = 1, sd = sigma_pro) # intial condition
for (i in 2:N_obs) {
  true_states_i[i] <- true_states_i[i-1] + rnorm(n = 1, sd = sigma_pro)
}

observations_i <- true_states_i + rnorm(n = length(true_states_i), sd = sigma_obs)

df_all <-
  data.frame(tru = true_states_i,
             obs = observations_i,
             Year = 1:N_obs)


# Fit random walk state-space model
fit <- stan(file = "sd_investigate_stan.stan", 
            data = list(N_obs = N_obs,
                        obs = df_all$obs),
            iter = 1e4)

df_fit <- 
  summary(fit, pars = c("tru_est"), probs = c(0.05, 0.95))$summary %>%
  as.data.frame()

df_all$tru_est <- df_fit$mean
df_all$tru_est_lo <- df_fit$`5%`
df_all$tru_est_hi <- df_fit$`95%`

# Plot tru, obs, and est
df2plot <- tidyr::gather(df_all, variable, value, -Year)
ggplot(data = df2plot, aes(x = Year)) +
  geom_ribbon(data = dplyr::filter(df2plot, variable %in% c("tru_est_lo",
                                                            "tru_est_hi")) %>%
                tidyr::spread(variable, value),
              aes(ymin = tru_est_lo, ymax = tru_est_hi),
              fill = "light blue") +
  # geom_point(data = dplyr::filter(df2plot, variable == "obs"),
  #            aes(y = value, color = variable)) +
  geom_line(data = dplyr::filter(df2plot, variable %in% c("tru", "tru_est")),
            aes(y = value, color = variable)) +
  theme_bw()


# Plot histogram of estimated for each year
df_samples <-
  extract(fit, pars = "tru_est", permuted = TRUE) %>%
  as.data.frame %>%
  tidyr::gather(variable, value)

print(fit)

derived_sigma_pro <- sd(df_all$tru_est[2:N_obs] - df_all$tru_est[1:(N_obs-1)])

#shinystan::launch_shinystan(fit)


