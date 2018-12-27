# Why do sd estimates from process model not equal emprical estimates
set.seed(123) # for reproducibility

nT <- 1000 # length of time series
true_states_i <- vector(mode = "numeric", length = nT)
true_states_i[1] <- rnorm(n = 1) # intial condition
sd_pro <- 1
for (i in 2:nT) {
  true_states_i[i] <- true_states_i[i-1] + rnorm(n = 1, sd = sd_pro) # step forward with process error
}

sd_obs <- .01
observations_i <- true_states_i + rnorm(n = length(true_states_i), sd = sd_obs) # add observation error

# plot it 
plot(true_states_i, xlab = "Year", ylab = "Value", type = "l")
points(observations_i)

# Compile model
library(TMB)

Version <- "sd_investigate"
compile( paste0(Version,".cpp") )

# Setup data
Data <- list("observations_i" = observations_i)

# Setup parameters and starting values
Parameters <- list("estimates_i" = rep(0, times = length(observations_i)),
                   "log_sigma_pro" = 0,
                   "log_sigma_obs" = 0)

# Build object
dyn.load(dynlib(Version))
Obj <- MakeADFun(data = Data, 
                 parameters = Parameters,
                 random = "estimates_i")

# Optimize
Opt <- TMBhelper::Optimize(obj = Obj, newtonsteps = 1)

# Report out the ADREPORT vars
Report_sd <- TMB::sdreport(Obj)

# Plot estimated and observed
library(dplyr)
library(tidyr)

df2plot_estimates <-
  data.frame(Estimate = Report_sd$value,
             estimateUpperCI = Report_sd$value + 1.96 * Report_sd$sd,
             estimateLowerCI = Report_sd$value - 1.96 * Report_sd$sd,
             variable = names(Report_sd$value)) %>%
  dplyr::filter(variable == "estimates_i") %>%
  dplyr::mutate(`True state` = true_states_i,
                Observations = observations_i,
                Year = 1:length(observations_i)) %>%
  dplyr::select(-variable) %>%
  tidyr::gather(variable, Value, -Year, -estimateUpperCI, -estimateLowerCI)


library(ggplot2)
(ggplot(data = df2plot_estimates, aes(x = Year, y = Value)) +
    geom_point(data = dplyr::filter(df2plot_estimates, variable == "Observations"), 
               shape = "o", size = 3.5) +
    geom_line(data = dplyr::filter(df2plot_estimates, variable != "Observations"),
              aes(color = variable)) +
    geom_ribbon(aes(ymin = estimateLowerCI, ymax = estimateUpperCI), 
                fill = "blue", color = NA, alpha = 0.2) +
    scale_color_manual(values = c("blue", "black")) +
    theme_bw() +
    theme(legend.title = element_blank()))

# Plot error sd estimates
library(magrittr)
estimates_i <- 
  df2plot_estimates %>% 
  dplyr::filter(variable == "Estimate") %$%
  Value

# Calculated derived sds (from fit)
derived_sigma_pro <- sd(estimates_i[2:nT] - estimates_i[1:(nT-1)])
derived_sigma_obs <- sd(estimates_i - observations_i)

df2plot <-
  data.frame(variable = "Process error sd",
             tru = sd_pro,
             sim = sd(true_states_i[2:nT] - true_states_i[1:(nT-1)]), #simulated sd
             der = derived_sigma_pro,
             est = Report_sd$value["sigma_pro"],
             sd  = Report_sd$sd[1]) %>%
  rbind(data.frame(variable = "Observation error sd",
                   tru = sd_obs,
                   sim = sd(true_states_i - observations_i), #simulated
                   der = derived_sigma_obs,
                   est = Report_sd$value["sigma_obs"],
                   sd  = Report_sd$sd[2]))


ggplot(df2plot, aes(y = variable)) +
  geom_point( aes(x = tru), color = "red", size  = 3) +
  geom_point( aes(x = est), color = "blue", size = 3) +
  geom_point( aes(x = der), color = "black", size = 3) +
  #geom_point( aes(x = sim), color = "purple", size = 3) +
  geom_errorbarh(aes(xmin = est - 1.96*sd,
                     xmax = est + 1.96*sd), height = 0.1, color = "blue") +
  theme_bw() +
  ylab("") +
  xlab("Value") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))# +
  #xlim(c(0,2))

# Do you get a better likelihood if you use the derived sd?
nll_derived <- Obj$fn(c(log(derived_sigma_pro), log(derived_sigma_obs)))
nll_mle     <- Obj$fn(Opt$par) # using MLE's

(nll_mle < nll_derived)
# MLE nll is lower than the derived nll, so it makes sense.

