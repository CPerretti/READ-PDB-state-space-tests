# Why do sd estimates from process model not equal emprical estimates
# < Change to Kalman filter
n_t <- 1000
sigma <- 1
y <- numeric(n_t)
y[1] <- 0

for (i in 2:n_t) {
  y[i] <- y[i-1] + rnorm(n = 1, sd = sigma)
}
  
plot(y, t = "l")

sigma_derived <- sd(y[2:n_t] - y[1:(n_t - 1)])

# Compile model
library(TMB)

compile("sd_investigate.cpp")

# Setup data
Data <- list("y" = y)

# Setup parameters and starting values
Parameters <- list("log_sigma" = 0)

# Build object
dyn.load(dynlib("sd_investigate"))
Obj <- MakeADFun(data = Data, 
                 parameters = Parameters)

# Optimize
Opt <- TMBhelper::Optimize(obj = Obj, newtonsteps = 1)

# Report out the ADREPORT vars
Report_sd <- TMB::sdreport(Obj)



# Plot predicted vs observed vs true
library(dplyr)
library(tidyr)
library(ggplot2)
# Plot estimated parameters vs true parameters
df2plot <-
  data.frame(variable = "sigma",
             tru = sigma,
             der = sigma_derived,
             est = Report_sd$value["sigma"],
             sd  = Report_sd$sd[names(Report_sd$value) == "sigma"])


ggplot(df2plot, aes(y = variable)) +
  geom_point( aes(x = tru), color = "red", size  = 3) +
  geom_point( aes(x = der), color = "black", size  = 3) +
  geom_point( aes(x = est), color = "blue", size = 3) +
  geom_errorbarh(aes(xmin = est - 1.96*sd,
                     xmax = est + 1.96*sd), height = 0.1, color = "blue") +
  theme_bw() +
  ylab("") +
  xlab("Value") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
