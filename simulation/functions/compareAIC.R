compareAIC <- function(sim_labelAccept, 
                       fitSim_noneAccept, 
                       fitSim_fixedAccept, 
                       fitSim_randomAccept,
                       err) {
  
  aic_noneAccept   <- data.frame(model = "no misreporting",
                                  sim_labelAccept,
                                  aic = sapply(fitSim_noneAccept, AIC))
  
  aic_fixedAccept  <- data.frame(model = "fixed",
                                  sim_labelAccept,
                                  aic = sapply(fitSim_fixedAccept, AIC))
  
  aic_randomAccept <- data.frame(model = "random walk",
                                  sim_labelAccept,
                                  aic =sapply(fitSim_randomAccept, AIC))
  
  aicAccept <-
    rbind(aic_noneAccept,
          aic_fixedAccept,
          aic_randomAccept) %>%
    mutate(scenario  = as.factor(paste(scenario, "scenario")),
           scenario  = factor(scenario, levels = c("no misreporting scenario",
                                                   "fixed scenario",
                                                   "random walk scenario",
                                                   "uniform random scenario")))
  

  
  # Plot boxplot of AICs
  colors2use <- RColorBrewer::brewer.pal(3, "Dark2")
  
  p <-
    ggplot(aicAccept, aes(x = model, y = aic, color = model)) +
    geom_boxplot() +
    facet_wrap(~scenario, scales = "free") +
    scale_color_manual(values = colors2use[1:3]) +
    theme_bw() +
    theme(legend.position = "none")
  
  print(p)
  
  # Plot histogram of AICs
    ggplot(aicAccept, aes(x = aic, color = model, fill = model)) +
    geom_histogram(alpha = 0.5) +
    facet_wrap(~scenario) +
    scale_color_manual(values = colors2use[1:3]) +
    scale_fill_manual(values = colors2use[1:3]) +
    theme_bw()
  
  
  # Plot histogram of differences between RW AIC and fixed AIC
  p <-
    ggplot(aicAccept %>% 
           spread(model, aic) %>% 
           mutate(`random walk aic - fixed aic` = `random walk` - fixed) %>%
           select(replicate, scenario, `random walk aic - fixed aic`),
         aes(x = `random walk aic - fixed aic`)) +
    geom_histogram() +
    facet_wrap(~scenario) +
    #scale_color_manual(values = colors2use[1:3]) +
    #scale_fill_manual(values = colors2use[1:3]) +
    theme_bw()
  
  print(p)
  
  
  # Compare error in each variable vs AIC
  meanErrWithAIC <- 
    err %>%
    group_by(scenario, replicate, 
             variable, model) %>%
    summarise(mape  = mean(abs_error_pc, na.rm = T)) %>%
    spread(model, mape) %>%
    mutate(`random walk mape - fixed mape` = `random walk` - fixed) %>%
    left_join(aicAccept %>% 
                spread(model, aic) %>% 
                mutate(`random walk aic - fixed aic` = `random walk` - fixed) %>%
                select(replicate, scenario, `random walk aic - fixed aic`))
  
  ggplot(meanErrWithAIC %>% 
           filter(variable != "catch_observed",
                  `random walk mape - fixed mape` < 20),
         aes(x = `random walk mape - fixed mape`,
             y = `random walk aic - fixed aic`)) +
    geom_point() +
    facet_grid(scenario~variable, scales = "free") +
    theme_bw()
  
  # Print proportion of replicates that each model is chosen best.
    propAicAccept <- 
      aicAccept %>%
      group_by(replicate, scenario) %>%
      mutate(aicMin = ifelse(aic == min(aic), 1, 0)) %>%
      group_by(scenario, model) %>%
      summarise(propMinAic = round(sum(aicMin) / n(), 2)) %>%
      spread(model, propMinAic)
  
  
}
