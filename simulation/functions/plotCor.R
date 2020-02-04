plotCor <- function(fit, sim_labelAccept, model) {
  corDfByRep <-
    purrr::pmap_dfr(list(fit = fit,
                         replicate = sim_labelAccept$replicate,
                         scenario = sim_labelAccept$scenario),
                    .f = function(fit, replicate, scenario) {
                      corDfByRep <-
                        fit$sdrep$cov.fixed %>% 
                        cov2cor %>%
                        as.data.frame()
                      
                      colnames(corDfByRep) <- rownames(corDfByRep)
                      
                      corDfByRep %>%
                        mutate(param1 = rownames(.)) %>%
                        gather(param2, correlation, -param1) %>%
                        cbind(replicate, scenario)
                    } )
  corDf <-
    corDfByRep %>%
    mutate(scenario = paste(scenario, "scenario"),
           scenario = factor(scenario, 
                             levels = c("no misreporting scenario", 
                                        "fixed scenario",
                                        "random walk scenario",
                                        "uniform random scenario"))) %>%
    group_by(scenario, param1, param2) %>%
    summarise(correlation = mean(correlation))
  
  p <-
    ggplot(corDf, aes(x = param1, y = param2)) +
    geom_tile(aes(fill=correlation)) +
    geom_text(aes(label = round(correlation,1)), size = 3) +
    facet_wrap(~scenario) +
    scale_fill_gradientn(limits = c(-1,1),
                         colours=c("navyblue", "darkmagenta", "darkorange1")) +
    xlab("") +
    ylab("") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    ggtitle(model)
  
  ggsave(p, file = paste0("./figures/corrMat", model, ".png"),
         height = 7, width = 8)
  
  print(p)
  
}