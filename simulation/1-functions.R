# Functions for simulations

## Plot N-at-age simulated vs fit #######################
plotN <- function(N, fit) {
  # Setup N-at-age to plot
  df2plotN <- 
    N %>%
    t() %>%
    as.data.frame() %>%
    cbind(data.frame(fit = fit$pl$logN %>% exp %>% t)) %>%
    dplyr::mutate(year = fit$data$years) %>%
    tidyr::gather(variable, N, -year) %>%
    tidyr::separate(variable, c("source", "age"))
  
  # Plot N-at-age (all ages should match exactly when using errPro_exact)
  ggplot(data = df2plotN,
         aes(x = year, y = N, color = source)) +
    geom_line() +
    facet_wrap(~age, scales = "free") +
    ylab("Abundance (1000's)")
}


## Plot C-at-age simulated vs fit #######################
plotC <- function(C, fit) {
  # Setup Catch to plot
  df_fitCatch <- 
    data.frame(variable = names(fit$sdrep$value),
               value = fit$sdrep$value) %>%
    dplyr::filter(variable == "logCatch") %>%
    dplyr::rename(logCatch = value) %>%
    dplyr::mutate(fit.total = exp(logCatch),
                  year = fit$data$years) %>%
    dplyr::select(year, fit.total)
  
  df2plotC <- 
    C %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(simulated.total = rowSums(.),
                  year = fit$data$years) %>%
    dplyr::left_join(df_fitCatch) %>%
    tidyr::gather(variable, Catch, -year) %>%
    tidyr::separate(variable, c("source", "age"))
  
  # Plot Catch (should exactly match in total subplot when using errPro_exact)
  ggplot(data = df2plotC,
           aes(x = year, y = Catch, color = source)) +
    geom_line() +
    facet_wrap(~age, scales = "free") +
    ylab("Catch (mt)")

}

## Plot survey data vs original fit #######################
plotS <- function(S, fit) {
  # Set up Survey data to plot
  df_simS <- # convert 3-d array to long data.frame
    as.data.frame(S) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(year_fleet = row.names(.)) %>%
    tidyr::separate(col = year_fleet, 
                    into = c("year", "fleetNames"), 
                    extra = "merge") %>%
    tidyr::gather(variable, value, -year, -fleetNames) %>%
    tidyr::separate(variable, c("source", "age")) %>%
    dplyr::mutate(year = as.numeric(year),
                  age = as.numeric(age),
                  observed = NA)
  
  
  # Organize SAM survey fits
  fleets <- unique(fit$data$aux[,"fleet"])
  idx <- fit$data$aux[,"fleet"]%in%fleets 
  df_fitS <-
    data.frame(year = fit$data$aux[idx,"year"], 
               observed = exp(fit$data$logobs[idx]), 
               value = exp(fit$obj$report(c(fit$sdrep$par.fixed,
                                            fit$sdrep$par.random))$predObs[idx]), 
               age = fit$data$aux[idx,"age"], 
               fleetNames = attr(fit$data,"fleetNames")[fit$data$aux[idx,"fleet"]],
               source = "fit")
  
  df2plotS <- # combine simulated Survey and fit Survey
    bind_rows(df_simS, df_fitS)
  
  # Plot Survey (should match exactly when using errPro_exact)
  ggplot(data = df2plotS %>% dplyr::filter(age>1, fleetNames != "Residual catch"),
           aes(x = year)) +
    geom_point(aes(y = observed)) +
    geom_line(aes(y = value, color = source)) +
    facet_wrap(age~fleetNames, scales = "free_y", nrow = 7)
}

## Plot simulate.sam() replicates ##########################
# SAM simulate feature
plotSimSAM <- function(fit, nsim = 1, seed = NULL) {
  set.seed(seed) # for reproduciblilty and to match with full.data = TRUE
  simOut <- stockassessment:::simulate.sam(fit, nsim = nsim, full.data = FALSE)
  
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
  ggplot(df2plotSimOut %>% dplyr::filter(variable %in% c("N")),
         aes(x = year, y = value, group = replicate)) +
    geom_line() + 
    facet_wrap(~age, scales = "free") +
    ylab("Abundance (1000's)")
}


## Prep simulation data to be fit by SAM ##################
prepSimData <- function(S, fit, C) {
  # Convert survey data to SAM input format
  surveys <-
    lapply(X = plyr::alply(S[, , surveyIndex], 3, .dims = TRUE), t) %>%
    lapply(function(x) {colnames(x) <- sub(colnames(x), # Change colnames
                                           pattern = "simulated.", 
                                           replacement = ""); x}) %>%
    # Remove ages not in survey
    lapply(function(x) x[, which(colSums(x, na.rm = TRUE) != 0)]) 
  
  for(i in 1:length(surveys)) { # set survey times to match real data
    attr(surveys[[i]], "time", fit$data$sampleTimes[surveyIndex[i]] + c(-0.25, 0.25))
  }
  
  # Export the simulated surveys to a text file, then import it
  # using read.ices().
  
  # First set survey file header
  write.table(rbind("US Atlantic Herring Survey Data", 100 + length(surveys)), 
              file = "./sim_data/surveys.dat", sep = " ", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  # Then append surveys
  for (i in 1:length(surveys)) { # loop over surveys
    header_survey <-
      rbind(names(surveys)[i], 
            paste(min(as.numeric(rownames(surveys[[i]]))), 
                  max(as.numeric(rownames(surveys[[i]])))),
            paste(1, 1, fit$data$sampleTimes[surveyIndex[i]] - 0.25, 
                  fit$data$sampleTimes[surveyIndex[i]] + 0.25),
            paste(min(colnames(surveys[[i]])), max(colnames(surveys[[i]]))))
    a_survey <- cbind(1, as.data.frame(surveys[[i]]))
    
    # Save it to file so it can be read in
    write.table(header_survey, file = "./sim_data/surveys.dat", sep = " ", 
                row.names = FALSE, col.names = FALSE, quote = FALSE,
                append = TRUE)
    
    write.table(a_survey, file = "./sim_data/surveys.dat", sep = " ", 
                row.names = FALSE, col.names = FALSE, quote = FALSE,
                append = TRUE)
  }
  
  
  # Next, manipulate catch so it also can be read in by read.ices()
  catch <- t(C) 
  colnames(catch) <- sub(colnames(catch), 
                         pattern = "simulated.", 
                         replacement = "")
  catch <- catch[, which(colSums(catch, na.rm = TRUE) != 0)]
  header_catch <-
    rbind("US Atlantic Herring Total Catch Numbers at age (000s; combines all gear types)", 
          paste(1, 2),
          paste(min(as.numeric(rownames(catch))),
                max(as.numeric(rownames(catch)))),
          paste(min(colnames(catch)), max(colnames(catch))),
          "1")
  
  # Save it to file so it can be read in
  write.table(header_catch, file = "./sim_data/catch.dat", sep = " ", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  write.table(catch, file = "./sim_data/catch.dat", sep = " ", 
              row.names = FALSE, col.names = FALSE, quote = FALSE,
              append = TRUE)
}

