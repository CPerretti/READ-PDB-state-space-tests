# Functions for simulations

## Simulation model #######################################
sim <- function(fit) {
  nA <- ncol(fit$data$propF) # number of age-classes
  nT <- fit$data$noYears # length of time series
  
  # Set F (need to replicate some elements to match ModelConf)
  #f <- exp(fit$pl$logF[(fit$conf$keyLogFsta[1,] + 1),])
  # Calculate a new relization of f errors with sd's from the fit
  errF <- matrix(data = NA,
                 nrow = nA, 
                 ncol = nT - 1)
  
  # Set F sd  (notice the reduction of sd the natural scale)
  fit$pl$logSdLogFsta <- 
    (c(0.1, rep(0.33, length(fit$pl$logSdLogFsta)-1)) * exp(fit$pl$logSdLogFsta)) %>% 
    log
  sdLogF <- exp(fit$pl$logSdLogFsta)
  for (i in 1:(nT-1)) { # Create F error
    errF[, i] <- rnorm(n = nA,
                       sd = sdLogF[fit$conf$keyVarF[1, ] + 1])[fit$conf$keyLogFsta[1, ] + 1]
  }
  
  # Simulate F
  logF <- matrix(data = NA, # F container
                 nrow = nA, 
                 ncol = nT, 
                 dimnames = list(paste0("tru.", c(1:nA)), fit$data$years)) 
  
  logF[, 1] <- fit$pl$logF[, 1][fit$conf$keyLogFsta[1, ] + 1] # initial F from fitted model
  
  for (i in 2:nT) logF[, i] <- logF[, i - 1] + errF[, i - 1]
  
  f <- exp(logF)
  
  # Set M
  m <- t(fit$data$natMor)
  
  # Calcuate total mortality
  z <- f + m
  
  # Set up matrix to record N-at-age
  logN <- matrix(data = NA,
                 nrow = nA, 
                 ncol = nT, 
                 dimnames = list(paste0("tru.", c(1:nA)), fit$data$years)) 
  
  logN[, 1] <- fit$pl$logN[,1] # initial condition from fitted model
  
  # Calculate the process errors that were estimated in the fit so we can 
  # exactly replicate the fit
  # errPro_exact <- matrix(data = NA,
  #                        nrow = nA, 
  #                        ncol = nT-1)
  # 
  # errPro_exact[1, ] <- fit$pl$logN[1, 2:nT] - fit$pl$logN[1, 1:(nT-1)]
  # errPro_exact[-c(1, nA), ] <-  fit$pl$logN[-c(1, nA), 2:nT] -
  #   (fit$pl$logN[-c(nA-1, nA), 1:(nT-1)] -
  #      z[-c(nA-1, nA), 1:(nT-1)])
  # errPro_exact[nA, ] <- fit$pl$logN[nA, 2:nT] -
  #   log(exp(fit$pl$logN[nA-1, 1:(nT-1)]) *
  #         exp(-z[nA-1, 1:(nT-1)]) +
  #         exp(fit$pl$logN[nA, 1:(nT-1)]) *
  #         exp(-z[nA, 1:(nT-1)]))
  # Calculate a new relization of process errors with sd's from the fit
  # Set N process sd
  errPro <- matrix(data = NA,
                   nrow = nA, 
                   ncol = nT-1)
  # Lower process error (!)
  fit$pl$logSdLogN[(fit$conf$keyVarLogN + 1)] <-
    fit$pl$logSdLogN[(fit$conf$keyVarLogN + 1)] %>%
    exp %>% "*"(0.33) %>% log # NOTICE THE 0.033 factor reduction on the natural scale
  sdLogN <- exp(fit$pl$logSdLogN[(fit$conf$keyVarLogN + 1)])
  for (i in 1:(nT-1)) { # Create process error (N-at-age)
    errPro[, i] <-  rnorm(n = nA, sd = sdLogN)
  }
  
  #errPro <- errPro_exact # Use if you want the herring 13 fit N
  
  ## Simulate process model #################################
  # N fit
  # Simulate N-at-age
  for (i in 2:nT) {
    logN[1, i] <- logN[1, i-1] + errPro[1, i-1]
    logN[-c(1, nA), i] <- logN[-c(nA-1, nA), i-1] - 
      z[-c(nA-1, nA), i-1] + 
      errPro[-c(1, nA), i-1]
    logN[nA, i] <- log(exp(logN[nA-1, i-1]) * exp(-z[nA-1, i-1]) +
                         exp(logN[nA, i-1]) * 
                         exp(-z[nA, i-1])) + errPro[nA, i-1]
  }
  
  N <- exp(logN)
  
  ## Simulate observation model #############################
  # Observation error (for catch and surveys)
  errObs <- array(data = NA, # error container (3-d: age x year x fleet)
                  dim = c(nA, nT, fit$data$noFleets),
                  dimnames = list(paste0("error.", c(1:nA)), 
                                  fit$data$years, 
                                  attr(fit$data,"fleetNames")))
  # Need to replicate some sd's to match config file
  index <- as.vector(t(fit$conf$keyVarObs + 1))
  index[index == 0] <- NA
  sdLogObs <- exp(fit$pl$logSdLogObs[index])
  
  # Make observation error (can only do uncorrelated error right now)
  for (j in 1:nT) { # all surveys in a year
    errObs[, j, ] <- rnorm(n = length(sdLogObs), sd = sdLogObs)
  }
  
  # Check that sds are correctly assigned to each age x survey combo
  # sdInput <- c()
  # for (i in 1:fit$data$noFleets) {
  #   sdInput <- c(sdInput, apply(errObs[, , i], 1, sd))  
  # }
  # 
  # cbind(sdInput, sdLogObs)
  
  
  # Simulate Catch
  
  # Calculate catch on N-scale (1000s)
  logCtru_N <- log(f / z * (1 - exp(-z)) * N)
  rownames(logCtru_N) <- 1:nA
  logCobs_N <- logCtru_N + errObs[, , "Residual catch"]
  
  Ctru_N <- exp(logCtru_N)
  Cobs_N <- exp(logCobs_N)
  
  # Convert to MT (1000s * kg = mt)
  Ctru_mt <- Ctru_N * t(fit$data$catchMeanWeight)
  Cobs_mt <- Cobs_N * t(fit$data$catchMeanWeight)
  
  
  # Simulate Survey (SIMULATIONS HAVE DATA IN ALL YEARS RIGHT NOW)
  logStru_N <- array(data = NA, # survey container (3-d: age x year x survey)
                     dim = c(nA, nT, fit$data$noFleets),
                     dimnames = list(c(1:nA), 
                                     fit$data$years, 
                                     attr(fit$data,"fleetNames")))
  
  logSobs_N <- logStru_N
  
  logSq <- matrix(data = NA, # survey q-at-age matrix
                  nrow = nrow(fit$conf$keyLogFpar), 
                  ncol = ncol(fit$conf$keyLogFpar))
  logSq[which(fit$conf$keyLogFpar != -1)] <- # fill with herring fit values
    fit$pl$logFpar[fit$conf$keyLogFpar[fit$conf$keyLogFpar != -1] + 1]
  Sq <- exp(logSq)
  
  surveyIndex <- # some fleets are fishermen not surveys
    (1:fit$data$noFleets)[fit$data$fleetTypes == 2] 
  for (i in surveyIndex) {
    logStru_N[, , i] <- log(Sq[i,] * exp(-z * fit$data$sampleTimes[i]) * N) 
    logSobs_N[, , i] <- logStru_N[, , i] + errObs[, , i]
  }
  
  Stru_N <- exp(logStru_N)
  Sobs_N <- exp(logSobs_N)
  
  
  trueParams <- list(sdrep = fit$sdrep, pl = fit$pl)
  trueParams$pl$logN <- log(N)
  dimnames(f) <- list(paste0("tru.", c(1:nA)), fit$data$years)
  trueParams$pl$logF <- log(f)
  return(list(trueParams = trueParams,
              N = N,
              Cobs_mt = Cobs_mt, Cobs_N = Cobs_N, 
              Ctru_mt = Ctru_mt, Ctru_N = Ctru_N, 
              Sobs_N = Sobs_N, Stru_N = Stru_N))
  
}

# Setup data and params for sam model #####################
setupModel <- function() {
  
  # Read in data
  cn <- read.ices("./sim_data/catch.dat") # catch abundace-at-age
  cw <- read.ices("../atlherring_example/data/Herrcw.dat") # catch mean weight-at-age
  dw <- cw # discards mean weight-at-age (using catch mean weight-at-age for now)
  lw <- cw # landings mean weight-at-age (using catch mean weight-at-age for now)
  pf <- read.ices("../atlherring_example/data/Herrpf.dat") # proportion of f before spawning
  lf <- pf; lf[,] <- 1 # fraction of catch that is landed (set to 1 for now)
  mo <- read.ices("../atlherring_example/data/Herrmo.dat") # maturity-at-age ogive
  nm <- read.ices("../atlherring_example/data/Herrnm.dat") # natural mortality-at-age
  pm <- read.ices("../atlherring_example/data/Herrpm.dat") # proportion of m before spawning
  sw <- read.ices("../atlherring_example/data/Herrsw.dat") # stock weight-at-age (kg)
  surveys <- read.ices("./sim_data/surveys.dat") #surveys
  #surveys <- read.ices("../atlherring_example/data/Herrsurvey_BigSep_NoAcoust.dat") #surveys
  
  # setup the data as needed for SAM
  dat <- setup.sam.data(surveys = surveys,
                        residual.fleet = cn,
                        prop.mature = mo,
                        stock.mean.weight = sw,
                        dis.mean.weight = dw,
                        land.mean.weight = lw,
                        land.frac = lf,
                        prop.f = pf,
                        prop.m = pm,
                        natural.mortality = nm,
                        catch.mean.weight = cw)
  
  # Load model configuration file
  conf <- loadConf(dat = dat, file = "../atlherring_example/ModelConf_original.txt")
  
  par <- defpar(dat, conf) # some default starting values
  
  return(list(dat = dat, conf = conf, par = par))
}
## Plot F-at-age simulated vs fit #######################
plotF <- function(simOut, fit) {
  # Setup F-at-age to plot
  df2plotF <-
    simOut$trueParams$pl$logF %>%
    exp %>%
    t() %>%
    as.data.frame() %>%
    cbind(data.frame(fit = fit$pl$logF[fit$conf$keyLogFsta[1, ] + 1, ] %>% exp %>% t)) %>%
    dplyr::mutate(year = fit$data$years) %>%
    tidyr::gather(variable, f, -year) %>%
    tidyr::separate(variable, c("source", "age")) %>%
    dplyr::mutate(age = paste0("age-", age))
  
  # Plot N-at-age (all ages should match exactly when using errPro_exact)
  ggplot(data = df2plotF,
         aes(x = year, y = f, color = source)) +
    geom_line() +
    facet_wrap(~age, scales = "free") +
    ylab("F") +
    ggtitle("F-at-age") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 12))
}

## Plot N-at-age simulated vs fit #######################
plotN <- function(simOut, fit) {
  N <- simOut$N
  
  # Setup N-at-age to plot
  df2plotN <- 
    N %>%
    t() %>%
    as.data.frame() %>%
    cbind(data.frame(fit = fit$pl$logN %>% exp %>% t)) %>%
    dplyr::mutate(year = fit$data$years) %>%
    tidyr::gather(variable, N, -year) %>%
    tidyr::separate(variable, c("source", "age")) %>%
    dplyr::mutate(age = paste0("age-", age))
  
  # Plot N-at-age (all ages should match exactly when using errPro_exact)
  ggplot(data = df2plotN,
         aes(x = year, y = N, color = source)) +
    geom_line() +
    facet_wrap(~age, scales = "free") +
    ylab("Abundance (1000's)") +
    ggtitle("Population abundance-at-age") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 12))
}

## Plot C-at-age simulated vs fit #######################
plotC <- function(simOut, fit) {
  Cobs_mt <- simOut$Cobs_mt
  Ctru_mt <- simOut$Ctru_mt
  
  # Setup Catch to plot
  df_Cfit_mt <- 
    data.frame(variable = names(fit$sdrep$value),
               value = fit$sdrep$value) %>%
    dplyr::filter(variable == "logCatch") %>%
    dplyr::rename(logCatch = value) %>%
    dplyr::mutate(Catch_mt = exp(logCatch),
                  year = fit$data$years,
                  age = "total",
                  source = "fit") %>%
    dplyr::select(year, age, Catch_mt, source)
  
  # Pull out simulated observations
  df_Cobs_mt <- 
    Cobs_mt %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(total = rowSums(.),
                  year = fit$data$years) %>%
    tidyr::gather(age, Catch_mt, -year) %>%
    dplyr::mutate(source = "observed")
  
  # Pull out simulated truth
  df_Ctru_mt <- 
    Ctru_mt %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(total = rowSums(.),
                  year = fit$data$years) %>%
    tidyr::gather(age, Catch_mt, -year) %>%
    dplyr::mutate(source = "tru")  
  
  # Pull out simulated true
  df2plot <- 
    bind_rows(df_Cfit_mt, 
              df_Cobs_mt,
              df_Ctru_mt) %>%
    dplyr::mutate(age = paste0("age-", age))
  
  # Plot Catch (should exactly match in total subplot when using errPro_exact)
  ggplot(data = df2plot,
           aes(x = year, y = Catch_mt, color = source)) +
    geom_line() +
    facet_wrap(~age, scales = "free") +
    ylab("Catch (MT)") +
    ggtitle("Fishery catch-at-age") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 12))

}

## Plot survey data vs original fit #######################
plotS <- function(simOut, fit) {
  Sobs_N <- simOut$Sobs_N
  Stru_N <- simOut$Stru_N
  
  # Set up Survey data to plot
  df_Sobs <- # convert 3-d array to long data.frame
    as.data.frame(Sobs_N) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(year_fleet = row.names(.)) %>%
    tidyr::separate(col = year_fleet, 
                    into = c("year", "fleetNames"), 
                    extra = "merge") %>%
    tidyr::gather(age, value, -year, -fleetNames) %>%
    dplyr::mutate(year = as.numeric(year),
                  age = as.numeric(age),
                  source = "observed")
  
  df_Stru <- # convert 3-d array to long data.frame
    as.data.frame(Stru_N) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(year_fleet = row.names(.)) %>%
    tidyr::separate(col = year_fleet, 
                    into = c("year", "fleetNames"), 
                    extra = "merge") %>%
    tidyr::gather(age, value, -year, -fleetNames) %>%
    dplyr::mutate(year = as.numeric(year),
                  age = as.numeric(age),
                  source = "tru")
  
  
  # Organize SAM survey fits
  fleets <- unique(fit$data$aux[,"fleet"])
  idx <- fit$data$aux[,"fleet"]%in%fleets 
  df_Sfit <-
    data.frame(year = fit$data$aux[idx,"year"], 
               #observed = exp(fit$data$logobs[idx]), 
               value = exp(fit$obj$report(c(fit$sdrep$par.fixed,
                                            fit$sdrep$par.random))$predObs[idx]), 
               age = fit$data$aux[idx,"age"], 
               fleetNames = attr(fit$data,"fleetNames")[fit$data$aux[idx,"fleet"]],
               source = "fit")
  
  df2plot <- # combine simulated Survey and fit Survey
    dplyr::bind_rows(df_Sobs, 
                     df_Stru, 
                     df_Sfit) %>%
    dplyr::mutate(age = paste0("age-", age))
  
  # Plot Survey (should match exactly when using errPro_exact)
  ggplot(data = df2plot %>% 
                dplyr::filter(age > 1, # for atl herring only
                              fleetNames != "Residual catch"),
           aes(x = year, y = value, color = source)) +
    geom_line() +
    facet_wrap(age ~ fleetNames, scales = "free", ncol = 5) +
    ggtitle("Survey catch-at-age") +
    ylab("Survey catch (1000s)") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 12))
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
prepSimData <- function(Sobs_N, fit, Cobs_N) {
  # Convert survey data to SAM input format
  surveyIndex <- # some fleets are fishermen not surveys
    (1:fit$data$noFleets)[fit$data$fleetTypes == 2] 
  surveys <-
    lapply(X = plyr::alply(Sobs_N[, , surveyIndex], 3, .dims = TRUE), t) %>%
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
  catch <- t(Cobs_N) 
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


## Plot parameters fit vs true ############################
plotPars <- function(fitSim, simOut) {
  parsFixed <- which(names(fitSim[[1]]$pl) %in% names(fitSim[[1]]$obj$par))
  df_parsOut <- data.frame()
  for (h in parsFixed) {
    for (i in 1:nRep) {
      df_parsOut <-
        rbind(df_parsOut,
              data.frame(variable = paste(names(fitSim[[1]]$pl)[h], 
                                          1:length(fitSim[[1]]$pl[[h]]), sep = "."),
                         tru = simOut[[i]]$trueParams$pl[[h]],
                         est = fitSim[[i]]$pl[[h]],
                         sd  = fitSim[[i]]$plsd[[h]],
                         replicate = i))
    }  
  }
  
  
  df2plot <-
    df_parsOut %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(tru = unique(tru),
                     est_mean = mean(est),
                     est_se = sd(est)/sqrt(nRep))
  
  ggplot(df2plot, aes(y = variable)) +
    geom_point( aes(x = tru), color = "red", size  = 3) +
    geom_point( aes(x = est_mean), color = "blue", size = 3) +
    geom_errorbarh(aes(xmin = est_mean - 1.96*est_se,
                       xmax = est_mean + 1.96*est_se), color = "blue") +
    theme_bw() +
    ylab("") +
    xlab("Value") +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14))
}

## Calculate timeseries error #############################
calcTsError <- function(fitSim, simOut) {
  indNF <- which(names(fitSim[[1]]$pl) %in% c("logN", "logF"))
  err_logNF <- data.frame()
  for (h in indNF) {
    for (i in 1:nRep) {
      rownames(fitSim[[i]]$pl[[h]]) <- paste0("fit.", 1:nrow(fitSim[[i]]$pl[[h]]))
      Sd <- exp(fitSim[[i]]$plsd[[h]])
      rownames(Sd) <- paste0("Sd.", 1:nrow(fitSim[[i]]$pl[[h]]))
      err_logNF <-
        rbind(err_logNF,
              fitSim[[i]]$pl[[h]] %>%
                t() %>%
                exp %>%
                as.data.frame() %>%
                cbind(Sd %>%
                        t() %>%
                        as.data.frame()) %>%
                cbind(data.frame(simOut[[i]]$trueParams$pl[[h]] %>% t %>% exp)) %>%
                dplyr::mutate(year = as.numeric(fitSim[[i]]$data$years)) %>%
                tidyr::gather(variable, N, -year) %>%
                dplyr::mutate(variable = gsub(x = variable, 
                                              pattern = "X", 
                                              replacement = "tru.")) %>%
                tidyr::separate(variable, c("source", "age")) %>%
                tidyr::spread(source, N) %>%
                dplyr::mutate(age = as.numeric(age),
                              error = (fit - tru),
                              decile = ceiling(10 * pnorm(q    = log(tru), 
                                                          mean = log(fit), 
                                                          sd   = log(Sd))),
                              replicate = i,
                              variable = ifelse(names(fitSim[[i]]$pl[h]) == "logN", "N",
                                                if(names(fitSim[[i]]$pl[h]) == "logF") "F")))
    }
  }
  return(err_logNF)
}

## Plot timeseries fit vs true ############################
plotTsError <- function(err_logNF) {

  
  # Calculate median error
  err_logNFannual <-
    err_logNF %>%
    dplyr::group_by(variable, age, year) %>%
    dplyr::summarise(median_error = median(error))
  
  # Plot median error for N and F in each year for each age
  p <-
    ggplot(err_logNFannual %>% 
             dplyr::filter(variable == "N") %>%
             dplyr::mutate(as.numeric(year)),
           aes(x = year, y = median_error)) +
      geom_line() +
      geom_hline(yintercept = 0) +
      theme_bw() +
      ylab("Median error (fit - true)") +
      xlab("Year") +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14)) +
      facet_wrap(~age) +
      ggtitle("N estimation error-at-age")
  print(p)
  
  p <-
    ggplot(err_logNFannual %>% 
             dplyr::filter(variable == "F") %>%
             dplyr::mutate(as.numeric(year)),
           aes(x = year, y = median_error)) +
      geom_line() +
      geom_hline(yintercept = 0) +
      theme_bw() +
      ylab("Median error (fit - true)") +
      xlab("Year") +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14)) +
      facet_wrap(~age) +
      ggtitle("F estimation error-at-age")
  print(p)
  
  # Plot a few example fit vs tru time series
  p <-
    ggplot(err_logNF %>%
             dplyr::select(-Sd, -decile) %>%
             dplyr::filter(variable == "F",
                           age == 4,
                           replicate %in% 1:10) %>%
             tidyr::gather(source, value, 
                           -year, -age, -replicate, -variable),
           aes(x = year, y = value, color = source)) +
      geom_line() +
      facet_wrap(~replicate) + 
      ggtitle("Age-4 F estimation error examples")
  print(p)
  
  # Plot decile coverage for each age and variable
  p <- 
    ggplot(err_logNF,
           aes(x = decile)) +
      geom_histogram(binwidth=1, colour="white") +
      geom_hline(yintercept = ncol(fitSim[[1]]$pl$logN) * nRep / 10, 
                 color = "dark grey") +
      theme_bw() +
      facet_grid(variable~age) +
      ggtitle("Confidence interval coverage for each age")
  print(p)
}
