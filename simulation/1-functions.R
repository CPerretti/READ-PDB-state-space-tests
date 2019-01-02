# Functions for simulations

## Simulation model #######################################
sim <- function(fit) {
  nA <- ncol(fit$data$propF) # number of age-classes
  years2use <- fit$data$years #<<< DO SOMETHING with fit$data$aux to filter years
  nT <- fit$data$noYears # length of time series
  
  # Set F (need to replicate some elements to match ModelConf)
  #f <- exp(fit$pl$logF[(fit$conf$keyLogFsta[1,] + 1),])
  # Calculate a new relization of f errors with sd's from the fit
  errF <- matrix(data = NA,
                 nrow = nA, 
                 ncol = nT - 1)
  
  # Set F sd  (notice the reduction of sd on the natural scale)
  fit$pl$logSdLogFsta <- 
    (c(0.1, rep(0.33, length(fit$pl$logSdLogFsta)-1)) * exp(fit$pl$logSdLogFsta)) %>% 
    #(c(1, rep(1, length(fit$pl$logSdLogFsta)-1)) * exp(fit$pl$logSdLogFsta)) %>% 
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
    exp %>% "*"(0.5) %>% log # NOTICE the reduction on the natural scale
    #exp %>% "*"(1) %>% log
  sdLogN <- exp(fit$pl$logSdLogN[(fit$conf$keyVarLogN + 1)])
  for (i in 1:(nT-1)) { # Create process error (N-at-age)
    errPro[, i] <-  rnorm(n = nA, sd = sdLogN)
  }
  
  #errPro <- errPro_exact # Use if you want the herring 13 fit N
  
  ## Simulate population model #################################
  # N fit
  # Simulate N-at-age
  for (i in 2:nT) {
    logN[1, i] <- logN[1, i-1] + errPro[1, i-1]
    logN[-c(1, nA), i] <- logN[-c(nA-1, nA), i-1] - 
                           z[-c(nA-1, nA), i-1] + 
                           errPro[-c(1, nA), i-1]
    logN[nA, i] <- log(exp(logN[nA-1, i-1] - z[nA-1, i-1]) +
                         exp(logN[nA, i-1] - z[nA, i-1])) + errPro[nA, i-1]
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
  fit$pl$logSdLogObs <- fit$pl$logSdLogObs %>% exp %>% "*"(1) %>% log # REDUCTION!!!!!!
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
  #logCtru_N <- log(f / z * (1 - exp(-z)) * N)
  logCtru_N <- logN - log(z) + log(1 - exp(-z)) + logF
  rownames(logCtru_N) <- 1:nA
  logCobs_N <- logCtru_N + errObs[, , "Residual catch"]
  
  Ctru_N <- exp(logCtru_N)
  Cobs_N <- exp(logCobs_N)
  
  # Convert to MT (1000s * kg = mt)
  Ctru_mt <- Ctru_N * t(fit$data$catchMeanWeight)
  Cobs_mt <- Cobs_N * t(fit$data$catchMeanWeight)
  
  
  # Simulate Survey
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
    #logStru_N[, , i] <- log(Sq[i,] * exp(-z * fit$data$sampleTimes[i]) * N)
    logStru_N[, , i] <- logN - z * fit$data$sampleTimes[i] + logSq[i,]
    logSobs_N[, , i] <- logStru_N[, , i] + errObs[, , i]
  }
  
  # Remove survey observations where they are missing in the real data
  for (i in surveyIndex) {
    years2include <- unique(fit$data$aux[, "year"][fit$data$aux[, "fleet"] == i])
    logStru_N[, !(colnames(logStru_N[,,i]) %in% years2include), i] <- NA
    logSobs_N[, !(colnames(logStru_N[,,i]) %in% years2include), i] <- NA
  }
  
  Stru_N <- exp(logStru_N)
  Sobs_N <- exp(logSobs_N)
  

  
  
  trueParams <- list(sdrep = fit$sdrep, pl = fit$pl)
  trueParams$pl$logN <- logN
  dimnames(f) <- list(paste0("tru.", c(1:nA)), fit$data$years)
  trueParams$pl$logF <- logF
  return(list(trueParams = trueParams,
              N = N, 
              logN = logN,
              logCobs_N = logCobs_N,
              logCtru_N = logCtru_N, 
              logStru_N = logStru_N, 
              logSobs_N = logSobs_N,
              Cobs_mt = Cobs_mt, 
              Cobs_N = Cobs_N, 
              Ctru_mt = Ctru_mt, 
              Ctru_N = Ctru_N, 
              Sobs_N = Sobs_N, 
              Stru_N = Stru_N))
  
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
    lapply(function(x) x[, which(colSums(x, na.rm = TRUE) != 0)]) %>%
    # Remove years (rows) without any data (¡caution if using data other than herring!)
    lapply(function(x) x[which(rowSums(x, na.rm = TRUE) != 0), ])
    
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

# Setup data and params for sam model #####################
setupModel <- function(conf) {
  
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
  conf <- conf
  
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
  if (!is.null(simOut$N)) {
    N <- simOut$N  
  } else { N <- exp(simOut$logN); rownames(N) <- paste0("tru.", 1:nrow(N))}
  
  
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
  
  # Put them together
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


## Plot parameters fit vs true ############################
plotPars <- function(fitSim, simOut) {
  nRepAccept <- length(simOut)
  parsFixed <- which(names(fitSim[[1]]$pl) %in% c("logSdLogN", "logSdLogFsta"))#names(fitSim[[1]]$obj$par))
  df_parsOut <- data.frame()
  for (h in parsFixed) {
    for (i in 1:nRepAccept) {
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
                     est_se = sd(est)/sqrt(nRepAccept))
  
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

## Calculate random effect timeseries error #######################
calcNFTsError <- function(fitSim, simOut) {
  indNF <- which(names(fitSim[[1]]$pl) %in% c("logN", "logF"))
  errRe <- data.frame()
  for (h in indNF) {
    for (i in 1:nRepAccept) {
      rownames(fitSim[[i]]$pl[[h]]) <- paste0("fit.", 1:nrow(fitSim[[i]]$pl[[h]]))
      sdLog <- fitSim[[i]]$plsd[[h]]
      rownames(sdLog) <- paste0("sdLog.", 1:nrow(fitSim[[i]]$pl[[h]]))
      errRe <-
        rbind(errRe,
              fitSim[[i]]$pl[[h]] %>%
                t() %>%
                exp %>%
                as.data.frame() %>%
                cbind(sdLog %>%
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
                              error_pc = 100 * (fit - tru) / tru,
                              decile = ceiling(10 * pnorm(q    = log(tru), 
                                                          mean = log(fit), 
                                                          sd   = sdLog)),
                              replicate = i,
                              variable = ifelse(names(fitSim[[i]]$pl[h]) == "logN", "N",
                                                if(names(fitSim[[i]]$pl[h]) == "logF") "F")))
    }
  }
  return(errRe)
}

## Calculate random effect timeseries error for simulate.sam output #######################
calcNFTsErrorSAM <- function(fitSimSAM, simOutSAM4error) {
  indNF <- which(names(fitSimSAM[[1]]$pl) %in% c("logN", "logF"))
  indNFsimOut <- which(names(simOutSAM4error[[1]]) %in% c("logN", "logF"))
  errRe <- data.frame()
  for (h in 1:length(indNF)) {
    for (i in 1:nRepSAMAccept) {
      rownames(fitSimSAM[[i]]$pl[[indNF[h]]]) <- paste0("fit.", 1:nrow(fitSimSAM[[i]]$pl[[indNF[h]]]))
      sdLog <- fitSimSAM[[i]]$plsd[[indNF[h]]]
      rownames(sdLog) <- paste0("sdLog.", 1:nrow(fitSimSAM[[i]]$pl[[indNF[h]]]))
      errRe <-
        rbind(errRe,
              fitSimSAM[[i]]$pl[[indNF[h]]] %>%
                t() %>%
                exp %>%
                as.data.frame() %>%
                cbind(sdLog %>%
                        t() %>%
                        as.data.frame()) %>%
                cbind(data.frame(simOutSAM4error[[i]][[indNFsimOut[h]]] %>% t %>% exp)) %>%
                dplyr::mutate(year = as.numeric(fitSimSAM[[i]]$data$years)) %>%
                tidyr::gather(variable, N, -year) %>%
                dplyr::mutate(variable = gsub(x = variable, 
                                              pattern = "X", 
                                              replacement = "tru.")) %>%
                tidyr::separate(variable, c("source", "age")) %>%
                tidyr::spread(source, N) %>%
                dplyr::mutate(age = as.numeric(age),
                              error = (fit - tru),
                              error_pc = 100 * (fit - tru) / tru,
                              decile = ceiling(10 * pnorm(q    = log(tru), 
                                                          mean = log(fit), 
                                                          sd   = sdLog)),
                              replicate = i,
                              variable = ifelse(names(fitSimSAM[[i]]$pl[indNF[h]]) == "logN", "N",
                                                if(names(fitSimSAM[[i]]$pl[indNF[h]]) == "logF") "F")))
    }
  }
  return(errRe)
}

## Calculate observed time series error ###################
calcCatchError <- function(fitSim, simOut) {
  
  errC <- data.frame()
  for (i in 1:length(fitSim)) {
    # Fit catch
    df_Cfit_mt <- 
      data.frame(variable = names(fitSim[[i]]$sdrep$value),
                 value = fitSim[[i]]$sdrep$value,
                 sd    = fitSim[[i]]$sdrep$sd) %>%
      dplyr::filter(variable == "logCatch") %>%
      dplyr::rename(logCatch = value,
                    sdLog = sd) %>%
      dplyr::mutate(fit = exp(logCatch),
                    year = fitSim[[i]]$data$years,
                    age = "total") %>%
      dplyr::select(year, age, fit, sdLog)
    
    # True catch
    df_Ctru_mt <- 
      simOut[[i]]$Ctru_mt %>%
      t() %>%
      as.data.frame() %>%
      dplyr::mutate(total = rowSums(.),
                    year = fitSim[[1]]$data$years) %>%
      tidyr::gather(age, Catch_mt, -year) %>%
      dplyr::rename(tru = Catch_mt) %>%
      dplyr::filter(age == "total")
    
    
    errC <-
      rbind(errC,
            df_Cfit_mt %>%
              dplyr::left_join(df_Ctru_mt) %>%
              dplyr::mutate(error = fit - tru,
                            error_pc = 100 * (fit - tru) / tru,
                            decile = ceiling(10 * pnorm(q    = log(tru),
                                                        mean = log(fit),
                                                        sd   = sdLog)),
                            replicate = i,
                            variable = "catch"))
  }  
  
  
  return(errC)
}

## Plot timeseries error ##################################
plotTsError <- function(err) {

  nRepAccept <- length(unique(err$replicate))
  
  # # Calculate median error
  # errAnnual <-
  #   err %>%
  #   dplyr::group_by(variable, age, year) %>%
  #   dplyr::summarise(median_error = median(error))
  
  # Plot median error for N and F in each year for each age
  # p <-
  #   ggplot(errAnnual %>% 
  #            dplyr::filter(variable == "N") %>%
  #            dplyr::mutate(as.numeric(year)),
  #          aes(x = year, y = median_error)) +
  #     geom_line() +
  #     geom_hline(yintercept = 0) +
  #     theme_bw() +
  #     ylab("Median raw error (fit - true)") +
  #     xlab("Year") +
  #     theme(axis.title = element_text(size = 16),
  #           axis.text = element_text(size = 14)) +
  #     facet_wrap(~age) +
  #     ggtitle("N estimation error-at-age")
  # print(p)
  
  # p <-
  #   ggplot(errAnnual %>% 
  #            dplyr::filter(variable == "F") %>%
  #            dplyr::mutate(as.numeric(year)),
  #          aes(x = year, y = median_error)) +
  #     geom_line() +
  #     geom_hline(yintercept = 0) +
  #     theme_bw() +
  #     ylab("Median raw error (fit - true)") +
  #     xlab("Year") +
  #     theme(axis.title = element_text(size = 16),
  #           axis.text = element_text(size = 14)) +
  #     facet_wrap(~age) +
  #     ggtitle("F estimation error-at-age")
  # print(p)
  
  # if(any(errAnnual$variable == "catch")) {
  #   p <-
  #     ggplot(errAnnual %>% 
  #              dplyr::filter(variable == "catch") %>%
  #              dplyr::mutate(as.numeric(year)),
  #            aes(x = year, y = median_error)) +
  #     geom_line() +
  #     geom_hline(yintercept = 0) +
  #     theme_bw() +
  #     ylab("Median raw error (fit - true)") +
  #     xlab("Year") +
  #     theme(axis.title = element_text(size = 16),
  #           axis.text = element_text(size = 14)) +
  #     facet_wrap(~age) +
  #     ggtitle("Catch estimation error")
  #   print(p)
  # }

  
  # Plot a few example fit vs tru time series
  p <-
    ggplot(err %>%
             dplyr::select(-sdLog, -decile, -error_pc) %>%
             dplyr::filter(variable == "F",
                           age == 2,
                           replicate %in% 1:10) %>%
             tidyr::gather(source, value, 
                           -year, -age, -replicate, -variable),
           aes(x = year, y = value, color = source)) +
      geom_line() +
      facet_wrap(~replicate) + 
      ggtitle("Age-2 F estimation error examples")
  print(p)
  
  # Plot decile coverage for N and F
  p <- 
    ggplot(err %>% dplyr::filter(variable %in% c("N", "F")),
           aes(x = decile)) +
      geom_histogram(binwidth=1, colour="white") +
      geom_hline(yintercept = ncol(fitSim[[1]]$pl$logN) * nRepAccept / 10, 
                 color = "dark grey") +
      theme_bw() +
      facet_grid(variable~age) +
      ggtitle("Confidence interval coverage for each age")
  print(p)
  
  # Plot decile coverage for Catch
  if(any(err$variable == "catch")) {
    p <- 
      ggplot(err %>% dplyr::filter(variable == "catch"),
             aes(x = decile)) +
      geom_histogram(binwidth=1, colour="white") +
      geom_hline(yintercept = ncol(fitSim[[1]]$pl$logN) * nRepAccept / 10, 
                 color = "dark grey") +
      theme_bw() +
      facet_wrap(~age) +
      ggtitle("Confidence interval coverage for catch")
    print(p)
  }
  
  # Plot relationship between percent error and true value for N and F
  # p <-
  #   ggplot(err  %>% dplyr::filter(variable %in% c("N", "F")),
  #          aes(x = tru, y = error_pc)) +
  #     geom_point(alpha = 0.2) +
  #     theme_bw() +
  #     facet_wrap(variable~age, scales = "free", nrow = 2) +
  #     xlab("True value") +
  #     ylab("Percent error") +
  #     ggtitle("Percent error vs true value")
  # print(p)

    
  err2plot <-
    err %>%  
    dplyr::filter(variable %in% c("N")) %>%
    dplyr::filter(year != min(year)) %>% # exclude initial condition which is identical across replicates
    dplyr::select(year, age, replicate, tru, fit) %>%
    dplyr::rename(N_tru = tru,
                  N_fit = fit) %>%
    dplyr::left_join({err %>%  
                      dplyr::filter(year != min(year)) %>%
                      dplyr::filter(variable %in% c("F")) %>%
                      dplyr::select(year, age, replicate, error_pc) %>%
                      dplyr::rename(F_error_pc = error_pc)}) %>%
    dplyr::left_join({err %>%  
                      dplyr::filter(year != min(year)) %>%
                      dplyr::filter(variable %in% c("N")) %>%
                      dplyr::select(year, age, replicate, error_pc) %>%
                      dplyr::rename(N_error_pc = error_pc)})
  
  
  # N
  # p <-
  #   ggplot(err2plot %>% dplyr::filter(N_tru < quantile(N_tru, .9)),
  #          aes(x = N_tru)) +
  #   geom_point(aes(y = N_error_pc), alpha = 0.2) +
  #   geom_hline(yintercept = 0) +
  #   facet_wrap(~age) +
  #   theme_bw() +
  #   xlab("N (1000's)") +
  #   ylab("Percent error of N estimate") +
  #   ggtitle("Percent error in N vs N value (<90th percentile)")
  # print(p)
  
  # F
  # p <-
  #   ggplot(err2plot %>% dplyr::filter(N_tru < quantile(N_tru, .9)),
  #          aes(x = N_tru)) +
  #     geom_point(aes(y = F_error_pc), alpha = 0.2) +
  #     geom_hline(yintercept = 0) +
  #     facet_wrap(~age) +
  #     theme_bw() +
  #     xlab("N (1000's)") +
  #     ylab("Percent error of F estimate") +
  #   ggtitle("Percent error in F vs N value (<90th percentile)")
  # print(p)
    
  
  # Plot mean percent error vs percentile of N_tru
  err2plot_percentile <-
    err2plot %>%
    dplyr::group_by(age) %>%
    dplyr::mutate(N_percentile = ntile(N_tru, 100)) %>%
    dplyr::group_by(age, N_percentile) %>%
    dplyr::summarise(F_error_pc_mean = mean(F_error_pc),
                     N_error_pc_mean = mean(N_error_pc),
                     nObsF = length(F_error_pc),
                     nObsN = length(N_error_pc),
                     F_error_pc_hi   = F_error_pc_mean + 1.96 * sd(F_error_pc)/sqrt(nObsF),
                     F_error_pc_lo  = F_error_pc_mean - 1.96 * sd(F_error_pc)/sqrt(nObsF),
                     N_error_pc_hi  = N_error_pc_mean + 1.96 * sd(N_error_pc)/sqrt(nObsN),
                     N_error_pc_lo  = N_error_pc_mean - 1.96 * sd(N_error_pc)/sqrt(nObsN),
                     N_tru_value     = mean(N_tru),
                     N_fit_max       = quantile(N_fit, .95),
                     N_fit_min       = quantile(N_fit, .05))

  
  # N
  p <-
    ggplot(err2plot_percentile,
           aes(x = N_percentile)) +
    geom_line(aes(y = N_error_pc_mean)) +
    geom_ribbon(aes(ymin = N_error_pc_lo, ymax = N_error_pc_hi), alpha = 0.3) +
    geom_hline(yintercept = 0) +
    facet_wrap(~age) +
    theme_bw() +
    xlab("Percentile of N") +
    ylab("Percent error of N estimate") +
    ggtitle("Mean percent error of N vs Percentile of N")
  print(p)
  
  # F
  p <-
    ggplot(err2plot_percentile,
           aes(x = N_percentile)) +
    geom_line(aes(y = F_error_pc_mean)) +
    geom_ribbon(aes(ymin = F_error_pc_lo, ymax = F_error_pc_hi), alpha = 0.3) +
    geom_hline(yintercept = 0) +
    facet_wrap(~age) +
    theme_bw() +
    xlab("Percentile of N") +
    ylab("Percent error of F estimate") +
    ggtitle("Mean percent error of F vs Percentile of N")
  print(p)
  
  # Plot percentile of N vs value of N
  p <-
    ggplot(err2plot_percentile %>% dplyr::filter(N_percentile < 10),
           aes(x = N_percentile)) +
    geom_ribbon(aes(ymin = N_fit_min, ymax = N_fit_max), fill = "blue", alpha = 0.3) +
    facet_wrap(~age, scales = "free_y") +
    xlab("N percentile") +
    ylab("Range of N fit (1000's)")
  print(p)
  
  
  
  # Plot relationship between percent error and true value for N and F
  # if(any(errAnnual$variable == "catch")) {
  #   p <-
  #     ggplot(err  %>% dplyr::filter(variable == "catch"),
  #            aes(x = tru, y = error_pc)) +
  #     geom_point(alpha = 0.2) +
  #     theme_bw() +
  #     facet_wrap(~age, scales = "free", nrow = 2) +
  #     xlab("True value") +
  #     ylab("Percent error") +
  #     ggtitle("Percent error vs true value for Catch")
  #   print(p)
  # }
    
  
}

## Plot mean error over all replicates ####################
plotTsMeanError <- function(err) {
  errMean <-
    err %>%
    dplyr::group_by(variable, age) %>%
    dplyr::summarise(error_mean = mean(error),
                     error_mean_se   = sd(error) / sqrt(nRepAccept),
                     error_pc_mean = mean(error_pc),
                     error_pc_mean_se = sd(error_pc) / sqrt(nRepAccept))
  
  # Plot raw error
  p <-
    ggplot(errMean, aes(x = age)) +
      geom_hline(aes(yintercept = 0), color = "black") +
      geom_point(aes(y = error_mean), color = "blue") +
      geom_errorbar(aes(ymin = error_mean - 1.96 * error_mean_se,
                        ymax = error_mean + 1.96 * error_mean_se),
                    width = 0.2,
                    color = "blue") +
      facet_wrap(~variable, scales = "free", nrow = 2) +
      ylab("Mean raw error (fit - true)") +
      xlab("Age") +
      ggtitle("Mean raw error over all replicates")
  print(p)
  
  # Plot percent error
  p <-
    ggplot(errMean, aes(x = age)) +
      geom_hline(aes(yintercept = 0), color = "black") +
      geom_point(aes(y = error_pc_mean), color = "blue") +
      geom_errorbar(aes(ymin = error_pc_mean - 1.96 * error_pc_mean_se,
                        ymax = error_pc_mean + 1.96 * error_pc_mean_se),
                    width = 0.2,
                    color = "blue") +
      facet_wrap(~variable, scales = "free", nrow = 2) +
      ylab("Mean percent error (100 * (fit - true) / true)") +
      xlab("Age") +
      ggtitle("Mean percent error over all replicates")
  print(p)
}



## Customized setup.sam.data ##############################
# to avoid logging data when it is alread logged
cp_setup.sam.data <- function (fleets = NULL, surveys = NULL, residual.fleet = NULL, 
                               prop.mature = NULL, stock.mean.weight = NULL, catch.mean.weight = NULL, 
                               dis.mean.weight = NULL, land.mean.weight = NULL, natural.mortality = NULL, 
                               prop.f = NULL, prop.m = NULL, land.frac = NULL, recapture = NULL) 
{
  fleet.idx <- 0
  type <- NULL
  time <- NULL
  name <- NULL
  corList <- list()
  idxCor <- matrix(NA, nrow = length(fleets) + length(surveys) + 
                     1, ncol = nrow(natural.mortality))
  colnames(idxCor) <- rownames(natural.mortality)
  dat <- data.frame(year = NA, fleet = NA, age = NA, aux = NA)
  weight <- NULL
  doone <- function(m) {
    year <- rownames(m)[row(m)]
    fleet.idx <<- fleet.idx + 1
    fleet <- rep(fleet.idx, length(year))
    age <- as.integer(colnames(m)[col(m)])
    aux <- as.vector(m)
    dat <<- rbind(dat, data.frame(year, fleet, age, aux))
    if ("weight" %in% names(attributes(m))) {
      weight <<- c(weight, as.vector(attr(m, "weight")))
    }
    else {
      if ("cov" %in% names(attributes(m))) {
        weight <<- c(weight, unlist(lapply(attr(m, "cov"), 
                                           diag)))
      }
      else {
        if ("cov-weight" %in% names(attributes(m))) {
          weight <<- c(weight, 1/unlist(lapply(attr(m, 
                                                    "cov-weight"), diag)))
        }
        else {
          weight <<- c(weight, rep(NA, length(year)))
        }
      }
    }
    if ("cov" %in% names(attributes(m))) {
      attr(m, "cor") <- lapply(attr(m, "cov"), cov2cor)
    }
    if ("cov-weight" %in% names(attributes(m))) {
      attr(m, "cor") <- lapply(attr(m, "cov-weight"), 
                               cov2cor)
    }
    if ("cor" %in% names(attributes(m))) {
      thisCorList <- attr(m, "cor")
      whichCorOK <- which(unlist(lapply(thisCorList, function(x) !any(is.na(x)))))
      thisCorList <- thisCorList[whichCorOK]
      corList <<- c(corList, thisCorList)
      nextIdx <- if (all(is.na(idxCor))) {
        0
      }
      else {
        max(idxCor, na.rm = TRUE)
      }
      idxCor[fleet.idx, colnames(idxCor) %in% rownames(m)][whichCorOK] <<- nextIdx:(nextIdx + 
                                                                                      length(thisCorList) - 1)
    }
  }
  if (!is.null(residual.fleet)) {
    doone(residual.fleet)
    type <- c(type, 0)
    time <- c(time, 0)
    name <- c(name, "Residual catch")
  }
  if (!is.null(fleets)) {
    if (is.data.frame(fleets) | is.matrix(fleets)) {
      doone(fleets)
      type <- c(type, 1)
      time <- c(time, 0)
      name <- c(name, "Comm fleet")
    }
    else {
      dummy <- lapply(fleets, doone)
      type <- c(type, rep(1, length(fleets)))
      time <- c(time, rep(0, length(fleets)))
      name <- c(name, strtrim(gsub("\\s", "", names(dummy)), 
                              50))
    }
  }
  if (!is.null(surveys)) {
    if (is.data.frame(surveys) | is.matrix(surveys)) {
      doone(surveys)
      thistype <- ifelse(min(as.integer(colnames(surveys))) < 
                           (-0.5), 3, 2)
      type <- c(type, thistype)
      time <- c(time, mean(attr(surveys, "time")))
      name <- c(name, "Survey fleet")
    }
    else {
      dummy <- lapply(surveys, doone)
      type <- c(type, unlist(lapply(surveys, function(x) ifelse(min(as.integer(colnames(x))) < 
                                                                  (-0.5), 3, 2))))
      time <- c(time, unlist(lapply(surveys, function(x) mean(attr(x, 
                                                                   "time")))))
      name <- c(name, strtrim(gsub("\\s", "", names(dummy)), 
                              50))
    }
  }
  if (is.null(land.frac)) {
    land.frac <- matrix(1, nrow = nrow(residual.fleet), 
                        ncol = ncol(residual.fleet))
  }
  if (is.null(dis.mean.weight)) {
    dis.mean.weight <- catch.mean.weight
  }
  if (is.null(land.mean.weight)) {
    land.mean.weight <- catch.mean.weight
  }
  if (is.null(prop.f)) {
    prop.f <- matrix(0, nrow = nrow(residual.fleet), ncol = ncol(residual.fleet))
  }
  if (is.null(prop.m)) {
    prop.m <- matrix(0, nrow = nrow(residual.fleet), ncol = ncol(residual.fleet))
  }
  dat$aux[which(dat$aux <= 0)] <- NA
  dat <- dat[!is.na(dat$year), ]
  if (!is.null(recapture)) {
    tag <- data.frame(year = recapture$ReleaseY)
    fleet.idx <- fleet.idx + 1
    tag$fleet <- fleet.idx
    tag$age <- recapture$ReleaseY - recapture$Yearclass
    tag$aux <- exp(recapture$r)
    tag <- cbind(tag, recapture[, c("RecaptureY", "Yearclass", 
                                    "Nscan", "R", "Type")])
    dat[names(tag)[!names(tag) %in% names(dat)]] <- NA
    dat <- rbind(dat, tag)
    weight <- c(weight, rep(NA, nrow(tag)))
    type <- c(type, 5)
    time <- c(time, 0)
    name <- c(name, "Recaptures")
  }
  dat <- dat[complete.cases(dat[, 1:3]), ]
  o <- order(as.numeric(dat$year), as.numeric(dat$fleet), 
             as.numeric(dat$age))
  attr(dat, "type") <- type
  names(time) <- NULL
  attr(dat, "time") <- time
  names(name) <- NULL
  attr(dat, "name") <- name
  dat <- dat[o, ]
  weight <- weight[o]
  newyear <- min(as.numeric(dat$year)):max(as.numeric(dat$year))
  newfleet <- min(as.numeric(dat$fleet)):max(as.numeric(dat$fleet))
  mmfun <- function(f, y, ff) {
    idx <- which(dat$year == y & dat$fleet == f)
    ifelse(length(idx) == 0, NA, ff(idx) - 1)
  }
  idx1 <- outer(newfleet, newyear, Vectorize(mmfun, c("f", 
                                                      "y")), ff = min)
  idx2 <- outer(newfleet, newyear, Vectorize(mmfun, c("f", 
                                                      "y")), ff = max)
  attr(dat, "idx1") <- idx1
  attr(dat, "idx2") <- idx2
  attr(dat, "minAgePerFleet") <- tapply(as.integer(dat[, "age"]), 
                                        INDEX = dat[, "fleet"], FUN = min)
  attr(dat, "maxAgePerFleet") <- tapply(as.integer(dat[, "age"]), 
                                        INDEX = dat[, "fleet"], FUN = max)
  attr(dat, "year") <- newyear
  attr(dat, "nyear") <- max(as.numeric(dat$year)) - min(as.numeric(dat$year)) + 
    1
  cutY <- function(x) x[rownames(x) %in% newyear, ]
  attr(dat, "prop.mature") <- cutY(prop.mature)
  attr(dat, "stock.mean.weight") <- cutY(stock.mean.weight)
  attr(dat, "catch.mean.weight") <- cutY(catch.mean.weight)
  attr(dat, "dis.mean.weight") <- cutY(dis.mean.weight)
  attr(dat, "land.mean.weight") <- cutY(land.mean.weight)
  attr(dat, "natural.mortality") <- cutY(natural.mortality)
  attr(dat, "prop.f") <- cutY(prop.f)
  attr(dat, "prop.m") <- cutY(prop.m)
  attr(dat, "land.frac") <- cutY(land.frac)
  ret <- list(noFleets = length(attr(dat, "type")), 
              fleetTypes = as.integer(attr(dat, "type")),
              sampleTimes = attr(dat, "time"), 
              noYears = attr(dat, "nyear"), 
              years = attr(dat, "year"), 
              minAgePerFleet = attr(dat, "minAgePerFleet"), 
              maxAgePerFleet = attr(dat, "maxAgePerFleet"), 
              nobs = nrow(dat), idx1 = attr(dat, "idx1"), 
              idx2 = attr(dat,"idx2"), idxCor = idxCor, 
              aux = data.matrix(dat[,-4]), 
              logobs = dat[, 4], #log(dat[, 4]), CHANGED TO ACCOMODATE LOGGED OBSERVATIONS
              weight = as.numeric(weight), 
              propMat = attr(dat, "prop.mature"), 
              stockMeanWeight = attr(dat, "stock.mean.weight"), 
              catchMeanWeight = attr(dat, "catch.mean.weight"), 
              natMor = attr(dat, "natural.mortality"), 
              landFrac = attr(dat, "land.frac"), 
              disMeanWeight = attr(dat, "dis.mean.weight"), 
              landMeanWeight = attr(dat, "land.mean.weight"), 
              propF = attr(dat, "prop.f"), propM = attr(dat, "prop.m"), 
              corList = corList)
  attr(ret, "fleetNames") <- attr(dat, "name")
  return(ret)
}


