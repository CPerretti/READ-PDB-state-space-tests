# Functions for simulations

#### Adjust plot properties ###############################


## Simulation model #######################################
sim <- function(fit, keyLogScale, noScaledYears, container_i) {
  nA <- ncol(fit$data$propF) # number of age-classes
  nT <- fit$data$noYears # length of time series
  
  
  # Setup keyLogScale to have unique logScale for each age that is fished
  nAs <- sum(keyLogScale[1,] > -1)
  
  switch(container_i$scenario,
         `uniform random` = {
           logScale <- matrix(data = log(runif(nAs * noScaledYears, 1.5, 10)),
                              nrow = nAs, ncol = noScaledYears)
         },
         `random walk` = { # RW from 1 with reflecting boundary at 1
           logSdLogScale <- log(0.2)
           rw_logScale_mat <- matrix(data = NA, nrow = nAs, ncol = noScaledYears)
           errS <- matrix(data = rnorm(nAs * noScaledYears, 0, exp(logSdLogScale)),
                          nrow = nAs, ncol = noScaledYears) #uncorrelated error
           #errS[errS[,1] < 0, 1] <- -errS[errS[,1] < 0, 1] #reflecting boundary
           rw_logScale_mat[,1] <- 1 + errS[,1]
           for(i in 2:noScaledYears){
             rw_logScale_mat[,i] <- rw_logScale_mat[,i-1] + errS[,i]
             #ind <- rw_logScale_mat[,i] < 1 # index to reflect
             #rw_logScale_mat[ind, i] <- rw_logScale_mat[ind, i-1] - errS[ind,i]
           }
           
           logScale <- matrix(data = rw_logScale_mat, nrow = nAs, ncol = noScaledYears)
         },
         fixed = { 
           # Misreporting on all ages
           logScale <- matrix(data = rep(log(runif(1, 1.5, 10)),
                                         times = nAs * noScaledYears),
                              nrow = nAs, ncol = noScaledYears)
           # Misreporting only on ages 1-3
           # logScale <- matrix(data = rep(c(log(runif(1, 1.5, 10)), 0), 
           #                               each = nAs * noScaledYears / 2),
           #                    nrow = nAs, ncol = noScaledYears, byrow = T)
         },
         `no misreporting` = {
           logScale <- matrix(data = log(1), nrow = nAs, ncol = noScaledYears)
         }
  ) 
  
  # Set F (need to replicate some elements to match ModelConf)
  #f <- exp(fit$pl$logF[(fit$conf$keyLogFsta[1,] + 1),])
  # Calculate a new relization of f errors with sd's from the fit
  errF <- matrix(data = NA,
                 nrow = nA, 
                 ncol = nT - 1)
  
  # Set F sd
  fit$pl$logSdLogFsta <- 
    (c(1, rep(1, length(fit$pl$logSdLogFsta)-1)) * exp(fit$pl$logSdLogFsta)) %>%
    "*"(1) %>%
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
  
  logF[, 1] <- rnorm(n = nA, 
                     mean = fit$pl$logF[, 1][fit$conf$keyLogFsta[1, ] + 1],
                     sd = fit$plsd$logF[, 1][fit$conf$keyLogFsta[1, ] + 1]) # initial F
  
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
  
  logN[, 1] <- rnorm(n = nA, 
                     mean = fit$pl$logN[, 1],
                     sd = fit$plsd$logN[, 1])# initial N
  
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
  # Possibly lower N process error (!)
  fit$pl$logSdLogN[(fit$conf$keyVarLogN + 1)] <-
    fit$pl$logSdLogN[(fit$conf$keyVarLogN + 1)] %>%
    exp %>% "*"(1) %>% log
  
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
  
  SSB <- (exp(logN)*exp(-exp(logF)*t(fit$data$propF)-t(fit$data$natMor)*t(fit$data$propM))*
         t(fit$data$propMat)*t(fit$data$stockMeanWeight)) %>% colSums
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
  fit$pl$logSdLogObs <- fit$pl$logSdLogObs %>% exp %>% "*"(1) %>% log
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
  
  logCobs_N0 <- logCtru_N + errObs[, , "Residual catch"]
  logCobs_N <- logCobs_N0
  logScale_full <- cbind(matrix(data = 0, nrow = nrow(logScale), ncol = nT - noScaledYears),
                         logScale)
  
  for (a in 1:nA){
    if (fit$conf$keyLogFsta[1, a] > -1) { # If there is fishing on an age
      logCobs_N[a,] <- logCobs_N0[a,] - logScale_full[keyLogScale[1,a] + 1,] # Misreport catch
    }
  }
    
  
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
  logSq[which(fit$conf$keyLogFpar != -1)] <- # fill with real fit values
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
  trueParams$pl$logScale <- logScale
  if (container_i$scenario == "random walk") trueParams$pl$logSdLogScale <- logSdLogScale
  trueParams$pl$logN <- logN
  dimnames(f) <- list(paste0("tru.", c(1:nA)), fit$data$years)
  trueParams$pl$logF <- logF
  return(list(trueParams = trueParams,
              N = N, 
              logN = logN,
              SSB = SSB,
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
    attr(surveys[[i]], "time", rep(fit$data$sampleTimes[surveyIndex[i]],2))
  }
  
  # Export the simulated surveys to a text file, then import it
  # using read.ices().
  
  # First set survey file header
  write.table(rbind("Fake Fish Survey Data", 100 + length(surveys)), 
              file = "./sim_data/surveys.dat", sep = " ", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  # Then append surveys
  for (i in 1:length(surveys)) { # loop over surveys
    header_survey <-
      rbind(names(surveys)[i], 
            paste(min(as.numeric(rownames(surveys[[i]]))),
                  max(as.numeric(rownames(surveys[[i]])))),
            paste(1, 1, fit$data$sampleTimes[surveyIndex[i]], 
                  fit$data$sampleTimes[surveyIndex[i]]),
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
    rbind("Fake Fish Total Catch Numbers at age (000s; combines all gear types)", 
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
setupModel <- function(conf, example_dir, noScaledYears, confLogScale) {
  
  # Read in data
  cn <- read.ices("./sim_data/catch.dat") # catch abundace-at-age
  surveys <- read.ices("./sim_data/surveys.dat") #surveys
  
  cw <- read.ices(paste0("../", example_dir, "/data/cw.dat")) # catch mean weight-at-age
  #dw <- cw # discards mean weight-at-age (using catch mean weight-at-age for now)
  #lw <- cw # landings mean weight-at-age (using catch mean weight-at-age for now)
  dw <- read.ices(paste0("../", example_dir, "/data/dw.dat")) 
  lw <- read.ices(paste0("../", example_dir, "/data/lw.dat"))
  pf <- read.ices(paste0("../", example_dir, "/data/pf.dat")) # proportion of f before spawning
  #lf <- pf; lf[,] <- 1 # fraction of catch that is landed (set to 1 for now)
  lf <- read.ices(paste0("../", example_dir, "/data/lf.dat"))
  mo <- read.ices(paste0("../", example_dir, "/data/mo.dat")) # maturity-at-age ogive
  nm <- read.ices(paste0("../", example_dir, "/data/nm.dat")) # natural mortality-at-age
  pm <- read.ices(paste0("../", example_dir, "/data/pm.dat")) # proportion of m before spawning
  sw <- read.ices(paste0("../", example_dir, "/data/sw.dat")) # stock weight-at-age (kg)
  
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
  
  # Try SAM misreported catch code
  switch(confLogScale$logScaleType,
         `random walk` = {
           conf$noScaledYears  <- confLogScale$noScaledYears
           conf$keyLogScale    <- confLogScale$keyLogScale
           conf$keyVarLogScale <- rep(0, sum(confLogScale$keyLogScale[1,] > -1))
           conf$keyScaledYears <- (max(dat$years) - confLogScale$noScaledYears + 1) : max(dat$years)
         },
         fixed = {
           conf$noScaledYears  <- confLogScale$noScaledYears
           conf$keyLogScale    <- confLogScale$keyLogScale
           conf$keyScaledYears <- (max(dat$years) - confLogScale$noScaledYears + 1) : 
                                   max(dat$years)
           conf$keyParScaledYA <- confLogScale$keyParScaledYA
           conf$constRecBreaks <- numeric(0)
         },
         `no misreporting` = {
           conf$constRecBreaks <- numeric(0)
         }
         )
  
  if (confLogScale$logScaleType == "random walk") {
    par <- stockassessment2::defpar(dat, conf) # some default starting values  
  } else {par <- stockassessment::defpar(dat, conf)}
  
  
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
    dplyr::mutate(age = paste0("age-", age)) #%>%
    #dplyr::mutate(source = ifelse(source == "fit", "NS cod fit", "Simulated"))
  
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

## Plot N-at-age simulated vs fit #######################
plotN_quants <- function(simOut, fit){
  # if (!is.null(simOut$N)) {
  #   N <- simOut$N  
  # } else { N <- exp(simOut$logN); rownames(N) <- paste0("tru.", 1:nrow(N))}
  # 
  
  df2plot_fit <- 
    data.frame(fit = fit$pl$logN %>% exp %>% t) %>%
    dplyr::mutate(year = fit$data$years) %>%
    tidyr::gather(variable, N, -year) %>%
    tidyr::separate(variable, c("source", "age")) %>%
    dplyr::mutate(age = paste0("age-", age))
  
  df2plot0 <- data.frame()
  for (i in 1:length(simOut)) {
    N <- simOut[[i]]$N
    df2plot0 <-
      rbind(df2plot0,
            {N %>%
             t() %>%
             as.data.frame() %>%
             #cbind(data.frame(fit = fit$pl$logN %>% exp %>% t)) %>%
             dplyr::mutate(year = fit$data$years) %>%
             tidyr::gather(variable, N, -year) %>%
             tidyr::separate(variable, c("source", "age")) %>%
             dplyr::mutate(age = paste0("age-", age),
                           replicate = i)})
  }

  df2plot <-
    df2plot0 %>%
    dplyr::group_by(year, age) %>%
    dplyr::summarise(#pc05 = quantile(N, probs = 0.05),
                     #pc10 = quantile(N, probs = 0.10),
                     #pc25 = quantile(N, probs = 0.25),
                     #`50th percetnile` = quantile(N, probs = 0.50)#,
                     `Mean of simulations` = mean(N)#,
                     #pc75 = quantile(N, probs = 0.75),
                     #pc90 = quantile(N, probs = 0.90),
                     #pc95 = quantile(N, probs = 0.95)
                     ) %>%
    tidyr::gather(variable, N, -year, -age)
    
  
  # Plot N-at-age 
  ggplot(data = df2plot, aes(x = year, y = N)) +
    geom_line(aes(color = variable)) +
    geom_line(data = df2plot_fit, aes(color = "NS cod fit")) +
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
    dplyr::mutate(source = "Fit") %>%
    dplyr::select(year, age, Catch_mt, source)
  
  # Pull out simulated observations
  df_Cobs_mt <- 
    Cobs_mt %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(total = rowSums(.),
                  year = fit$data$years) %>%
    tidyr::gather(age, Catch_mt, -year) %>%
    #dplyr::mutate(source = "observed")
    dplyr::mutate(source = "Observed")
  
  # Pull out simulated truth
  df_Ctru_mt <- 
    Ctru_mt %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(total = rowSums(.),
                  year = fit$data$years) %>%
    tidyr::gather(age, Catch_mt, -year) %>%
    #dplyr::mutate(source = "tru")  
    dplyr::mutate(source = "True")
  
  # Put them together
  df2plot <- 
    bind_rows(df_Cfit_mt, 
              df_Cobs_mt,
              df_Ctru_mt) %>%
    dplyr::mutate(age = paste0("age-", age))
  
  # Plot Catch (should exactly match in total subplot when using errPro_exact)
  ggplot(data = df2plot %>% 
           dplyr::filter(age == "age-total",
                         source != "Fit"),
         aes(x = year, y = Catch_mt, color = source)) +
    geom_line() +
    #facet_wrap(~age, scales = "free") +
    ylab("Catch (MT)") +
    ggtitle("Fishery catch") + #-at-age") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 12))
  
}

## Plot C_N-at-age simulated vs fit ####################### 
# Numbers at age *not* MT.
# fitSim[[1]]$data$aux #<-- Use to tell which is commerical catch, and which age
# fitSim[[1]]$data$logobs #<-- Use to pull out observations
# fitSim[[1]]$rep #<-- Use to pull out model fit


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
               value = exp(fit$rep$predObs[idx]), 
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
plotPars <- function(container, models2plot) {
colors2use <- RColorBrewer::brewer.pal(3, "Dark2")
  
df_parsOut <- data.frame()
  for (model in models2plot) {
    switch(model,
           `random walk` = {
             fitSim <- container$fitSim_random
           },
           fixed = {
             fitSim <- container$fitSim_fixed
           },
           `no misreporting` = {
             fitSim <- container$fitSim_none
           }
    )
    pars2plot <- which(names(fitSim[[1]]$pl) %in% names(fitSim[[1]]$obj$par))
    
    for (h in pars2plot) {
      for (i in 1:nrow(container)) {
        h_tru <- which(names(container$simOut[[i]]$trueParams$pl) == names(fitSim[[i]]$pl[h]))
        df_parsOut <-
          rbind(df_parsOut,
                data.frame(variable = paste(names(fitSim[[1]]$pl)[h], 
                                            1:length(fitSim[[1]]$pl[[h]]), sep = "."),
                           # Don't plot tru if it isn't the same length as estimated
                           # becuase it means the values don't match up.
                           tru = (if ((length(h_tru) >0) && 
                                      length(container$simOut[[i]]$trueParams$pl[[h_tru]]) == 
                                      length(fitSim[[i]]$pl[[h]])) {
                             container$simOut[[i]]$trueParams$pl[[h_tru]] 
                           } else {NA}),
                           est = fitSim[[i]]$pl[[h]],
                           sd  = fitSim[[i]]$plsd[[h]],
                           replicate = container$replicate[[i]],
                           scenario  = paste(container$scenario[[i]], "scenario"),
                           model = model))
      }  
    }
  }  
    
  
  df2plot <-
    df_parsOut %>%
    dplyr::mutate(model  = factor(model, levels = c("no misreporting", 
                                                    "fixed",
                                                    "random walk"))) %>%
    dplyr::group_by(model, scenario, variable) %>%
    dplyr::summarise(tru = unique(tru),
                     est_mean = mean(est),
                     est_se = sd(est)/sqrt(length(unique(replicate))))
  
  scenarios2plot <- unique(df2plot$scenario)
  for (i in 1:length(scenarios2plot)) {
    p <-
      ggplot(df2plot %>% 
               dplyr::filter(scenario == scenarios2plot[i]),
             aes(y = variable)) +
        geom_point( aes(x = tru), shape = 124, 
                    color = "black", size  = 12) +
        geom_point( aes(x = est_mean, color = model), size = 3) +
        geom_errorbarh(aes(xmin = est_mean - 1.96*est_se,
                           xmax = est_mean + 1.96*est_se,
                           color = model)) +
        theme_bw() +
        ylab("") +
        xlab("Value") +
        theme(axis.title = element_text(size = 16),
              axis.text = element_text(size = 14)) +
        scale_color_manual(values = colors2use) +
        ggtitle(scenarios2plot[i])
    print(p)
  }
}

## Calculate random effect timeseries error #######################
calcReTsError <- function(fitSim, simOut, confLogScale) {
  indRe <- which(names(fitSim$pl) %in% c("logN", "logF"))
  if (confLogScale$logScaleType %in% c("random walk", "fixed")) {
    indRe <- c(indRe, which(names(fitSim$pl) == "logScale"))
  }
  errRe <- data.frame()

  for (h in indRe) {
    varName <- substr(names(fitSim$pl[h]), 4, 999)
    year <- if (varName == "Scale") {
      as.numeric(fitSim$conf$keyScaledYears)
    } else {as.numeric(fitSim$data$years)}

    #for (i in 1:nRepAccept) {
    if (confLogScale$logScaleType == "fixed" & varName == "Scale") {
      fitSim$pl[[h]] <- matrix(fitSim$pl[[h]][confLogScale$keyParScaledYA + 1],
                               nrow = nrow(confLogScale$keyParScaledYA)) %>% t()
      fitSim$plsd[[h]] <- matrix(fitSim$plsd[[h]][confLogScale$keyParScaledYA + 1],
                                 nrow = nrow(confLogScale$keyParScaledYA)) %>% t()
    }
      rownames(fitSim$pl[[h]]) <- paste0("fit.", 1:nrow(fitSim$pl[[h]]))
      sdLog <- fitSim$plsd[[h]]
      rownames(sdLog) <- paste0("sdLog.", 1:nrow(fitSim$pl[[h]]))
      h_tru <- which(names(simOut$trueParams$pl) == names(fitSim$pl[h]))
      
      errRe <-
        rbind(errRe,
              fitSim$pl[[h]] %>%
              t() %>%
              exp %>%
              as.data.frame() %>%
              cbind(sdLog %>%
                      t() %>%
                      as.data.frame()) %>%
              cbind(data.frame(simOut$trueParams$pl[[h_tru]] %>% t %>% exp)) %>%
              dplyr::mutate(year = year) %>%
              tidyr::gather(variable, N, -year) %>%
              dplyr::mutate(variable = gsub(x = variable, 
                                            pattern = "X", 
                                            replacement = "tru.")) %>%
              tidyr::separate(variable, c("source", "age")) %>%
              tidyr::spread(source, N) %>%
              dplyr::mutate(age = as.numeric(age),
                            error = (fit - tru),
                            abs_error = abs(error),
                            error_pc = 100 * (fit - tru) / tru,
                            abs_error_pc = abs(error_pc),
                            decile = ceiling(10 * pnorm(q    = log(tru), 
                                                        mean = log(fit), 
                                                        sd   = sdLog)),
                            variable = varName))
    #}
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
                              abs_error_pc = abs(error_pc),
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
calcCSSBError <- function(fitSim, simOut) {
  
  #errCSSB <- data.frame()
  #for (i in 1:length(fitSim)) {
    # Fit catch and SSB
    df_CSSBfit <- 
      data.frame(variable = names(fitSim$sdrep$value),
                 value = fitSim$sdrep$value,
                 sd    = fitSim$sdrep$sd) %>%
      dplyr::filter(variable %in% c("logCatch", "logssb")) %>%
      dplyr::rename(sdLog = sd) %>%
      dplyr::mutate(fit = exp(value),
                    year = rep(fitSim$data$years, 2),
                    age = "total",
                    variable = as.character(variable),
                    variable = ifelse(variable == "logCatch", "catch", variable),
                    variable = ifelse(variable == "logssb", "ssb", variable)) %>%
      dplyr::select(year, age, fit, sdLog, variable)
    
    # True catch and SSB
    df_CSSBtru <- 
      simOut$Ctru_mt %>%
      t() %>%
      as.data.frame() %>%
      dplyr::mutate(total = rowSums(.),
                    year = fitSim$data$years) %>%
      tidyr::gather(age, Catch_mt, -year) %>%
      dplyr::rename(tru = Catch_mt) %>%
      dplyr::filter(age == "total") %>%
      dplyr::mutate(variable = "catch") %>%
      rbind(data.frame(tru = simOut$SSB,
                       year = names(simOut$SSB),
                       age = "total",
                       variable = "ssb")) %>%
      dplyr::mutate(year = as.integer(year))
     
    
    errCSSB <-
            df_CSSBfit %>%
              dplyr::left_join(df_CSSBtru) %>%
              dplyr::mutate(error = fit - tru,
                            abs_error = abs(error),
                            error_pc = 100 * (fit - tru) / tru,
                            abs_error_pc = abs(error_pc),
                            decile = ceiling(10 * pnorm(q    = log(tru),
                                                        mean = log(fit),
                                                        sd   = sdLog)))
  #}  
  
  
  return(errCSSB)
}

## Plot timeseries error ##################################
plotTsError <- function(container) {
  
  colors2use <- RColorBrewer::brewer.pal(3, "Dark2")
    
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
  # p <-
  #   ggplot(err %>%
  #            dplyr::select(-sdLog, -decile, -error_pc) %>%
  #            dplyr::filter(variable == "F",
  #                          age == 2,
  #                          replicate %in% 1:10) %>%
  #            tidyr::gather(source, value, 
  #                          -year, -age, -replicate, -variable),
  #          aes(x = year, y = value, color = source)) +
  #   geom_line() +
  #   facet_wrap(~replicate) + 
  #   ggtitle("Age-2 F estimation error examples")
  # print(p)
  
  # Plot decile coverage for N and F
  # p <- 
  #   ggplot(err %>% dplyr::filter(variable %in% c("N", "F")),
  #          aes(x = decile)) +
  #   geom_histogram(binwidth=1, colour="white") +
  #   geom_hline(yintercept = noYears * nRepAccept / 10, 
  #              color = "dark grey") +
  #   theme_bw() +
  #   facet_grid(variable~age) +
  #   ggtitle("Confidence interval coverage for each age")
  # print(p)
  
  # Plot decile coverage for Catch
  # if(any(err$variable == "catch")) {
  #   p <- 
  #     ggplot(err %>% dplyr::filter(variable == "catch"),
  #            aes(x = decile)) +
  #     geom_histogram(binwidth=1, colour="white") +
  #     geom_hline(yintercept = noYears * nRepAccept / 10, 
  #                color = "dark grey") +
  #     theme_bw() +
  #     facet_wrap(~age) +
  #     ggtitle("Confidence interval coverage for catch")
  #   print(p)
  # }
  
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
    
  err <- data.frame()
  for (i in 1:nrow(container)) {
    err <-
      rbind(err, 
            {rbind(data.frame(container$err_random[[i]], model = "random walk"),
                   data.frame(container$err_fixed[[i]],  model = "fixed"),
                   data.frame(container$err_none[[i]],   model = "no misreporting")) %>%
                dplyr::mutate(replicate = container$replicate[i],
                              scenario  = as.factor(paste(container$scenario[i], "scenario")),
                              scenario  = factor(scenario, levels = c("no misreporting scenario",
                                                                      "fixed scenario",
                                                                      "random walk scenario",
                                                                      "uniform random scenario")),
                              model     = factor(model, levels = c("no misreporting", 
                                                                   "fixed",
                                                                   "random walk")))})
  }
  
  # err2plot <-
  #   err %>%  
  #   dplyr::filter(variable %in% c("N")) %>%
  #   dplyr::select(model, scenario, year, age, replicate, tru, fit) %>%
  #   dplyr::rename(N_tru = tru,
  #                 N_fit = fit) %>%
  #   dplyr::left_join({err %>%  
  #       dplyr::filter(variable %in% c("F")) %>%
  #       dplyr::select(model, scenario, year, age, replicate, error_pc, abs_error) %>%
  #       dplyr::rename(F_error_pc = error_pc,
  #                     F_abs_error = abs_error)}) %>%
  #   dplyr::left_join({err %>%  
  #       dplyr::filter(variable %in% c("N")) %>%
  #       dplyr::select(model, scenario, year, age, replicate, error_pc, abs_error) %>%
  #       dplyr::rename(N_error_pc = error_pc,
  #                     N_abs_error = abs_error)}) %>%
  #   dplyr::left_join({err %>%  
  #       dplyr::filter(variable %in% c("Scale")) %>%
  #       dplyr::select(model, scenario, year, age, replicate, error_pc, abs_error) %>%
  #       dplyr::rename(Scale_error_pc = error_pc,
  #                     Scale_abs_error = abs_error)})
  
  
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
  # err2plot_percentile <-
  #   err2plot %>%
  #   dplyr::group_by(model, scenario, age) %>%
  #   dplyr::mutate(N_percentile = ntile(N_tru, 100)) %>%
  #   dplyr::group_by(model, scenario, age, N_percentile) %>%
  #   dplyr::summarise(F_error_pc_mean = mean(F_error_pc),
  #                    N_error_pc_mean = mean(N_error_pc),
  #                    nObsF = length(F_error_pc),
  #                    nObsN = length(N_error_pc),
  #                    F_error_pc_hi   = F_error_pc_mean + 1.96 * sd(F_error_pc)/sqrt(nObsF),
  #                    F_error_pc_lo  = F_error_pc_mean - 1.96 * sd(F_error_pc)/sqrt(nObsF),
  #                    N_error_pc_hi  = N_error_pc_mean + 1.96 * sd(N_error_pc)/sqrt(nObsN),
  #                    N_error_pc_lo  = N_error_pc_mean - 1.96 * sd(N_error_pc)/sqrt(nObsN),
  #                    N_tru_value     = mean(N_tru),
  #                    N_fit_max       = quantile(N_fit, .95),
  #                    N_fit_min       = quantile(N_fit, .05))
  
  
  # N
  # p <-
  #   ggplot(err2plot_percentile,
  #          aes(x = N_percentile)) +
  #   geom_line(aes(y = N_error_pc_mean)) +
  #   geom_ribbon(aes(ymin = N_error_pc_lo, ymax = N_error_pc_hi), alpha = 0.3) +
  #   geom_hline(yintercept = 0) +
  #   facet_wrap(~age) +
  #   theme_bw() +
  #   xlab("Percentile of N") +
  #   ylab("Percent error of N estimate") +
  #   ggtitle("Mean percent error of N vs Percentile of N")
  # print(p)
  
  # F
  # p <-
  #   ggplot(err2plot_percentile,
  #          aes(x = N_percentile)) +
  #   geom_line(aes(y = F_error_pc_mean)) +
  #   geom_ribbon(aes(ymin = F_error_pc_lo, ymax = F_error_pc_hi), alpha = 0.3) +
  #   geom_hline(yintercept = 0) +
  #   facet_wrap(~age) +
  #   theme_bw() +
  #   xlab("Percentile of N") +
  #   ylab("Percent error of F estimate") +
  #   ggtitle("Mean percent error of F vs Percentile of N")
  # print(p)
  
  # Plot percentile of N vs value of N
  # p <-
  #   ggplot(err2plot_percentile %>% dplyr::filter(N_percentile < 10),
  #          aes(x = N_percentile)) +
  #   geom_ribbon(aes(ymin = N_fit_min, ymax = N_fit_max), fill = "blue", alpha = 0.3) +
  #   facet_wrap(~age, scales = "free_y") +
  #   xlab("N percentile") +
  #   ylab("Range of N fit (1000's)")
  # print(p)
  
  
  # Plot mean error in catch and ssb vs year
  err2plot_CSSB <-
    err %>%  
    dplyr::filter(variable %in% c("catch", "ssb", "F", "N")) %>%
    dplyr::select(model, scenario, year, variable, age, replicate, error_pc, abs_error_pc) %>%
    dplyr::group_by(model, scenario, year, variable) %>%
    dplyr::summarise(error_pc_mean = mean(error_pc, na.rm = T),
                     mape  = mean(abs_error_pc, na.rm = T),
                     nObs = length(abs_error_pc),
                     error_pc_hi   = error_pc_mean + 1.96 * sd(error_pc, na.rm = T)/sqrt(nObs),
                     error_pc_lo  = error_pc_mean - 1.96 * sd(error_pc, na.rm = T)/sqrt(nObs),
                     mape_hi   = mape + 1.96 * sd(abs_error_pc, na.rm = T)/sqrt(nObs),
                     mape_lo   = mape - 1.96 * sd(abs_error_pc, na.rm = T)/sqrt(nObs))
  p <-
    ggplot(err2plot_CSSB %>% 
             dplyr::filter(model %in% c("fixed", "random walk")),#,
                           #scenario %in% c("no misreporting scenario", "fixed scenario"),
                           #variable == "catch"
            #               ),
           aes(x = year, color = model, fill = model)) +
    geom_line(aes(y = error_pc_mean)) +
    geom_ribbon(aes(ymin = error_pc_lo, ymax = error_pc_hi), color = NA, alpha = 0.3) +
    geom_hline(yintercept = 0) +
    facet_grid(scenario~variable, scales = "free_y") +
    theme_bw() +
    xlab("Year") +
    ylab("Mean percent error") +
    scale_color_manual(values = colors2use[2:3]) +
    scale_fill_manual(values = colors2use[2:3]) +
    ggtitle("Estimation error")  +
    theme(axis.title   = element_text(size = 14),
          plot.title   = element_text(size = 16),
          strip.text   = element_text(size = 9))
  
  print(p)
  
  # Plot error of Scale estimates
  if (any(err$variable == "Scale")) {
    
    # Plot average error for each scenario
    err2plot_Scale_overall <-
      err %>%  
      dplyr::filter(variable %in% c("Scale")) %>%
      dplyr::select(model, scenario, year, variable, age, replicate, 
                    error_pc, abs_error_pc, abs_error) %>%
      dplyr::group_by(model, scenario, year, variable) %>%
      dplyr::summarise(error_pc_mean = mean(error_pc, na.rm = T),
                       mape  = mean(abs_error_pc, na.rm = T),
                       mae   = mean(abs_error),
                       nObs  = length(abs_error_pc),
                       error_pc_hi = error_pc_mean + 1.96 * sd(error_pc, na.rm = T)/sqrt(nObs),
                       error_pc_lo = error_pc_mean - 1.96 * sd(error_pc, na.rm = T)/sqrt(nObs),
                       mape_hi   = mape + 1.96 * sd(abs_error_pc, na.rm = T)/sqrt(nObs),
                       mape_lo   = mape - 1.96 * sd(abs_error_pc, na.rm = T)/sqrt(nObs),
                       mae_hi   = mae + 1.96 * sd(mae, na.rm = T)/sqrt(nObs),
                       mae_lo   = mae - 1.96 * sd(mae, na.rm = T)/sqrt(nObs))
    
    p <-
      ggplot(err2plot_Scale_overall,
             aes(x = year, color = model, fill = model)) +
      geom_line(aes(y = mape)) +
      geom_ribbon(aes(ymin = mape_lo, ymax = mape_hi), color = NA, alpha = 0.3) +
      geom_hline(yintercept = 0) +
      facet_grid(~scenario) +
      theme_bw() +
      xlab("Year") +
      ylab("Mean absolute percent error") +
      scale_color_manual(values = colors2use[2:3]) +
      scale_fill_manual(values = colors2use[2:3]) +
      ggtitle("Scale estimation error") 
    
    print(p)
    
    
    # Plot example fit vs true scale
    err2plot_Scale <-
      err %>%
      dplyr::mutate(fit_975 = exp(log(fit) + 1.96 * sdLog),
                    fit_025 = exp(log(fit) - 1.96 * sdLog),
                    replicate = paste("replicate", replicate),
                    age = paste("age-", age)) %>%
      dplyr::filter(variable %in% c("Scale"))
    scenarios2plot <- unique(err2plot_Scale$scenario)
    for (i in 1:length(scenarios2plot)) {
      p <-
        ggplot(err2plot_Scale %>% 
                 dplyr::filter(scenario == scenarios2plot[i]) %>%
                 dplyr::filter(replicate %in% unique(replicate)[1:4]), 
               aes(x = year)) +
        geom_line(aes(y = fit, color = model)) +
        geom_line(aes(y = tru), color = "black") +
        geom_ribbon(aes(ymin = fit_975, ymax = fit_025, fill = model), 
                    alpha = 0.3, color = NA) +
        facet_grid(replicate ~ age, scales = "free_y") +
        theme_bw() +
        xlab("Year") +
        ylab("Estimated and True scale parameter value") +
        scale_color_manual(values = colors2use[2:3]) +
        scale_fill_manual(values = colors2use[2:3]) +
        ggtitle(paste0("Scale parameter error (", 
                       scenarios2plot[i], ")"))
      print(p)
    }
    
    # Example scale time series
      p <-
        ggplot(err2plot_Scale %>%
                 dplyr::filter(replicate %in% unique(replicate)[1]), 
               aes(x = year)) +
        geom_line(aes(y = fit, color = model)) +
        geom_hline(yintercept = 1, color = "dark grey") +
        geom_line(aes(y = tru), color = "black") +
        geom_ribbon(aes(ymin = fit_975, ymax = fit_025, fill = model), 
                    alpha = 0.3, color = NA) +
        facet_grid(scenario~age, scales = "free_y") +
        theme_bw() +
        xlab("Year") +
        ylab("Scale parameter value") +
        #scale_y_continuous(breaks=seq(1,11,2)) +
        scale_x_continuous(breaks = seq(2005, 2014, 4)) +
        scale_color_manual(values = colors2use[2:3]) +
        scale_fill_manual(values = colors2use[2:3]) +
        ggtitle("An example of each misreporting scenario") +
        theme(axis.title   = element_text(size = 14),
              plot.title   = element_text(size = 16),
              strip.text   = element_text(size = 12))
      print(p)

  }
  
  

  
  
  # Plot average true catch vs average estimated catch
  # err2plot_C_Av <-
  #   err %>%  
  #   dplyr::filter(variable %in% c("catch")) %>%
  #   dplyr::select(year, age, replicate, tru, fit) %>%
  #   dplyr::group_by(year) %>%
  #   dplyr::summarise(mean_tru = mean(tru),
  #                    mean_fit = mean(fit))
  
  # p <-
  #   ggplot(err2plot_C_Av, aes(x = year)) +
  #   geom_line(aes(y = mean_tru), col = "black") +
  #   geom_line(aes(y = mean_fit), col = "red") +
  #   #geom_ribbon(aes(ymin = error_pc_lo, ymax = error_pc_hi), alpha = 0.3) +
  #   theme_bw() +
  #   xlab("Year") +
  #   ylab("Mean catch (MT)")
  # print(p)
  
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
  
  
  # Plot pairs of F error vs Scale error for uniform scenario
  errPairs <- 
    err %>%
    dplyr::filter(year %in% container$setupMod_fixed[[1]]$conf$keyScaledYears) %>%
    dplyr::group_by(model, scenario, variable) %>%
    dplyr::summarise(mape = mean(abs_error_pc, na.rm = T),
                     mape_se = sd(abs_error_pc, na.rm = T)/
                               sqrt(length(!is.na(abs_error_pc)))) 
  
  ggplot(errPairs %>%
           dplyr::select(-mape_se) %>%
           tidyr::spread(variable, mape) %>%
           tidyr::gather(variable, mape, -model, -scenario, -Scale) %>%
           dplyr::filter(model != "no misreporting",
                         scenario == "uniform random scenario"),
         aes(x = Scale, y = mape, color = model, group = variable, shape = variable)) +
    geom_point(size = 4) +
    geom_line(color = "black") +
    scale_color_manual(values = colors2use[2:3]) +
    xlab("Estimation error of Scale parameter (MAPE)") +
    ylab("Estimation error of variable (MAPE)") +
    ggtitle("Uniform random scenario")
  
}

## Plot mean error over all replicates ####################
plotTsMeanError <- function(err, nRepAccept) {
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



## Customized sam.fit() used to fit random effect logScale #######
sam.fit_cp <-
function (data, conf, parameters, newtonsteps = 3, rm.unidentified = FALSE, 
          run = TRUE, 
          lower = stockassessment2:::getLowerBounds(parameters), 
          upper = stockassessment2:::getUpperBounds(parameters), 
          sim.condRE = TRUE, ignore.parm.uncertainty = FALSE, rel.tol = 1e-10,
          ...) 
{
  definit <- stockassessment2::defpar(data, conf)
  
  # if (!identical(parameters, relist(unlist(parameters), skeleton = definit))) {
  #   warning("Initial values are not consistent, so running with default init values from defpar()")
  #   parameters <- definit
  # }
  data <- stockassessment2:::clean.void.catches(data, conf)
  tmball <- c(data, conf, simFlag = as.numeric(sim.condRE))
  if (is.null(tmball$resFlag)) {
    tmball$resFlag <- 0
  }
  nmissing <- sum(is.na(data$logobs))
  parameters$missing <- numeric(nmissing)
  ran <- c("logN", "logF", "logScale", "missing")
  
  
  obj <- TMB::MakeADFun(tmball, parameters, random = ran, DLL = "stockassessment2", 
                   ...)
  obj$env$tracepar <- TRUE
  if (rm.unidentified) {
    gr <- obj$gr()
    safemap <- obj$env$parList(gr)
    safemap <- safemap[!names(safemap) %in% ran]
    safemap <- lapply(safemap, function(x) factor(ifelse(abs(x) > 
                                                           1e-15, 1:length(x), NA)))
    ddd <- list(...)
    if (!is.null(ddd$map)) {
      safemap <- c(ddd$map, safemap)
      ddd$map <- safemap
      ddd$data <- tmball
      ddd$parameters <- parameters
      ddd$random <- ran
      ddd$DLL <- "stockassessment2"
      obj <- do.call(TMB::MakeADFun, ddd)
    }
    else {
      obj <- TMB::MakeADFun(tmball, parameters, random = ran, 
                            map = safemap, DLL = "stockassessment2", ...)
    }
  }
  lower2 <- rep(-Inf, length(obj$par))
  upper2 <- rep(Inf, length(obj$par))
  for (nn in names(lower)) lower2[names(obj$par) == nn] = lower[[nn]]
  for (nn in names(upper)) upper2[names(obj$par) == nn] = upper[[nn]]

  if (!run) 
    return(list(sdrep = NA, pl = parameters, plsd = NA, 
                data = data, conf = conf, opt = NA, obj = obj))
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 1, 
                                                        eval.max = 2000,
                                                        iter.max = 1000, 
                                                        rel.tol = rel.tol), 
                lower = lower2, upper = upper2)
  for (i in seq_len(newtonsteps)) {
    g <- as.numeric(obj$gr(opt$par))
    h <- optimHess(opt$par, obj$fn, obj$gr)
    opt$par <- opt$par - solve(h, g)
    opt$objective <- obj$fn(opt$par)
  }
  rep <- obj$report()
  sdrep <- TMB::sdreport(obj, opt$par, ignore.parm.uncertainty = ignore.parm.uncertainty)
  idx <- c(which(names(sdrep$value) == "lastLogN"), which(names(sdrep$value) == 
                                                            "lastLogF"))
  sdrep$estY <- sdrep$value[idx]
  sdrep$covY <- sdrep$cov[idx, idx]
  idx <- c(which(names(sdrep$value) == "beforeLastLogN"), 
           which(names(sdrep$value) == "beforeLastLogF"))
  sdrep$estYm1 <- sdrep$value[idx]
  sdrep$covYm1 <- sdrep$cov[idx, idx]
  pl <- as.list(sdrep, "Est")
  plsd <- as.list(sdrep, "Std")
  sdrep$cov <- NULL
  ret <- list(sdrep = sdrep, pl = pl, plsd = plsd, data = data, 
              conf = conf, opt = opt, obj = obj, rep = rep, low = lower, 
              hig = upper)
  attr(ret, "RemoteSha") <- substr(packageDescription("stockassessment2")$RemoteSha, 
                                   1, 12)
  attr(ret, "Version") <- packageDescription("stockassessment2")$Version
  class(ret) <- "sam"
  return(ret)
}

#### Altered version of retro function ####
retro_cp <- function (fit, year = NULL, ncores = detectCores(), ...) 
{
  data <- fit$data
  y <- fit$data$aux[, "year"]
  f <- fit$data$aux[, "fleet"]
  suf <- sort(unique(f))
  maxy <- sapply(suf, function(ff) max(y[f == ff]))
  if (length(year) == 1) {
    mat <- sapply(suf, function(ff) {
      my <- maxy[ff]
      my:(my - year + 1)
    })
    if (year == 1) 
      mat <- matrix(mat, nrow = 1)
  }
  if (is.vector(year) & length(year) > 1) {
    mat <- sapply(suf, function(ff) year)
  }
  if (is.matrix(year)) {
    mat <- year
  }
  if (nrow(mat) > length(unique(y))) 
    stop("The number of retro runs exceeds number of years")
  if (ncol(mat) != length(suf)) 
    stop("Number of retro fleets does not match")
  setup <- lapply(1:nrow(mat), function(i) do.call(rbind, 
                                                   lapply(suf, function(ff) if (mat[i, ff] <= maxy[ff]) {
                                                     cbind(mat[i, ff]:maxy[ff], ff)
                                                   })))
  if (ncores > 1) {
    cl <- makeCluster(ncores)
    on.exit(stopCluster(cl))
    clusterExport(cl, varlist = "fit", envir = environment())
    clusterEvalQ(cl, {source("1-functions.R")}) 
    runs <- parLapply(cl, setup, function(s) runwithout_cp(fit, year = s[, 1], fleet = s[, 2], ...))
  }
  else {
    runs <- lapply(setup, function(s) runwithout_cp(fit, year = s[, 1], fleet = s[, 2], ...))
  }
  converg <- unlist(lapply(runs, function(x) x$opt$conv))
  if (any(converg != 0)) 
    warning(paste0("retro run(s) ", paste0(which(converg != 
                                                   0), collapse = ","), " did not converge."))
  attr(runs, "fit") <- fit
  class(runs) <- "samset"
  runs
}

#### Another customized function to run retros ####
runwithout_cp <- function (fit, year = NULL, fleet = NULL, map = fit$obj$env$map, 
          ...) 
{
  data <- reduce_cp(fit$data, year = year, fleet = fleet, conf = fit$conf)
  conf <- attr(data, "conf")
  par <- stockassessment2::defpar(data, conf)
  par[!names(par) %in% c("logN", "logF", "logScale")] <- fit$pl[!names(fit$pl) %in%
                                                      c("missing", "logN", "logF", "logScale")]
  ret <- sam.fit_cp(data, conf, par, rm.unidentified = TRUE,
                    map = map, lower = fit$low, upper = fit$hig, ...)
  return(ret)
}


### And another need to run retros ####
reduce_cp <- function (data, year = NULL, fleet = NULL, age = NULL, conf = NULL) 
{
  nam <- c("year", "fleet", "age")[c(length(year) > 0, length(fleet) > 
                                       0, length(age) > 0)]
  if ((length(year) == 0) & (length(fleet) == 0) & (length(age) == 
                                                    0)) {
    idx <- rep(TRUE, nrow(data$aux))
  }
  else {
    idx <- !do.call(paste, as.data.frame(data$aux[, nam, 
                                                  drop = FALSE])) %in% do.call(paste, as.data.frame(cbind(year = year, 
                                                                                                          fleet = fleet, age = age)))
  }
  data$aux <- data$aux[idx, ]
  data$logobs <- data$logobs[idx]
  data$weight <- data$weight[idx]
  suf <- sort(unique(data$aux[, "fleet"]))
  data$noFleets <- length(suf)
  data$fleetTypes <- data$fleetTypes[suf]
  data$sampleTimes <- data$sampleTimes[suf]
  data$years <- min(as.numeric(data$aux[, "year"])):max(as.numeric(data$aux[, 
                                                                            "year"]))
  ages <- min(as.numeric(data$aux[, "age"])):max(as.numeric(data$aux[, 
                                                                     "age"]))
  data$noYears <- length(data$years)
  mmfun <- function(f, y, ff) {
    idx <- which(data$aux[, "year"] == y & data$aux[, "fleet"] == 
                   f)
    ifelse(length(idx) == 0, NA, ff(idx) - 1)
  }
  data$idx1 <- outer(suf, data$years, Vectorize(mmfun, c("f", 
                                                         "y")), ff = min)
  data$idx2 <- outer(suf, data$years, Vectorize(mmfun, c("f", 
                                                         "y")), ff = max)
  data$nobs <- length(data$logobs[idx])
  data$propMat <- data$propMat[rownames(data$propMat) %in% 
                                 data$years, colnames(data$propMat) %in% ages]
  data$stockMeanWeight <- data$stockMeanWeight[rownames(data$stockMeanWeight) %in% 
                                                 data$years, colnames(data$stockMeanWeight) %in% ages]
  data$natMor <- data$natMor[rownames(data$natMor) %in% data$years, 
                             colnames(data$natMor) %in% ages]
  data$propF <- data$propF[rownames(data$propF) %in% data$years, 
                           colnames(data$propF) %in% ages]
  data$propM <- data$propM[rownames(data$propM) %in% data$years, 
                           colnames(data$propM) %in% ages]
  data$landFrac <- data$landFrac[rownames(data$landFrac) %in% 
                                   data$years, colnames(data$landFrac) %in% ages]
  data$catchMeanWeight <- data$catchMeanWeight[rownames(data$catchMeanWeight) %in% 
                                                 data$years, colnames(data$catchMeanWeight) %in% ages]
  data$disMeanWeight <- data$disMeanWeight[rownames(data$disMeanWeight) %in% 
                                             data$years, colnames(data$disMeanWeight) %in% ages]
  data$landMeanWeight <- data$landMeanWeight[rownames(data$landMeanWeight) %in% 
                                               data$years, colnames(data$landMeanWeight) %in% ages]
  data$aux[, "fleet"] <- match(data$aux[, "fleet"], suf)
  data$minAgePerFleet <- tapply(as.integer(data$aux[, "age"]), 
                                INDEX = data$aux[, "fleet"], FUN = min)
  data$maxAgePerFleet <- tapply(as.integer(data$aux[, "age"]), 
                                INDEX = data$aux[, "fleet"], FUN = max)
  attr(data, "fleetNames") <- attr(data, "fleetNames")[suf]
  if (!missing(conf)) {
    .reidx <- function(x) {
      if (any(x >= (-0.5))) {
        xx <- x[x >= (-0.5)]
        x[x >= (-0.5)] <- match(xx, sort(unique(xx))) - 
          1
      }
      x
    }
    conf$keyLogFsta <- .reidx(conf$keyLogFsta[suf, , drop = FALSE])
    conf$keyLogFpar <- .reidx(conf$keyLogFpar[suf, , drop = FALSE])
    conf$keyQpow <- .reidx(conf$keyQpow[suf, , drop = FALSE])
    conf$keyVarF <- .reidx(conf$keyVarF[suf, , drop = FALSE])
    conf$keyVarObs <- .reidx(conf$keyVarObs[suf, , drop = FALSE])
    conf$obsCorStruct <- conf$obsCorStruct[suf]
    conf$keyCorObs <- conf$keyCorObs[suf, , drop = FALSE]
    yidx <- conf$keyScaledYears %in% data$aux[data$aux[, 
                                                       "fleet"] == 1, "year"]
    if (length(conf$keyScaledYears) > 0) {
      conf$noScaledYears <- sum(yidx)
      conf$keyScaledYears <- as.vector(conf$keyScaledYears)[yidx]
      # conf$keyParScaledYA <- .reidx(conf$keyParScaledYA[yidx, 
      #                                                   , drop = FALSE])
    }
    conf$keyBiomassTreat <- conf$keyBiomassTreat[suf]
    conf$obsLikelihoodFlag <- conf$obsLikelihoodFlag[suf]
    attr(data, "conf") <- conf
  }
  data
}


calcCatchAdviceError <- function(fit, simOut, confLogScale){

  idxno <- which(fit$data$years==max(fit$data$years))
  aveYears <- 5
  ave.sw<-colMeans(fit$data$stockMeanWeight[(idxno-aveYears+1):idxno,,drop=FALSE])
  ave.cw<-colMeans(fit$data$catchMeanWeight[(idxno-aveYears+1):(idxno-1),,drop=FALSE])
  ave.pm<-colMeans(fit$data$propMat[(idxno-aveYears+1):idxno,,drop=FALSE])
  ave.nm<-colMeans(fit$data$natMor[(idxno-aveYears+1):idxno,,drop=FALSE])
  ave.sel_est <- colMeans(faytable(fit)[(idxno-aveYears+1):idxno,,drop=FALSE]/
                            apply(faytable(fit)[(idxno-aveYears+1):idxno,,drop=FALSE],1,max))
  ave.sel_tru <- colMeans(t(exp(simOut$trueParams$pl$logF))[(idxno-aveYears+1):idxno,,drop=FALSE]/
                            apply(t(exp(simOut$trueParams$pl$logF))[(idxno-aveYears+1):idxno,,drop=FALSE],1,max) )
  
  naa_est <- vector("numeric", length(fitReal$conf$minAge:fitReal$conf$maxAge))
  naa_est[1] <- 1
  naa_tru <- vector("numeric", length(fitReal$conf$minAge:fitReal$conf$maxAge))
  naa_tru[1] <- 1
  full_f <- seq(0, 10, 0.01)
  maxFage_est <- which.max(ave.sel_est)
  maxFage_tru <- which.max(ave.sel_tru)
  ssbr_est <- vector("numeric", length = length(full_f))
  ssbr_tru <- vector("numeric", length = length(full_f))
  for (h in 1:length(full_f)) {
    for (i in 2:fit$conf$maxAge)  {
      naa_est[i] <- naa_est[i-1] * exp(-(full_f[h] * ave.sel_est[i-1] + ave.nm[i-1]))
      naa_tru[i] <- naa_tru[i-1] * exp(-(full_f[h] * ave.sel_tru[i-1] + ave.nm[i-1]))
      if (i == fit$conf$maxAge & fit$conf$maxAgePlusGroup) { # if plus group
        naa_est[i] <- naa_est[i] * 1/(1-exp(-(full_f[h] * ave.sel_est[maxFage_est] + ave.nm[maxFage_est])))
        naa_tru[i] <- naa_tru[i] * 1/(1-exp(-(full_f[h] * ave.sel_tru[maxFage_tru] + ave.nm[maxFage_tru])))
      }
    }
    ssbr_est[h] <- sum(naa_est * ave.sw * ave.pm)
    ssbr_tru[h] <- sum(naa_tru * ave.sw * ave.pm)
  }
  
  # plot(full_f[1:400], ssbr_tru[1:400]/ssbr_tru[1], type = "l", col ="black", 
  #      xlab = "F", ylab = expression(SSB[F] / SSB[F == 0]), cex.lab = 1.1)
  # lines(full_f[1:400], ssbr_est[1:400]/ssbr_est[1], type = "l", col ="blue")
  # legend(2.7, 0.8, legend=c("True", "Estimated"),
  #        col=c("black", "blue"), lty = 1, cex = 1.2)
  f40ind_tru <- which.min(abs(ssbr_tru/ssbr_tru[1] - .4))
  f40ind_est <- which.min(abs(ssbr_est/ssbr_est[1] - .4))
  f40_tru <- full_f[f40ind_tru] * ave.sel_tru
  f40_est <- full_f[f40ind_est] * ave.sel_est
  z40_tru <- f40_tru + ave.nm
  z40_est <- f40_est + ave.nm
  
  catch40_mt_tru <- exp(simOut$trueParams$pl$logN[,idxno]) * 
                      (1 - exp(-z40_tru)) * (f40_tru/z40_tru) * ave.cw
  catch40_mt_est <- exp(fit$pl$logN[,idxno]) * 
                      (1 - exp(-z40_est)) * (f40_est/z40_est) * ave.cw
  

  return(data.frame(model = confLogScale$logScaleType,
                    variable = "catch_at_F40_mt",
                    tru = sum(catch40_mt_tru),
                    fit = sum(catch40_mt_est)) %>% #,
                    #age = substr(names(catch40_mt_tru), 5,6), 
                    #row.names = NULL) %>%
           dplyr::mutate(error = (fit - tru),
                         abs_error = abs(error),
                         error_pc = 100 * (fit - tru) / tru,
                         abs_error_pc = abs(error_pc)))

}
