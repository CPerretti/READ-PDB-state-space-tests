# expand sd_investigate approach to multiple realizations and range of process and observation errors

library(TMB)
library(ggplot2)

# compile cpp code and load dll
compile("srw.cpp")
dyn.load(dynlib("srw"))

# simulation settings
nloops <- 100
n.obs <- 35
true.process.error.vals <- c(0.1,0.3,0.5)
true.obs.error.vals <- c(0.1,0.3,0.5)
ncases <- length(true.process.error.vals)*length(true.obs.error.vals)
res <- matrix(NA, nrow=ncases*nloops, ncol=10)

# loop through the errors
icount <- 0
case <- 0

for (ii in 1:length(true.process.error.vals)){
  true.process.error <- true.process.error.vals[ii]
  for (jj in 1:length(true.obs.error.vals)){
    true.obs.error <- true.obs.error.vals[jj]
    case <- case + 1
    for (iloop in 1:nloops){
      icount <- icount + 1
      true.population <- rep(0,n.obs)
      true.population[1] <- rnorm(1,0,true.process.error)
      for (i in 2:n.obs){
        true.population[i] <- true.population[i-1] + rnorm(1,0,true.process.error) 
      }
      observed <- round(true.population + rnorm(n.obs,0,true.obs.error), 6)
      
      dat <- list(
        observed=observed
      )
      parameters <- list(
        population=rep(0,n.obs),
        log_process_error=0,
        log_obs_error=0
      )
      
      obj <- MakeADFun(dat,parameters,DLL="srw", random=c("population"), silent=TRUE)
      opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(iter.max=1000,eval.max=1000))
      rep <- sdreport(obj)
      srep <- summary(sdreport(obj))
      est.prc.err <- as.vector(round(srep[rownames(srep) == "process_error",],5))
      est.obs.err <- as.vector(round(srep[rownames(srep) == "obs_error",],5))
      popest <- srep[rownames(srep) == "population", 1]
      derived.prc.err <- sd(popest[2:n.obs] - popest[1:(n.obs - 1)])
      derived.obs.err <- sd(popest - observed)
      res[icount,] <- c(case,true.process.error,true.obs.error,iloop,est.prc.err,est.obs.err,derived.prc.err,derived.obs.err)
    }
  }
}

test.in <- function(vec){
  within <- 1
  if (is.na(sum(vec[5:8]))) within <- 0
  else{
    if (vec[2] < (vec[5] - 2 * vec[6])) within <- 0
    if (vec[2] > (vec[5] + 2 * vec[6])) within <- 0
    if (vec[3] < (vec[7] - 2 * vec[8])) within <- 0
    if (vec[3] > (vec[7] + 2 * vec[8])) within <- 0
  }
  return(within)
}

test.fail <- function(vec){
  fail <- 0
  if (vec[5] <= 0.001) fail <- 1
  if (vec[7] <= 0.001) fail <- 1
  return(fail)
}

pdf("srw_sim.pdf", onefile = TRUE)

my.col <- rainbow(ncases)
my.x <- range(res[,5], na.rm=T)
my.y <- range(res[,7], na.rm=T)
for (icase in 1:ncases){
  cres <- res[res[,1] == icase,]
  plot(cres[,5],cres[,7],col=my.col[icase],xlab="Process Error",ylab="Observation Error",xlim=my.x,ylim=my.y)
  abline(v=cres[1,2],lty=2)
  abline(h=cres[1,3],lty=3)
  prop.in <- 100 * sum(apply(cres, 1, FUN=test.in)) / nloops
  prop.fail <- 100 * sum(apply(cres, 1, FUN=test.fail)) / nloops
  title(main=paste("Case ",icase," (",prop.in,"% true within CI, ",prop.fail,"% fail)", sep=""),outer=F)
}

# now compare true, estimated, and derived process and observation errors
resdf <- data.frame(case = as.factor(res[,1]),
                    true.process.error = res[,2],
                    true.obs.error = res[,3],
                    iloop = res[,4],
                    est.prc.err = res[,5],
                    est.prc.err.sd = res[,6],
                    est.obs.err = res[,7],
                    est.obs.err.sd = res[,8],
                    derived.prc.err = res[,9],
                    derived.obs.err = res[,10])

g1 <- ggplot(resdf, aes(x=est.prc.err, y=derived.prc.err, color=case)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  ggtitle("Process Error") +
  theme_bw()
print(g1)

g2 <- ggplot(resdf, aes(x=est.obs.err, y=derived.obs.err, color=case)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  ggtitle("Observation Error") +
  theme_bw()
print(g2)

g3 <- ggplot(resdf, aes(x=true.process.error, y=est.prc.err, color=case)) +
  geom_boxplot() +
  theme_bw()
print(g3)

g4 <- ggplot(resdf, aes(x=true.process.error, y=derived.prc.err, color=case)) +
  geom_boxplot() +
  theme_bw()
print(g4)

g5 <- ggplot(resdf, aes(x=true.obs.error, y=est.obs.err, color=case)) +
  geom_boxplot() +
  theme_bw()
print(g5)

g6 <- ggplot(resdf, aes(x=true.obs.error, y=derived.obs.err, color=case)) +
  geom_boxplot() +
  theme_bw()
print(g6)

dev.off() # close pdf
