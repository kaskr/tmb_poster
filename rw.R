library(MASS)

simdata <- function(){
  set.seed(1);
  local({
    rho=0.9
    sds=seq(0.5,2,length=stateDim)
    sdObs=rep(1,stateDim);
    corrMat=matrix(0.0,stateDim,stateDim)
    for(i in 1:stateDim){
      for(j in 1:stateDim){
        corrMat[i,j] = rho^abs(i-j)
      }
    }
    Sigma=corrMat*(sds %o% sds)
    d=matrix(NA,timeSteps,stateDim)
    obs=d;
    ##init state
    d[1,] = rnorm(stateDim);
    i=1;
    obs[i,] = d[i,] + rnorm(stateDim,rep(0,stateDim),sdObs)
    for(i in 2:timeSteps){
      d[i,] = d[i-1,] + mvrnorm(1,rep(0,stateDim),Sigma=Sigma)
      obs[i,] = d[i,] + rnorm(stateDim,rep(0,stateDim),sdObs)
    }
    matplot(d,type="l",las=1,lty=1,xlab="Time",ylab="Simulated states and measurements")
    matpoints(obs,pch=1);
  },.GlobalEnv)
}
stateDim=3
timeSteps=100
simdata()
library(TMB)
compile("rw.cpp")
dyn.load("rw.so")
data <- list(obs=t(obs))
parameters <- list(
  X=data$obs*0,
  transf_rho=0.1,
  logsds=sds*0,
  logsdObs=sdObs*0
  )
obj <- MakeADFun(data,parameters,random="X",DLL="rw")
newtonOption(obj, smartsearch=FALSE)

obj$fn()
obj$gr()
system.time(opt <- do.call("optim",obj))
pl <- obj$env$parList() ## <-- List of predicted random effects
#matpoints(t(pl$X),type="l",col="blue",lwd=2)
#summary(as.mle(opt))
sdreport(obj)
