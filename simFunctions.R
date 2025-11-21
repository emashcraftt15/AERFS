###############################
# CGRFS simulation functions
#
# Subodh Selukar
# 2025-11-21
###############################


### data generation functions
require(mvtnorm)
#install.packages('mvtnorm')
#library(mvtnorm)
## helper function for inverse cdf method
# Inputs: (details under master function)
# 1. dist
# 2. paramVec=vector of parameters consistent with dist
#     NOTE: THERE IS NO ERROR HANDLING! 
# 3. vecIn=vector of standard normal distributed random variables to be transformed to dist(paramVec) distribution
makeMargVec <- function(dist,
                        paramVec,
                        vecIn){
  if (dist=="exp"){
    vecOut <- qexp(
      p=pnorm(vecIn),
      rate=paramVec[1]
    )
  } else if (dist=="wei"){
    vecOut <- qweibull(
      p=pnorm(vecIn),
      shape=paramVec[1],
      scale=paramVec[2]
    )
  } else if (dist=="gam"){
    vecOut <- qgamma(
      p=pnorm(vecIn),
      shape=paramVec[1],
      rate=paramVec[2]
    )
  } else if (dist=="logn"){
    vecOut <- qlnorm(
      p=pnorm(vecIn),
      meanlog=paramVec[1],
      sdlog=paramVec[2]
    )
  } # only these 4 options are currently supported
  
  return(vecOut)
}

## master function: uses "NORTA" method for converting correlated normal variables to specified correlated marginals
# Inputs: 
#   1. n=sample size ofeach dataset
#   2. numSims=number of simulations with this set of parameters
#   3. seedIn=simulation seed
#   4. corMat=matrix of correlations, needs to be 5x5 to agree with vector of 5 event times
#     NOTE: diagonal elements must be 1, i.e., marginal distributions need to be standardized to variance 1
#   5. params=list of length 6(=5 event times+1 censoring)
#     NOTE: these parameters are for the interval between events, NOT time from HCT
#     each list element should be a numeric vector with length equal to number of parameters for the respective distribution (uniform )
#     the list must be ordered as CG1,CG1.end,CG2,CG2.end,RD,C
#   6. dist=one of exp,wei,gam,logn for each event time distribution ** currently only supports one distribution family at a time
# Output: data frame with n*numSims rows and 7 columns (CG1,CG1.end,CG2,CG2.end,RD,C,simIdx)

makeSimDat <- function(n=100,numSims=1,
                   seedIn=2023,
                   corMat=NULL, # if null, assumed to all be independent
                   params,
                   dist){
  
  set.seed(seedIn)
  
  if (is.null(corMat)){
    corMat <- diag(nrow=5,ncol=5)
  }
  
  # first, generate correlated normal data
  normDat <- rmvnorm(n*numSims,
                     mean=rep(0,nrow(corMat)),
                     sigma=corMat) # marginally standard normal, with off-diagonal elements corresponding to (approximately) desired marginal correlations
  
  # next, apply inverse cdf method from helper function to get simulation data
  # for CG1 to CG2.end, final data needs to be time from transplant, not interval since last CG event; need to add the respective interval to the previous time from HCT
  CG1 <- makeMargVec(dist=dist,
                     paramVec=params[[1]],
                     vecIn=normDat[,1])
  CG1.end <- CG1+
    makeMargVec(dist=dist,
                     paramVec=params[[2]],
                     vecIn=normDat[,2])
  CG2 <- CG1.end+
    makeMargVec(dist=dist,
                paramVec=params[[3]],
                vecIn=normDat[,3])
  CG2.end <- CG2+
    makeMargVec(dist=dist,
                paramVec=params[[4]],
                vecIn=normDat[,4])
  RD <- makeMargVec(dist=dist,
                     paramVec=params[[5]],
                     vecIn=normDat[,5])
  
  # assuming independent, uniform censoring
  C <- runif(n*numSims,params[[6]][1],params[[6]][2])
  
  out <- data.frame(CG1,CG1.end,
             CG2,CG2.end,
             RD,
             C,
             simIdx=rep(1:numSims,each=n) # index the simulations (each of nrow=n) for ease later
  )
  
  out$id <- 1:nrow(out)
  

  return( out )
}

### Calculate empirical CGRFS using simulated data: this is an oracle approach (data that would be censored in practice are observed here)
calcEmpCGRFS <- function(simDatIn,
                         timesIn=NULL){
  if (is.null(timesIn)) timesIn <- seq(0,max(simDatIn),length.out=10000)
  funCGRFS <- Vectorize(function(t) mean(
    (
      t < simDatIn$CG1 | # before CG1 
        (t > simDatIn$CG1.end & t < simDatIn$CG2) | # after CG1.end
        (t > simDatIn$CG2.end) # after CG2.end
    ) & (t < simDatIn$RD) # always need to be before RD
  ))
  
  return(funCGRFS(timesIn))
}

### Calculate empirical pepeQ using simulated data: this is an oracle approach (data that would be censored in practice are observed here)
calcEmpQ <- function(simDatIn,
                         timesIn=NULL){
  if (is.null(timesIn)) timesIn <- seq(0,max(simDatIn),length.out=1000)
  
  # Q is conditional probability of Pr( CG1 OR CG2|not RD)
  
  empCDF <- Vectorize(function(t) mean(
    ((t > simDatIn$CG1 & t < simDatIn$CG1.end) | # in CG1 at t
      (t > simDatIn$CG2 & t < simDatIn$CG2.end)) [ # in CG2 at t
      (t <= pmin(simDatIn$R,simDatIn$D)) # NOT (died or relapsed) before t
  ])
  )
  
  empCDF(timesIn)
}

### Calculate true CGRFS using simulation parameters
calcTruCGRFS <- function(simParamsIn, # assumes use of exponential distributions as used in makeSimDat and identical simulation parameter input
                        timesIn=NULL){
  
  if (is.null(timesIn)) timesIn <- seq(0,simParamsIn[6],length.out=1000)
  
  int1 <- function(t) integrate(function(b) (pexp(t,simParamsIn[1])-pexp(pmax(t-b,0),simParamsIn[1]))*dexp(b,simParamsIn[2]),
                                lower=0,upper=Inf)$value
  int2 <- function(t) integrate(function(b) (pgamma(t,sum(simParamsIn[1:3]),1)-pgamma(pmax(t-b,0),sum(simParamsIn[1:3]),1))*dexp(b,simParamsIn[4]),
                                lower=0,upper=Inf)$value
  
  truCDF <- Vectorize(function(t){ # union of 3 mutually exclusive sets = CDF
    int1(t)*(1-pexp(t,sum(simParamsIn[5])))+ # Pr(CG1 AND not RD)
      int2(t)*(1-pexp(t,sum(simParamsIn[5])))+ # Pr(CG2 AND not RD)
      pexp(t,sum(simParamsIn[5])) # Pr(RD)
  }
  )
  
  1-truCDF(timesIn) # return survival curve
}

### Calculate true prevalence using simulation parameters
calcTruQ <- function(simParamsIn, # assumes use of exponential distributions as used in makeSimDat and identical simulation parameter input
                         timesIn=NULL){
  
  # Q is conditional probability of Pr( CG1 OR CG2|not RD)
  
  if (is.null(timesIn)) timesIn <- seq(0,simParamsIn[8],length.out=1000)
  
  int1 <- function(t) integrate(function(b) (pexp(t,simParamsIn[1])-pexp(pmax(t-b,0),simParamsIn[1]))*dexp(b,simParamsIn[2]),
                                lower=0,upper=Inf)$value
  int2 <- function(t) integrate(function(b) (pgamma(t,sum(simParamsIn[1:3]),1)-pgamma(pmax(t-b,0),sum(simParamsIn[1:3]),1))*dexp(b,simParamsIn[4]),
                                lower=0,upper=Inf)$value
  
  truCDF <- Vectorize(function(t){ # Pr( (CG1 or CG2) AND (not RD))/Pr(not RD)
    (int1(t)*(1-pexp(t,sum(simParamsIn[6:7])))+ # Pr(CG1 AND not RD)
      int2(t)*(1-pexp(t,sum(simParamsIn[6:7]))))/ # Pr(CG2 AND not RD)
      pexp(t,sum(simParamsIn[6:7]),lower.tail=FALSE) # Pr(not RD)
  }
  )
  
  truCDF(timesIn) # return curve
}

### CGRFS adjusted to match boundaries; in simulation, just use CGRFS output, don't recalculate.
## Output: data.frame(times,est(cGRFS.adj),s1,s2,s3,s4,s5)
## Input: output from calcCGRFS
calcCGRFS.Adj.Sim <- function(cgrfsOut){ 
  
  
  out.CGRFSAdj <- pmin(
    pmax(
      cgrfsOut$surv1,cgrfsOut$est
    ),
    cgrfsOut$surv4
  )
  
  return(
    data.frame(
      time=cgrfsOut$time,
      estAdj=out.CGRFSAdj
    )
  )
}
