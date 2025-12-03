###############################
# analysis functions
#
# Subodh Selukar
# 2025-11-21
###############################

library(survival)


##### Helper functions

### Helper function to return data frame of times and estimates using survival package
## Input: 
# 1. simpleDat = data.frame(Y,S) where Y is observed time for respective TTE and S is event indicator
# 2. timesIn (as described above)
wrapSurv <- function(simpleDat,timesIn=NULL){
  
  if (is.null(timesIn)) timesIn <- simpleDat$Y
  
  survOut <- with(simpleDat,survfit(Surv(Y,S)~1))
  summSurv <- summary(survOut,times=timesIn,
                      extend=TRUE) # extra option needed in case timesIn go beyond max(Y)
  
  return(
    data.frame(time=summSurv$time,est=summSurv$surv)
  )
}

### Helper function to create simpleDat, used in wrapSurv
## Input:
# 1. eventMat: matrix/data frame containing all of the times of "event" components; variable number of columns, but nrow(eventMat)=n
# variable number allows for arbitrary number of components to "event" definition
# 2. censTime: vector of censoring times with compatible length (censTime should never be NA)
makeSimpleDat <- function(eventMat,censTime){
  eventTime <- apply(eventMat,1,min,na.rm=TRUE) # in real data, event times may be NA because of censoring/competing risks
  
  Y <- pmin(eventTime,censTime,na.rm=TRUE) # observed time
  S <- (eventTime <= censTime & !is.na(eventTime))  # indicator of event status
  
  return(
    data.frame(Y,S)
  )
}

##### Functions to calculate different survival estimators

### datIn: (common input to all functions)
## Rows are independent subjects
## Columns are time to the following events: (NOT dates)
# 1. First moderate-severe cGVHD, CG1
# 2. End of first mod-sev cGVHD, CG1.end
# 3. Second mod-sev cGVHD, CG2
# 4. End of second mod-sev cGVHD, CG2.end
# 5. First grade 3-4 aGVHD, AG
# 6. Relapse or disease progression, R
# 7. All-cause death, D
# NOTE: any of the above can be recorded as missing in real data if they are not observed
# 8. Last follow-up, C 
# if C is missing, need to assign C = max(CG1,CG2,AG,R,D)

### timesIn: times at which to calculate estimates (default NULL uses all observed times, common to all functions)

### OS (time to death)
## Output: data.frame(times,est(DFS))
calcOS <- function(datIn,timesIn=NULL){
  osMat <- data.frame(datIn$D)
  
  simpleDat <- makeSimpleDat(osMat,datIn$C)
  
  return(
    wrapSurv(simpleDat,timesIn)
  )
}

### DFS (time to relapse/progression or death)
## Output: data.frame(times,est(DFS))
calcDFS <- function(datIn,timesIn=NULL){
  dfsMat <- data.frame(datIn$R,datIn$D)
  
  simpleDat <- makeSimpleDat(dfsMat,datIn$C)
  
  return(
    wrapSurv(simpleDat,timesIn)
  )
}

### anycGVHD.DFS (time to first mod-sev cGVHD, relapse/progression or death)
## Output: data.frame(times,est(DFS))
calcAnycGVHD.DFS <- function(datIn,timesIn=NULL){
  cghvd.dfsMat <- data.frame(datIn$CG1,datIn$R,datIn$D) # CG1 must be < CG2 by construction so we don't need to address
  
  simpleDat <- makeSimpleDat(cghvd.dfsMat,datIn$C)
  
  return(
    wrapSurv(simpleDat,timesIn)
  )
}

### GRFS (time to first gr3-4 aGVHD, mod-sev cGVHD, relapse/progression or death)
## Output: data.frame(times,est(DFS))
calcGRFS <- function(datIn,timesIn=NULL){
  grfsMat <- data.frame(datIn$AG,datIn$CG1,datIn$R,datIn$D) # CG1 must be < CG2 by construction so we don't need to address
  
  simpleDat <- makeSimpleDat(grfsMat,datIn$C)
  
  return(
    wrapSurv(simpleDat,timesIn)
  )
}

### Helper function to create s1Dat
## Output: data.frame(CG1,RD,C)
makeS1Dat <- function(datIn,CG1,RD,C){
 return(makeSimpleDat(data.frame(datIn$CG1,datIn$RD),datIn$C))
 
}

### Helper function to create s2Dat
## Output: data.frame(CG2,RD,C)
makeS2Dat <- function(datIn,CG2,RD,C){
   return(makeSimpleDat(data.frame(datIn$CG2,datIn$RD),datIn$C))
 
}

### Helper function to create s3Dat
## Output: data.frame(CG1.end, RD, C)
makeS3Dat <- function(datIn,CG1.end,RD,C){
  return( makeSimpleDat(data.frame(datIn$CG1.end,datIn$RD),datIn$C))
 
}


### Helper function to create s4Dat
## Output: data.frame(RD,C)
makeS4Dat <- function(datIn,RD,C){
  return( makeSimpleDat(data.frame(datIn$RD),datIn$C))

}

### Helper function to create s5Dat
## Output: data.frame(CG2.end,RD,C)
makeS5Dat <- function(datIn,CG2.end,RD,C){
return(
  makeSimpleDat(
    data.frame(datIn$CG2.end,datIn$RD),datIn$C)
  ) 

}



### Helper function to create survival est at parameter of interest
## Output: data.frame(est,time)
makeSurvEst <- function(datIn,timesIn){
  return(
    wrapSurv(
      datIn,timesIn
    )
  )
}






### CGRFS (Solomon et al. 2017)
## Output: data.frame(times,est(cGRFS),s1,s2,s3,s4,s5)
## Input: 
# (standard inputs, datIn and timesIn)
# timesLength: length of the times vector when timesIn is not provided (NA=use observed times, default to 10,000 steps)


calcCGRFS <- function(datIn,timesIn=NULL,timesLength=10000){ 
  
  
  s1Dat <- makeSimpleDat(data.frame(datIn$CG1,datIn$R,datIn$D),datIn$C)
  s2Dat <- makeSimpleDat(data.frame(datIn$CG2,datIn$R,datIn$D),datIn$C)
  s3Dat <- makeSimpleDat(data.frame(datIn$CG1.end,datIn$R,datIn$D),datIn$C)
  s4Dat <- makeSimpleDat(data.frame(datIn$R,datIn$D),datIn$C)
  s5Dat <- makeSimpleDat(data.frame(datIn$CG2.end,datIn$R,datIn$D),datIn$C)
  
  if (is.null(timesIn) & is.na(timesLength)) {
    timesIn <- sort(unique(c(s1Dat$Y,s2Dat$Y,s3Dat$Y,s4Dat$Y,s5Dat$Y))) # calculate curves across all observed times for all endpoints
  } else if (is.null(timesIn) & !is.na(timesLength)) {
    timesIn <- seq(0,max(datIn[,-8],na.rm=TRUE),length.out=timesLength) # make a time grid of timesLength points from 0 to max observed event time
  }
  
  s1 <- wrapSurv(
    s1Dat,timesIn
  )
  s2 <- wrapSurv(
    s2Dat,timesIn
  )
  s3 <- wrapSurv(
    s3Dat,timesIn
  )
  s4 <- wrapSurv(
    s4Dat,timesIn
  )
  s5 <- wrapSurv(
    s5Dat,timesIn
  )
  
  out.CGRFS <- s1$est+(s2$est-s3$est)+(s4$est-s5$est)
  
  return(
    data.frame(
      time=timesIn,
      est=out.CGRFS,
      surv1=s1$est,
      surv2=s2$est,
      surv3=s3$est,
      surv4=s4$est,
      surv5=s5$est
    )
  )
}


### CGRFS (Solomon et al. 2017); second version of CGRFS for simulations to accept estimates (of compatible length) of s1-s5
## Output: data.frame(time,est(cGRFS),surv1,surv2,surv3,surv4,surv5)
## Input: 
# (makesurvdat output s1 - s5, timesIn)
# timesLength: length of the times vector when timesIn is not provided (NA=use observed times, default to 10,000 steps)

calcCGRFS2 <- function(s1,s2,s3,s4,s5,timesIn){ 
 
  
  CGRFS <- s1$est+(s2$est-s3$est)+(s4$est-s5$est)
  
  return(data.frame(
    time=timesIn,
    est=CGRFS,
    surv1=s1$est,
    surv2=s2$est,
    surv3=s3$est,
    surv4=s4$est,
    surv5=s5$est
  )
  )
  

}



### CGRFS adjusted to match boundaries
## Output: data.frame(times,est(cGRFS.adj))
## Input: 
# (standard inputs, datIn and timesIn)
# timesLength: length of the times vector when timesIn is not provided (NA=use observed times, default to 10,000 steps)
calcCGRFS.Adj <- function(datIn,timesIn=NULL,timesLength=10000){ 
  if (is.null(timesIn) & is.na(timesLength)) {
    timesIn <- sort(unique(c(s1Dat$Y,s2Dat$Y,s3Dat$Y,s4Dat$Y,s5Dat$Y))) # calculate curves across all observed times for all endpoints
    # invalid to use times > last observed event?
  } else if (is.null(timesIn) & !is.na(timesLength)) {
    timesIn <- seq(0,max(datIn[,-8],na.rm=TRUE),length.out=timesLength) # make a time grid of timesLength points from 0 to max observed event time
  }
  
  out.CGRFS <- calcCGRFS(datIn=datIn,timesIn=timesIn)
  
  out.CGRFSAdj <- pmin(
    pmax(
      out.CGRFS$surv1,out.CGRFS$est
    ),
    out.CGRFS$surv4
  )
  
  return(
    data.frame(
      time=timesIn,
      estAdj=out.CGRFSAdj
    )
  )
}

### CGRFS adjusted to match boundaries 2
## Output: data.frame(times,est(cGRFS.adj))
## Input: 
# (standard inputs, datIn(cgrfsDat) and timesIn)
# timesLength: length of the times vector when timesIn is not provided (NA=use observed times, default to 10,000 steps)

calcCGRFS.Adj2 <- function(timesIn,cgrfsDat){ 
  
  
  out.CGRFSAdj <- pmin(
    pmax(
      cgrfsDat$surv1,cgrfsDat$est
    ),
    cgrfsDat$surv4
  )
  
  return(
    data.frame(
      time=timesIn,
      estAdj=out.CGRFSAdj
    )
  )
  
}



### Pepe prevalence, Q (1991)
## Output: data.frame(times,est(pepeQ),s1,s2,s3,s4,s5)
## Input: 
# (standard inputs, datIn and timesIn)
# timesLength: length of the times vector when timesIn is not provided (NA=use observed times, default to 10,000 steps)
pepeQ <- function(datIn,timesIn=NULL,timesLength=10000){ 
  
  # use the same notation as with CGRFS for consistency
  s1Dat <- makeSimpleDat(data.frame(datIn$CG1,datIn$R,datIn$D),datIn$C)
  s2Dat <- makeSimpleDat(data.frame(datIn$CG2,datIn$R,datIn$D),datIn$C)
  s3Dat <- makeSimpleDat(data.frame(datIn$CG1.end,datIn$R,datIn$D),datIn$C)
  s4Dat <- makeSimpleDat(data.frame(datIn$R,datIn$D),datIn$C)
  s5Dat <- makeSimpleDat(data.frame(datIn$CG2.end,datIn$R,datIn$D),datIn$C)
  
  if (is.null(timesIn) & is.na(timesLength)) {
    timesIn <- sort(unique(c(s1Dat$Y,s2Dat$Y,s3Dat$Y,s4Dat$Y,s5Dat$Y))) # calculate curves across all observed times for all endpoints
  } else if (is.null(timesIn) & !is.na(timesLength)) {
    timesIn <- seq(0,max(datIn[,-8],na.rm=TRUE),length.out=timesLength) # make a time grid of timesLength points from 0 to max observed event time
  }
  
  s1 <- wrapSurv(
    s1Dat,timesIn
  )
  s2 <- wrapSurv(
    s2Dat,timesIn
  )
  s3 <- wrapSurv(
    s3Dat,timesIn
  )
  s4 <- wrapSurv(
    s4Dat,timesIn
  )
  s5 <- wrapSurv(
    s5Dat,timesIn
  )
  
  out.Q <- (
    (s3$est-s1$est)+(s5$est-s2$est)
  )/s4$est # conditional on alive and disease-free
  
  return(
    data.frame(
      time=timesIn,
      est=out.Q,
      surv1=s1$est,
      surv2=s2$est,
      surv3=s3$est,
      surv4=s4$est,
      surv5=s5$est
    )
  )
}

### AERFS 

## internal function: calculate number in CG who are alive relapse-free
# Input: datIn used in AERFS and timesVec
# Output: data.frame(Y,S) as in makeSimpleDat 
calcM.AERFS <- function(datIn,timesIn){
  

  calcM.functionOfTime <- Vectorize(function(t) sum( # compute M(t) for one t
    (datIn$C > t & ( is.na(datIn$RD) | (!is.na(datIn$RD) & datIn$RD > t) )) & # uncensored, alive, relapse-free at t AND...
      (
        (!is.na(datIn$CG1) & # check CG1
           (
             (is.na(datIn$CG1.end) & t >= datIn$CG1) | # no CG1.end and in CG1
               (!is.na(datIn$CG1.end) & t >= datIn$CG1 & t < datIn$CG1.end) # yes CG1.end and in CG1 ** by convention, CG1.end counts as OFF cGVHD
             
           )
        ) | 
          (!is.na(datIn$CG2) & # check CG2
             (
               (is.na(datIn$CG2.end) & t >= datIn$CG2) | # no CG2.end and in CG2
                 (!is.na(datIn$CG2.end) & t >= datIn$CG2 & t < datIn$CG2.end) # yes CG2.end and in CG2 ** by convention, CG2.end counts as OFF cGVHD
               
             )
          )
        ) 
    ) 
  )
  
  calcM.functionOfTime(timesIn) # calculates at all times in timesIn
}

## Main AERFS function
## Output: data.frame(times,est(AERFS),s4,m,n)
## Input: 
# (standard inputs, datIn and timesIn)
# timesLength: length of the times vector when timesIn is not provided (NA=use observed times, default to 10,000 steps)
calcAERFS <- function(datIn,timesIn=NULL,timesLength=10000){ 
  
  s4Dat <- makeSimpleDat(data.frame(datIn$RD),datIn$C)
  
  
  if (is.null(timesIn) & !is.na(timesLength)) {
    timesIn <- seq(0,max(datIn[,-8],na.rm=TRUE),length.out=timesLength) # make a time grid of timesLength points from 0 to max observed event time
  }
  
  s4 <- summary(survfit(Surv(s4Dat$Y,s4Dat$S)~1),
                times=timesIn,
                extend=TRUE) # extra option needed in case timesIn go beyond max(Y)
  
  M <- calcM.AERFS(datIn,timesIn)
  
  out.AERFS <- s4$surv*(1-(M/s4$n.risk))
  
  return(
    data.frame(
      time=timesIn,
      est=out.AERFS,
      s4=s4$surv,
      M=M,
      N=s4$n.risk
    )
  )
}
### Helper function to create full summary(survfit(surv))
## Output: summary at timesIn of sXDat
summOfSxDat <- function(datIn,timesIn){
  return(summary(survfit(Surv(datIn$Y,datIn$S)~1),
                 times=timesIn,
                 extend=TRUE) # extra option needed in case timesIn go beyond max(Y)
         ) 
}


## Main AERFS function 2: this version addresses the rare case when the number at risk decreases to 0 before the end of follow-up times
## Output: data.frame(times,est(AERFS),s4,m,n)
## Input: 
# (standard inputs, datIn and timesIn, s4summDat)
calcAERFS2 <- function(datIn,timesIn,s4summDat){ 
  
  M <- calcM.AERFS(datIn,timesIn)
  
  out.AERFS <- s4summDat$surv*(1-(M/s4summDat$n.risk))
  
  when0 <- min(which(s4summDat$n.risk == 0 ))
  if (!is.infinite(when0)) {
    out.AERFS[when0:length(out.AERFS)] <-s4summDat$surv[when0]
  } 
  
  return(
    data.frame(
      time=timesIn,
      est=out.AERFS,
      s4=s4summDat$surv,
      M=M,
      N=s4summDat$n.risk
    )
  )
}



## Mean squared difference function
## Calculates and saves mean absolute difference of estimates
## Output: vector of outcomes for ETM, AERFS, CGRFS, and TrueCGRFS
## Input: timesIn, ETM estimates at timesIn, AERFS estimates at timesIn,
# CGRFS estimates at timesIn, EmpCGRFS estimates at timesIn

calcMeanAbsDiff <- function(timesIn,EtmSumFitTimes,AERFSsumFitTimes,CGRFSsumFitTimes,empCGRFSsumFitTimes) {
  compareDf <- cbind(timesIn,ETM=EtmSumFitTimes,AERFS=AERFSsumFitTimes,CGRFS=CGRFSsumFitTimes,empCGRFS=empCGRFSsumFitTimes)
  

  return(
    vec.abs.means <- c(mean(abs(compareDf[,2]-compareDf[,5])),
                       mean(abs(compareDf[,3]-compareDf[,5])),
                       mean(abs(compareDf[,4]-compareDf[,5])))
  )
 
}

## Mean squared difference function
## Calculates and saves mean squared difference of estimates
## Output: vector of outcomes for ETM, AERFS, CGRFS, and TrueCGRFS
## Input: timesIn, ETM estimates at timesIn, AERFS estimates at timesIn,
# CGRFS estimates at timesIn, EmpCGRFS estimates at timesIn
calcMeanSquaredDiff <- function(timesIn,EtmSumFitTimes,AERFSsumFitTimes,CGRFSsumFitTimes, empCGRFSsumFitTimes) {
  compareDf <- cbind(timesIn,ETM=EtmSumFitTimes,AERFS=AERFSsumFitTimes,CGRFS=CGRFSsumFitTimes,empCGRFS=empCGRFSsumFitTimes)
  
  return(
    vec.means.squ <- c(mean(compareDf[,2]-compareDf[,5])**2,
                       mean(compareDf[,3]-compareDf[,5])**2,
                       mean(compareDf[,4]-compareDf[,5])**2
      
    )
  )
  
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


