###############################
# master function
#
# Subodh Selukar
# 2025-11-21
###############################

### Function to run all helper functions to conduct simulation

## Output: CSV files of simulation output
## Input:
# N (sample size), numSims (number of simulations),
# seedIn (seed), corMat (correlation matrix),
# params (parameter estimates), dist (sampeling distribution)

##default time grid is 1000
##must enter/change csv path in function to local drive
MasterFunction <- function(
                           sim.out.n,
                           sim.out.numSims,
                           sim.out.seedIn,
                           sim.out.corMat,
                           sim.out.params,
                           sim.out.dist
                           ){
  ### 1. Simulate all data
  sim.out <- makeSimDat(n = sim.out.n, numSims= sim.out.numSims, seedIn = sim.out.seedIn,
                        corMat = sim.out.corMat, params = sim.out.params, dist = sim.out.dist)
  
  ###   a. indexing censoring time
        x <- sim.out.params[[6]][[2]]
        y <- sim.out.params[[5]]
        
  ###   b. creating main.timesIn
         Main.timesIn <-seq(0,x, by=(qexp(0.1,rate=y)/1000))
       
         
  ### 2. Calculate Emperical CGRFS
  
  EmpCGRFS <- calcEmpCGRFS(simDatIn = sim.out, timesIn = Main.timesIn)

  
  
  ### 3. Create empty matrix of size numSims x 6
  MeanDifferences <-matrix(,sim.out.numSims,6)
  
  ### 3.a Create empty matrix of size numSims x 6 to hold logical values
  Boundries <- matrix(,sim.out.numSims,6)
  
  ### 3.b create empty matrix to hold 3 different estimates for each simulation
  EstimatesETM <- matrix(,length(Main.timesIn),sim.out.numSims)
  EstimatesAERFS <- matrix(,length(Main.timesIn),sim.out.numSims)
  EstimatesCGRFS <- matrix(,length(Main.timesIn),sim.out.numSims)
  
  ### 4. Start loop over number of datasets
    
    for(i in 1:sim.out.numSims) {
  ###   a. temp data of ith sim.out
          temp <- sim.out |>
            subset(simIdx== i)
          
  ###   b. calculate analysis functions
          ### CGRFS
            #Creating S1 - S5 Datasets
            s1Dat <- makeS1Dat(temp,CG1,RD,C)
            s2Dat <- makeS2Dat(temp,CG2,RD,C)
            s3Dat <- makeS3Dat(temp,CG1.end,RD,C)
            s4Dat <- makeS4Dat(temp,RD,C)
            s5Dat <- makeS5Dat(temp,CG2.end,RD,C)
          
            #Creating S1-S5 survival estimates to calculate CGRFS
            s1 <- makeSurvEst(datIn=s1Dat,timesIn=Main.timesIn)
            s2 <- makeSurvEst(datIn=s2Dat,timesIn=Main.timesIn)
            s3 <- makeSurvEst(datIn=s3Dat,timesIn=Main.timesIn)
            s4 <- makeSurvEst(datIn=s4Dat,timesIn=Main.timesIn)
            s5 <- makeSurvEst(datIn=s5Dat,timesIn=Main.timesIn)
          
            #output CGRFS estimates at Main.timesIn
            cgrfs_est <- calcCGRFS2(s1,s2,s3,s4,s5,timesIn=Main.timesIn)
          
          ### ETM
            # Transform simulation data to fit ETM function
            EtmDatOut <- makeEtmDat(datIn=temp)
            
            # ETM summary function
            EtmSum <- makeEtmSum(EtmDatOut, timeYears = Main.timesIn)
            
            
          ### AERFS
            #output AERFS estimates at main.timesIn
            
            #summary of S4dat for aerfs estimates
            s4.aerfs <- summOfSxDat(datIn=s4Dat,timesIn = Main.timesIn)
            
            
            aerfs_est <- calcAERFS2(datIn = temp,timesIn = Main.timesIn,s4summDat = s4.aerfs)
            
  ###   c. calculate Mean differences and save in vector
          absmeans<-  calcMeanAbsDiff(timesIn = Main.timesIn,
                            EtmSumFitTimes = EtmSum$`0 0`$P,
                            AERFSsumFitTimes = aerfs_est$est,
                            CGRFSsumFitTimes = cgrfs_est$est,
                            empCGRFSsumFitTimes = EmpCGRFS)
            
          meansqu <- calcMeanSquaredDiff(timesIn = Main.timesIn,
                                EtmSumFitTimes = EtmSum$`0 0`$P,
                                AERFSsumFitTimes = aerfs_est$est,
                                CGRFSsumFitTimes = cgrfs_est$est,
                                empCGRFSsumFitTimes = EmpCGRFS)
          
  ###   c.1 determine if ETM, AERFS, and CGRFS are in bounds of S4 and S1
          ETM_S4<-any(ifelse(EtmSum$`0 0`$P > s4$est,"Yes","No") =="Yes")
          AERFS_S4<-any(ifelse(aerfs_est$est > s4$est,"Yes","No") =="Yes")
          CGRFS_S4<-any(ifelse(cgrfs_est$est > s4$est,"Yes","No") =="Yes")
          
          ETM_S1<-any(ifelse(EtmSum$`0 0`$P < s1$est,"Yes","No") =="Yes")
          AERFS_S1<-any(ifelse(aerfs_est$est < s1$est,"Yes","No") =="Yes")
          CGRFS_S1<-any(ifelse(cgrfs_est$est < s1$est,"Yes","No") =="Yes")
         
          
  ###   d. save vectors to matrix
             MeanDifferences[i,1:3] <- absmeans
             MeanDifferences[i,4:6] <- meansqu
             
  ###   d.1 save into matrix           
             Boundries[i,1] <-ETM_S4
             Boundries[i,2] <-AERFS_S4
             Boundries[i,3] <-CGRFS_S4
             Boundries[i,4] <-ETM_S1
             Boundries[i,5] <-AERFS_S1
             Boundries[i,6] <-CGRFS_S1
             
  ###   d.2 save estimates into matrix
             EstimatesETM[,i] <-EtmSum$`0 0`$P
             EstimatesAERFS[,i] <-aerfs_est$est
             EstimatesCGRFS[,i] <-cgrfs_est$est
             
            
    }
        
            
  return(list(MeanDifferences,Boundries,EmpCGRFS,EstimatesETM,EstimatesAERFS,EstimatesCGRFS))
  
        
}
