###############################
# simulation settings
#
# Subodh Selukar
# 2025-11-21
###############################

### This program runs simulation and calls masterFunction.R

### This script holds calls to 
# 1. ETM transition matrix
# 2. Correlation matrices
# 3. Parameter estimates
# 4. Output simulation results (need to update file paths to your local computer)

###Creating taskid - will be equal to the row in sim.settings
taskid <- as.numeric(Sys.getenv("LSB_JOBINDEX"))


###Create vectors of simulation settings for each variable
v.params <- c("Set1","Set2","Set3","Set4","Set5")
v.ss <- c(50,75,100,250,1000)
v.cens <- c(0.50,0.10)
v.sims <- c(10000)
v.seed <- c(2023)
v.cor <-c("NULL","Mat1", "Mat2")
v.d <- c("exp")


###Create grid of combinations of simulation settings
sim.settings <- expand.grid(v.params,v.ss, v.cens, v.sims, v.seed, v.cor, v.d, stringsAsFactors=FALSE)
#row corresponds to taskid


###Assign row to ThisJob to
ThisJob <-sim.settings[taskid,]



###call transition matrix for etmfunctions
tra <-matrix(FALSE,3,3)
tra[1,2:3] <-FALSE
tra[2,3] <-TRUE
tra[2,1] <-TRUE
tra[3,1] <-TRUE
tra[3,2] <-TRUE
tra


###call HAPNK1 corr matrix
r1 <- c(1,	-0.1822128,	0.5,	0,	0.5694575)
r2 <- c(-0.1822128,	1,	0,	0.7,	-0.3373307)
r3 <- c(0.5,	0,	1,	-0.1822128,	0.5694575)
r4 <- c(0,	0.7,	-0.1822128,	1,	-0.3373307)
r5 <- c(0.5694575,	-0.3373307,	0.5694575,	-0.3373307,	1)
HAPNK1corr <- rbind(r1,r2,r3,r4,r5)


### make hypothetical matrices for simulations
colCG1 <- c(1,-0.2,0.1,-0.1,0.5)
colCG1end <- c(-0.2,1,-0.1,0.1,-0.4)
colCG2 <- c(0.1,-0.1,1,-0.2,0.5)
colCG2end <- c(-0.1,0.1,-0.2,1,-0.4)
colRD <- c(0.5,-0.4,0.5,-0.4,1)

hypMat1 <- cbind(colCG1,colCG1end,colCG2,colCG2end,colRD)

hypMat2 <- hypMat1
hypMat2[5,3:4] <- hypMat2[5,3:4]/2
hypMat2[3:4,5] <- hypMat2[3:4,5]/2


###Calling function programs
source("/etmFunctions.R")
source("/anaFunctions.R")
source("/simFunctions.R")
source("/masterFunction.R")

print(ThisJob)
### Assign parameter lists and correlation matrices

if (ThisJob$Var1 == "Set1") {
  params <- list(1,1,1,1,1e-2,c(0,qexp(ThisJob$Var3, rate=1e-2)))
} else if (ThisJob$Var1 == "Set2") {
  params <- list(1,1e-1,1,1e-1,1e-2,c(0,qexp(ThisJob$Var3, rate=1e-2)))
} else if (ThisJob$Var1 == "Set3") 
  {params <- list(1,1,1e-1,1,1e-2,c(0,qexp(ThisJob$Var3, rate=1e-2)))
} else if (ThisJob$Var1 == "Set4") {
  params <- list(1,1,1,1e-1,1e-2,c(0,qexp(ThisJob$Var3, rate=1e-2)))
} else if (ThisJob$Var1 == "Set5") {
  params <- list(1,1e-1,1,1,1e-2,c(0,qexp(ThisJob$Var3, rate=1e-2)))
} 

if (ThisJob$Var6 == "NULL") {
  cormtx <- NULL 
} else if (ThisJob$Var6 == "Mat1") {
  cormtx <- hypMat1
} else if (ThisJob$Var6 == "Mat2") {
  cormtx <- hypMat2
}

options(warn=-1)

SimOut <- MasterFunction(sim.out.n = ThisJob$Var2, sim.out.numSims = ThisJob$Var4, sim.out.seedIn = ThisJob$V5,
                         sim.out.corMat = cormtx, 
                         sim.out.params = params,
                         sim.out.dist = ThisJob$Var7)


###   print matrix of differences          
write.csv(SimOut[[1]],
          file=paste("/SimOutput/MeanDifs",ThisJob$Var1,ThisJob$Var3,ThisJob$Var2,ThisJob$Var6,ThisJob$Var7,".csv", sep = "_"), 
          row.names=TRUE )

write.csv(SimOut[[2]],
          file=paste("/SimOutput/Boundries",ThisJob$Var1,ThisJob$Var3,ThisJob$Var2,ThisJob$Var6,ThisJob$Var7,".csv",sep = "_"), 
          row.names=TRUE )


###    print True CGRFS
write.csv(SimOut[[3]],
          file=paste("/SimOutput/EmpCGRFS",ThisJob$Var1,ThisJob$Var3,ThisJob$Var2,ThisJob$Var6,ThisJob$Var7,".csv",sep = "_"),
          row.names =TRUE)


###   print Estimates
write.csv(SimOut[[4]],
          file=paste("/SimOutput/EstimatesETM",ThisJob$Var1,ThisJob$Var3,ThisJob$Var2,ThisJob$Var6,ThisJob$Var7,".csv",sep = "_"),
          row.names = TRUE)

write.csv(SimOut[[5]],
          file=paste("/SimOutput/EstimatesAERFS",ThisJob$Var1,ThisJob$Var3,ThisJob$Var2,ThisJob$Var6,ThisJob$Var7,".csv",sep = "_"),
          row.names = TRUE)

write.csv(SimOut[[6]],
          file=paste("/SimOutput/EstimatesCGRFS",ThisJob$Var1,ThisJob$Var3,ThisJob$Var2,ThisJob$Var6,ThisJob$Var7,".csv",sep = "_"),
          row.names = TRUE)

