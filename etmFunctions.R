###############################
# etm functions
#
# Subodh Selukar
# 2025-11-21
###############################

### Transforms simulation data to format acceptable for etm()
### Runs etm() on transformed data 
## Output: summary(etm())
## Input: datIn

options("install.lock"=FALSE)
#install.packages("etm")
library(etm)
library(dplyr)
library(Rcpp)



makeEtmDat <- function(datIn
                      ){

######## transition from state 0 to state 1 (GVHD 1)
### CG1 < RD & CG1 < C
### exit time = CG1

Trans_0_1 <- datIn |>
  filter(CG1 < RD & CG1 < C)|>
  mutate(from = 0,
         to = 1,
         entry = 0,
         exit = (CG1), 
         id = id) |>
  select(id, entry, exit, from, to)

######## transition from state 1 to state 0 (GVHD1 recovery)
### 
Trans_1_0 <- datIn |>
  filter(CG1.end < RD & CG1.end < C)|> 
  mutate(from = 1,
         to = 0,
         entry = (CG1),
         exit = (CG1.end),
         id = id) |>
  select(id, entry, exit, from, to)


######## transition from state 0 to state 1 (GVHD2)

Trans_0_1_2 <- datIn |>
  filter(CG2 < RD & CG2 < C)|>
  mutate(from = 0,
         to = 1,
         entry = (CG1.end),
         exit = (CG2),
         id = id) |>
  select(id, entry, exit, from, to)

######## transition from state 1 to state 0 (GVHD2 recovery)

Trans_1_0_2 <- datIn |>
  filter(CG2.end < RD & CG2.end < C)|> ## SRS
  mutate(from = 1,
         to = 0,
         entry = (CG2),
         exit = (CG2.end),
         id = id) |>
  select(id, entry, exit, from, to)


######## transition from state 1 (GVHD1) to state 2 (relapse/death)
### CG1.end > RD & RD < C & CG1 < RD

Trans_1_2 <- datIn |>
  filter(CG1.end > RD & RD < C & CG1 < RD)|>
  mutate(from = 1,
         to = 2,
         entry = (CG1),
         exit = (RD),
         id = id) |>
  select(id, entry, exit, from, to)


######## transition from state 1 (GVHD2) to state 2 (relapse/death)
### CG2.end > RD & RD < C & CG2 < RD

Trans_1_2.2 <- datIn |>
  filter(CG2.end > RD & RD < C & CG2 < RD)|>
  mutate(from = 1,
         to = 2,
         entry = (CG2),
         exit = (RD),
         id = id) |>
  select(id, entry, exit, from, to)



######## transition from state 0 to state 2 (relapse/death) #no gvhd
### CG1 > RD & RD < C 

Trans_0_2 <- datIn |>
  filter(CG1 > RD & RD < C ) |>
  mutate(from = 0,
         to = 2,
         entry = 0,
         exit = RD,
         id = id) |>
  select(id, entry, exit, from, to)

######## transition from state 0 to state 2 (relapse/death) #After recovery of CG1
### CG1.end < RD & RD > CG2 & RD < C 

Trans_0_2.2 <- datIn |>
  filter(CG1.end < RD & RD > CG2 & RD < C  ) |>
  mutate(from = 0,
         to = 2,
         entry = CG1.end,
         exit = RD,
         id = id
  ) |>
  select(id, entry, exit, from, to)

######## transition from state 0 to state 2 (relapse/death) #After recovery of CG2
### CG2.end < RD  & RD < C 

Trans_0_2.3 <- datIn |>
  filter(CG2.end < RD  & RD < C  ) |>
  mutate(from = 0,
         to = 2,
         entry = CG2.end,
         exit = RD,
         id = id
  ) |>
  select(id, entry, exit, from, to)


######## transition from state 1 to state cens 
###CG1.end > C & C < RD & CG1 < C
Trans_1_cens <- datIn |>
  filter(CG1.end > C & C < RD & CG1 < C)|>
  mutate(from = 1,
         to = 'cens',
         entry = CG1,
         exit = C,
         id = id) |>
  select(id, entry, exit, from, to)


######## transition from state 1 (CG2) to state cens 
###CG2.end > C & C < RD & CG2 < C
Trans_1_cens.2 <- datIn |>
  filter(CG2.end > C & C < RD & CG2 < C)|>
  mutate(from = 1,
         to = 'cens',
         entry = CG2,
         exit = C,
         id = id) |>
  select(id, entry, exit, from, to)



######## transition from state 0 to state cens  # no cvhd
### CG1 > C & RD > C 

Trans_0_cens <- datIn |>
  filter(CG1 > C & RD > C )|>
  mutate(from = 0,
         to = 'cens',
         entry = 0,
         exit = (C),
         id = id) |>
  select(id, entry, exit, from, to)


######## transition from state 0 to state cens  # from recovery 1
### ### CG1.end < C & C < CG2 & RD > C 

Trans_0_cens_1 <- datIn |>
  filter(CG1.end < C & C < CG2 & RD > C  )|> ## SRS
  mutate(from = 0,
         to = 'cens',
         entry = CG1.end,
         exit = (C),
         id = id) |>
  select(id, entry, exit, from, to)


######## transition from state 0 to state cens (LFU from recovery 2)
### CG2.end < C  & RD > C 
Trans_0_cens_2 <- datIn |>
  filter(CG2.end < C  & RD > C ) |>
  mutate(from = 0,
         to = "cens",
         entry = CG2.end,
         exit = C,
         id = id) |>
  select(id, entry, exit, from, to)



######## Rbind Trans_state_state datasets

DatOut <- rbind(Trans_0_1, Trans_1_0,
              Trans_0_1_2, Trans_1_0_2,
              Trans_0_2, Trans_0_2.2, Trans_0_2.3,
              Trans_1_2, Trans_1_2.2, 
              Trans_1_cens, Trans_1_cens.2,
              Trans_0_cens,
              Trans_0_cens_1,
              Trans_0_cens_2)

return(DatOut)
}


### Function to make ETM summary data
## Output: summary of ETM at timeYears
## Input:
# standard input ( datIn, timeYears), transition matrix tra 
makeEtmSum <- function(datIn , #makeEtmDat output
                       timeYears ,  #times to print
                       transMatrix = tra #transition matrix
                       
    ){


fit <-etm(datIn, state.names = c("2","1","0") , transMatrix, cens.name = "cens", s = 0)

sum.fit <-summary(fit,times = timeYears)

return(sum.fit)
}




