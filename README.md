AERFS README



This repository contains R scripts to analyze the CGRFS endpoint using different candidate estimators. The R scripts support a user to estimate CGRFS with real data and to replicate the simulation study in Selukar et al.



anaFunctions.R holds functions to calculate and output CGRFS and AERFS estimates. The script includes helper functions related to those estimators, including data structures, and includes code for calculating Pepe's Q function from Pepe et al. (1991) Stat Med.

etmFunctions.R holds functions to transform simulated data to appropriate format for etm() and runs etm() on transformed data. The multi-state model for this etm call corresponds to the 3-state model in Selukar et al.

simFunctions.R holds functions to simulate data. The script also includes a function calculate the oracle estimates (Monte Carlo-based "true" functions) of CGRFS and Pepe Q by pooling all simulated datasets.

masterFunction.R holds a function to run all helper functions in anaFunctions.R, etmFunctions.R, simFunctions.R and outputs simulation results.

simSettings.R creates a simulation grid, calls all of the scripts, and runs the simulation for a selected element of the simulation grid.
