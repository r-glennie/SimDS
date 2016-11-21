# SimDs
Code library to simulate line transect distance sampling with moving animals, using R and C++. 

Requires R packages: 
Rcpp
RcppArmadillo

Requires C++ Armadillo library. 

Packaged with MCDS engine from Program Distance (http://distancesampling.org/). Runs on Linux machine using WINE. 

PERFORMING A SIMULATION 
--------------------------------------------

1. Set parameters for simulation in simds.r file. 

2. If needed, alter following functions in simds.cc file: 
  FitModel(): alter to custom model fitting code, uses MCDS under WINE currently. 
  Clean(): function moves and renames simulated data and MCDS output files to a results directory

3. If using MCDS, set options in mcds/command.txt, see Program Distance user manual for description. 

4. Once simulation engine is complete. If using MCDS, set options in outputreader.py and run to create CSV file of parameter estimates for analysis in R. 
