# Author: Richard Glennie
# 2016
# 
# Running a distance sampling simulation where animals move

library(Rcpp)
library(RcppArmadillo) 
sourceCpp("simds.cc")

# set simulation parameters 
population.size <- 100
region.size <- c(1000, 1000)
# simulation assumes all lines are randomly placed 
# and have same length/width
line.size <- c(60, 1000) 
# behaviour = 0 (no movement), 1 (linear movement constant speed, random direction),
# 2 (OU home-range movement) 
behaviour <- 1
observer.speed <- 1
num.transects <- 50
# parameter contains (c, d) detection parameters for continuous hazard 
# (see p Advanced Distance Sampling) 
# Afterward, parameter contains movement parameters: 
# behaviour = 0 (no parameters added), 1 (add animal speed), 2 (add c, tau parameters (see Gillespie et al (1996) for parameterisation) 
parameter <- c(100, 3, 1) 

# run simulation 
sim.dat <- c(population.size, region.size, line.size, behaviour, observer.speed, num.transects) 

num.simulations <- 100
dt <- 1

Simulate(num.simulations, parameter, sim.dat, dt) 
