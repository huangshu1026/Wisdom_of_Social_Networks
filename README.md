# Wisdom of Social Networks 

This repository contains R code for replicating simulations, tables and figures 
described in the working paper "Wisdom of Social Networks" 
by Shu Huang, Stephen B. Broomell, and Russell Golman. 

## How to use
1. Clone this repository to your local computer
2. For each script to run, be sure to change the working directory in each script (searching "setwd()")

## Algorithm
1. The file "CommunityDetection.R" contains the replication of Newman (2006)'s community detection algorithm that maximizes the modularity, and the algorithm we compute the optimal weights given community detection results. 
2. The file "BayesianNetModel.R" contains the Bayesian network algorithm. We used a modified MCMC algorithm to obtain the posterior probabilities of the given number of possible type configurations. Also, this document also includes the algorithm that we use to compute the optimal weights given Bayesian types and their probabilities. 
3. The file "ExternalEstimation.R" contains the functions we use to simulate multiple social networks given a specific type configuration and then estimate the within-type correlation and across-type correlation. 
4. The file "GlobalFunctions.R" contains multiple functions for general use, e.g., convert the type vector into type matrix; generate networks given a type matrix and connection probabilities; find the nearest positive-definite matrix for a non-singular matrix. 

## Data
1. There is no data file.

## Simulation
1. The file "simulation_grid.R" contains the sample code for simulation in the working paper. 
