
# Being Bayesian about 
# Learning Bayesian networks
# from ordinal data

# R code for the Bayesian approach for modelling the dependencies
# among ordinal data in form of latent Bayesian networks (BoB method)

# paper submitted to the
# International Journal of Approximate Reasoning

# author: Marco Grzegorczyk (Groningen University, NL)
# Email to: m.a.grzegorczyk@rug.nl

# Specify your working directory (where the files with R codes have been saved)
# strMCMC_BoB.R
# strMCMC_BoB_mixed.R
# strMCMC_BGe.R
# MORE_FUNCTIONS.R
# STOLEN_CODES.R
# compute_AUCs.R

#################################################################

# For example use directory X:/BoB
setwd("X:/BoB")

# Read two of the provided files with the source commands:
source("strMCMC_BoB.R")
source("MORE_FUNCTIONS.R") 

#################################################################

# This R code shows how to reproduce the BoB results from Figure 1.

# It makes use of R code provided by 
# Xiang Ge Luo, Giusi Moffa, Jack Kuipers
# Learning Bayesian networks from ordinal data
# https://github.com/xgluo/OSEM/blob/master/R/ordinalScore.R

# Install the R package 'pcalg' from CRAN and load it
# https://www.rdocumentation.org/packages/pcalg/versions/2.7-11/topics/randDAG
# https://cran.r-project.org/web//packages/pcalg/pcalg.pdf

install.packages("pcalg") 
library(pcalg)

# Load the R code from Xiang Ge Luo, Giusi Moffa, Jack Kuipers

source("STOLEN_CODES.R")

############################################################

# Generate data set with n=20 nodes and N=500 observations

n <- 20
N <- 500

# Generate a regular DAG with 20 nodes and with 4 neighbors (cf. Luo et al)
trueDAG <- randDAG(n = n, d = 4, method = "er", wFUN = list(mywFUN))

# Generate the ordinal dataset

ordinal_data <- generateOrdinal(N, n, trueDAG, exp_levels = 4, concent_param = 2)

# To adjust for the BoB method transpose and shift the values by 1:

data = t(ordinal_data) + 1

############################################################
# Analyse the data with the BoB method

# Total number of MCMC iterations (reduced)
iterations = 3000 

# We used 
# iterations = 3000000

# Thin out by the factor
thin_factor  =   15 

# We used 
# thin_factor = 1500

# Probability with which the latent values are re-sampled in each iteration
q_resample   = 0.05


# Run the MCMC simulation for the BoB method
out_BoB = strMCMC_BoB(data,iterations,thin_factor,q_resample)

results_BoB = cpdag_list(out_BoB,100)

BoB_SCORES = results_BoB[[3]]

BoB_SCORES

############################################################

# For the AUROC and AUPR computation install the R package 'pracma'
# https://cran.r-project.org/web/packages/pracma/index.html

install.packages("pracma") 
library(pracma)

# Load the source code for computing the AUROC and the AUPR score
source("compute_AUCs.R") 

# Transform trueDAG into adjacency matrix TRUE_DAG
TRUE_DAG = as(trueDAG,"matrix")
TRUE_DAG[which(TRUE_DAG != 0)] <- 1

# COMPUTE AUROC
res = compute_AUROC_NEW(TRUE_DAG, BoB_SCORES)
res[[1]]

# COMPUTE AUPR
res = compute_AUPR_NEW(TRUE_DAG, BoB_SCORES)
res[[1]]


# The AUCs are rather low.
# The number of MCMC iterations must be drastically increased!

# The BoB AUPR results reported Figure 1 can re reproduced by 
# analysing 10 data sets for each of the nine combinations of n and N
# We considered the network sizes n=10,20,30 
# and the numbers of observations N=100,200,500

# The OSEM results have been obtained with the R code provided by
# Xiang Ge Luo, Giusi Moffa, Jack Kuipers
# Learning Bayesian networks from ordinal data
# https://github.com/xgluo/OSEM




