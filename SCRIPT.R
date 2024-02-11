
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

#################################################################

# For example use directory X:/BoB
setwd("X:/BoB")

# Read the provided files with the source command
source("strMCMC_BoB.R")
source("strMCMC_BoB_mixed.R") 
source("strMCMC_BGe.R") 
source("MORE_FUNCTIONS.R") 


#################################################################

# For illustration purposes we use the data from Geiger and Heckerman (1994)
# Load the continuous and the discretized version of the data 
# cf. Table 1 in Section S1 of the Supplement
# There are n=3 variables and N=20 observations (cases)

load("data_GH.Rdata")
load("data_GH_dis.Rdata")

# continuous (latent) data
data_con = data_GH

# discretized (observed) data 
data_dis = data_GH_dis

######################################################################
# Interestingly, we will find that the DAG X1 -> X2 <- X3
# has the highest posterior probability

# Hence, we suppose that Figure 1 in Geiger and Heckerman (1994) 
# should show the DAG X1 -> X2 <- X3 rather than X1 -> X3 <- X2
#######################################################################

# Optionally load your own continuous and/or ordinal data set 
# 'data_con' and/or 'data_dis'

# Both must be n-by-N data matrices, where 
# data_con[i,j] is the j-th continuous value of node X_i
# data_dis[i,j] is the j-th discrete value of node X_i


# The discrete levels are supposed to be integers, 
# and must start with the integer 1

##################################################################

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

# iterations/thin_factor (here: 3000/15 = 200) DAGS will be sampled 

set.seed(1994)

####################################################################
# BoB method on ordinal data
####################################################################

# First, we run the MCMC simulation for the BoB method
out_BoB = strMCMC_BoB(data_dis,iterations,thin_factor,q_resample)

# out_BoB[[1]][[1]] is the adjacency matrix of the initial (=empty) DAG
# out_BoB[[2]][[1]] is the score        of the initial (=empty) DAG

# out_BoB[[1]][[k+1]] is the adjacency matrix of the k-th sampled DAG
# out_BoB[[2]][[k+1]] is the score            of the k-th sampled DAG

# Process the results by computing CPDAGs and averaging across them
# Skip the first 100 DAGs (to take burn in phase into account)

results_BoB = cpdag_list(out_BoB,100)

BoB_SCORES = results_BoB[[3]]

# Edge scores of the BoB method:

BoB_SCORES

####################################################################
# BGe method on continuous data
####################################################################

# Second, we run the MCMC simulation for the BGe method
out_BGe = strMCMC_BGe(data_con,iterations,thin_factor)

# out_BGe[[1]][[1]] is the adjacency matrix of the initial (=empty) DAG
# out_BGe[[2]][[1]] is the score        of the initial (=empty) DAG

# out_BGe[[1]][[k+1]] is the adjacency matrix of the k-th sampled DAG
# out_BGe[[2]][[k+1]] is the score            of the k-th sampled DAG

# Process the results by computing CPDAGs and averaging across them
# Skip the first 100 DAGs (to take burn in phase into account)

results_BGe = cpdag_list(out_BGe,100)

BGe_SCORES = results_BGe[[3]]

# Edge scores of the BGe method:

BGe_SCORES


####################################################################
# BoB method on mixed data (discrete and continuous mixed)
####################################################################

# Third, we run the MCMC simulation for the BOB method on mixed data

# Let's assume that the 1st variable is continuous
# and that the 2nd and 3rd variable are ordinal


# We mix rows of data_con and data_dis correspondingly:

data_mixed = rbind(data_con[1,],data_dis[2:3,])

# and we define an indicator vector that indicates the ordinal variables

INDICATOR_VEC = c(1,1,1)

# Which variables are ordinal?
# which(INDICATOR_VECTOR)

out = strMCMC_BoB_mixed(data_dis,iterations,thin_factor,q_resample, INDICATOR_VEC)

# out[[1]][[1]] is the adjacency matrix of the initial (=empty) DAG
# out[[2]][[1]] is the score        of the initial (=empty) DAG

# out[[1]][[k+1]] is the adjacency matrix of the k-th sampled DAG
# out[[2]][[k+1]] is the score            of the k-th sampled DAG

# Process the results by computing CPDAGs and averaging across them
# Skip the first 100 DAGs (to take burn in phase into account)

results_BoB_mixed = cpdag_list(out,100)

BoB_mixed_SCORES = results_BoB_mixed[[3]]

# Edge scores of the BGe method:

BoB_mixed_SCORES


#########################################################################

# JUST A NOTE 
# 'strMCMC_BoB_mixed' can also be used for pure ordinal data
# That is, 'strMCMC_BoB' is a special case of 'strMCMC_BoB_mixed'

# EXAMPLE

INDICATOR_VEC = c(1,1,1)

out_BoB = strMCMC_BoB_mixed(data_dis,iterations,thin_factor,q_resample, INDICATOR_VEC)





