# Authors: Yannick Vandendijck and Giancarlo Camarda (updated by Oswaldo Gressani)
# File last updated on 01/03/2022
#_____________________________________________________________________
# Simulation of Poisson data without kink


rm(list=ls(all=TRUE))

sim.N = 100
seed.in = 2017
kink=FALSE

## R-Libraries ----
library(fields)
library(Matrix)
library(clusterPower)
library(graphics)
library(Rcpp)
library(RcppArmadillo)
library(inline)
library(mvtnorm)


source('SourceCode_ReadInData.R')
source('SourceCode_CreatePenalty.R')
source('SourceCode_PIRWLS_Transformations.R')
source('SourceCode_PIRWLS_NoTransformation.R')
source('SourceCode_GridSearchPenalties.R')
source('SourceCode_TransformedData.R')
source('SourceCode_MultivariateNormal.R')


## Social contact data & demographic data ----
##____________________________________________
in_data = Load_Social_Contact_Data()
str(in_data)

## Preparing some matrices and vectors
max_age = 76
lower_index = 1
upper_index = max_age + 1

m = n = upper_index
mn = m*n
One = matrix(1, ncol=m)
E = in_data$tilde_e[lower_index:upper_index] %*% One
p = in_data$P[lower_index:upper_index]
P = p %*% One

### load true GAMMA matrix
load("Gamma_matrix.RData")
GAMMA.true = est.gamma.matrix
range(GAMMA.true)

plot(diag(GAMMA.true))
plot(GAMMA.true[10,])

par(mfrow = c(1,1))

# plot of the matrix
plot.x1 = (lower_index:upper_index) - 1
plot.x2 = (lower_index:upper_index) - 1
image.plot(plot.x1, plot.x2, GAMMA.true, col=topo.colors(19), zlim=c(0.00, 0.85),
           xlab="Age of the respondent", ylab="Age of the contact",
           cex.lab=1.6, cex.axis=1.6)


plot.x1 = (lower_index:upper_index) - 1
plot.x2 = (lower_index:upper_index) - 1
image.plot(plot.x1, plot.x2, log(GAMMA.true), col=topo.colors(19), zlim=c(-4.20, 0.00),
           xlab="Age of the respondent", ylab="Age of the contact",
           cex.lab=1.6, cex.axis=1.6)





