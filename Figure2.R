# Authors: Yannick Vandendijck and Giancarlo Camarda (updated by Oswaldo Gressani)
# File last updated on 01/03/2022
#_____________________________________________________________________

rm(list=ls(all=TRUE))

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


## Social contact data & demographic data ----
##____________________________________________
in_data = Load_Social_Contact_Data()
dev.off()
#str(in_data)


## Preparing some matrices and vectors
max_age = 76
lower_index = 1
upper_index = max_age + 1

m = n = upper_index
One = matrix(1, ncol=m)
Y = in_data$mat_cont[lower_index:upper_index, lower_index:upper_index]
E = in_data$tilde_e[lower_index:upper_index] %*% One
p = in_data$P[lower_index:upper_index]
P = p %*% One

GAMMA = Y/E
ETA = log(GAMMA)
Z = GAMMA*P

par(mfrow=c(1,1))

barplot(in_data$tilde_e[lower_index:upper_index], 
        xlab="Age of the respondent",
        ylab="Number of respondents", 
        ylim=c(0,25), axes=FALSE, cex.lab=1.6, space=0.0, xaxt="n")
axis(side=2, cex.axis=1.6)
axis(side=1, cex.axis=0.8, at=1:77-0.5, tick=TRUE, labels=0:76)


rm(list=ls(all=TRUE))
source('SourceCode_ReadInData_with_weights.R')

## Social contact data & demographic data ----
##____________________________________________
in_data = Load_Social_Contact_Data()
# str(in_data)


## Preparing some matrices and vectors
max_age = 76
lower_index = 1
upper_index = max_age + 1

m = n = upper_index
One = matrix(1, ncol=m)
Y = in_data$mat_cont[lower_index:upper_index, lower_index:upper_index]
E = in_data$tilde_e[lower_index:upper_index] %*% One
p = in_data$P[lower_index:upper_index]
P = p %*% One

GAMMA = Y/E
ETA = log(GAMMA)
Z = GAMMA*P

range(ETA,finite=TRUE)
plot.x1 = (lower_index:upper_index) -1
plot.x2 = (lower_index:upper_index) -1
image.plot(plot.x1, plot.x2, ETA, col=topo.colors(20), zlim=c(-3.80,2.40),
           xlab="Age of the respondent", ylab="Age of the contact",
           cex.lab=1.6, cex.axis=1.6)


