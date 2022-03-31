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
library(svcm)

source('SourceCode_ReadInData_with_weights.R')
source('SourceCode_CreatePenalty.R')
source('SourceCode_PIRWLS_Transformations.R')
source('SourceCode_PIRWLS_NoTransformation.R')
source('SourceCode_TransformedData.R')
source('SourceCode_MultivariateNormal.R')
source('SourceCode_GridSearchPenalties.R')

## Social contact data & demographic data ----
##____________________________________________
in_data = Load_Social_Contact_Data()
str(in_data)


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

range(GAMMA,finite=TRUE)
image.plot(plot.x1, plot.x2, GAMMA, col=topo.colors(20), zlim=c(0,11),
	xlab="Age of the respondent", ylab="Age of the contact",
	cex.lab=1.6, cex.axis=1.6)


barplot(in_data$tilde_e[lower_index:upper_index], 
        xlab="Age of the respondent",
        ylab="Number of respondents", 
        ylim=c(0,25), axes=FALSE, cex.lab=1.6, space=0.0, xaxt="n")
axis(side=2, cex.axis=1.6)
axis(side=1, cex.axis=0.8, at=1:77-0.5, tick=TRUE, labels=0:76)


## All analysis ----

pen_method2_no_kink = Penalty_Transformed_Dataset(m = m, n = n, kink = FALSE)
pen_method2_kink = Penalty_Transformed_Dataset(m = m, n = n, kink = TRUE, max.kink.age=30)
pen_method3_no_kink = Penalty_Transformed_Penalty(m = m, n = n, kink = FALSE)
pen_method3_kink = Penalty_Transformed_Penalty(m = m, n = n, kink = TRUE, max.kink.age=30)


## Analsis with no smoothing over cohorts - Poisson distribution
##______________________________________________________________
# t0 <- Sys.time()  
# output.cs.1 = cleverSearch_No_Cohort_Smoothing(m=m, n=n, Y=Y, E=E, p=p, dist="Poisson", lambda.x=0.5, lambda.y=0.5, length.grid=51)
# t1 <- Sys.time()
# time <- t1-t0
# time
# output.cs.1

# optimal model - lambda1=0.4560054, lambda2=0.4560054, AIC=24110.68 (weights)
fit0_poisson_method1 = No_Cohort_Smoothing_No_Symmetry(m=m, n=n, Y=Y, E=E, lambda1=0.4560054, lambda2=0.4560054, 
                                                       max.iter = 25, dist="Poisson")
fit1_poisson_method1 = No_Cohort_Smoothing_Symmetry(m=m, n=n, Y=Y, E=E, p=p, lambda1=0.4560054, lambda2=0.4560054, 
                                                    start.values=fit0_poisson_method1$eta.start,
                                                    max.iter=100, get.V.matrix=FALSE, dist="Poisson")
str(fit1_poisson_method1)


## Analsis with transformed data set - Poisson distribution
##_________________________________________________________
# t0 <- Sys.time()  
# output.cs.2 = cleverSearch_Cohort_Smoothing(m=m, n=n, Y=Y, E=E, p=p, PEN=pen_method2_no_kink, dist="Poisson", 
#                                            lambda.x=0.5, lambda.y=0.5, length.grid=51)
# t1 <- Sys.time()
# time <- t1-t0
# time
# output.cs.2

	# optimal model: lambda1=0.50, lambda2=0.3154787 AIC=24489.58 (weights)
fit0.poisson.method2 = Cohort_Smoothing_No_Symmetry(m=m, n=n, Y=Y, E=E, input.Pen=pen_method2_no_kink, 
	                                                  lambda1=0.50, lambda2=0.3154787, max.iter=25, dist="Poisson")
fit1.poisson.method2 = Cohort_Smoothing_Symmetry(m=m, n=n, Y=Y, E=E, p=p, input.Pen=pen_method2_no_kink, 
	                                              lambda1=0.50, lambda2=0.3154787, start.values=fit0.poisson.method2$eta.start, 
	                                              max.iter=100, get.V.matrix=FALSE, dist="Poisson")
str(fit1.poisson.method2)

## Analsis with transformed penalty - Poisson distribution
##________________________________________________________
# t0 <- Sys.time()  
# output.cs.3 = cleverSearch_Cohort_Smoothing(m=m, n=n, Y=Y, E=E, p=p, PEN=pen_method3_no_kink, dist="Poisson", 
#                                             lambda.x=0.5, lambda.y=0.5, length.grid=51)
# t1 <- Sys.time()
# time <- t1-t0
# time
# output.cs.3

	# optimal model: lambda1=0.5, lambda2=0.3154787, AIC=24537.65 (weights)
fit0.poisson.method3 = Cohort_Smoothing_No_Symmetry(m=m, n=n, Y=Y, E=E, input.Pen=pen_method3_no_kink,
	                                                  lambda1=0.5, lambda2=0.3154787, max.iter=25, dist="Poisson")
fit1.poisson.method3 = Cohort_Smoothing_Symmetry(m=m, n=n, Y=Y, E=E, p=p, input.Pen=pen_method3_no_kink,
	                                                  lambda1=0.5, lambda2=0.3154787, start.values=fit0.poisson.method3$eta.start, max.iter=100, get.V.matrix=FALSE, dist="Poisson")
str(fit1.poisson.method3)


#________________________________________________________________________________________________________________#
#________________________________________________________________________________________________________________#
#________________________________________________________________________________________________________________#
#________________________________________________________________________________________________________________#


## Analsis with no smoothing over cohorts - Negative Binomial 1 distribution
##__________________________________________________________________________
# t0 <- Sys.time()  
# output.cs.4 = cleverSearch_No_Cohort_Smoothing(m=m, n=n, Y=Y, E=E, p=p, dist="Neg.Bin1", lambda.x=20, lambda.y=20, length.grid=51)
# t1 <- Sys.time()
# time <- t1-t0
# time
# output.cs.4

	# optimal model - lambda1=15.17155, lambda=16.63528, AIC=21097.68 (weights)
fit0.neg.bin1.method1 = No_Cohort_Smoothing_No_Symmetry(m=m, n=n, Y=Y, E=E, lambda1=15.17155, lambda2=16.63528,  
                                                      max.iter=25, max.iter.phi=25, dist="Neg.Bin1")
fit1.neg.bin1.method1 = No_Cohort_Smoothing_Symmetry(m=m, n=n, Y=Y, E=E, p=p,	lambda1=15.17155, lambda2=16.63528, 
                                                     start.values=fit0.neg.bin1.method1$eta.start, max.iter=100, max.iter.phi=100,
                                                     dist="Neg.Bin1", phi.input=fit0.neg.bin1.method1$phi, 
                                                     get.V.matrix=FALSE, fix.phi=FALSE)
str(fit1.neg.bin1.method1)


## Analsis with transformed data set - Negative Binomial 1 distribution
##_____________________________________________________________________
# t0 <- Sys.time()  
# output.cs.5 = cleverSearch_Cohort_Smoothing(m=m, n=n, Y=Y, E=E, p=p, PEN=pen_method2_no_kink, dist="Neg.Bin1", 
#                                                lambda.x=30, lambda.y=900, length.grid=51)
# t1 <- Sys.time()
# time <- t1-t0
# time
# output.cs.5

	# optimal model - lambda1= 22.75733, lambda2=1714.91465, AIC=21102.05 (weights)
fit0.neg.bin1.method2 = Cohort_Smoothing_No_Symmetry(m=m, n=n, Y=Y, E=E, input.Pen=pen_method2_no_kink,
                                                     lambda1= 22.75733, lambda2=1714.91465, max.iter=25, max.iter.phi=25, dist="Neg.Bin1")
fit1.neg.bin1.method2 = Cohort_Smoothing_Symmetry(m=m, n=n, Y=Y, E=E, p=p, input.Pen=pen_method2_no_kink,
                                                  lambda1= 22.75733, lambda2=1714.91465, start.values=fit0.neg.bin1.method2$eta.start, max.iter=100, max.iter.phi=100,
                                                  dist="Neg.Bin1", phi.input=fit0.neg.bin1.method2$phi, 
                                                  get.V.matrix=FALSE, fix.phi=FALSE)
str(fit1.neg.bin1.method2)

## Analsis with transformed penalty - Negative Binomial 1 distribution
##____________________________________________________________________
# t0 <- Sys.time()  
# output.cs.6 = cleverSearch_Cohort_Smoothing(m=m, n=n, Y=Y, E=E, p=p, PEN=pen_method3_no_kink, dist="Neg.Bin1", 
#                                             lambda.x=30, lambda.y=900, length.grid=51)
# t1 <- Sys.time()
# time <- t1-t0
# time
# output.cs.6

	# optimal model: lambda1=27.36033, lambda2=1564.02075, AIC=21115.15
fit0.neg.bin1.method3 = Cohort_Smoothing_No_Symmetry(m=m, n=n, Y=Y, E=E, input.Pen=pen_method3_no_kink,
                                                     lambda1=27.36033, lambda2=1564.02075, max.iter=25, max.iter.phi=25, dist="Neg.Bin1")
fit1.neg.bin1.method3 = Cohort_Smoothing_Symmetry(m=m, n=n, Y=Y, E=E, p=p, input.Pen=pen_method3_no_kink,
                                                  lambda1=27.36033, lambda2=1564.02075, start.values=fit0.neg.bin1.method3$eta.start, max.iter=100, max.iter.phi=100,
                                                  dist="Neg.Bin1", phi.input=fit0.neg.bin1.method3$phi, get.V.matrix=FALSE, fix.phi=FALSE)
str(fit1.neg.bin1.method3)


#________________________________________________________________________________________________________________#
#________________________________________________________________________________________________________________#
#________________________________________________________________________________________________________________#
#________________________________________________________________________________________________________________#


## Analsis with transformed data set - Negative Binomial 1 distribution - with Kink
##_________________________________________________________________________________
# t0 <- Sys.time()  
# output.cs.7 = cleverSearch_Cohort_Smoothing(m=m, n=n, Y=Y, E=E, p=p, PEN=pen_method2_kink, dist="Neg.Bin1", 
#                                             lambda.x=30, lambda.y=1000, length.grid=51)
# t1 <- Sys.time()
# time <- t1-t0
# time
# output.cs.7

# optimal model - lambda1=30, lambda=1584.893, AIC=21078.19
fit0.neg.bin1.method2.Kink = Cohort_Smoothing_No_Symmetry(m=m, n=n, Y=Y, E=E, input.Pen=pen_method2_kink,
                                                     lambda1=30, lambda2=1584.893, max.iter=25, max.iter.phi=25, dist="Neg.Bin1")
fit1.neg.bin1.method2.Kink = Cohort_Smoothing_Symmetry(m=m, n=n, Y=Y, E=E, p=p, input.Pen=pen_method2_kink,
                                                  lambda1=30, lambda2=1584.893, start.values=fit0.neg.bin1.method2.Kink$eta.start, 
                                                  max.iter=100, max.iter.phi=100, dist="Neg.Bin1", 
                                                  phi.input=fit0.neg.bin1.method2.Kink$phi, get.V.matrix=FALSE, fix.phi=FALSE)
str(fit1.neg.bin1.method2.Kink)

## Analsis with transformed penalty - Negative Binomial 1 distribution - with Kink
##________________________________________________________________________________
# t0 <- Sys.time()  
# output.cs.8 = cleverSearch_Cohort_Smoothing(m=m, n=n, Y=Y, E=E, p=p, PEN=pen_method3_kink, dist="Neg.Bin1", 
#                                             lambda.x=40, lambda.y=1000, length.grid=51)
# t1 <- Sys.time()
# time <- t1-t0
# time
# output.cs.8

# optimal model: lambda1=40, lambda2=1584.893, AIC=21086.02
fit0.neg.bin1.method3.Kink = Cohort_Smoothing_No_Symmetry(m=m, n=n, Y=Y, E=E, input.Pen=pen_method3_kink,
                                                     lambda1=40, lambda2=1584.893, max.iter=25, max.iter.phi=25, dist="Neg.Bin1")
fit1.neg.bin1.method3.Kink = Cohort_Smoothing_Symmetry(m=m, n=n, Y=Y, E=E, p=p, input.Pen=pen_method3_kink,
                                                  lambda1=40, lambda2=1584.893, start.values=fit0.neg.bin1.method3.Kink$eta.start, 
                                                  max.iter=100, max.iter.phi=100, dist="Neg.Bin1",
                                                  phi.input=fit0.neg.bin1.method3.Kink$phi, get.V.matrix=FALSE, fix.phi=FALSE)
str(fit1.neg.bin1.method3.Kink)


#-------------------------- Summarize results in a Table (Table 3 in main manuscript)

Table3 <- matrix(0,nrow = 8, ncol = 7)
colnames(Table3) <- c("Distribution", "lambda1", "lambda2", "ED", "-2logL","AIC", "phi")
rownames(Table3) <- c("M0-Poiss", "M0-NB", "M1-Poiss", "M1-NB", "M2-Poiss", "M2-NB", "M1-NB-Kink", "M2-NB-Kink")
Table3[,1] <- c("Poisson", "NegBin","Poisson","NegBin","Poisson","NegBin","NegBin","NegBin")
Table3[,2] <- c(0.46,15.17,0.50,22.76,0.50,27.36,30.00,40.00)
Table3[,3] <- c(0.46,16.64,0.32,1714.91,0.32,1564.02,1584.89,1584.89)
Table3[,4] <- round(c(fit1_poisson_method1$edf,
                fit1.neg.bin1.method1$edf,
                fit1.poisson.method2$edf,
                fit1.neg.bin1.method2$edf,
                fit1.poisson.method3$edf,
                fit1.neg.bin1.method3$edf,
                fit1.neg.bin1.method2.Kink$edf,
                fit1.neg.bin1.method3.Kink$edf
                ),2)
Table3[,5] <- round(c(fit1_poisson_method1$ll,
                      fit1.neg.bin1.method1$ll,
                      fit1.poisson.method2$ll,
                      fit1.neg.bin1.method2$ll,
                      fit1.poisson.method3$ll,
                      fit1.neg.bin1.method3$ll,
                      fit1.neg.bin1.method2.Kink$ll,
                      fit1.neg.bin1.method3.Kink$ll
                      ),2)
Table3[,6] <- round(c(fit1_poisson_method1$aic,
                      fit1.neg.bin1.method1$aic,
                      fit1.poisson.method2$aic,
                      fit1.neg.bin1.method2$aic,
                      fit1.poisson.method3$aic,
                      fit1.neg.bin1.method3$aic,
                      fit1.neg.bin1.method2.Kink$aic,
                      fit1.neg.bin1.method3.Kink$aic
                      ),2)
Table3[,7] <- round(c(NA,fit1.neg.bin1.method1$phi,NA,fit1.neg.bin1.method2$phi,NA,
                fit1.neg.bin1.method3$phi,fit1.neg.bin1.method2.Kink$phi,
                fit1.neg.bin1.method3.Kink$phi),2)


## Figure 4 ----------------------------------------------------------------

plot.x1 = (lower_index:upper_index) - 1
plot.x2 = (lower_index:upper_index) - 1

range(fit1.neg.bin1.method1$est.ETA.S)
range(fit1.neg.bin1.method2$est.ETA.S)
range(fit1.neg.bin1.method3$est.ETA.S)

range(fit1.neg.bin1.method1$est.GAMMA.S * P)
br <- quantile(fit1.neg.bin1.method1$est.GAMMA.S * P, seq(0, 1, length=20))
br[1] <- 1450

# For M0 (ETA)
par(mfrow=c(2,3))
image.plot(plot.x1, plot.x2, fit1.neg.bin1.method1$est.ETA.S, col=topo.colors(19), zlim=c(-4.30, 1.10),
           xlab="Age of the respondent", ylab="Age of the contact",
           cex.lab=1.6, cex.axis=1.6)

# For M1 (ETA)
image.plot(plot.x1, plot.x2, fit1.neg.bin1.method2$est.ETA.S, col=topo.colors(19), zlim=c(-4.30, 1.10),
           xlab="Age of the respondent", ylab="Age of the contact",
           cex.lab=1.6, cex.axis=1.6)

# For M2 (ETA)
image.plot(plot.x1, plot.x2, fit1.neg.bin1.method3$est.ETA.S, col=topo.colors(19), zlim=c(-4.30, 1.10),
           xlab="Age of the respondent", ylab="Age of the contact",
           cex.lab=1.6, cex.axis=1.6)

# For M0 (Gamma * P)
image.plot(plot.x1-1, plot.x2-1, fit1.neg.bin1.method1$est.GAMMA.S * P, col=topo.colors(19), breaks = br,
           xlab="Age of the respondent", ylab="Age of the contact",
           cex.lab=1.6, cex.axis=1.6)

# For M1 (Gamma * P)
image.plot(plot.x1-1, plot.x2-1, fit1.neg.bin1.method2$est.GAMMA.S * P, col=topo.colors(19), breaks = br,
           xlab="Age of the respondent", ylab="Age of the contact",
           cex.lab=1.6, cex.axis=1.6)

# For M2 (Gamma *P)
image.plot(plot.x1-1, plot.x2-1, fit1.neg.bin1.method3$est.GAMMA.S * P, col=topo.colors(19), breaks = br,
           xlab="Age of the respondent", ylab="Age of the contact",
           cex.lab=1.6, cex.axis=1.6)
# dev.off()



# Figure 5 ------------------------------------------------------------

par(mfrow = c(1,1))
# For M2
range(fit1.neg.bin1.method3.Kink$est.ETA.S)
br <- quantile(fit1.neg.bin1.method2.Kink$est.GAMMA.S * P, seq(0, 1, length=20))
br[1] = 1450

image.plot(plot.x1, plot.x2, fit1.neg.bin1.method3.Kink$est.ETA.S, col=topo.colors(19),
           xlab="Age of the respondent", ylab="Age of the contact",  zlim=c(-4.30, 1.10),
           cex.lab=1.6, cex.axis=1.6)

image.plot(plot.x1-1, plot.x2-1, fit1.neg.bin1.method3.Kink$est.GAMMA.S * P, col=topo.colors(19),
           xlab="Age of the respondent", ylab="Age of the contact", breaks=br,
           cex.lab=1.6, cex.axis=1.6)

plot(diag(fit1.neg.bin1.method3$est.ETA.S), type="l", lwd=3.5, ylim=c(-3.2, 2.0), 
    xlab="Index main diagonal element", ylab="Estimated log contact rate", cex.lab=1.6, cex.axis=1.6)
points(diag(ETA), cex=diag(E)/10, pch=16)
lines(diag(fit1.neg.bin1.method3.Kink$est.ETA.S), lwd=3.5, col=2)
legend(x=35, y=2.4, legend=c("w/o kink","w kink"), bty="n", col=1:2, lty=1, cex=1.6)



as.data.frame(Table3)


