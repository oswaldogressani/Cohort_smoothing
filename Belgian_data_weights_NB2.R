# Author: Yannick Vandendijck (based on programs of Giancarlo Camarda)
# Paper: 
# File last updated on 31/10/2022
#_____________________________________________________________________

rm(list=ls(all=TRUE))

## R-Libraries ----
.libPaths(c("C:/Users/YVanden2/R-4.1.2", .libPaths())) # directory where libraries are located
library(fields)
library(Matrix)
library(clusterPower)
library(graphics)
library(inline)
library(mvtnorm)
library(svcm)
library(Rcpp)
library(RcppArmadillo)



dir = "C:/Users/YVanden2/OneDrive - JNJ/UHasselt/Social Contact Rates/SCR_paper1/Article R code_cleversearch/"
source(paste(dir,'SourceCode_ReadInData_with_weights.R',sep=""))
source(paste(dir,'SourceCode_CreatePenalty.R',sep=""))
source(paste(dir,'SourceCode_PIRWLS_Transformations.R',sep=""))
source(paste(dir,'SourceCode_PIRWLS_NoTransformation.R',sep=""))
source(paste(dir,'SourceCode_TransformedData.R',sep=""))
source(paste(dir,'SourceCode_GridSearchPenalties.R',sep=""))


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


barplot(in_data$tilde_e[lower_index:upper_index], xlab="Age of the respondent", ylab="Number of respondents", 
        ylim=c(0,25), axes=FALSE, cex.lab=1.6, space=0.0, )
axis(side=2, cex.axis=1.6)
axis(side=1, cex.axis=0.8, at=1:77-0.5, tick=TRUE, labels=0:76)


#________________________________________________________________________________________________________________#
#________________________________________________________________________________________________________________#
#________________________________________________________________________________________________________________#
#________________________________________________________________________________________________________________#
## Penalty matrices ----

pen_method2_no_kink = Penalty_Transformed_Dataset(m = m, n = n, kink = FALSE)
pen_method2_kink = Penalty_Transformed_Dataset(m = m, n = n, kink = TRUE, max.kink.age=30)
pen_method3_no_kink = Penalty_Transformed_Penalty(m = m, n = n, kink = FALSE)
pen_method3_kink = Penalty_Transformed_Penalty(m = m, n = n, kink = TRUE, max.kink.age=30)


#________________________________________________________________________________________________________________#
#________________________________________________________________________________________________________________#
#________________________________________________________________________________________________________________#
#________________________________________________________________________________________________________________#


## Model M0 - NB2 distribution ----
##__________________________________________________________________________

## grid search
# grid.l1 = grid.l2 = c(1, 5, 10, 50, 100, 500, 1000, 5000, 10000)
# grid.l1 = grid.l2 = c(1, 5, 10, 50, 100, 500)
# grid.search.neg.bin2.M0 = Grid_Search_No_Cohort_Smoothing(m=m, n=n, Y=Y, E=E, p=p, grid.l1=grid.l1, grid.l2=grid.l2,
#                                                                dist="Neg.Bin2")
# grid.search.neg.bin2.M0$out
#write.table(grid.search.neg.bin2.M0$out, "gridsearch1_NB2_M0.txt", row.names = F)


# output.cs.1 = cleverSearch_No_Cohort_Smoothing(m=m, n=n, Y=Y, E=E, p=p, dist="Neg.Bin2", 
#                                                lambda.x=5, lambda.y=5, length.grid=101)
# output.cs.1


# optimal model
t0 <- Sys.time()
fit0.neg.bin2.M0 = No_Cohort_Smoothing_No_Symmetry(m=m, n=n, Y=Y, E=E,
                                                   lambda1=4.560054, lambda2=4.560054, dist="Neg.Bin2",
                                                   max.iter=25, max.iter.phi=25,
                                                   fix.phi=FALSE, phi.input=1.0)
fit1.neg.bin2.M0 = No_Cohort_Smoothing_Symmetry(m=m, n=n, Y=Y, E=E, p=p,
                                                lambda1=4.560054, lambda2=4.560054, dist="Neg.Bin2",
                                                start.values=fit0.neg.bin2.M0$eta.start,
                                                max.iter=100, max.iter.phi=100,
                                                phi.input=1.0, fix.phi=FALSE, get.V.matrix=FALSE)
t1 <- Sys.time()
time <- t1-t0
time*60


str(fit1.neg.bin2.M0)
fit1.neg.bin2.M0$edf
fit1.neg.bin2.M0$ll
fit1.neg.bin2.M0$aic
fit1.neg.bin2.M0$bic



## Model M2 - NB2 distribution ----
##____________________________________________________________________

## grid search
# grid.l1 = grid.l2 = c(1, 5, 10, 50, 100, 500, 1000, 5000, 10000)
# grid.search.neg.bin2.M2 = Grid_Search_Cohort_Smoothing(m=m, n=n, Y=Y, E=E, p=p, PEN=pen_method3_no_kink,
#                                                        grid.l1=grid.l1,grid.l2=grid.l2, dist="Neg.Bin2")
# grid.search.neg.bin2.M2$out
#write.table(grid.search.neg.bin2.M2$out, "gridsearch1_NB2_M2.txt", row.names = F)


# output.cs.1 = cleverSearch_Cohort_Smoothing(m=m, n=n, Y=Y, E=E, p=p,
#                                             PEN=pen_method3_no_kink, dist="Neg.Bin2",
#                                             lambda.x=5, lambda.y=5, length.grid=101)
# output.cs.1


# optimal model: lambda1=4.158819, lambda2=5.482391, AIC=21107.61, BIC=23235.0

t0 <- Sys.time()
fit0.neg.bin2.M2 = Cohort_Smoothing_No_Symmetry(m=m, n=n, Y=Y, E=E, input.Pen = pen_method3_no_kink,
                                                lambda1=4.158819, lambda2=5.482391, dist="Neg.Bin2",
                                                max.iter=25, max.iter.phi=25,
                                                fix.phi=FALSE, phi.input=1.0)
fit1.neg.bin2.M2 = Cohort_Smoothing_Symmetry(m=m, n=n, Y=Y, E=E, p=p, input.Pen = pen_method3_no_kink,
                                             lambda1=4.158819, lambda2=5.482391, dist="Neg.Bin2",
                                             start.values=fit0.neg.bin2.M2$eta.start,
                                             max.iter=100, max.iter.phi=100,
                                             phi.input=1.0, fix.phi=FALSE, get.V.matrix=FALSE)

t1 <- Sys.time()
time <- t1-t0
time*60

str(fit1.neg.bin2.M2)
fit1.neg.bin2.M2$edf
fit1.neg.bin2.M2$ll
fit1.neg.bin2.M2$aic
fit1.neg.bin2.M2$bic



## Model M2 kink - NB2 distribution ----
##____________________________________________________________________

## grid search
# grid.l1 = grid.l2 = c(1, 5, 10, 50, 100, 500, 1000, 5000, 10000)
# grid.search.neg.bin2.M2.kink = Grid_Search_Cohort_Smoothing(m=m, n=n, Y=Y, E=E, p=p, PEN=pen_method3_kink,
#                                                        grid.l1=grid.l1,grid.l2=grid.l2, dist="Neg.Bin2")
# grid.search.neg.bin2.M2.kink$out
#write.table(grid.search.neg.bin2.M2.kink$out, "gridsearch1_NB2_M2_kink.txt", row.names = F)


# output.cs.1 = cleverSearch_Cohort_Smoothing(m=m, n=n, Y=Y, E=E, p=p, PEN=pen_method3_kink, dist="Neg.Bin2",
#                                             lambda.x=5, lambda.y=5, length.grid=101)
# output.cs.1


# optimal model: lambda1=4.560054, lambda2=5.482391, AIC=21101.2, BIC=
t0 <- Sys.time()
fit0.neg.bin2.M2.kink = Cohort_Smoothing_No_Symmetry(m=m, n=n, Y=Y, E=E, input.Pen = pen_method3_kink,
                                                     lambda1=4.560054, lambda2=5.482391, dist="Neg.Bin2",
                                                     max.iter=25, max.iter.phi=25,
                                                     fix.phi=FALSE, phi.input=1.0)
fit1.neg.bin2.M2.kink = Cohort_Smoothing_Symmetry(m=m, n=n, Y=Y, E=E, p=p, input.Pen = pen_method3_kink,
                                                  lambda1=4.560054, lambda2=5.482391, dist="Neg.Bin2",
                                                  start.values=fit0.neg.bin2.M2.kink$eta.start,
                                                  max.iter=100, max.iter.phi=100,
                                                  phi.input=1.0, fix.phi=FALSE, get.V.matrix=FALSE)
t1 <- Sys.time()
time <- t1-t0
time*60

str(fit1.neg.bin2.M2.kink)
fit1.neg.bin2.M2.kink$edf
fit1.neg.bin2.M2.kink$ll
fit1.neg.bin2.M2.kink$aic
fit1.neg.bin2.M2.kink$bic
fit1.neg.bin2.M2.kink$phi



## Model M1 kink - NB2 distribution ----
##____________________________________________________________________

# ## grid search
# grid.l1 = grid.l2 = c(1, 5, 10, 50, 100, 500, 1000, 5000, 10000)
# grid.search.neg.bin2.M1.kink = Grid_Search_Cohort_Smoothing(m=m, n=n, Y=Y, E=E, p=p, PEN=pen_method2_kink, 
#                                                             grid.l1=grid.l1, grid.l2=grid.l2, dist="Neg.Bin2")
# 
# grid.search.neg.bin2.M1.kink$out
# #write.table(grid.search.neg.bin2.M1.kink$out, "gridsearch1_NB2_M1_kink.txt", row.names = F)


# output.cs.1 = cleverSearch_Cohort_Smoothing(m=m, n=n, Y=Y, E=E, p=p, PEN=pen_method2_kink, dist="Neg.Bin2",
#                                             lambda.x=10, lambda.y=10, length.grid=101)
# output.cs.1


# optimal model: lambda1=4.365158, lambda2=5.248075, AIC=21081.79, BIC=
t0 <- Sys.time()
fit0.neg.bin2.M1.kink = Cohort_Smoothing_No_Symmetry(m=m, n=n, Y=Y, E=E, input.Pen=pen_method2_kink,
                                                     lambda1=4.365158, lambda2=5.248075, dist="Neg.Bin2", 
                                                     max.iter=25, max.iter.phi=25)
fit1.neg.bin2.M1.kink = Cohort_Smoothing_Symmetry(m=m, n=n, Y=Y, E=E, p=p, input.Pen=pen_method2_kink,
                                                  lambda1=4.365158, lambda2=5.248075, 
                                                  start.values=fit0.neg.bin2.M1.kink$eta.start, phi.input=fit0.neg.bin2.M1.kink$phi,
                                                  max.iter=100, max.iter.phi=100, dist="Neg.Bin2",
                                                  get.V.matrix=FALSE, fix.phi=FALSE)
t1 <- Sys.time()
time <- t1-t0
time*60


str(fit1.neg.bin2.M1.kink)
fit1.neg.bin2.M1.kink$edf
fit1.neg.bin2.M1.kink$ll
fit1.neg.bin2.M1.kink$aic
fit1.neg.bin2.M1.kink$bic
fit1.neg.bin2.M1.kink$phi




## Model M1 - NB2 distribution ----
##____________________________________________________________________

# ## grid search
# grid.l1 = grid.l2 = c(1, 5, 10, 50, 100, 500, 1000, 5000, 10000)
# grid.search.neg.bin2.M1.kink = Grid_Search_Cohort_Smoothing(m=m, n=n, Y=Y, E=E, p=p, PEN=pen_method2_kink, 
#                                                             grid.l1=grid.l1, grid.l2=grid.l2, dist="Neg.Bin2")
# 
# grid.search.neg.bin2.M1.kink$out
# #write.table(grid.search.neg.bin2.M1.kink$out, "gridsearch1_NB2_M1_kink.txt", row.names = F)


# output.cs.1 = cleverSearch_Cohort_Smoothing(m=m, n=n, Y=Y, E=E, p=p, PEN=pen_method2_no_kink, 
#                                             dist="Neg.Bin2",
#                                             lambda.x=5, lambda.y=5, length.grid=101)
# output.cs.1


 
# optimal model: lambda1=4.158819, lambda2=, AIC=21087.??, BIC=
t0 <- Sys.time()
fit0.neg.bin2.M1 = Cohort_Smoothing_No_Symmetry(m=m, n=n, Y=Y, E=E, input.Pen=pen_method2_no_kink,
                                                lambda1=4.158819, lambda2=5.000000, dist="Neg.Bin2", 
                                                max.iter=25, max.iter.phi=25)
fit1.neg.bin2.M1 = Cohort_Smoothing_Symmetry(m=m, n=n, Y=Y, E=E, p=p, input.Pen=pen_method2_no_kink,
                                             lambda1=4.158819, lambda2=5.000000,
                                             start.values=fit0.neg.bin2.M1$eta.start, 
                                             phi.input=fit0.neg.bin2.M1$phi,max.iter=100, 
                                             max.iter.phi=100, dist="Neg.Bin2",
                                             get.V.matrix=FALSE, fix.phi=FALSE)

t1 <- Sys.time()
time <- t1-t0
time*60

str(fit1.neg.bin2.M1)
fit1.neg.bin2.M1$edf
fit1.neg.bin2.M1$ll
fit1.neg.bin2.M1$aic
fit1.neg.bin2.M1$bic
fit1.neg.bin2.M1$phi






#________________________________________________________________________________________________________________#
#________________________________________________________________________________________________________________#
#________________________________________________________________________________________________________________#
#________________________________________________________________________________________________________________#
#________________________________________________________________________________________________________________#
#________________________________________________________________________________________________________________#

    
M0 = fit1.neg.bin2.M0
M2 = fit1.neg.bin2.M2.kink
str(M0)
str(M2)


### Matrix plot ----
plot.x1 = (lower_index:upper_index) - 1
plot.x2 = (lower_index:upper_index) - 1

image.plot(plot.x1, plot.x2, M0$est.ETA.S, col=topo.colors(19), zlim=c(-4.30, 1.10),
           xlab="Age of the respondent", ylab="Age of the contact",
           cex.lab=1.6, cex.axis=1.6)

image.plot(plot.x1, plot.x2, M2$est.ETA.S, col=topo.colors(19), zlim=c(-4.30, 1.10),
           xlab="Age of the respondent", ylab="Age of the contact",
           cex.lab=1.6, cex.axis=1.6)


### Matrix plot 2 ----
br = quantile(M2$est.GAMMA.S * P, seq(0, 1, length=20))
br[1] = 1000

image.plot(plot.x1-1, plot.x2-1, M0$est.GAMMA.S * P, col=topo.colors(19),
           xlab="Age of the respondent", ylab="Age of the contact", breaks=br,
           cex.lab=1.6, cex.axis=1.6)

image.plot(plot.x1-1, plot.x2-1, M2$est.GAMMA.S * P, col=topo.colors(19),
           xlab="Age of the respondent", ylab="Age of the contact", breaks=br,
           cex.lab=1.6, cex.axis=1.6)


### Diagonal plot ----
plot(diag(M0$est.ETA.S), type="l", lwd=3.5, ylim=c(-3.2, 2.0), 
     xlab="Index main diagonal element", ylab="Estimated log contact rate", cex.lab=1.6, cex.axis=1.6)
points(diag(ETA), cex=diag(E)/10, pch=16)
lines(diag(M2$est.ETA.S), lwd=3.5, col=2)
legend(x=45, y=2.4, legend=c("w/o kink","w kink"), bty="n", col=1:2, lty=1, cex=1.6)


### Other ages ----
plot(M0$est.ETA.S[17,], type="l", lwd=3.5, ylim=c(-4.0, 1.20),
     main="Respondent of age 16 years",
     xlab="Age of the contact", ylab="Estimated log contact rate", cex.lab=1.6, cex.axis=1.6)
points(ETA[17,], cex=diag(E)/10, pch=16)
lines(M2$est.ETA.S[17,], lwd=2.5, col=2, lty=2)
legend(x=45, y=1.2, legend=c("w/o kink","w kink"), bty="n", col=1:2, lty=1:2, cex=1.6, lwd=3)


plot(M0$est.ETA.S[,17], type="l", lwd=3.5, ylim=c(-4.0, 1.20),
     main="Contact of age 16 years",
     xlab="Age of the contact", ylab="Estimated log contact rate", cex.lab=1.6, cex.axis=1.6)
points(ETA[,17], cex=diag(E)/10, pch=16)
lines(M2$est.ETA.S[,17], lwd=2.5, col=2, lty=2)
legend(x=45, y=1.2, legend=c("w/o kink","w kink"), bty="n", col=1:2, lty=1:2, cex=1.6, lwd=3)


plot(M0$est.ETA.S[51,], type="l", lwd=3.5, ylim=c(-4.0, 1.20),
     main="Respondent of age 50 years",
     xlab="Age of the contact", ylab="Estimated log contact rate", cex.lab=1.6, cex.axis=1.6)
points(ETA[51,], cex=diag(E)/10, pch=16)
lines(M2$est.ETA.S[51,], lwd=2.5, col=2, lty=2)
legend(x=45, y=1.2, legend=c("w/o kink","w kink"), bty="n", col=1:2, lty=1:2, cex=1.6, lwd=3)
