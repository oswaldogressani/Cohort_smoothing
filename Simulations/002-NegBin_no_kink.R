# Authors: Yannick Vandendijck and Giancarlo Camarda (updated by Oswaldo Gressani)
# File last updated on 01/03/2022
#_____________________________________________________________________
# Simulation of NegBin data without kink


rm(list=ls(all=TRUE))

sim.N = 100
seed.in = 2017 + (3)*sim.N
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


#_________________________________________________________________________________________________________________#
#_________________________________________________________________________________________________________________#
###---- START SIMULATION ----

pen_method2_no_kink = Penalty_Transformed_Dataset(m = m, n = n, kink = FALSE)
pen_method2_kink = Penalty_Transformed_Dataset(m = m, n = n, kink = TRUE, max.kink.age = 30)
pen_method3_no_kink = Penalty_Transformed_Penalty(m = m, n = n, kink = FALSE)
pen_method3_kink = Penalty_Transformed_Penalty(m = m, n = n, kink = TRUE, max.kink.age = 30)

matrix.conv.ind1 = matrix(NA, nrow=sim.N, ncol=5)
matrix.conv.ind2 = matrix(NA, nrow=sim.N, ncol=5)
matrix.edf = matrix(NA, nrow=sim.N, ncol=5)
matrix.ll = matrix(NA, nrow=sim.N, ncol=5)
matrix.aic = matrix(NA, nrow=sim.N, ncol=5)
matrix.bic = matrix(NA, nrow=sim.N, ncol=5)
matrix.phi = matrix(NA, nrow=sim.N, ncol=5)
estimated.ETAmatrix1 = estimated.ETAmatrix2 = estimated.ETAmatrix3 = estimated.ETAmatrix4 = estimated.ETAmatrix5 = list()
estimated.COVmatrix1 = estimated.COVmatrix2 = estimated.COVmatrix3 = estimated.COVmatrix4 = estimated.COVmatrix5 = list()



for (k in 1:sim.N){
  set.seed(seed.in + k)
	### simulate data
	mu = as.matrix(E*GAMMA.true)
	mu.vec = as.vector(mu)
	sim.contacts.vector = rnbinom(mn, mu=mu.vec, size=mu.vec*(1/2.0))
	sim.mat.cont = matrix(sim.contacts.vector, nrow=m, ncol=n)
	Y <- sim.mat.cont
	y <- as.vector(t(Y))

	### actual rates, log-rates and crude total population contacts
	GAMMA <- Y/E
	ETA <- log(GAMMA)
	Z <- GAMMA*P

	### grid search
	# grid.l1 = grid.l2 = c(50, 100,200,400)
	# grid.l1 = c(20, 50, 75, 100)
	# grid.l2 = c(500,750, 1000, 1500)
	# grid.search = Grid_Search_No_Cohort_Smoothing(m=m, n=n, Y=Y, E=E, p=p, grid.l1=grid.l1, grid.l2=grid.l2,
	#                                              dist="Neg.Bin1", fix.phi=TRUE, phi.input=2.0)
	# grid.search = Grid_Search_Cohort_Smoothing(m=m, n=n, Y=Y, E=E, p=p, grid.l1=grid.l1, grid.l2=grid.l2,fix.phi=TRUE, phi.input=2.0,
	#                                               dist="Neg.Bin1", PEN=pen_method3_kink)
	# z = t(matrix(grid.search$out[,5],5,5))
	# range(z)

	
	## Analsis with no smoothing over cohorts - NegBin1 distribution 100 - 100
	fit0.NegBin1.method1 = No_Cohort_Smoothing_No_Symmetry(m=m, n=n, Y=Y, E=E, 
	                                                       lambda1 = 100, lambda2=100, max.iter = 25, max.iter.phi=25, dist="Neg.Bin1")
	fit1.NegBin1.method1 = No_Cohort_Smoothing_Symmetry(m=m, n=n, Y=Y, E=E, p=p, 
	                                                    lambda1=100, lambda2=100, start.values=fit0.NegBin1.method1$eta.start, phi.input=fit0.NegBin1.method1$phi, 
	                                                    max.iter=100, max.iter.phi=100, fix.phi=FALSE, get.V.matrix=TRUE, dist="Neg.Bin1")
	
	## Analsis with transformed data set - NegBin1 distribution - 75 - 750
	fit0.NegBin1.method2 = Cohort_Smoothing_No_Symmetry(m=m, n=n, Y=Y, E=E, input.Pen=pen_method2_no_kink, 
	                                                    lambda1=75, lambda2=750, max.iter=25, max.iter.phi=25, dist="Neg.Bin1")
	fit1.NegBin1.method2 = Cohort_Smoothing_Symmetry(m=m, n=n, Y=Y, E=E, p=p, input.Pen=pen_method2_no_kink, 
	                                                 lambda1=75, lambda2=750, start.values=fit0.NegBin1.method2$eta.start, phi.input=fit0.NegBin1.method2$phi, 
	                                                 max.iter=100, max.iter.phi=100, fix.phi=FALSE, get.V.matrix=TRUE, dist="Neg.Bin1")

	## Analsis with transformed penalty - NegBin1 distribution - 75- 750
	fit0.NegBin1.method3 = Cohort_Smoothing_No_Symmetry(m=m, n=n, Y=Y, E=E, input.Pen=pen_method3_no_kink,
	                                                    lambda1=75, lambda2=750, max.iter=25, max.iter.phi=25, dist="Neg.Bin1")
	fit1.NegBin1.method3 = Cohort_Smoothing_Symmetry(m=m, n=n, Y=Y, E=E, p=p, input.Pen=pen_method3_no_kink,
	                                                 lambda1=75, lambda2=750, start.values=fit0.NegBin1.method3$eta.start, phi.input=fit0.NegBin1.method3$phi,
	                                                 max.iter=100, max.iter.phi=100, fix.phi=FALSE, get.V.matrix=TRUE, dist="Neg.Bin1")

		## Analsis with transformed data set - NegBin1 distribution - with kink 75- 750
	fit1.NegBin1.method4 = Cohort_Smoothing_Symmetry(m=m, n=n, Y=Y, E=E, p=p, input.Pen=pen_method2_kink, 
	                                                 lambda1=75, lambda2=750, start.values=fit0.NegBin1.method2$eta.start, phi.input=fit0.NegBin1.method2$phi,
	                                                 max.iter=100, max.iter.phi=100, fix.phi=FALSE, get.V.matrix=TRUE, dist="Neg.Bin1")
	
	## Analsis with transformed penalty - NegBin1 distribution  - with kink - 75- 750
	fit1.NegBin1.method5 = Cohort_Smoothing_Symmetry(m=m, n=n, Y=Y, E=E, p=p, input.Pen=pen_method3_kink,
	                                                 lambda1=75, lambda2=750, start.values=fit0.NegBin1.method3$eta.start, phi.input=fit0.NegBin1.method3$phi, 
	                                                 max.iter=100, max.iter.phi=100, fix.phi=FALSE, get.V.matrix=TRUE, dist="Neg.Bin1")
	

	var_estimates_method1 = t(matrix(diag(fit1.NegBin1.method1$V), m, n))
	ll_matrix_method1 = fit1.NegBin1.method1$est.ETA.S - 1.96*sqrt(var_estimates_method1)
	ul_matrix_method1 = fit1.NegBin1.method1$est.ETA.S + 1.96*sqrt(var_estimates_method1)
	ind_matrix_method1 = (ll_matrix_method1 <= log(GAMMA.true) & log(GAMMA.true) <= ul_matrix_method1)

	var_estimates_method2 = t(matrix(diag(fit1.NegBin1.method2$V), m, n))
	ll_matrix_method2 = fit1.NegBin1.method2$est.ETA.S - 1.96*sqrt(var_estimates_method2)
	ul_matrix_method2 = fit1.NegBin1.method2$est.ETA.S + 1.96*sqrt(var_estimates_method2)
	ind_matrix_method2 = (ll_matrix_method2 <= log(GAMMA.true) & log(GAMMA.true) <= ul_matrix_method2)
	
	var_estimates_method3 = t(matrix(diag(fit1.NegBin1.method3$V), m, n))
	ll_matrix_method3 = fit1.NegBin1.method3$est.ETA.S - 1.96*sqrt(var_estimates_method3)
	ul_matrix_method3 = fit1.NegBin1.method3$est.ETA.S + 1.96*sqrt(var_estimates_method3)
	ind_matrix_method3 = (ll_matrix_method3 <= log(GAMMA.true) & log(GAMMA.true) <= ul_matrix_method3)
	
	var_estimates_method4 = t(matrix(diag(fit1.NegBin1.method4$V), m, n))
	ll_matrix_method4 = fit1.NegBin1.method4$est.ETA.S - 1.96*sqrt(var_estimates_method4)
	ul_matrix_method4 = fit1.NegBin1.method4$est.ETA.S + 1.96*sqrt(var_estimates_method4)
	ind_matrix_method4 = (ll_matrix_method4 <= log(GAMMA.true) & log(GAMMA.true) <= ul_matrix_method4)
	
	var_estimates_method5 = t(matrix(diag(fit1.NegBin1.method5$V), m, n))
	ll_matrix_method5 = fit1.NegBin1.method5$est.ETA.S - 1.96*sqrt(var_estimates_method5)
	ul_matrix_method5 = fit1.NegBin1.method5$est.ETA.S + 1.96*sqrt(var_estimates_method5)
	ind_matrix_method5 = (ll_matrix_method5 <= log(GAMMA.true) & log(GAMMA.true) <= ul_matrix_method5)

	### output and save all the important stuff
  matrix.conv.ind1[k,1] = fit1.NegBin1.method1$conv.ind[1]
	matrix.conv.ind1[k,2] = fit1.NegBin1.method2$conv.ind[1]
	matrix.conv.ind1[k,3] = fit1.NegBin1.method3$conv.ind[1]
	matrix.conv.ind1[k,4] = fit1.NegBin1.method4$conv.ind[1]
	matrix.conv.ind1[k,5] = fit1.NegBin1.method5$conv.ind[1]
	matrix.conv.ind2[k,1] = fit1.NegBin1.method1$conv.ind[2]
	matrix.conv.ind2[k,2] = fit1.NegBin1.method2$conv.ind[2]
	matrix.conv.ind2[k,3] = fit1.NegBin1.method3$conv.ind[2]
	matrix.conv.ind2[k,4] = fit1.NegBin1.method4$conv.ind[2]
	matrix.conv.ind2[k,5] = fit1.NegBin1.method5$conv.ind[2]
	matrix.edf[k,1] = fit1.NegBin1.method1$edf
	matrix.edf[k,2] = fit1.NegBin1.method2$edf
	matrix.edf[k,3] = fit1.NegBin1.method3$edf
	matrix.edf[k,4] = fit1.NegBin1.method4$edf
	matrix.edf[k,5] = fit1.NegBin1.method5$edf
	matrix.ll[k,1] = fit1.NegBin1.method1$ll
	matrix.ll[k,2] = fit1.NegBin1.method2$ll
	matrix.ll[k,3] = fit1.NegBin1.method3$ll
	matrix.ll[k,4] = fit1.NegBin1.method4$ll
	matrix.ll[k,5] = fit1.NegBin1.method5$ll
	matrix.aic[k,1] = fit1.NegBin1.method1$aic
	matrix.aic[k,2] = fit1.NegBin1.method2$aic
	matrix.aic[k,3] = fit1.NegBin1.method3$aic
	matrix.aic[k,4] = fit1.NegBin1.method4$aic
	matrix.aic[k,5] = fit1.NegBin1.method5$aic
	matrix.bic[k,1] = fit1.NegBin1.method1$bic
	matrix.bic[k,2] = fit1.NegBin1.method2$bic
	matrix.bic[k,3] = fit1.NegBin1.method3$bic
	matrix.bic[k,4] = fit1.NegBin1.method4$bic
	matrix.bic[k,5] = fit1.NegBin1.method5$bic
	matrix.phi[k,1] = fit1.NegBin1.method1$phi
	matrix.phi[k,2] = fit1.NegBin1.method2$phi
	matrix.phi[k,3] = fit1.NegBin1.method3$phi
	matrix.phi[k,4] = fit1.NegBin1.method4$phi
	matrix.phi[k,5] = fit1.NegBin1.method5$phi
	
	estimated.ETAmatrix1[[k]] = fit1.NegBin1.method1$est.ETA.S
	estimated.ETAmatrix2[[k]] = fit1.NegBin1.method2$est.ETA.S
	estimated.ETAmatrix3[[k]] = fit1.NegBin1.method3$est.ETA.S
	estimated.ETAmatrix4[[k]] = fit1.NegBin1.method4$est.ETA.S
	estimated.ETAmatrix5[[k]] = fit1.NegBin1.method5$est.ETA.S
	estimated.COVmatrix1[[k]] = ind_matrix_method1
	estimated.COVmatrix2[[k]] = ind_matrix_method2
	estimated.COVmatrix3[[k]] = ind_matrix_method3
	estimated.COVmatrix4[[k]] = ind_matrix_method4
	estimated.COVmatrix5[[k]] = ind_matrix_method5

	print("-----------------------------------------")
	print(paste("SIMULATION RUN ",k," of ",sim.N,sep=""))
	print("-----------------------------------------")
}


# save results
save(matrix.conv.ind1, matrix.conv.ind2, matrix.edf, matrix.ll, matrix.aic, matrix.bic, matrix.phi,
     estimated.ETAmatrix1, estimated.ETAmatrix2, estimated.ETAmatrix3, estimated.ETAmatrix4, estimated.ETAmatrix5,
     estimated.COVmatrix1, estimated.COVmatrix2, estimated.COVmatrix3, estimated.COVmatrix4, estimated.COVmatrix5,
     file="newoutput/NegBin_no_kink.RData")


#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#

load("newoutput/NegBin_no_kink.RData")



### SOME PLOTS
#-------------
x1 = 1:77-1
x2 = 1:77-1

plot.index=73
par(mfrow=c(2,2))
# plots
range(log(GAMMA.true))
image.plot(x1-1, x2-1, log(GAMMA.true), col=topo.colors(19), zlim=c(-4.5,0.0),
           xlab="Age of the respondent", ylab="Age of the contact",
           cex.lab=1.6, cex.axis=1.6)
image.plot(x1-1, x2-1,  estimated.ETAmatrix1[[plot.index]], col=topo.colors(19),
           xlab="Age of the respondent", ylab="Age of the contact", zlim=c(-4.5,0.0),
           cex.lab=1.6, cex.axis=1.6)
image.plot(x1-1, x2-1,  estimated.ETAmatrix2[[plot.index]], col=topo.colors(19),
           xlab="Age of the respondent", ylab="Age of the contact", zlim=c(-4.5,0.0),
           cex.lab=1.6, cex.axis=1.6)
image.plot(x1-1, x2-1,  estimated.ETAmatrix3[[plot.index]], col=topo.colors(19),
           xlab="Age of the respondent", ylab="Age of the contact", zlim=c(-4.5,0.0),
           cex.lab=1.6, cex.axis=1.6)

par(mfrow=c(2,2))
# plots
range(GAMMA.true*P)
image.plot(x1-1, x2-1, GAMMA.true*P, col=topo.colors(19), zlim=c(0.0,115900.0),
           xlab="Age of the respondent", ylab="Age of the contact",
           cex.lab=1.6, cex.axis=1.6)
image.plot(x1-1, x2-1,  exp(estimated.ETAmatrix1[[plot.index]])*P, col=topo.colors(19),
           xlab="Age of the respondent", ylab="Age of the contact", zlim=c(0.0,115900.0),
           cex.lab=1.6, cex.axis=1.6)
image.plot(x1-1, x2-1,  exp(estimated.ETAmatrix2[[plot.index]])*P, col=topo.colors(19),
           xlab="Age of the respondent", ylab="Age of the contact", zlim=c(0.0,115900.0),
           cex.lab=1.6, cex.axis=1.6)
image.plot(x1-1, x2-1,  exp(estimated.ETAmatrix3[[plot.index]])*P, col=topo.colors(19),
           xlab="Age of the respondent", ylab="Age of the contact", zlim=c(0.0,115900.0),
           cex.lab=1.6, cex.axis=1.6)


### BIAS ETA MATRIX
#------------------
ETA.TRUE= log(GAMMA.true)
bias.ETA1 = matrix(NA,nrow=n,ncol=m)
bias.ETA2 = matrix(NA,nrow=n,ncol=m)
bias.ETA3 = matrix(NA,nrow=n,ncol=m)
bias.ETA4 = matrix(NA,nrow=n,ncol=m)
bias.ETA5 = matrix(NA,nrow=n,ncol=m)
for (xx in 1:m){
  for (yy in 1:n){
    bias.ETA1[xx,yy] = mean( ETA.TRUE[xx,yy]-sapply(estimated.ETAmatrix1, function(mm) mm[xx,yy]) )
    bias.ETA2[xx,yy] = mean( ETA.TRUE[xx,yy]-sapply(estimated.ETAmatrix2, function(mm) mm[xx,yy]) )
    bias.ETA3[xx,yy] = mean( ETA.TRUE[xx,yy]-sapply(estimated.ETAmatrix3, function(mm) mm[xx,yy]) )
    bias.ETA4[xx,yy] = mean( ETA.TRUE[xx,yy]-sapply(estimated.ETAmatrix4, function(mm) mm[xx,yy]) )
    bias.ETA5[xx,yy] = mean( ETA.TRUE[xx,yy]-sapply(estimated.ETAmatrix5, function(mm) mm[xx,yy]) )
  }
}

range(bias.ETA1);range(bias.ETA2);range(bias.ETA3);range(bias.ETA4);range(bias.ETA5)

par(mfrow=c(3,2))
par(mar=c(5.1,5.1,2.1,1.5))
image.plot(x1-1,x2-1,ETA.TRUE,zlim=c(-4.2,1.0),col=topo.colors(20), main="True surface", 
           xlab="Age of the respondent", ylab="Age of the contact", cex.lab=1.6, cex.axis=1.6)
image.plot(x1-1,x2-1,bias.ETA5,zlim=c(-0.60,0.60),col=topo.colors(20), main="Model M0 without kink",
           xlab="Age of the respondent", ylab="Age of the contact", cex.lab=1.6, cex.axis=1.6)
image.plot(x1-1,x2-1,bias.ETA1,zlim=c(-0.60,0.60),col=topo.colors(20), main="Model M1 without kink",
           xlab="Age of the respondent", ylab="Age of the contact", cex.lab=1.6, cex.axis=1.6)
image.plot(x1-1,x2-1,bias.ETA2,zlim=c(-0.60,0.60),col=topo.colors(20), main="Model M2 without kink",
           xlab="Age of the respondent", ylab="Age of the contact", cex.lab=1.6, cex.axis=1.6)
image.plot(x1-1,x2-1,bias.ETA3,zlim=c(-0.60,0.60),col=topo.colors(20), main="Model M1 with kink",
           xlab="Age of the respondent", ylab="Age of the contact", cex.lab=1.6, cex.axis=1.6)
image.plot(x1-1,x2-1,bias.ETA4,zlim=c(-0.60,0.60),col=topo.colors(20), main="Model M2 with kink",
           xlab="Age of the respondent", ylab="Age of the contact", cex.lab=1.6, cex.axis=1.6)


sum(bias.ETA1^2)
sum(bias.ETA2^2)
sum(bias.ETA3^2)
sum(bias.ETA4^2)
sum(bias.ETA5^2)


diag.ETA1 = diag.ETA2 = diag.ETA3 = diag.ETA4 = diag.ETA5 = vector()
for (xx in 1:m){
  diag.ETA1[xx] = mean( sapply(estimated.ETAmatrix1, function(mm) mm[xx,xx]) )
  diag.ETA2[xx] = mean( sapply(estimated.ETAmatrix2, function(mm) mm[xx,xx]) )
  diag.ETA3[xx] = mean( sapply(estimated.ETAmatrix3, function(mm) mm[xx,xx]) )
  diag.ETA4[xx] = mean( sapply(estimated.ETAmatrix4, function(mm) mm[xx,xx]) )
  diag.ETA5[xx] = mean( sapply(estimated.ETAmatrix5, function(mm) mm[xx,xx]) )
}


par(mfrow=c(1,2))
par(mar=c(5.1,5.1,2.1,1.5))
plot(diag(ETA.TRUE), xlab="Index main diagonal element", ylab="Log contact rate", cex.axis=1.6, cex.lab=1.4, ylim=c(-2.2, 0.00))
lines(diag.ETA1, col=1, lwd=2.8)
lines(diag.ETA2, col=2, lwd=2.8)
lines(diag.ETA3, col=3, lwd=2.8, lty=2)
legend(30, 0.1, c("M0 w/o kink","M1 w/o kink","M2 w/o kink"), col=1:3, lty=c(1,1,2), lwd=2.8, bty="n")

plot(diag(ETA.TRUE), xlab="Index main diagonal element", ylab="Log contact rate", cex.axis=1.6, cex.lab=1.4, ylim=c(-2.2, 0.00))
lines(diag.ETA4, col=2, lwd=2.8)
lines(diag.ETA5, col=3, lwd=2.8, lty=2)
legend(30, 0.1, c("M1 w kink","M2 w kink"), col=2:3, lty=c(1,2), lwd=2.8, bty="n")



### BIAS GAMMA MATRIX
#--------------------
GAMMA.TRUE= GAMMA.true
bias.GAMMA1 = matrix(NA,nrow=n,ncol=m)
bias.GAMMA2 = matrix(NA,nrow=n,ncol=m)
bias.GAMMA3 = matrix(NA,nrow=n,ncol=m)
bias.GAMMA4 = matrix(NA,nrow=n,ncol=m)
bias.GAMMA5 = matrix(NA,nrow=n,ncol=m)
for (xx in 1:m){
  for (yy in 1:n){
    bias.GAMMA1[xx,yy] = mean( GAMMA.TRUE[xx,yy]-sapply(estimated.ETAmatrix1, function(mm) exp(mm[xx,yy])) )
    bias.GAMMA2[xx,yy] = mean( GAMMA.TRUE[xx,yy]-sapply(estimated.ETAmatrix2, function(mm) exp(mm[xx,yy])) )
    bias.GAMMA3[xx,yy] = mean( GAMMA.TRUE[xx,yy]-sapply(estimated.ETAmatrix3, function(mm) exp(mm[xx,yy])) )
    bias.GAMMA4[xx,yy] = mean( GAMMA.TRUE[xx,yy]-sapply(estimated.ETAmatrix4, function(mm) exp(mm[xx,yy])) )
    bias.GAMMA5[xx,yy] = mean( GAMMA.TRUE[xx,yy]-sapply(estimated.ETAmatrix5, function(mm) exp(mm[xx,yy])) )
  }
}

sum(bias.GAMMA1^2)
sum(bias.GAMMA2^2)
sum(bias.GAMMA3^2)
sum(bias.GAMMA4^2)
sum(bias.GAMMA5^2)

range(bias.GAMMA1);range(bias.GAMMA2);range(bias.GAMMA3);range(bias.GAMMA4);range(bias.GAMMA5)

par(mfrow=c(3,2))
par(mar=c(5.1,5.1,2.1,1.5))
image.plot(x1-1,x2-1,GAMMA.TRUE,zlim=c(0.0,1.30),col=topo.colors(20), main="True surface", 
           xlab="Age of the respondent", ylab="Age of the contact", cex.lab=1.6, cex.axis=1.6)
image.plot(x1-1,x2-1,bias.GAMMA5,zlim=c(-0.15,0.15),col=topo.colors(20), main="Model M0 without kink",
           xlab="Age of the respondent", ylab="Age of the contact", cex.lab=1.6, cex.axis=1.6)
image.plot(x1-1,x2-1,bias.GAMMA1,zlim=c(-0.15,0.15),col=topo.colors(20), main="Model M1 without kink",
           xlab="Age of the respondent", ylab="Age of the contact", cex.lab=1.6, cex.axis=1.6)
image.plot(x1-1,x2-1,bias.GAMMA2,zlim=c(-0.15,0.15),col=topo.colors(20), main="Model M2 without kink",
           xlab="Age of the respondent", ylab="Age of the contact", cex.lab=1.6, cex.axis=1.6)
image.plot(x1-1,x2-1,bias.GAMMA3,zlim=c(-0.15,0.15),col=topo.colors(20), main="Model M1 with kink",
           xlab="Age of the respondent", ylab="Age of the contact", cex.lab=1.6, cex.axis=1.6)
image.plot(x1-1,x2-1,bias.GAMMA4,zlim=c(-0.15,0.15),col=topo.colors(20), main="Model M2 with kink",
           xlab="Age of the respondent", ylab="Age of the contact", cex.lab=1.6, cex.axis=1.6)


diag.GAMMA1 = diag.GAMMA2 = diag.GAMMA3 = diag.GAMMA4 = diag.GAMMA5 = vector()
for (xx in 1:m){
  diag.GAMMA1[xx] = mean( sapply(estimated.ETAmatrix1, function(mm) exp(mm[xx,xx])) )
  diag.GAMMA2[xx] = mean( sapply(estimated.ETAmatrix2, function(mm) exp(mm[xx,xx])) )
  diag.GAMMA3[xx] = mean( sapply(estimated.ETAmatrix3, function(mm) exp(mm[xx,xx])) )
  diag.GAMMA4[xx] = mean( sapply(estimated.ETAmatrix4, function(mm) exp(mm[xx,xx])) )
  diag.GAMMA5[xx] = mean( sapply(estimated.ETAmatrix5, function(mm) exp(mm[xx,xx])) )
}

par(mfrow=c(1,2))
par(mar=c(5.1,5.1,2.1,1.5))
plot(diag(GAMMA.TRUE), xlab="Index main diagonal element", ylab="Contact rate", cex.axis=1.6, cex.lab=1.4, ylim=c(0.0, 1.00))
lines(diag.GAMMA1, col=1, lwd=2.8)
lines(diag.GAMMA2, col=2, lwd=2.8)
lines(diag.GAMMA3, col=3, lwd=2.8, lty=2)
legend(30, 1.00, c("M0 w/o kink","M1 w/o kink","M2 w/o kink"), col=1:3, lty=c(1,1,2), lwd=2.8, bty="n")

plot(diag(GAMMA.TRUE), xlab="Index main diagonal element", ylab="Contact rate", cex.axis=1.6, cex.lab=1.4, ylim=c(0.0, 1.00))
lines(diag.GAMMA4, col=2, lwd=2.8)
lines(diag.GAMMA5, col=3, lwd=2.8, lty=2)
legend(30, 1.00, c("M1 w kink","M2 w kink"), col=2:3, lty=c(1,2), lwd=2.8, bty="n")



### MSE ETA MATRIX
mse.ETA1 = matrix(NA,nrow=n,ncol=m)
mse.ETA2 = matrix(NA,nrow=n,ncol=m)
mse.ETA3 = matrix(NA,nrow=n,ncol=m)
mse.ETA4 = matrix(NA,nrow=n,ncol=m)
mse.ETA5 = matrix(NA,nrow=n,ncol=m)
for (xx in 1:m){
  for (yy in 1:n){
    mse.ETA1[xx,yy] = mean( (ETA.TRUE[xx,yy]-sapply(estimated.ETAmatrix1, function(mm) mm[xx,yy]))^2 )
    mse.ETA2[xx,yy] = mean( (ETA.TRUE[xx,yy]-sapply(estimated.ETAmatrix2, function(mm) mm[xx,yy]))^2 )
    mse.ETA3[xx,yy] = mean( (ETA.TRUE[xx,yy]-sapply(estimated.ETAmatrix3, function(mm) mm[xx,yy]))^2 )
    mse.ETA4[xx,yy] = mean( (ETA.TRUE[xx,yy]-sapply(estimated.ETAmatrix4, function(mm) mm[xx,yy]))^2 )
    mse.ETA5[xx,yy] = mean( (ETA.TRUE[xx,yy]-sapply(estimated.ETAmatrix5, function(mm) mm[xx,yy]))^2 )
  }
}
sum(mse.ETA1)
sum(mse.ETA2)
sum(mse.ETA3)
sum(mse.ETA4)
sum(mse.ETA5)

range(mse.ETA3)
par(mar=c(5.1,6.1,2.1,1.5))
image.plot(x1-1,x2-1,mse.ETA3,zlim=c(0.00,0.31),col=topo.colors(20), 
           xlab="Age of the respondent", ylab="Age of the contact", cex.lab=1.6, cex.axis=1.6)



### MSE GAMMA MATRIX
mse.GAMMA1 = matrix(NA,nrow=n,ncol=m)
mse.GAMMA2 = matrix(NA,nrow=n,ncol=m)
mse.GAMMA3 = matrix(NA,nrow=n,ncol=m)
mse.GAMMA4 = matrix(NA,nrow=n,ncol=m)
mse.GAMMA5 = matrix(NA,nrow=n,ncol=m)
for (xx in 1:m){
  for (yy in 1:n){
    mse.GAMMA1[xx,yy] = mean( (GAMMA.TRUE[xx,yy]-sapply(estimated.ETAmatrix1, function(mm) exp(mm[xx,yy])))^2 )
    mse.GAMMA2[xx,yy] = mean( (GAMMA.TRUE[xx,yy]-sapply(estimated.ETAmatrix2, function(mm) exp(mm[xx,yy])))^2 )
    mse.GAMMA3[xx,yy] = mean( (GAMMA.TRUE[xx,yy]-sapply(estimated.ETAmatrix3, function(mm) exp(mm[xx,yy])))^2 )
    mse.GAMMA4[xx,yy] = mean( (GAMMA.TRUE[xx,yy]-sapply(estimated.ETAmatrix4, function(mm) exp(mm[xx,yy])))^2 )
    mse.GAMMA5[xx,yy] = mean( (GAMMA.TRUE[xx,yy]-sapply(estimated.ETAmatrix5, function(mm) exp(mm[xx,yy])))^2 )
  }
}
sum(mse.GAMMA1)
sum(mse.GAMMA2)
sum(mse.GAMMA3)
sum(mse.GAMMA4)
sum(mse.GAMMA5)

range(mse.GAMMA3)
par(mar=c(5.1,6.1,2.1,1.5))
image.plot(x1-1, x2-1, mse.GAMMA3, zlim=c(0.00,0.012), col=topo.colors(20), 
           xlab="Age of the respondent", ylab="Age of the contact", cex.lab=1.6, cex.axis=1.6)



### COVERAGE
cov1 = cov2 = cov3 = cov4 = cov5 = vector()
for (i in 1:sim.N){
  cov1[i] = mean(estimated.COVmatrix1[[i]])
  cov2[i] = mean(estimated.COVmatrix2[[i]])
  cov3[i] = mean(estimated.COVmatrix3[[i]])
  cov4[i] = mean(estimated.COVmatrix4[[i]])
  cov5[i] = mean(estimated.COVmatrix5[[i]])
}
cbind(cov1, cov2, cov3, cov4, cov5)

mean(cov1)*100
mean(cov2)*100
mean(cov3)*100
mean(cov4)*100
mean(cov5)*100


cov.matrix = matrix(NA,nrow=n,ncol=m)
for (xx in 1:m){
  for (yy in 1:n){
    cov.matrix[xx,yy] = mean( sapply(estimated.COVmatrix3, function(mm) mm[xx,yy]) )
  }
}

range(cov.matrix)
par(mar=c(5.1,6.1,2.1,1.5))
image.plot(x1-1, x2-1, cov.matrix, zlim=c(0.00,1.00), col=topo.colors(20), 
           xlab="Age of the respondent", ylab="Age of the contact", cex.lab=1.6, cex.axis=1.6)


mean(matrix.phi[,1]);quantile(matrix.phi[,1],c(0.025, 0.975))
mean(matrix.phi[,2]);quantile(matrix.phi[,2],c(0.025, 0.975))
mean(matrix.phi[,3]);quantile(matrix.phi[,3],c(0.025, 0.975))
mean(matrix.phi[,4]);quantile(matrix.phi[,4],c(0.025, 0.975))
mean(matrix.phi[,5]);quantile(matrix.phi[,5],c(0.025, 0.975))








