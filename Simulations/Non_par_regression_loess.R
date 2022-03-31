# Authors: Yannick Vandendijck and Giancarlo Camarda (updated by Oswaldo Gressani)
# File last updated on 01/03/2022
#_____________________________________________________________________

rm(list=ls(all=TRUE))

## R-Libraries
##############

library(fields)
library(Matrix)
library(clusterPower)
library(graphics)

dir = "C:/Users/lucp2486/Documents/Onderzoek/Social Contact Rates/SCR_paper1/Article R code Simulation/"
source(paste(dir,'SourceCode_TransformedData.R',sep=""))
source(paste(dir,'SourceCode_ReadInData.R',sep=""))
source(paste(dir,'SourceCode_CreatePenalty.R',sep=""))
source(paste(dir,'SourceCode_PIRWLS_Transformations.R',sep=""))
source(paste(dir,'SourceCode_PIRWLS_NoTransformation.R',sep=""))
source(paste(dir,'SourceCode_GridSearchPenalties.R',sep=""))


## Social contact data & demographic data
##---------------------------------------

in.data = Load_Social_Contact_Data()
str(in.data)


## Preparing some matrices and vectors
##------------------------------------

x1 <- 1:77					## age of the contacts (ages in which we have at least 1 respondent)
x2 <- x1					## age of the respondent
m = n = length(x1)
mn <- m*n
One <- matrix(1, ncol=m)

## observed contact matrix and vector
Y <- in.data$mat_cont[1:m,1:n]
y <- as.vector(t(Y))

## exposures (from the survey, # of respondent by age)
E <- in.data$tilde_e%*%One
e <- as.vector(t(E))

## known population
p <- in.data$P[1:m]
P <- p%*%One

## actual rates, log-rates and crude total population contacts
GAMMA <- Y/E
ETA <- log(GAMMA)
Z <- GAMMA*P

range(ETA,finite=TRUE)
image.plot(x1-1, x2-1, ETA, col=topo.colors(20), zlim=c(-3.15,2.40),
	xlab="Age of the respondent", ylab="Age of the contact",
	cex.lab=1.6, cex.axis=1.6)

range(GAMMA,finite=TRUE)
image.plot(x1, x2, GAMMA, col=topo.colors(20), zlim=c(0,11),
	xlab="Age of the respondent", ylab="Age of the contact",
	cex.lab=1.6, cex.axis=1.6)


#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#


index.matrix = matrix(NA, nrow=n, ncol=n)
for (i in 1:n){
	for (j in 1:n){
	  index.matrix[i,j] = i-j
	}
}
edit(index.matrix)

range(GAMMA)
GAMMA.vec = as.vector(GAMMA)
GAMMA.grid = as.data.frame(expand.grid(1:n,1:n))
loess.data = data.frame(y=GAMMA.vec, x1=GAMMA.grid[,1], x2=as.vector(index.matrix))


k=1
trace.hat.vector = vector()
span.grid = seq(0.01,0.15,0.001)
for (i in span.grid){
	loess.fit = loess(y~x1:x2, data=loess.data, span=i, degree=1, surface="direct", family="gaussian")
	trace.hat.vector[k] = loess.fit$trace.hat
	print(k); k=k+1
}

plot(span.grid, trace.hat.vector)
abline(h=50, col="red", lwd=2)
cbind(span.grid, trace.hat.vector)



# final fit (span=0.064 for 50 - span=0.113 for 30)
loess.fit = loess(y~x1:x2, data=loess.data, span=0.113, degree=1, surface="direct", family="gaussian")
summary(loess.fit)


pred.matrix = matrix(loess.fit$fitted, nrow=77, ncol=77)
image(pred.matrix)
range(pred.matrix[pred.matrix>0])
plot(diag(pred.matrix))

pred.matrix[pred.matrix<0] = range(pred.matrix[pred.matrix>0])[1]/2
pop.matrix.contacts = pred.matrix * P
sym.pop.matrix.contacts = matrix(NA, nrow=77, ncol=77)
for (i in 1:n){
	for (j in 1:n){
	sym.pop.matrix.contacts[i,j] = (pop.matrix.contacts[i,j]+pop.matrix.contacts[j,i])/2
	}
}
image(sym.pop.matrix.contacts, col=topo.colors(19))
est.gamma.matrix = sym.pop.matrix.contacts/P
image(est.gamma.matrix, col=topo.colors(19))

plot(diag(est.gamma.matrix))
plot(est.gamma.matrix[5,])


# output est.gamma.matrix such that it can be used as an input for the simulation study
file.name = "C:/Users/lucp2486/Documents/Onderzoek/Social Contact Rates/SCR_paper1/Article R code Simulation/Gamma_matrix_30.RData"
save(est.gamma.matrix, file=file.name)






#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#


ETA.vec = as.vector(ETA)
ETA.vec[ETA.vec < -100] = -3.14
ETA.grid = as.data.frame(expand.grid(1:n,1:n))
loess.data = data.frame(y=ETA.vec, x1=ETA.grid[,1], x2=ETA.grid[,2])

# final fit
loess.fit = loess(y~x1+x2, data=loess.data, span=0.10, degree=1, surface="direct", family="gaussian")
summary(loess.fit)

loess.fit = predict(loess.fit, newdata=data.frame(x1=ETA.grid[,1], x2=ETA.grid[,2]))
range(loess.fit)
pred.matrix = matrix(loess.fit, nrow=77, ncol=77)
image(pred.matrix)

plot(diag(ETA))

y.test = diag(GAMMA)
x.test = 1:77
loess.test = loess(y.test ~ x.test, span=0.40)
plot(log(loess.test$fitted), ylim=c(-3.5, 0.9) )

## Other option include smooth.2d
## smooth.2d( data.ETA$y, x=data.ETA$x, theta=2.5, nrow=n, ncol=n )
## Other option include smooth.2d
## smooth.2d( data.ETA$y, x=data.ETA$x, theta=2.5, nrow=n, ncol=n )

## calculate Mallow's Cp: does not work for bandwith selection in this case
t1 = sum((GAMMA.vec - loess.fit$fitted)^2)
t2 = 2 * loess.fit$s^2 * loess.fit$trace.hat
t3 = -loess.fit$s^2
Mallows.Cp = t1+t2+t3

lambda.grid = seq(0.0005,0.01,0.0001)
k=1
Mallows.Cp.vector = vector()
for (i in lambda.grid){
	loess.fit = loess(y~x1:x2, data=loess.data, span=i, degree=1, cell=0.2)
	t1 = sum((GAMMA.vec - loess.fit$fitted)^2)
	t2 = 2 * loess.fit$s^2 * loess.fit$trace.hat
	t3 = -loess.fit$s^2
	Mallows.Cp = t1+t2+t3
	Mallows.Cp.vector[k] = Mallows.Cp
	print(k)
	k=k+1
}
cbind(lambda.grid, Mallows.Cp.vector)
plot(lambda.grid, Mallows.Cp.vector)

## perform 10k-cross validation: Does not work for bandwith selection in this case
CV.smoothing.parameter.loess = function(k, param, seeding){
  set.seed(seeding)
  del.points = sample(1:n^2,n^2, replace=FALSE)
  del.indices = round(seq(1,n^2,length=k+1))
  del.indices[1] = 0
  CV.score = vector()
  current.grid = expand.grid(1:n, 1:n)
  # perform cross validations
  for (i in 1:k){
    current.del.points = del.points[(del.indices[i]+1) : del.indices[i+1]]
    current.loess.data = loess.data[-current.del.points,]
    current.loess.fit = loess(y~x1:x2, data=current.loess.data, span=param, degree=1, cell=0.2)
    current.pred = predict(current.loess.fit,newdata=data.frame(x1=current.grid[,1] , x2=current.grid[,2]))
    CV.score[i] = sum((current.pred[current.del.points] - loess.data$y[current.del.points])^2)
  }
  return(sum(CV.score))
}

CV.smoothing.parameter.loess(k=10, param=0.01, seeding=2017)

