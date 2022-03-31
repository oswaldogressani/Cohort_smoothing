# R-function to perform Penalized-Iterative ReWeighted Least Squares (P-IRWLS)
# 
# Author: Yannick Vandendijck (based on programs of Giancarlo Camarda)
# File created on 04/11/2016
#-----------------------------------------------------------------------------

No_Cohort_Smoothing_No_Symmetry = function(m, n, Y, E, lambda1, lambda2, 
                                           max.iter=20, max.iter.phi=20, 	
                                           dist="Poisson", fix.phi=FALSE,
                                           phi.input=1.0){

	y <- as.vector(t(Y))
	e <- as.vector(t(E))
	
	### penalty stuff
	Dx1 <- diff(Diagonal(m), diff=2)
	Dx2 <- diff(Diagonal(n), diff=2)
	Px1 <- kronecker(Diagonal(m), t(Dx1)%*%Dx1)
	Px2 <- kronecker(t(Dx2)%*%Dx2, Diagonal(n))
	Pen <- lambda1*Px1 + lambda2*Px2

	### starting values
	eta <- Matrix(0, m*n, 1)
	eta = log((y+1)/(e+1))


#### POISSON ####
#################
	### Iteration Steps P-IRWLS - Poisson
	if (dist=="Poisson"){
	phi=0.00
	for (it in 1:max.iter){ 
  	mu.vector = e * exp(eta)
  	mu <- Matrix(mu.vector, m*n, 1)
  	Psi <- Diagonal(x=as.vector(mu))
		pseudo <- eta + 1/mu*(y-mu)
  	r <- Psi %*% pseudo
  	WpP <- Psi + Pen
  	etanew <- solve(WpP, r)
  	deta <- max(abs(etanew - eta))
  	eta <- etanew
  	cat(it, deta, "\n")
  	if(deta<=10^-4 & it>=3) break
	}}


#### NEG BIN 1 ####
###################
	### Iteration Steps P-IRWLS - negative binomial 1
	if (dist=="Neg.Bin1"){
	phi=1.00
	if (fix.phi==TRUE) {phi=phi.input}
	for (it1 in 1:max.iter.phi){
	  for (it in 1:max.iter){ 
  		mu.vector = e * exp(eta)
  		mu <- Matrix(mu.vector, m*n, 1)
  		Psi <- Diagonal( x=as.vector(mu)/(1+phi) )
		  pseudo <- eta + 1/mu*(y-mu)
  		r <- Psi %*% pseudo
  		WpP <- Psi + Pen
  		etanew <- solve(WpP, r)
  		deta <- max(abs(etanew - eta))
  		eta <- etanew
  		cat(it, deta, "\n")
  		if(deta<=10^-4 & it>=3) break
	  }
	  
	  # Effective degrees of freedom
	  H = solve(WpP,Psi)
	  edf = sum(diag(H))

	  # estimate disperion parameter phi
	  phi.old = phi
	  mu.fit = exp(eta) * e
	  phi.new = sum(1/mu.fit * (y-mu.fit)^2)/(m*n - edf) - 1
	  if (fix.phi==TRUE) {phi.new=phi.input}
	  print(paste("disperion parameter:",phi.new))
	  diff.phi = abs(phi.new - phi.old)
	  phi = phi.new
	  if (diff.phi<=10^-4) break
	}

	# final parameter estimation of eta
	  for (it in 1:max.iter){ 
  		mu.vector = e * exp(eta)
  		mu <- Matrix(mu.vector, m*n, 1)
  		Psi <- Diagonal( x=as.vector(mu)/(1+phi) )
		  pseudo <- eta + 1/mu*(y-mu)
  		r <- Psi %*% pseudo
  		WpP <- Psi + Pen
  		etanew <- solve(WpP, r)
  		deta <- max(abs(etanew - eta))
  		eta <- etanew
  		cat(it, deta, "\n")
  		if(deta<=10^-4 & it>=3) break
	  }
	}


#### NEG BIN 2 ####
###################
	### Iteration Steps P-IRWLS - negative binomial 2
	if (dist=="Neg.Bin2"){
	phi=1.00
	if (fix.phi==TRUE) {phi=phi.input}
	for (it1 in 1:max.iter.phi){
	  for (it in 1:max.iter){ 
  		mu.vector = e * exp(eta)
  		mu <- Matrix(mu.vector, m*n, 1)
  		Psi <- Diagonal( x=as.vector(mu)/(1+phi*as.vector(mu)) )
		  pseudo <- eta + 1/mu*(y-mu)
  		r <- Psi %*% pseudo
  		WpP <- Psi + Pen
  		etanew <- solve(WpP, r)
  		deta <- max(abs(etanew - eta))
  		eta <- etanew
  		cat(it, deta, "\n")
  		if(deta<=10^-4 & it>=3) break
	  }
	  
	  # Effective degrees of freedom
	  H = solve(WpP,Psi)
	  edf = sum(diag(H))

	  # estimate disperion parameter phi
	  phi.old = phi
	  mu.fit = exp(eta) * e
	  f.d1 <- function (aa) { sum(((y-mu.fit)^2)/(mu.fit*(1+aa*mu.fit))) - n*m + edf }
	  phi.new = uniroot(f.d1,c(0,10),tol=10^-5)$root
	  if (fix.phi==TRUE) {phi.new=phi.input}
	  print(paste("disperion parameter:",phi.new))
	  diff.phi = abs(phi.new - phi.old)
	  phi = phi.new
	  if (diff.phi<=10^-4) break
	}

	# final parameter estimation of eta
	  for (it in 1:max.iter){ 
  		mu.vector = e * exp(eta)
  		mu <- Matrix(mu.vector, m*n, 1)
  		Psi <- Diagonal( x=as.vector(mu)/(1+phi*as.vector(mu)) )
		  pseudo <- eta + 1/mu*(y-mu)
  		r <- Psi %*% pseudo
  		WpP <- Psi + Pen
  		etanew <- solve(WpP, r)
  		deta <- max(abs(etanew - eta))
  		eta <- etanew
  		cat(it, deta, "\n")
  		if(deta<=10^-4 & it>=3) break
	  }
	}


#### GENERAL STUFF ####
#######################
	eta.start = eta
	ETA.NS <- t(matrix(eta, n, n))
	GAMMA.NS <- exp(ETA.NS)
	
return(list(est.ETA.NS=ETA.NS, est.GAMMA.NS=GAMMA.NS, eta.start=eta.start,
            phi=phi))
}





No_Cohort_Smoothing_Symmetry= function(m, n, Y, E, p, start.values, lambda1, 
                                       lambda2, max.iter=10, 	max.iter.phi=20,
                                       dist="Poisson", phi.input=1.00, 
                                       get.V.matrix=FALSE, fix.phi=FALSE){

	y <- as.vector(t(Y))
	e <- as.vector(t(E))
	input.M = Create_Help_Data(m=m ,n=n, E=E, Y=Y)

	### penalty stuff
	Dx1 <- diff(Diagonal(m), diff=2)
	Dx2 <- diff(Diagonal(n), diff=2)
	Px1 <- kronecker(Diagonal(m), t(Dx1)%*%Dx1)
	Px2 <- kronecker(t(Dx2)%*%Dx2, Diagonal(n))
	Pen <- lambda1*Px1 + lambda2*Px2

	### Construct Constraints Matrix
	KAPPA0 <- matrix(log(p), n, n)
	KAPPA1 <- t(KAPPA0)
	KAPPA01 <- KAPPA0-KAPPA1
	kappa <- c(KAPPA01[which(input.M$posZ<0)])
	nc <- length(kappa) ## number of constraints

	# constraints matrix
	L <- Matrix(0, nc, n*n)
	pos0 <- rep(1:n, n)
	pos1 <- rep(1:n, each=n)
	#colnames(L) <- paste(pos1, pos0, sep=".")
	# placing 1
	cc1a <- outer(seq(0, by=n, length=n), seq(1, n), "+")
	cc1 <- cc1a[row(cc1a) < col(cc1a)]
	rr1 <- 1:nc
	L[cbind(rr1,sort(cc1))] <- 1
	# placing -1
	cc2a <- outer(seq(0, by=n, length=n), seq(1, n), "+")
	cc2 <- cc2a[row(cc2a) > col(cc2a)]
	rr2 <- 1:nc
	L[cbind(rr2,cc2)] <- -1

	eta <- as.vector(start.values)
	conv.ind = c(0,0)


#### POISSON ####
#################
	### Iteration Steps P-IRWLS - Poisson
	if (dist=="Poisson"){
	phi = 0.00
	for (it in 1:max.iter){
 	 	mu.vector <- e * exp(eta)
  	mu <- Matrix(mu.vector, n*n, 1)
  	Psi <- Diagonal(x=as.vector(mu))
		pseudo <- eta + 1/mu*(y-mu)
		r <- Psi %*% pseudo
  	WpP <- Psi + Pen
  	LHS <- Matrix(0, n*n + nc, n*n + nc)
  	LHS[1:(n*n),1:(n*n)] <- WpP
  	LHS[1:(n*n),1:nc + n*n] <- t(L)
  	LHS[1:nc + n*n, 1:(n*n)] <- L
  	RHS <- Matrix(0, (n*n)+nc, 1)
  	RHS[1:(n*n)] <- r
  	RHS[1:nc+(n*n)] <- kappa
  
  	eta.ome <- solve(LHS, RHS)
  	etanew <- eta.ome[1:(n*n)]
  	deta <- max(abs(eta - etanew))
  	eta <- etanew
  	cat(it, deta, '\n')
		if(deta<=10^-4 & it>=3) conv.ind[1] = conv.ind[2] =1 
		if(deta<=10^-4 & it>=3) break
	}
	# -2 loglikelihood
	mu.ll = exp(eta) * e
	ll = -2*sum( y*log(mu.ll) - mu.ll - lfactorial(y) )
	}



#### NEG BIN 1 ####
###################
	### Iteration Steps P-IRWLS - negative binomial 1
	if (dist=="Neg.Bin1"){
	phi = phi.input
	for (it1 in 1:max.iter.phi){
	
	  for (it in 1:max.iter){
 	 	mu.vector <- e * exp(eta)
  		mu <- Matrix(mu.vector, n*n, 1)
  		Psi <- Diagonal( x=as.vector(mu)/(1+phi) )
		  pseudo <- eta + 1/mu*(y-mu)
		  r <- Psi %*% pseudo
  		WpP <- Psi + Pen
  		LHS <- Matrix(0, n*n + nc, n*n + nc)
  		LHS[1:(n*n),1:(n*n)] <- WpP
  		LHS[1:(n*n),1:nc + n*n] <- t(L)
  		LHS[1:nc + n*n, 1:(n*n)] <- L
  		RHS <- Matrix(0, (n*n)+nc, 1)
  		RHS[1:(n*n)] <- r
  		RHS[1:nc+(n*n)] <- kappa
  
  		eta.ome <- solve(LHS, RHS)
  		etanew <- eta.ome[1:(n*n)]
  		deta <- max(abs(eta - etanew))
  		eta <- etanew
  		cat(it, deta, '\n')
		if(deta<=10^-4 & it>=3) conv.ind[1] = 1
		if(deta<=10^-4 & it>=3) break
	  }

	  # Effective degrees of freedom
	  H = solve(WpP,Psi)
	  edf = sum(diag(H))

	  # estimate disperion parameter phi
	  phi.old = phi
	  mu.fit = exp(eta) * e
	  if (fix.phi==FALSE) {phi.new = sum(1/mu.fit * (y-mu.fit)^2)/(m*n - edf) - 1}
	  if (fix.phi==TRUE) {phi.new=phi}
	  print(paste("disperion parameter:",phi.new))
	  diff.phi = abs(phi.new - phi.old)
	  phi = phi.new
	  if (diff.phi<=10^-4) conv.ind[2]=1
	  if (diff.phi<=10^-4) break
	}

	# final estimation of the eta
	  for (it in 1:max.iter){
 	 	mu.vector <- e * exp(eta)
  		mu <- Matrix(mu.vector, n*n, 1)
  		Psi <- Diagonal( x=as.vector(mu)/(1+phi) )
		pseudo <- eta + 1/mu*(y-mu)
		r <- Psi %*% pseudo
  		WpP <- Psi + Pen
  		LHS <- Matrix(0, n*n + nc, n*n + nc)
  		LHS[1:(n*n),1:(n*n)] <- WpP
  		LHS[1:(n*n),1:nc + n*n] <- t(L)
  		LHS[1:nc + n*n, 1:(n*n)] <- L
  		RHS <- Matrix(0, (n*n)+nc, 1)
  		RHS[1:(n*n)] <- r
  		RHS[1:nc+(n*n)] <- kappa
  
  		eta.ome <- solve(LHS, RHS)
  		etanew <- eta.ome[1:(n*n)]
  		deta <- max(abs(eta - etanew))
  		eta <- etanew
  		cat(it, deta, '\n')
		if(deta<=10^-4 & it>=3) conv.ind[1] = 1
		if(deta<=10^-4 & it>=3) break
	  }

	# -2 loglikelihood
	mu.ll = exp(eta) * e
	ll = -2*sum( lfactorial(y + mu.ll/phi - 1) - 
	               lfactorial(y+1-1) - lfactorial(mu.ll/phi - 1) +
	               mu.ll/phi*log((mu.ll/phi)/(mu.ll/phi+mu.ll)) +
	               y*log((mu.ll)/(mu.ll/phi+mu.ll)) )
	}


#### NEG BIN 2 ####
###################
	### Iteration Steps P-IRWLS - negative binomial 2
	if (dist=="Neg.Bin2"){
	phi = phi.input
	for (it1 in 1:max.iter.phi){
	
	  for (it in 1:max.iter){
 	 	mu.vector <- e * exp(eta)
  		mu <- Matrix(mu.vector, n*n, 1)
  		Psi <- Diagonal( x=as.vector(mu)/(1+phi*as.vector(mu)) )
		pseudo <- eta + 1/mu*(y-mu)
		r <- Psi %*% pseudo
  		WpP <- Psi + Pen
  		LHS <- Matrix(0, n*n + nc, n*n + nc)
  		LHS[1:(n*n),1:(n*n)] <- WpP
  		LHS[1:(n*n),1:nc + n*n] <- t(L)
  		LHS[1:nc + n*n, 1:(n*n)] <- L
  		RHS <- Matrix(0, (n*n)+nc, 1)
  		RHS[1:(n*n)] <- r
  		RHS[1:nc+(n*n)] <- kappa
  
  		eta.ome <- solve(LHS, RHS)
  		etanew <- eta.ome[1:(n*n)]
  		deta <- max(abs(eta - etanew))
  		eta <- etanew
  		cat(it, deta, '\n')
		if(deta<=10^-4 & it>=3) conv.ind[1] = 1
		if(deta<=10^-4 & it>=3) break
	  }

	  # Effective degrees of freedom
	  H = solve(WpP,Psi)
	  edf = sum(diag(H))

	  # estimate disperion parameter phi
	  phi.old = phi
	  mu.fit = exp(eta) * e
	  f.d1 <- function (aa) { sum(((y-mu.fit)^2)/(mu.fit*(1+aa*mu.fit))) - 
	      n*m + edf }
	  if (fix.phi==FALSE) {phi.new = uniroot(f.d1,c(0,10),tol=10^-5)$root}
	  if (fix.phi==TRUE) {phi.new=phi}
	  print(paste("disperion parameter:",phi.new))
	  diff.phi = abs(phi.new - phi.old)
	  phi = phi.new
	  if (diff.phi<=10^-4) conv.ind[2]=1
	  if (diff.phi<=10^-4) break
	}

	# final estimation of the eta
	  for (it in 1:max.iter){
 	 	mu.vector <- e * exp(eta)
  		mu <- Matrix(mu.vector, n*n, 1)
  		Psi <- Diagonal( x=as.vector(mu)/(1+phi*as.vector(mu)) )
		pseudo <- eta + 1/mu*(y-mu)
		r <- Psi %*% pseudo
  		WpP <- Psi + Pen
  		LHS <- Matrix(0, n*n + nc, n*n + nc)
  		LHS[1:(n*n),1:(n*n)] <- WpP
  		LHS[1:(n*n),1:nc + n*n] <- t(L)
  		LHS[1:nc + n*n, 1:(n*n)] <- L
  		RHS <- Matrix(0, (n*n)+nc, 1)
  		RHS[1:(n*n)] <- r
  		RHS[1:nc+(n*n)] <- kappa
  
  		eta.ome <- solve(LHS, RHS)
  		etanew <- eta.ome[1:(n*n)]
  		deta <- max(abs(eta - etanew))
  		eta <- etanew
  		cat(it, deta, '\n')
		if(deta<=10^-4 & it>=3) conv.ind[1] = 1
		if(deta<=10^-4 & it>=3) break
	  }

	# -2 loglikelihood
	mu.ll = exp(eta) * e
	ll = -2*sum( lfactorial(y + 1/phi - 1) - lfactorial(y+1-1) - 
	               lfactorial(1/phi - 1) + 1/phi*log((1/phi)/(1/phi+mu.ll)) + 
	               y*log((mu.ll)/(1/phi+mu.ll)) )
	}


#### GENERAL STUFF ####
#######################
	ETA.S <- t(matrix(eta, n, n))
	GAMMA.S <- exp(ETA.S)

	# Effective degrees of freedom
	H = solve(WpP,Psi)
	edf = sum(diag(H))
	if (dist=="Poisson"){
		aic = ll + 2*(edf)
		bic = ll + log(m*n)*(edf)
	}
	if (dist=="Neg.Bin1" | dist=="Neg.Bin2"){
		aic = ll + 2*(edf+1)
		bic = ll + log(m*n)*(edf+1)
	}

	# Matrix-covariance matrix
	if (get.V.matrix==FALSE) {V=0}
	if (get.V.matrix==TRUE) {V = solve(WpP)}

	# check constraints
	boxplot(kappa-as.vector(L%*%eta))           # approx. zero

return(list(est.ETA.S=ETA.S, est.GAMMA.S=GAMMA.S, conv.ind=conv.ind, edf=edf, 
            ll=ll, aic=aic, bic=bic, V=V, phi=phi))
}

