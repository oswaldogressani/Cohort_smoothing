# R-function to perform Penalized-Iterative ReWeighted Least Squares (P-IRWLS)
# 
# Authors: Yannick Vandendijck and Giancarlo Camarda (updated by Oswaldo Gressani)
# File last updated on 01/03/2022
#_____________________________________________________________________

Cohort_Smoothing_No_Symmetry = function(m, n, Y, E, input.Pen, 	lambda1, lambda2, 
                                        max.iter=10, max.iter.phi=20, 
                                        dist="Poisson", fix.phi=FALSE, 
                                        phi.input=phi.input){

	y <- as.vector(t(Y))
	e <- as.vector(t(E))
	input.M = Create_Help_Data(m=m ,n=n, E=E, Y=Y)

	ind.transformed.data = as.numeric( nrow(input.Pen$Ph) == 2*m*m -m )
	ind.transformed.penalty = as.numeric( nrow(input.Pen$Ph) == m*m )

		### Starting values and Penalty
		if (ind.transformed.data==1){
			xx1 <- 0:(m+n-2)/(m+n-2)
			mm <- length(xx1)
			mmn <- mm*m	
			eta <- Matrix(0, mmn, 1)
			eta[input.M$pos] = log((input.M$breve.y+1)/(input.M$breve.e+1))[input.M$pos]

			Pen <- lambda1 * input.Pen$Pv + lambda2 * input.Pen$Ph
		}

		if (ind.transformed.penalty==1){
			eta <- Matrix(0, n*n, 1)
			eta = log((y+1)/(e+1))

			Pen <- lambda1 * input.Pen$Ph + lambda2 * input.Pen$Pd
		}


#### POISSON ####
#################

		### Iteration Steps P-IRWLS
		if (dist=="Poisson"){
		  phi=0.00
		
		  if (ind.transformed.data==1){
		    for (it in 1:max.iter){ 
		    	mu <- Matrix(0, mmn, 1)
    			mu[input.M$pos] <- input.M$breve.e[input.M$pos] * 
    			  exp(eta)[input.M$pos] * input.M$ww[input.M$pos]
  			  Psi <- Matrix(0, mmn,mmn)
			    Psi[cbind(input.M$pos, input.M$pos)] <- mu@x
			    pseudo <- Matrix(0, mmn, 1)
			    pseudo[input.M$pos] <- eta[input.M$pos] +
			      1/mu@x*(input.M$breve.y[input.M$pos] - mu@x)
  			  r <- Psi %*% pseudo
  			  WpP <- Psi + Pen
  			  etanew <- solve(WpP, r)
			    deta <- max(abs(etanew - eta))
  			  eta <- etanew
  			  cat(it, deta, "\n")
  			    if(deta<=10^-4 & it>=3) break
		      }
		      eta.start <- eta				
		      eta.vec.tr.NS <- rep(NA,mmn)		
		      eta.vec.tr.NS[input.M$pos] <- eta[input.M$pos]
		      eta.tr.NS <- matrix(eta.vec.tr.NS, mm, n)
		  
		      # re-structure as original data
		      ETA.NS <- matrix(NA, m, n)
		      for(i in 1:n){
		    	  ETA.NS[i,] <- eta.tr.NS[!is.na(eta.tr.NS[,i]), i]
		      }
		      GAMMA.NS <- exp(ETA.NS)
		    }

		  if (ind.transformed.penalty==1){
		    for (it in 1:max.iter){ 
  			mu.vector = e * exp(eta)
  			mu <- Matrix(mu.vector, n*n, 1)
  			Psi <- Diagonal(x=as.vector(mu))
			  pseudo <- eta + 1/mu*(y-mu)
  			r <- Psi %*% pseudo
  			WpP <- Psi + Pen
  			etanew <- solve(WpP, r)
  			deta <- max(abs(etanew - eta))
  			eta <- etanew
  			cat(it, deta, "\n")
  			if(deta<=10^-4 & it>=3) break
		    }
		    eta.start = eta
		    ETA.NS <- t(matrix(eta, n, n))
		    GAMMA.NS <- exp(ETA.NS)
		  }
		}



#### NEG BIN 1 ####
###################

		### Iteration Steps P-IRWLS
		if (dist=="Neg.Bin1"){
		  phi=1.00
		  if (fix.phi==TRUE) {phi=phi.input}
		  if (ind.transformed.data==1){
		    for (it1 in 1:max.iter.phi){
		      for (it in 1:max.iter){ 
		    	  mu <- Matrix(0, mmn, 1)
  			    mu[input.M$pos] <- input.M$breve.e[input.M$pos] * 
  			      exp(eta)[input.M$pos] * input.M$ww[input.M$pos]
  			    Psi <- Matrix(0, mmn,mmn)
			      Psi[cbind(input.M$pos, input.M$pos)] <- mu@x/(1+phi)
			      pseudo <- Matrix(0, mmn, 1)
			      pseudo[input.M$pos] <- eta[input.M$pos] + 1/mu@x*
			        (input.M$breve.y[input.M$pos] - mu@x)
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
		        mu.fit = exp(eta[input.M$pos]) * input.M$breve.e[input.M$pos]
		        y.fit = input.M$breve.y[input.M$pos]
		        phi.new = sum(1/mu.fit * (y.fit-mu.fit)^2)/(m*n - edf) - 1
		        if (fix.phi==TRUE) {phi.new=phi.input}
		        print(paste("disperion parameter:",phi.new))
		        diff.phi = abs(phi.new - phi.old)
		        phi = phi.new
		        if (diff.phi<=10^-4) break
		      }

		    # final parameter estimation of phi
		    for (it in 1:max.iter){ 
		    	mu <- Matrix(0, mmn, 1)
  			  mu[input.M$pos] <- input.M$breve.e[input.M$pos] * 
  			    exp(eta)[input.M$pos] * input.M$ww[input.M$pos]
  			  Psi <- Matrix(0, mmn,mmn)
			    Psi[cbind(input.M$pos, input.M$pos)] <- mu@x/(1+phi)
			    pseudo <- Matrix(0, mmn, 1)
			    pseudo[input.M$pos] <- eta[input.M$pos] + 
			      1/mu@x*(input.M$breve.y[input.M$pos] - mu@x)
  			  r <- Psi %*% pseudo
  			  WpP <- Psi + Pen
  			  etanew <- solve(WpP, r)
			    deta <- max(abs(etanew - eta))
  			  eta <- etanew
  			  cat(it, deta, "\n")
  			  if(deta<=10^-4 & it>=3) break
		    }

		    eta.start <- eta				
		    eta.vec.tr.NS <- rep(NA,mmn)					
		    eta.vec.tr.NS[input.M$pos] <- eta[input.M$pos]
		    eta.tr.NS <- matrix(eta.vec.tr.NS, mm, n)
		  
		    # re-structure as original data
		    ETA.NS <- matrix(NA, m, n)
		    for(i in 1:n){
		    	ETA.NS[i,] <- eta.tr.NS[!is.na(eta.tr.NS[,i]), i]
		    }
		    GAMMA.NS <- exp(ETA.NS)
		  }


		  if (ind.transformed.penalty==1){
		    for (it1 in 1:max.iter.phi){
		      for (it in 1:max.iter){ 
  			    mu.vector = e * exp(eta)
  			    mu <- Matrix(mu.vector, n*n, 1)
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
		      mu.fit = e*exp(eta)
		      phi.new = sum(1/mu.fit * (y-mu.fit)^2)/(m*n - edf) - 1
		      if (fix.phi==TRUE) {phi.new=phi.input}
		      print(paste("disperion parameter:",phi.new))
		      diff.phi = abs(phi.new - phi.old)
		      phi = phi.new
		      if (diff.phi<=10^-4) break
		    }

		    # final parameter estimation of phi
		    for (it in 1:max.iter){ 
  			mu.vector = e * exp(eta)
  			mu <- Matrix(mu.vector, n*n, 1)
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

		    eta.start = eta
		    ETA.NS <- t(matrix(eta, n, n))
		    GAMMA.NS <- exp(ETA.NS)
		  }
		}


#### NEG BIN 2 ####
###################

		### Iteration Steps P-IRWLS
		if (dist=="Neg.Bin2"){
		  phi=1.00
		  if (fix.phi==TRUE) {phi=phi.input}
		
		  if (ind.transformed.data==1){
		    for (it1 in 1:max.iter.phi){
		      for (it in 1:max.iter){ 
		    	mu <- Matrix(0, mmn, 1)
  			  mu[input.M$pos] <- input.M$breve.e[input.M$pos] * 
  			    exp(eta)[input.M$pos] * input.M$ww[input.M$pos]
  			  Psi <- Matrix(0, mmn,mmn)
			    Psi[cbind(input.M$pos, input.M$pos)] <- mu@x/(1+phi*mu@x)
			    pseudo <- Matrix(0, mmn, 1)
			    pseudo[input.M$pos] <- eta[input.M$pos] + 1/mu@x*
			      (input.M$breve.y[input.M$pos] - mu@x)
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
		      mu.fit = exp(eta[input.M$pos]) * input.M$breve.e[input.M$pos]
		      y.fit = input.M$breve.y[input.M$pos]
		      f.d1 <- function (aa) { sum(((y.fit-mu.fit)^2)/
		                                    (mu.fit*(1+aa*mu.fit))) - n*m + edf }
		      if (fix.phi==FALSE) {phi.new = uniroot(f.d1,c(0,5),tol=10^-5)$root}
		      if (fix.phi==TRUE) {phi.new=phi.input}
		      print(paste("disperion parameter:",phi.new))
		      diff.phi = abs(phi.new - phi.old)
		      phi = phi.new
		      if (diff.phi<=10^-4) break
		    }

		    # final parameter estimation of phi
		    for (it in 1:max.iter){ 
		    	mu <- Matrix(0, mmn, 1)
  			  mu[input.M$pos] <- input.M$breve.e[input.M$pos] *
  			    exp(eta)[input.M$pos] * input.M$ww[input.M$pos]
  			  Psi <- Matrix(0, mmn,mmn)
			    Psi[cbind(input.M$pos, input.M$pos)] <- mu@x/(1+phi*mu@x)
			    pseudo <- Matrix(0, mmn, 1)
			    pseudo[input.M$pos] <- eta[input.M$pos] + 
			      1/mu@x*(input.M$breve.y[input.M$pos] - mu@x)
  			  r <- Psi %*% pseudo
  			  WpP <- Psi + Pen
  			  etanew <- solve(WpP, r)
			    deta <- max(abs(etanew - eta))
  			  eta <- etanew
  			  cat(it, deta, "\n")
  			  if(deta<=10^-4 & it>=3) break
		    }

		    eta.start <- eta				
		    eta.vec.tr.NS <- rep(NA,mmn)					
		    eta.vec.tr.NS[input.M$pos] <- eta[input.M$pos]
		    eta.tr.NS <- matrix(eta.vec.tr.NS, mm, n)
		  
		    # re-structure as original data
		    ETA.NS <- matrix(NA, m, n)
		    for(i in 1:n){
		    	ETA.NS[i,] <- eta.tr.NS[!is.na(eta.tr.NS[,i]), i]
		    }
		    GAMMA.NS <- exp(ETA.NS)
		  }



		  if (ind.transformed.penalty==1){
		    for (it1 in 1:max.iter.phi){
		      for (it in 1:max.iter){ 
  			mu.vector = e * exp(eta)
  			mu <- Matrix(mu.vector, n*n, 1)
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
		      mu.fit = e*exp(eta)
		      f.d1 <- function (aa) { sum(((y-mu.fit)^2)/
		                                    (mu.fit*(1+aa*mu.fit))) - n*m + edf }
		      if (fix.phi==FALSE) {phi.new = uniroot(f.d1,c(0,5),tol=10^-5)$root}
		      if (fix.phi==TRUE) {phi.new=phi.input}
		      print(paste("disperion parameter:",phi.new))
		      diff.phi = abs(phi.new - phi.old)
		      phi = phi.new
		      if (diff.phi<=10^-4) break
		    }

		    # final parameter estimation of phi
		    for (it in 1:max.iter){ 
  			  mu.vector = e * exp(eta)
  			  mu <- Matrix(mu.vector, n*n, 1)
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

		    eta.start = eta
		    ETA.NS <- t(matrix(eta, n, n))
		    GAMMA.NS <- exp(ETA.NS)
		  }
		}

return(list(est.ETA.NS=ETA.NS, est.GAMMA.NS=GAMMA.NS, eta.start=eta.start, 
            phi=phi))
}





###############################################################################

Cohort_Smoothing_Symmetry = function(m, n, Y, E, p, 	input.Pen, 	start.values,
                                     lambda1, lambda2, max.iter=10, 
                                     max.iter.phi=20, dist="Poisson",	
                                     phi.input=0.00, 	get.V.matrix=FALSE, 
                                     fix.phi=FALSE){

	y <- as.vector(t(Y))
	e <- as.vector(t(E))
	input.M = Create_Help_Data(m=m ,n=n, E=E, Y=Y)


		ind.transformed.data = as.numeric( nrow(input.Pen$Ph) == 2*m*m - m )
		ind.transformed.penalty = as.numeric( nrow(input.Pen$Ph) == m*m )

		### Construct Constraints Matrix and penalty
		if (ind.transformed.data==1){
			Pen <- lambda1 * input.Pen$Pv + lambda2 * input.Pen$Ph

			KAPPA0 <- matrix(log(p), m, m)
			KAPPA1 <- t(KAPPA0)
			KAPPA01 <- KAPPA0-KAPPA1
			kappa <- c(KAPPA01[which(input.M$posZ<0)])
			nc <- length(kappa) ## number of constraints

			xx1 <- 0:(m+n-2)/(m+n-2)
			mm <- length(xx1)
			mmn <- mm*n
			L <- Matrix(0, nc, mmn)
			pos0 <- rep(1:mm, m)
			pos1 <- rep(1:m, each=mm)
			colnames(L) <- paste(pos0, pos1, sep=".")
			# placing 1
			cc1a <- outer(seq(m+1, mm), seq(0, by=mm, length=n-1), "+")
			cc1 <- cc1a[row(cc1a)+col(cc1a)<=m]
			rr1 <- 1:nc
			L[cbind(rr1,cc1)] <- 1
		
			# placing -1
			cc2a <- outer(seq(0, by=mm-1, length=n-1), seq(mm+m-1,
			                                               length=m-1,by=mm), "+")
			cc2 <- cc2a[row(cc2a)+col(cc2a)<=m]
			rr2 <- 1:nc
			L[cbind(rr2,cc2)] <- -1
		}

		if (ind.transformed.penalty==1){
			Pen <- lambda1 * input.Pen$Ph + lambda2 * input.Pen$Pd

			KAPPA0 <- matrix(log(p), n, n)
			KAPPA1 <- t(KAPPA0)
			KAPPA01 <- KAPPA0-KAPPA1
			kappa <- c(KAPPA01[which(input.M$posZ<0)])
			nc <- length(kappa) ## number of constraints

			## constraints matrix
			L <- Matrix(0, nc, n*n)
			pos0 <- rep(1:n, n)
			pos1 <- rep(1:n, each=n)
			colnames(L) <- paste(pos1, pos0, sep=".")
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
		}

		eta <- as.vector(start.values)
		conv.ind = c(0,0)


#### POISSON ####
#################
		### Iteration Steps P-IRWLS
		if (dist=="Poisson"){
		  phi=0.00

		  if (ind.transformed.data==1){
		      for (it in 1:max.iter){
  			mu <- Matrix(0, mmn, 1)
  			mu[input.M$pos] <- input.M$breve.e[input.M$pos] * exp(eta)[input.M$pos] 
  			* input.M$ww[input.M$pos]
  			Psi <- Matrix(0, mmn,mmn)
  			Psi[cbind(input.M$pos, input.M$pos)] <- mu@x
  			pseudo <- Matrix(0, mmn, 1)
			  pseudo[input.M$pos] <- eta[input.M$pos] +
			    1/mu@x*(input.M$breve.y[input.M$pos] - mu@x)
  			r <- Psi %*% pseudo
	 		  WpP <- Psi + Pen
 			  #
			  LHS <- Matrix(0, mmn+nc, mmn+nc)
  			LHS[1:mmn,1:mmn] <- WpP
  			LHS[1:mmn,1:nc+mmn] <- t(L)
  			LHS[1:nc+mmn, 1:mmn] <- L
  			RHS <- Matrix(0, mmn+nc, 1)
  			RHS[1:mmn] <- r
			  RHS[1:nc+mmn] <- kappa

  			eta.ome <- solve(LHS, RHS)
  			etanew <- eta.ome[1:mmn]
  			deta <- max(abs(eta - etanew))
  			eta <- etanew
  			cat(it, deta, '\n')
			  if(deta<=10^-4 & it>=3) conv.ind[1] = conv.ind[2] = 1
  			if(deta<=10^-4 & it>=3) break
		      }

		      eta.all.hat.vector <- eta
		      eta.hat.vector <- rep(NA,mmn)
		      eta.hat.vector[input.M$pos] <- eta[input.M$pos]
		      eta.hat.matrix <- matrix(eta.hat.vector, mm, n)

		      # re-structure as original data
		      ETA.S <- matrix(NA, m, n)
		      for(i in 1:n){
		    	ETA.S[,i] <- eta.hat.matrix[!is.na(eta.hat.matrix[,i]), i]
		      }
		      ETA.S <- t(ETA.S)
		      GAMMA.S <- exp(ETA.S)

		      # -2 log likelihood
		      y.ll = input.M$breve.y[input.M$pos]
		      mu.ll = exp(eta[input.M$pos]) * input.M$breve.e[input.M$pos]
		      ll = -2*sum( y.ll*log(mu.ll) - mu.ll - lfactorial(y.ll) )

		      # Matrix-covariance matrix
		      if (get.V.matrix==FALSE) {V=0}
		      if (get.V.matrix==TRUE) {
		        V.help1 = solve(WpP)
		        V.help2 = V.help1[input.M$pos,]
		        V.help3 = V.help2[,input.M$pos]
		        V=V.help3
		      }

		      # check constraints
		      boxplot(kappa - as.vector(L%*%eta.hat.vector))           
		  }


		  if (ind.transformed.penalty==1){
		    eta <- as.vector(start.values)
		    for (it in 1:max.iter){
 		 	  mu.vector <- e * exp(eta)
  			mu <- Matrix(mu.vector, n*n, 1)
  			Psi <- Diagonal(x=as.vector(mu))
			  pseudo <- eta + 1/mu*(y-mu)
  			r <- Psi %*% pseudo
  			WpP <- Psi + Pen
  			#
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
			  if(deta<=10^-4 & it>=3) conv.ind[1] = conv.ind[2] = 1
  			if(deta<=10^-4 & it>=3) break
		    }

		    ETA.S <- t(matrix(eta, n, n))
		    GAMMA.S <- exp(ETA.S)

		    # Deviance
		    mu.ll = exp(eta) * e
		    ll = -2*sum( y*log(mu.ll) - mu.ll - lfactorial(y) )

		    # Matrix-covariance matrix
		    if (get.V.matrix==FALSE) {V=0}
		    if (get.V.matrix==TRUE) { V=solve(WpP) }

		    # check constraints
		    boxplot(kappa-as.vector(L%*%eta))           # approx. zero
		  }
		}



#### NEG BIN 1 ####
###################
		### Iteration Steps P-IRWLS
		if (dist=="Neg.Bin1"){
		  phi=phi.input

		  if (ind.transformed.data==1){
		    for (it1 in 1:max.iter.phi){
		      for (it in 1:max.iter){
  			    mu <- Matrix(0, mmn, 1)
  			    mu[input.M$pos] <- input.M$breve.e[input.M$pos] * 
  			      exp(eta)[input.M$pos] * input.M$ww[input.M$pos]
  			    Psi <- Matrix(0, mmn,mmn)
  			    Psi[cbind(input.M$pos, input.M$pos)] <- mu@x/(1+phi)
  			    pseudo <- Matrix(0, mmn, 1)
			      pseudo[input.M$pos] <- eta[input.M$pos] + 
			        1/mu@x*(input.M$breve.y[input.M$pos] - mu@x)
  			    r <- Psi %*% pseudo
	 		      WpP <- Psi + Pen
 			      #
			      LHS <- Matrix(0, mmn+nc, mmn+nc)
  			    LHS[1:mmn,1:mmn] <- WpP
  			    LHS[1:mmn,1:nc+mmn] <- t(L)
  			    LHS[1:nc+mmn, 1:mmn] <- L
  			    RHS <- Matrix(0, mmn+nc, 1)
  			    RHS[1:mmn] <- r
			      RHS[1:nc+mmn] <- kappa

  			eta.ome <- solve(LHS, RHS)
  			etanew <- eta.ome[1:mmn]
  			deta <- max(abs(eta - etanew))
  			eta <- etanew
  			cat(it, deta, '\n')
			if(deta<=10^-4 & it>=3) conv.ind[2] = 1
  			if(deta<=10^-4 & it>=3) break
		      }
	  	    
		      # Effective degrees of freedom
		      H = solve(WpP,Psi)
		      edf = sum(diag(H))

		      # estimate disperion parameter phi
		      phi.old = phi
		      mu.fit = exp(eta[input.M$pos]) * input.M$breve.e[input.M$pos]
		      y.fit = input.M$breve.y[input.M$pos]
		      phi.new = sum(1/mu.fit * (y.fit-mu.fit)^2)/(m*n - edf) - 1
		      if (fix.phi==TRUE) {phi.new=phi}
		      print(paste("disperion parameter:",phi.new))
		      diff.phi = abs(phi.new - phi.old)
		      phi = phi.new
		      if (diff.phi<=10^-4) conv.ind[1] = 1
		      if (diff.phi<=10^-4) break
		    }

		    # final estimation of the eta
		    for (it in 1:max.iter){
  			mu <- Matrix(0, mmn, 1)
  			mu[input.M$pos] <- input.M$breve.e[input.M$pos] *
  			  exp(eta)[input.M$pos] * input.M$ww[input.M$pos]
  			Psi <- Matrix(0, mmn,mmn)
  			Psi[cbind(input.M$pos, input.M$pos)] <- mu@x/(1+phi)
  			pseudo <- Matrix(0, mmn, 1)
			  pseudo[input.M$pos] <- eta[input.M$pos] + 
			    1/mu@x*(input.M$breve.y[input.M$pos] - mu@x)
  			r <- Psi %*% pseudo
	 		  WpP <- Psi + Pen
 			  #
			  LHS <- Matrix(0, mmn+nc, mmn+nc)
  			LHS[1:mmn,1:mmn] <- WpP
  			LHS[1:mmn,1:nc+mmn] <- t(L)
  			LHS[1:nc+mmn, 1:mmn] <- L
  			RHS <- Matrix(0, mmn+nc, 1)
  			RHS[1:mmn] <- r
			  RHS[1:nc+mmn] <- kappa

  			eta.ome <- solve(LHS, RHS)
  			etanew <- eta.ome[1:mmn]
  			deta <- max(abs(eta - etanew))
  			eta <- etanew
  			cat(it, deta, '\n')
			  if(deta<=10^-4 & it>=3) conv.ind[2] = 1
  			if(deta<=10^-4 & it>=3) break
		    }

		    eta.all.hat.vector <- eta
		    eta.hat.vector <- rep(NA,mmn)
		    eta.hat.vector[input.M$pos] <- eta[input.M$pos]
		    eta.hat.matrix <- matrix(eta.hat.vector, mm, n)

		    # re-structure as original data
		    ETA.S <- matrix(NA, m, n)
		    for(i in 1:n){
		    	ETA.S[,i] <- eta.hat.matrix[!is.na(eta.hat.matrix[,i]), i]
		    }
		    ETA.S <- t(ETA.S)
		    GAMMA.S <- exp(ETA.S)

		    # -2 log likelihood
		    y.ll = input.M$breve.y[input.M$pos]
		    mu.ll = exp(eta[input.M$pos]) * input.M$breve.e[input.M$pos]
		    ll = -2*sum( lfactorial(y.ll + mu.ll/phi - 1) - lfactorial(y.ll+1-1) -
		                   lfactorial(mu.ll/phi - 1) + mu.ll/phi*log((mu.ll/phi)/
		              (mu.ll/phi+mu.ll)) + y.ll*log((mu.ll)/(mu.ll/phi+mu.ll)) )

		    # Matrix-covariance matrix
		    if (get.V.matrix==FALSE) {V=0}
		    if (get.V.matrix==TRUE) {
		      V.help1 = solve(WpP)
		      V.help2 = V.help1[input.M$pos,]
		      V.help3 = V.help2[,input.M$pos]
		      V=V.help3
		    }


		    # check constraints
		    boxplot(kappa - as.vector(L%*%eta.hat.vector))           
		  }



		  if (ind.transformed.penalty==1){
		    for (it1 in 1:max.iter.phi){
		      for (it in 1:max.iter){
 		 	mu.vector <- e * exp(eta)
  			mu <- Matrix(mu.vector, n*n, 1)
  			Psi <- Diagonal( x=as.vector(mu)/(1+phi) )
			pseudo <- eta + 1/mu*(y-mu)
  			r <- Psi %*% pseudo
  			WpP <- Psi + Pen
  			#
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
			if(deta<=10^-4 & it>=3) conv.ind[2] = 1
  			if(deta<=10^-4 & it>=3) break
		      }

	  	      # Effective degrees of freedom
		      H = solve(WpP,Psi)
		      edf = sum(diag(H))

		      # estimate disperion parameter phi
		      phi.old = phi
		      mu.fit = e*exp(eta)
		      phi.new = sum(1/mu.fit * (y-mu.fit)^2)/(m*n - edf) - 1
		      if (fix.phi==TRUE) {phi.new=phi}
		      print(paste("disperion parameter:",phi.new))
		      diff.phi = abs(phi.new - phi.old)
		      phi = phi.new
		      if (diff.phi<=10^-4) conv.ind[1] = 1
		      if (diff.phi<=10^-4) break
		    }

		    # final parameter estimation of phi
		    for (it in 1:max.iter){
 		 	mu.vector <- e * exp(eta)
  			mu <- Matrix(mu.vector, n*n, 1)
  			Psi <- Diagonal( x=as.vector(mu)/(1+phi) )
			pseudo <- eta + 1/mu*(y-mu)
  			r <- Psi %*% pseudo
  			WpP <- Psi + Pen
  			#
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
			if(deta<=10^-4 & it>=3) conv.ind[2] = 1
  			if(deta<=10^-4 & it>=3) break
		    }

		    ETA.S <- t(matrix(eta, n, n))
		    GAMMA.S <- exp(ETA.S)

		    # Deviance
		    mu.ll = exp(eta) * e
		    ll = -2*sum( lfactorial(y + mu.ll/phi - 1) - lfactorial(y+1-1) -
		                   lfactorial(mu.ll/phi - 1) + mu.ll/phi*log((mu.ll/phi)/
		                  (mu.ll/phi+mu.ll)) + y*log((mu.ll)/(mu.ll/phi+mu.ll)) )

		    # Matrix-covariance matrix
		    if (get.V.matrix==FALSE) {V=0}
		    if (get.V.matrix==TRUE) {V = solve(WpP)}

		    # check constraints
		    boxplot(kappa-as.vector(L%*%eta))           
		  }
		}



#### NEG BIN 2 ####
###################
		### Iteration Steps P-IRWLS
		if (dist=="Neg.Bin2"){
		  phi=phi.input

		  if (ind.transformed.data==1){
		    for (it1 in 1:max.iter.phi){
		      for (it in 1:max.iter){
  			mu <- Matrix(0, mmn, 1)
  			mu[input.M$pos] <- input.M$breve.e[input.M$pos] *
  			  exp(eta)[input.M$pos] * input.M$ww[input.M$pos]
  			Psi <- Matrix(0, mmn,mmn)
  			Psi[cbind(input.M$pos, input.M$pos)] <- mu@x/(1+phi*mu@x)
  			pseudo <- Matrix(0, mmn, 1)
			pseudo[input.M$pos] <- eta[input.M$pos] + 
			  1/mu@x*(input.M$breve.y[input.M$pos] - mu@x)
  			r <- Psi %*% pseudo
	 		WpP <- Psi + Pen
 			#
			LHS <- Matrix(0, mmn+nc, mmn+nc)
  			LHS[1:mmn,1:mmn] <- WpP
  			LHS[1:mmn,1:nc+mmn] <- t(L)
  			LHS[1:nc+mmn, 1:mmn] <- L
  			RHS <- Matrix(0, mmn+nc, 1)
  			RHS[1:mmn] <- r
			RHS[1:nc+mmn] <- kappa

  			eta.ome <- solve(LHS, RHS)
  			etanew <- eta.ome[1:mmn]
  			deta <- max(abs(eta - etanew))
  			eta <- etanew
  			cat(it, deta, '\n')
			if(deta<=10^-4 & it>=3) conv.ind[2] = 1
  			if(deta<=10^-4 & it>=3) break
		      }
	  	    
		      # Effective degrees of freedom
		      H = solve(WpP,Psi)
		      edf = sum(diag(H))

		      # estimate disperion parameter phi
		      phi.old = phi
		      mu.fit = exp(eta[input.M$pos]) * input.M$breve.e[input.M$pos]
		      y.fit = input.M$breve.y[input.M$pos]
		      f.d1 <- function (aa) { sum(((y.fit-mu.fit)^2)/(mu.fit*(1+aa*mu.fit)))
		        - n*m + edf }
		      if (fix.phi==FALSE) {phi.new = uniroot(f.d1,c(0,5),tol=10^-5)$root}
		      if (fix.phi==TRUE) {phi.new=phi}
		      print(paste("disperion parameter:",phi.new))
		      diff.phi = abs(phi.new - phi.old)
		      phi = phi.new
		      if (diff.phi<=10^-4) conv.ind[1] = 1
		      if (diff.phi<=10^-4) break
		    }

		    # final estimation of the eta
		    for (it in 1:max.iter){
  			mu <- Matrix(0, mmn, 1)
  			mu[input.M$pos] <- input.M$breve.e[input.M$pos] *
  			  exp(eta)[input.M$pos] * input.M$ww[input.M$pos]
  			Psi <- Matrix(0, mmn,mmn)
  			Psi[cbind(input.M$pos, input.M$pos)] <- mu@x/(1+phi*mu@x)
  			pseudo <- Matrix(0, mmn, 1)
			pseudo[input.M$pos] <- eta[input.M$pos] + 
			  1/mu@x*(input.M$breve.y[input.M$pos] - mu@x)
  			r <- Psi %*% pseudo
	 		WpP <- Psi + Pen
 			#
			LHS <- Matrix(0, mmn+nc, mmn+nc)
  			LHS[1:mmn,1:mmn] <- WpP
  			LHS[1:mmn,1:nc+mmn] <- t(L)
  			LHS[1:nc+mmn, 1:mmn] <- L
  			RHS <- Matrix(0, mmn+nc, 1)
  			RHS[1:mmn] <- r
			RHS[1:nc+mmn] <- kappa

  			eta.ome <- solve(LHS, RHS)
  			etanew <- eta.ome[1:mmn]
  			deta <- max(abs(eta - etanew))
  			eta <- etanew
  			cat(it, deta, '\n')
			if(deta<=10^-4 & it>=3) conv.ind[2] = 1
  			if(deta<=10^-4 & it>=3) break
		    }

		    eta.all.hat.vector <- eta
		    eta.hat.vector <- rep(NA,mmn)
		    eta.hat.vector[input.M$pos] <- eta[input.M$pos]
		    eta.hat.matrix <- matrix(eta.hat.vector, mm, n)

		    # re-structure as original data
		    ETA.S <- matrix(NA, m, n)
		    for(i in 1:n){
		    	ETA.S[,i] <- eta.hat.matrix[!is.na(eta.hat.matrix[,i]), i]
		    }
		    ETA.S <- t(ETA.S)
		    GAMMA.S <- exp(ETA.S)

		    # -2 log likelihood
		    y.ll = input.M$breve.y[input.M$pos]
		    mu.ll = exp(eta[input.M$pos]) * input.M$breve.e[input.M$pos]
		    ll = -2*sum( lfactorial(y.ll + 1/phi - 1) - 
		                   lfactorial(y.ll+1-1) - lfactorial(1/phi - 1) +
		                   1/phi*log((1/phi)/(1/phi+mu.ll)) + 
		                   y.ll*log((mu.ll)/(1/phi+mu.ll)) )

		    # Matrix-covariance matrix
		    if (get.V.matrix==FALSE) {V=0}
		    if (get.V.matrix==TRUE) {
		      V.help1 = solve(WpP)
		      V.help2 = V.help1[input.M$pos,]
		      V.help3 = V.help2[,input.M$pos]
		      V=V.help3
		    }


		    # check constraints
		    boxplot(kappa - as.vector(L%*%eta.hat.vector))           
		  }



		  if (ind.transformed.penalty==1){
		    for (it1 in 1:max.iter.phi){
		      for (it in 1:max.iter){
 		 	mu.vector <- e * exp(eta)
  			mu <- Matrix(mu.vector, n*n, 1)
  			Psi <- Diagonal( x=as.vector(mu)/(1+phi*as.vector(mu)) )
			pseudo <- eta + 1/mu*(y-mu)
  			r <- Psi %*% pseudo
  			WpP <- Psi + Pen
  			#
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
			if(deta<=10^-4 & it>=3) conv.ind[2] = 1
  			if(deta<=10^-4 & it>=3) break
		      }

	  	      # Effective degrees of freedom
		      H = solve(WpP,Psi)
		      edf = sum(diag(H))

		      # estimate disperion parameter phi
		      phi.old = phi
		      mu.fit = e*exp(eta)
		      f.d1 <- function (aa) { sum(((y-mu.fit)^2)/
		                                    (mu.fit*(1+aa*mu.fit))) - n*m + edf }
		      if (fix.phi==FALSE) {phi.new = uniroot(f.d1,c(0,5),tol=10^-5)$root}
		      if (fix.phi==TRUE) {phi.new=phi}
		      print(paste("disperion parameter:",phi.new))
		      diff.phi = abs(phi.new - phi.old)
		      phi = phi.new
		      if (diff.phi<=10^-4) conv.ind[1] = 1
		      if (diff.phi<=10^-4) break
		    }

		    # final parameter estimation of phi
		    for (it in 1:max.iter){
 		 	mu.vector <- e * exp(eta)
  			mu <- Matrix(mu.vector, n*n, 1)
  			Psi <- Diagonal( x=as.vector(mu)/(1+phi*as.vector(mu)) )
			pseudo <- eta + 1/mu*(y-mu)
  			r <- Psi %*% pseudo
  			WpP <- Psi + Pen
  			#
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
			if(deta<=10^-4 & it>=3) conv.ind[2] = 1
  			if(deta<=10^-4 & it>=3) break
		    }

		    ETA.S <- t(matrix(eta, n, n))
		    GAMMA.S <- exp(ETA.S)

		    # Deviance
		    mu.ll = exp(eta) * e
		    ll = -2*sum( lfactorial(y + 1/phi - 1) - lfactorial(y+1-1) -
		                   lfactorial(1/phi - 1) + 
		                   1/phi*log((1/phi)/(1/phi+mu.ll)) +
		                   y*log((mu.ll)/(1/phi+mu.ll)) )

		    # Matrix-covariance matrix
		    if (get.V.matrix==FALSE) {V=0}
		    if (get.V.matrix==TRUE) {V = solve(WpP)}

		    # check constraints
		    boxplot(kappa-as.vector(L%*%eta))           
		  }
		}


#### GENERAL STUFF ####
#######################
	# Effective degrees of freedom (CORRECT FORMULA???)
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

return(list(est.ETA.S=ETA.S, est.GAMMA.S=GAMMA.S, conv.ind=conv.ind, 
            edf=edf, ll=ll, aic=aic, bic=bic, V=V, phi=phi))
}
