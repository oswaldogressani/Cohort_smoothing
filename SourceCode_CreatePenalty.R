# Authors: Yannick Vandendijck and Giancarlo Camarda (updated by Oswaldo Gressani)
# File last updated on 01/03/2022
#_____________________________________________________________________

Penalty_Transformed_Dataset = function(m, n, kink=FALSE, max.kink.age=30){

	if(kink==FALSE){
	 ### penalty stuff
	 xx1 <- 0:(m+n-2)/(m+n-2)
	 mm <- length(xx1)
	 Dxx1 <- diff(Diagonal(mm), diff=2)
	 Dx2 <- diff(Diagonal(n), diff=2)
	 Pxx1 <- kronecker(Diagonal(n), t(Dxx1)%*%Dxx1)
	 Px2 <- kronecker(t(Dx2)%*%Dx2, Diagonal(mm))
	}

	if(kink==TRUE){
	 mka = max.kink.age + 1
	 ### penalty stuff
	 xx1 <- 0:(m+n-2)/(m+n-2)
	 mm <- length(xx1)
	 tilde.I = Diagonal(mm)
	 tilde.I[mm/2+0.5,mm/2+0.5] = 0
	 Dxx1 <- diff(tilde.I, diff=2)
	 Dxx1[(m-2),(m-1)]= -1
	 Dxx1[(m-1),(m+1)]= -1
	 Dxx1[m,(m+1)]= -1
	 Dxx1.nonadj = diff(Diagonal(mm), diff=2)
	 Dx2 <- diff(Diagonal(n), diff=2)
	 Pxx1 <- kronecker(Diagonal(x=c(rep(1,mka),rep(0,n-mka))), t(Dxx1)%*%Dxx1) + 
	   kronecker(Diagonal(x=c(rep(0,mka),rep(1,n-mka))),
	             t(Dxx1.nonadj)%*%Dxx1.nonadj)
	 Px2 <- kronecker(t(Dx2)%*%Dx2, Diagonal(mm))
	}

return(list(Pv=Pxx1,Ph=Px2))

}



Penalty_Transformed_Penalty = function(m, n, kink=FALSE, max.kink.age=30){

  if(kink==FALSE){
	  ## Horizontal Direction
	  Dh = diff(Diagonal(n), diff=2)
	  Ph = kronecker(Diagonal(n) , t(Dh) %*% Dh)
  }

	if(kink==TRUE){
	 mka = max.kink.age+1
	 Ph = Matrix(0,nrow=n*n,ncol=n*n)
	 matrix.kink.adjustment = t(Matrix(c(1,-1,0,0,0, 0,1,0,-1,0, 0,0,0,-1,1),5,3))
	 for (i in 1:n){
  	  Dh  = diff(Diagonal(n), diff=2)
  	  if (i==1 & i<mka){
    		Dh[1,1:3] = matrix.kink.adjustment[-c(1:2),-c(1:2)]
  	  }
  	  if (i==2 & i<mka){
    	  	Dh[1:2,1:4] = matrix.kink.adjustment[-c(1),-c(1)]
   	  }
	    if (i==(n-1) & i<mka){
    		Dh[(n-3):(n-2) , (n-3):n] = matrix.kink.adjustment[-c(3) , -c(5)]
  	  }
  	  if (i==n & i<mka){
    		Dh[(n-2) , (n-2):n] = matrix.kink.adjustment[-c(2:3) , -c(4:5)]
  	  }
  	  if (i %in% (3:(n-2))  & i<mka){
   		  Dh[((i-2):(i)) , ((i-2):(i+2))] = matrix.kink.adjustment
  	  }
  	  Ph[((n*(i-1)+1):(n*i)) , ((n*(i-1)+1):(n*i))] = t(Dh)%*%Dh
	 }
	}

	 # Diagonal Direction (direction of cohorts)
	 Pd = Matrix(0,nrow=n*n,ncol=n*n)
	 for (i in 1:n){
	  for (j in 1:n){
		dim.input = n-abs(j-i)
		if (dim.input > 2){
		  D <- diff(Diagonal(dim.input), diff=2)
		  DD <- t(D) %*% D
		}
		if (dim.input == 2){
		  D <- diff(Diagonal(2), diff=1)
		  DD <- t(D) %*% D
		}
		if (dim.input == 1){
		  DD <- matrix(0,nrow=1,ncol=1)
		}

		if ( i <= j ){
		 index1 <- n*(i-1) + j
		 index2 <- seq( ((j-i)+1) , (n*(n-(j-i))) , (n+1) )
		 Pd[index1,index2] <- DD[i,]
		}
		if ( i > j ){
		 index1 <- n*(i-1) + j
		 index2 <- seq( (1+n*(i-j)) , (n*n-(i-j)) , (n+1) )
		 Pd[index1,index2] <- DD[j,]
		}
  	  }
	}

return(list(Ph=Ph, Pd=Pd))

}