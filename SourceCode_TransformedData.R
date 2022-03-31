# R-function to transform dataset such that horizontal and vertical direction 
# penalties can be used
# 
# Authors: Yannick Vandendijck and Giancarlo Camarda (updated by Oswaldo Gressani)
# File last updated on 01/03/2022
#_____________________________________________________________________


Create_Help_Data = function(m ,n, E, Y){

	### new dimensions
	xx1 <- 0:(m+n-2)/(m+n-2)
	mm <- length(xx1)
	mmn <- mm*m

	### where to place the old values
	Pos <- cbind(c(row(Y)-col(Y)+m), c(col(Y))) ## row,col
	posZ <- col(Y) - row(Y) ##

	### new (empty) matrices
	Etr <- matrix(NA,mm,n)
	Ytr <- matrix(NA,mm,n)
	### filling them up
	Etr[Pos] <- c(t(E))
	Ytr[Pos] <- c(t(Y))
	### in vectors  
	breve.y <- c(Ytr)
	breve.e <- c(Etr)

	### weights for cells with zeros
	Wtr <- matrix(0, mm, n)
	Wtr[Pos] <- 1
	pos <- which(Wtr==1)
	ww <- Matrix(0, mmn, 1)
	ww[pos] <- 1
	breve.y[is.na(breve.y)] <- 9999
	breve.e[is.na(breve.e)] <- 9999
	#cbind( tail(breve.y,100) , tail(as.vector(ww),100) )

return(list(Pos=Pos, pos=pos, posZ=posZ, Etr=Etr, Ytr=Ytr, 
	breve.y=breve.y, breve.e=breve.e, Wtr=Wtr, ww=ww))
}

