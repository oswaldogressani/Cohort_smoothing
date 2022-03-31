# Authors: Yannick Vandendijck and Giancarlo Camarda (updated by Oswaldo Gressani)
# File last updated on 01/03/2022
#_____________________________________________________________________


Grid_Search_No_Cohort_Smoothing = function(m, n, Y, E, p, grid.l1, grid.l2, dist="Poisson", fix.phi=FALSE, phi.input=1.0){

	l1 = grid.l1
	l2 = grid.l2

	kk=1
	for (ii in l1){
	 for (jj in l2){
  	out.estimation = No_Cohort_Smoothing_No_Symmetry(m=m, n=n, Y=Y, E=E, lambda1=ii, lambda2=jj,
  	                                                   max.iter=25, max.iter.phi=25, dist=dist, fix.phi=fix.phi, phi.input=phi.input)
   	sv = out.estimation$eta.start
	  sv.phi = out.estimation$phi
	  if (fix.phi==TRUE) {sv.phi=phi.input}
  	final.est = No_Cohort_Smoothing_Symmetry(m=m, n=n, Y=Y, E=E, p=p, lambda1=ii, lambda2=jj, 
	                                             start.values=sv, max.iter=50, max.iter.phi=50,
	                                             dist=dist, phi.input=sv.phi, get.V.matrix=FALSE, fix.phi=fix.phi)
	  
  	if (kk==1){
      out.l1=ii; out.l2=jj; 
   		out.edf=final.est$edf; out.ll=final.est$ll;
   		out.aic=final.est$aic; out.bic=final.est$bic;
  	  }
  	if (kk>1){
   		out.l1=c(out.l1,ii); out.l2=c(out.l2,jj);
   		out.edf=c(out.edf,final.est$edf); out.ll=c(out.ll,final.est$ll);
   		out.aic=c(out.aic,final.est$aic); out.bic=c(out.bic,final.est$bic);
  	}
	  out.all = cbind(out.l1,out.l2,out.edf,out.ll,out.aic,out.bic,1:kk)
  	# print(paste("lambda1=",ii," & ","lambda2=",jj))
	  # print(paste("grid search ",kk," out of",length(l1)*length(l2)))
	  kk=kk+1
 	 }
	}

	index.min = which.min(out.all[,"out.aic"])
	r.l1 = out.all[index.min,"out.l1"]
	r.l2 = out.all[index.min,"out.l2"]

return(list(r.l1=r.l1, r.l2=r.l2, out=out.all))
}



Grid_Search_Cohort_Smoothing = function(PEN, m, n, Y, E, p, grid.l1, grid.l2, dist="Poisson", fix.phi=FALSE, phi.input=1.0){

	l1 = grid.l1
	l2 = grid.l2

	kk=1
	for (ii in l1){
	 for (jj in l2){
  	out.estimation = Cohort_Smoothing_No_Symmetry(m=m, n=n, Y=Y, E=E, input.Pen=PEN, 
  	                                              lambda1=ii, lambda2=jj, 
  	                                              max.iter=25, max.iter.phi=25, dist=dist, fix.phi=fix.phi, phi.input=phi.input)
 	  sv = out.estimation$eta.start
	  sv.phi = out.estimation$phi
	  if (fix.phi==TRUE) {sv.phi=phi.input}
  	final.est = Cohort_Smoothing_Symmetry(m=m, n=n, Y=Y, E=E, p=p, input.Pen=PEN,
  	                                      lambda1=ii, lambda2=jj, start.values=sv, max.iter=50, 
  	                                      max.iter.phi=50, dist=dist, phi.input=sv.phi, 
  	                                      get.V.matrix=FALSE, fix.phi=fix.phi)
	  if (kk==1){
      out.l1=ii; out.l2=jj; 
   		out.edf=final.est$edf; out.ll=final.est$ll;
   		out.aic=final.est$aic; out.bic=final.est$bic;
  	}
  	if (kk>1){
   		out.l1=c(out.l1,ii); out.l2=c(out.l2,jj);
   		out.edf=c(out.edf,final.est$edf); out.ll=c(out.ll,final.est$ll);
   		out.aic=c(out.aic,final.est$aic); out.bic=c(out.bic,final.est$bic);
  	}
	  out.all = cbind(out.l1,out.l2,out.edf,out.ll,out.aic,out.bic,1:kk)
  	# print(paste("lambda1=",ii," & ","lambda2=",jj))
	  # print(paste("grid search ",kk," out of",length(l1)*length(l2)))
  	kk=kk+1
 	 }
	}

	index.min = which.min(out.all[,"out.aic"])
	r.l1 = out.all[index.min,"out.l1"]
	r.l2 = out.all[index.min,"out.l2"]

return(list(r.l1=r.l1, r.l2=r.l2, out=out.all))
}
