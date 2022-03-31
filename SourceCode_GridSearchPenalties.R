# Authors: Yannick Vandendijck and Giancarlo Camarda (updated by Oswaldo Gressani)
# File last updated on 01/03/2022
#_____________________________________________________________________

cleverSearch_No_Cohort_Smoothing = function(m=m, n=n, Y=Y, E=E, p=p, 
                                            dist="Poisson", lambda.x, 
                                            lambda.y, length.grid=6){
  fix.phi=FALSE
  phi.input=1.0

  No_Cohort_Smoothing_Optimizer = function(lambdas){
	  l1 = lambdas[1]
	  l2 = lambdas[2]
  	out.estimation = No_Cohort_Smoothing_No_Symmetry(m=m, n=n, Y=Y, E=E, 
  	                                               lambda1=l1, lambda2=l2,
  	                                               max.iter=25, max.iter.phi=25,
  	                                               dist=dist, fix.phi=fix.phi, 
  	                                               phi.input=phi.input)
    sv = out.estimation$eta.start
	  sv.phi = out.estimation$phi
	  if (fix.phi==TRUE) {sv.phi=phi.input}
    final.est = No_Cohort_Smoothing_Symmetry(m=m, n=n, Y=Y, E=E, p=p, 
                                             lambda1=l1, lambda2=l2, 
	                                           start.values=sv, max.iter=50, 
                                             max.iter.phi=50,
	                                           dist=dist, phi.input=sv.phi, 
                                             get.V.matrix=FALSE, fix.phi=fix.phi)
	  out.aic=final.est$aic
    return(out.aic)
  }

  opt <- cleversearch(No_Cohort_Smoothing_Optimizer,
                    lower=c(log10(lambda.x)-1, log10(lambda.y)-1),
                    upper=c(log10(lambda.x)+1, log10(lambda.y)+1),
                    length.grid, logscale=TRUE, verbose=TRUE)
  return(opt)
}



cleverSearch_Cohort_Smoothing = function(m=m, n=n, Y=Y, E=E, p=p, PEN, 
                                         dist="Poisson", lambda.x, lambda.y,
                                         length.grid=6){
  fix.phi=FALSE
  phi.input=1.0
  
  Cohort_Smoothing_Optimizer = function(lambdas){
    l1 = lambdas[1]
    l2 = lambdas[2]
    out.estimation = Cohort_Smoothing_No_Symmetry(m=m, n=n, Y=Y, E=E,
                                                  input.Pen=PEN, lambda1=l1, 
                                                  lambda2=l2,
                                                  max.iter=25, max.iter.phi=25,
                                                  dist=dist, fix.phi=fix.phi, 
                                                  phi.input=phi.input)
    sv = out.estimation$eta.start
    sv.phi = out.estimation$phi
    if (fix.phi==TRUE) {sv.phi=phi.input}
    final.est = Cohort_Smoothing_Symmetry(m=m, n=n, Y=Y, E=E, p=p, 
                                          input.Pen=PEN, lambda1=l1,
                                          lambda2=l2, 
                                             start.values=sv, max.iter=50, 
                                          max.iter.phi=50,
                                             dist=dist, phi.input=sv.phi,
                                          get.V.matrix=FALSE, fix.phi=fix.phi)
    out.aic=final.est$aic
    return(out.aic)
  }
  
  opt <- cleversearch(Cohort_Smoothing_Optimizer,
                      lower=c(log10(lambda.x)-1, log10(lambda.y)-1),
                      upper=c(log10(lambda.x)+1, log10(lambda.y)+1),
                      length.grid, logscale=TRUE, verbose=TRUE)
  return(opt)
}
