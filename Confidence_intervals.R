
###**********************************************************************************************************************###
###**********************************************************************************************************************###

### SIMULATION ----

## Sample from the multivariate normal distribution for M2 kink
beta = as.vector(t(fit1.neg.bin1.M2.Kink$est.ETA.S))
Sigma = as.matrix(fit1.neg.bin1.M2.Kink$V)
dim(Sigma)

start.time <- Sys.time()
out.mvr = rmvnorm.rcpp(5000, beta, Sigma)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
out.mvr.M2kink = out.mvr


## Sample from the multivariate normal distribution for M2 kink
beta = as.vector(t(fit1.neg.bin1.M2$est.ETA.S))
Sigma = as.matrix(fit1.neg.bin1.M2$V)
dim(Sigma)

start.time <- Sys.time()
out.mvr = rmvnorm.rcpp(5000, beta, Sigma)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
out.mvr.M2 = out.mvr



# ## plot a simulated contact matrices
# simulated.contact.matrix1 = t(matrix( out.mvr[3,], 77, 77))
# simulated.contact.matrix2 = t(matrix( out.mvr[500,], 77, 77))
# simulated.contact.matrix3 = t(matrix( out.mvr[900,], 77, 77))
# 
# range(simulated.contact.matrix1)
# range(simulated.contact.matrix2)
# range(simulated.contact.matrix3)
# 
# # plot log(eta) values
# image.plot(plot.x1-1, plot.x2-1, simulated.contact.matrix3, col=topo.colors(20), zlim=c(-4.60,1.26),
#            xlab="Age of the respondent", ylab="Age of the contact",
#            cex.lab=1.6, cex.axis=1.6)
# 
# # plot on population level and make sure that it is symmetric
# br = quantile(c(exp(simulated.contact.matrix1)*P, exp(simulated.contact.matrix2)*P, exp(simulated.contact.matrix3)*P), seq(0, 1, length=20))
# test = exp(simulated.contact.matrix3) * P
# test2=matrix(NA,77,77)
# for (i in 1:m){
#   for (j in 1:m){
#     test2[i,j] = test2[j,i] = (test[i,j]+test[j,i])/2
#   }
# }
# image.plot(plot.x1-1, plot.x2-1, test2, col=topo.colors(19), breaks=br,
#            xlab="Age of the respondent", ylab="Age of the contact",
#            cex.lab=1.6, cex.axis=1.6)



## Get variability around diagonal ----
diag_matrix = matrix(NA, nrow=1000, ncol=77)
for (s in 1:5000){
  simulated.contact.matrix = t(matrix( out.mvr.M2kink[s,], 77, 77))
  diag_matrix[s,] = diag(simulated.contact.matrix)
}
ll.M2kink = apply(diag_matrix, 2, function(x) quantile(x, p=0.025))
ul.M2kink = apply(diag_matrix, 2, function(x) quantile(x, p=0.975))

diag_matrix = matrix(NA, nrow=5000, ncol=77)
for (s in 1:5000){
  simulated.contact.matrix = t(matrix( out.mvr.M2[s,], 77, 77))
  diag_matrix[s,] = diag(simulated.contact.matrix)
}
ll.M2 = apply(diag_matrix, 2, function(x) quantile(x, p=0.025))
ul.M2 = apply(diag_matrix, 2, function(x) quantile(x, p=0.975))



par(mfrow=c(1,1))
plot(plot.x-1, diag(fit1.neg.bin1.M2$est.ETA.S), type="l", lwd=3.5, ylim=c(-3.2, 2.0), 
     xlab="Index main diagonal element", ylab="Estimated log contact rate", cex.lab=1.6, cex.axis=1.6)
lines(plot.x-1, ll.M2, lty=2, lwd=2)
lines(plot.x-1, ul.M2, lty=2, lwd=2)
lines(diag(fit1.neg.bin1.M2.Kink$est.ETA.S), lwd=3.5, col=2)
lines(plot.x-1, ll.M2kink, lty=2, col=2, lwd=2)
lines(plot.x-1, ul.M2kink, lty=2, col=2, lwd=2)
points(diag(ETA), cex=diag(E)/10, pch=16)
legend(x=55, y=2.4, legend=c("w/o kink","w kink"), bty="n", col=1:2, lty=1, cex=1.6)

dev.off()






## 95% Confidence intervals matrix ----
list_out1 = list()
for (s in 1:2500){
  simulated.contact.matrix = t(matrix( out.mvr.M2kink[s,], 77, 77))
  list_out1[[s]] = simulated.contact.matrix
}
list_out2 = list()
for (s in 2501:5000){
  simulated.contact.matrix = t(matrix( out.mvr.M2kink[s,], 77, 77))
  list_out2[[s-2500]] = simulated.contact.matrix
}

list_out = c(list_out1, list_out2)

test = apply(simplify2array(list_out), 1:2, quantile, prob = c(0.025, 0.975))
M_LL = test[1,,]
M_UL = test[2,,]

## check whether it was done correctly
plot(plot.x-1, diag(fit1.neg.bin1.M2.Kink$est.ETA.S), type="l", lwd=3.5, ylim=c(-3.2, 2.0), 
     xlab="Index main diagonal element", ylab="Estimated log contact rate", cex.lab=1.6, cex.axis=1.6)
lines(plot.x-1, diag(M_LL), lty=2, lwd=2)
lines(plot.x-1, diag(M_UL), lty=2, lwd=2)


range(M_LL)
range(M_UL)
image.plot(plot.x1-1, plot.x2-1, M_LL, col=topo.colors(20), zlim=c(-7.20, 1.35),
           xlab="Age of the respondent", ylab="Age of the contact",
           cex.lab=1.6, cex.axis=1.6)

image.plot(plot.x1-1, plot.x2-1, fit1.neg.bin1.M2.Kink$est.ETA.S, col=topo.colors(20), zlim=c(-7.20, 1.35),
           xlab="Age of the respondent", ylab="Age of the contact",
           cex.lab=1.6, cex.axis=1.6)

image.plot(plot.x1-1, plot.x2-1, M_UL, col=topo.colors(20), zlim=c(-7.20, 1.35),
           xlab="Age of the respondent", ylab="Age of the contact",
           cex.lab=1.6, cex.axis=1.6)

# save all as 1000x800


