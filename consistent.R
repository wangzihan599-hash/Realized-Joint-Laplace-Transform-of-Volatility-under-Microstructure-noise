
#jump
stabfun = function(alpha, Times){   # produce jumps
  n = length(Times)
  Ti = Times
  gamrd = runif(n, -pi/2, pi/2)
  exprd = rexp(n)
  Tidif = diff(c(0, Ti))
  deltaJ = Tidif^(1/alpha)*sin(alpha*gamrd)/(cos(gamrd)^(1/alpha))*(cos((1
           - alpha)*gamrd)/exprd)^((1-alpha)/alpha)
  Jti = cumsum(deltaJ)
  return(Jti)
}
# Ti = 1:1000/1000
# res = stabfun(alpha=0.5, Times=Ti)
# plot(Ti, res, type="p", col=4, cex=0.6)
# res1 = stabfun(alpha=1, Times=Ti)
# plot(Ti, res1, type="p", col=4, cex=0.6)
# res2 = stabfun(alpha=1.9, Times=Ti)
# plot(Ti, res2, type="p", col=4, cex=0.6)

#X_1t  X_{2t} Y1 Y2
datagenfun = function(n, alpha1, alpha2, omega1, omega2,  gam ){
  Ti = 1:n/n
  J1 = stabfun(alpha=alpha1, Times=Ti)  # stable random variable
  J2 = stabfun(alpha=alpha2, Times=Ti)
  deltan = max(Ti)/n
  dW = rnorm(n, 0, sqrt(1/n))
  corrmatB = matrix(c(deltan, gam*deltan, gam*deltan, deltan), nrow=2)
  cholmatB = chol(corrmatB)
  db = matrix(rnorm(n*2, 0, 1), nrow=n)
  dB = db%*%t(cholmatB)
  
  tau1 = rep(NA, n)
  tau1[1] = dB[1, 1]/(1+0.025*deltan)
  for(i in 2:n){
    tau1[i] = (dB[i, 1] + tau1[i-1])/(1+0.025*deltan)
  }
  tau2 = rep(NA, n)
  tau2[1] = dB[1, 2]/(1+0.025*deltan)
  for(i in 2:n){
    tau2[i] = (dB[i, 2] + tau2[i-1])/(1+0.025*deltan)
  }
  
  volat1 = exp(-0.3125+0.125*tau1)
  volat2 = exp(-0.3125+0.125*tau2)
  X1 = 0.03*Ti - cumsum(0.3*volat1*(dB[ , 1]))+ cumsum(sqrt(0.91)*dW*volat1) + 0.1*J1  
  X2 = 0.03*Ti - cumsum(0.3*volat2*(dB[ , 2]))+ cumsum(sqrt(0.91)*dW*volat2) + 0.1*J2  
  
  noise1 = rnorm(n, 0, omega1)
  noise2 = rnorm(n, 0, omega2)
  Y1 = X1 + noise1
  Y2 = X2 + noise2
  
  list(T = Ti, Y1 = Y1, Y2=Y2, X1=X1, X2=X2, volat1=volat1, volat2=volat2)
}



preavefunV = function(u, v, y1, y2, theta, Ti, omega1, omega2){
  gfun = function(u){ pmin(u, 1-u) }
  ny = length(y1)
  kn = max(3, as.integer(theta*(ny^(1/2)))) 
  ln = max(3, as.integer(theta*(ny^(1/2+1/10))))
  weight = gfun(u=1:(kn-1)/kn)
  y1dif = diff(y1) 
  y2dif = diff(y2)
  y1dif[abs(y1dif)>(6*omega1)]<-0
  y2dif[abs(y2dif)>(6*omega2)]<-0
  y1preave = NULL
  y2preave = NULL
  omegahat12 = NULL
  omegahat22 = NULL
  for(i in 1:(ny-ln)) {
    y1preave = c(y1preave,  weight%*%y1dif[i:(i+kn-2)])
    omegahat12 = c(omegahat12, 1/(2*ln)*sum((y1dif[i:(i+ln-1)])^2))
  }
  for(i in 1:(ny-ln)) {
    y2preave = c(y2preave,  weight%*%y2dif[i:(i+kn-2)])
    omegahat22 = c(omegahat22, 1/(2*ln)*sum((y2dif[i:(i+ln-1)])^2))
  }
  deltan = max(Ti)/ny
  phi2 = sum(weight^2)/kn
  #weight1 = gfun(u=0:(kn)/kn)
  #subterm = NULL
  #for(i in 1:kn){
  # subterm = c(subterm, ((weight1[i]-weight1[i-1])/(1/kn))^2)
  #}
  # phi1 = mean(subterm)
  phi1 = 1  
  A = NULL 
  B = NULL
  for(i in 1:as.integer((ny-kn-ln)/(2*kn))){  #if i choose  i in 0:as.integer((nx-kn-ln)/(2*kn)  the last term has a little problem
    A = c(A, cos(sqrt(2*u)*y1preave[2*i*kn]/(sqrt(phi2*deltan*kn)))*cos(sqrt(2*v)*y2preave[2*i*kn+kn]/(sqrt(phi2*deltan*kn))))
    B = c(B, exp(u*phi1*omegahat12[2*i*kn]/(phi2*theta^2))*exp(v*phi1*omegahat22[2*i*kn+kn]/(phi2*theta^2)))
  }
  est = 2*kn*deltan*sum(A*B)
  list(y1preave=y1preave, y2preave=y2preave, est=est, omegahat12=omegahat12, omegahat22=omegahat22)
}



mainfun = function(n=23400, u=1.2, v=1.2, alpha1=0.5, alpha2=0.5,  omega1=0.001, omega2=0.0001, theta=0.06,  gam=0.3){
  turepar = NULL
  res1 = NULL
  preaveresult = NULL
  y1preave = NULL
  y2preave = NULL
  dat = datagenfun(n=n, alpha1=alpha1, alpha2=alpha2, omega1=omega1, omega2=omega2,  gam=gam)
  T = dat$T
  Y1 = dat$Y1
  Y2 = dat$Y2
  sig1 = dat$volat1
  sig2 = dat$volat2
  res1 = preavefunV(u=u, v=v, y1=Y1, y2=Y2,  theta=theta, Ti=T, omega1=omega1, omega2=omega2)   
  preaveresult = res1$est #number
  dT = diff(c(0, T))
  truepar = sum(exp(-u*sig1^2-v*sig2^2)*dT)  #number
  #standarddeviation = sqrt(covariance*kn*deltan)
  #gamma = sqrt(covariance)
  list(  truepar=truepar, preaveresult=preaveresult)
}

repfun = function(N, n=23400, u=1.2, v=1.2,  alpha1=0.5, alpha2=0.5, omega1=0.001, omega2=0.001, theta=0.6,  gam=0.3){
 
  truepar = NULL
  preaveresult = NULL
  for (k in 1:N){
    res = mainfun(n=n, u=u, v=v,  alpha1=alpha1, alpha2=alpha2,  omega1=omega1, omega2=omega2, theta=theta,  gam=gam) 
    
    preaveresult = c(preaveresult, res$preaveresult)
    truepar = c(truepar, res$truepar)
    cat("k = ", k, "\n")
  }
  #at alst I choose median of bias of Vn, sd of Vn, mse of Vn
  #bias of relative bias of  Vn, sd of relative bias of Vn, mse of relative of bias of Vn.three of them are lager
  
  biasmean = mean(preaveresult-truepar)
  sd = sd(preaveresult)
  mse = mean( (preaveresult-truepar)^2 )
  
  relbias.mean = mean( (preaveresult-truepar)/truepar )
  relbias.sd = sd((preaveresult-truepar)/truepar )
  relbias.mse = mean(((preaveresult-truepar)/truepar)^2)
  
  output = rbind( cbind(biasmean,  sd, mse, relbias.mean, relbias.sd, relbias.mse))
  colnames(output) = c("bias.mean", "sd", "mse", "relbias.mean" , "relbias.sd", "relbias.mse")
  list( output=output)
}







