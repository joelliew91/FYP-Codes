prior_svvg<-function(x){
	mu = dnorm(x$mu,mean = 0,sd = 10,log=T)
	gamma = dnorm(x$gamma,mean = 0,sd=10,log=T)
	sigma_j = dgamma(x$sigma_j,1,1,log=T)
	rho_bar = prior_rho(x)
	k_sigmaV_v = prior_k_sigmaV_v(x)
    #print(k_sigmaV_v)
	return(mu+gamma+sigma_j+rho_bar+k_sigmaV_v)
}

prior_k_sigmaV_v<-function(x){
	d_max = 2
	d = 4*x$k*x$v/x$sigma_v^2
    #print(d)
	I = 0
	if(0.5*d<d_max)
		if(1<0.5*d)
			I = 1
	val = log(x$k)*(1-1)-1*x$k+log(I)+dgamma(x$sigma_v,1,1,log=T)+dgamma(x$v,1,1,log=T)
	return(val)
}
prior_rho<-function(x){
	w = (x$rho+1)/2
	prior = dbeta(w,1/7,1,log=T)
	return(prior)

}
