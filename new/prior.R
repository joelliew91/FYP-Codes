prior_svvg<-function(x,prior){
	mu = dnorm(x$mu,mean = prior$mu_mn,sd = prior$mu_sd,log=T)
	gamma = dnorm(x$gamma,mean = prior$gamma_mn,sd=prior$gamma_sd,log=T)
	sigma_j = dgamma(x$sigma_j,prior$sigma_j_tau,prior$sigma_j_v,log=T)
	rho_bar = prior_rho(x,prior)
	k_sigmaV_v = prior_k_sigmaV_v(x,prior)
	return(mu+gamma+sigma+rho_bar+k_sigmaV_v)
}

prior_k_sigmaV_v<-function(x,prior){
	d_max = 2
	d = 4*x$k*x$v/x$sigma_v^2
	I = 0
	if(0.5*d<d_max)
		if(1<0.5*d)
			I = 1
	val = log(x$k)*(prior$k_tau-1)-prior$k_v*x$k+log(I)+dgamma(x$sigma_v,prior$sigma_v_tau,prior$sigma_v_v,log=T)+dgamma(x$v,prior$v_tau,prior$v_v,log=T)
	return(val)
}
prior_rho<-function(x,prior){
	w = (x$rho+1)/2
	prior = dbeta(w,prior$rhob_tau,prior$rhob_v,log=T)
	return(prior)

}
