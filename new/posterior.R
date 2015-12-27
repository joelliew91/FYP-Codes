library(Bessel)

source('~/Desktop/FYP/Codes/new/prior.R')

posterior<-function(yt,yu,x,delta,z,v,vt,vu,w,set){
    val1 = prior_svvg(x,prior)  # log form
    val2 = sum(variance_gamma(z,x,delta))    # log form
    val3 = sum(post_vt(vt,vu,x,delta))   # log form
    val4 = sum(log(aux_g(v,w,x,delta,vt,vu,set)))
    val5 = sum(likelihood(yt,yu,x,delta,z,v,vt,vu))
    return(val1+val2+val3+val4+val5)
}

likelihood<-function(yt,yu,x,delta,z,v,vt,vu){
    v_bar = (vt-vu-x$k*x$v*delta+x$k*v)/x$sigma_v
    mn = yu+x$mu*delta+z+x$rho*v_bar
    dev = (1-x$rho^2)*v
    val = dnorm(yt,mean = mn,sd = sqrt(dev),log=T)
}

variance_gamma<-function(z,x,delta){
	a = delta
	sigma = x$sigma_j*sqrt(delta)
	val1 = 2*exp(x$gamma*z^2/sigma)/(a^(delta/a)*sqrt(2*pi)*gamma(delta/a))
	val2 = (z^2/(x$gamma^2+2*sigma^2/a))^(delta/(2*a)-0.5)
	v = z*sqrt(x$gamma^2+2*sigma^2/a)/sigma^2
	nu = delta/a-0.5
	val3 = besselK(v,nu)
	return(log(val1*val2*val3)) 

}

post_vt<-function(vt,vu,x,delta){
	d = 4*x$k*x$v/x$sigma_v^2
	val1 = x$sigma_v^2*(1-exp(-x$k*delta))/(4*ex$k)
	tau = 4*x$k*exp(x$k*delta)*vu/(x$sigma_v^2*(1-exp(-x$k*delta)))
	return(dchisq(vt,df=d,ncp = tau,log=T)+log(val1))
}

aux_g<-function(v,w,x,delta,vt,vu,set){
    I_c = 1
    I = 0
    if(v<set$max_v)
        if(v>set$min_v)
            if(w<set$max_w)
                if(w>set$min_w){
                    I = 1
                    I_c = 0
                }
                
    temp = phi(w,x,delta,vt,vu)
    g = cos(v*w)*Re(temp)+sin(v*w)*Im(temp)
    p = abs(g)*I+abs(g)*exp(-v-w)*I_c
    return(p)
}

phi<-function(y,x,delta,vt,vu){
    d = 4*x$k*x$v/x$sigma_v^2
    nu = sqrt(x$k^2-2*x$sigma_v^2*y*complex(imaginary=1))
    val1 = nu*exp(-0.5*(nu-x$k)*delta)*(1-exp(-x$k*(delta)))/(x$k*(1-exp(-nu*delta)))
    val2 = exp((vu+vt)/x$sigma_v^2*(x$k*(1+exp(-x$k*delta))/(1-exp(-x$k*delta))-nu*(1+exp(-nu*delta))/(1-exp(-nu*delta))))
    v1 = sqrt(vt*vu)*4*nu*exp(-0.5*nu*delta)/(x$sigma_v^2*(1-exp(-nu*delta)))
    v2 = sqrt(vt*vu)*4*x$k*exp(-0.5*x$k*delta)/(x$sigma_v^2*(1-exp(-x$k*delta)))
    val3 = besselI(v1,nu=0.5*d-1)/besselI(v2,nu=0.5*d-1)
    return(val1*val2*val3)
}