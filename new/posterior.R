library(Bessel)

source('~/Desktop/FYP/Codes/new/prior.R')

posterior_latent_z<-function(x,yt,yu,delta,z,v,vt,vu,w,set){
    val1 = prior_svvg(x)
    val2 = likelihood(yt,yu,x,delta,z,v,vt,vu)
    val3 = variance_gamma(z,x,delta)
    return(val1+val2+val3)
}

posterior_latent_w<-function(x,yt,yu,delta,z,v,vt,vu,w,set){
    val1 = prior_svvg(x)
    val2 = aux_g(v,w,x,delta,vt,vu,set)
    return(val1+val2)
}
posterior_latent_vt<-function(x,yt,yu,delta,z,v,vt,vu,w,set){
    val1 = sum(post_vt(vt,vu,x,delta))
    val2 = sum(aux_g(v,w,x,delta,vt,vu,set))
    val3 = sum(likelihood(yt,yu,x,delta,z,v,vt,vu))
    return(val1+val2+val3)
}

posterior_latent_vs<-function(x,yt,yu,delta,z,v,vt,vu,w,set){
    val1 = aux_g(v,w,x,delta,vt,vu,set)
    val2 = likelihood(yt,yu,x,delta,z,v,vt,vu)
    return(val1+val2)
}



posterior<-function(x,yt,yu,delta,z,v,vt,vu,w,set){

    val1 = prior_svvg(x)  # log form
    #print(val1)
    val2 = sum(variance_gamma(z,x,delta))    # log form
    #print(val2)
    val3 = sum(post_vt(vt,vu,x,delta))   # log form
    #print(val3)
    val4 = sum(aux_g(v,w,x,delta,vt,vu,set))
    #print(val4)
    val5 = sum(likelihood(yt,yu,x,delta,z,v,vt,vu))
    #print(val5)
    return(val1+val2+val3+val4+val5)
}

likelihood<-function(yt,yu,x,delta,z,v,vt,vu){
    v_bar = (vt-vu-x$k*x$v*delta+x$k*v)/x$sigma_v
    mn = yu+x$mu*delta+z+x$rho*v_bar
    dev = (1-x$rho^2)*v

    val = dnorm(yt,mean = mn,sd = sqrt(dev),log=T)
    
    return(val)
}

variance_gamma<-function(z,x,delta){
	a = delta
	sigma = x$sigma_j*sqrt(delta)
	val1 = 2*exp(x$gamma*z^2/sigma)/(a^(delta/a)*sqrt(2*pi)*gamma(delta/a))
	val2 = (z^2/(x$gamma^2+2*sigma^2/a))^(delta/(2*a)-0.5)
	v = z*sqrt(x$gamma^2+2*sigma^2/a)/sigma^2
	nu = delta/a-0.5
	val3 = BesselK(v,nu)
	return(log(val1*val2*val3)) 

}

post_vt<-function(vt,vu,x,delta){
	d = 4*x$k*x$v/x$sigma_v^2
	val1 = x$sigma_v^2*(1-exp(-x$k*delta))/(4*x$k)
	tau = 4*x$k*exp(x$k*delta)*vu/(x$sigma_v^2*(1-exp(-x$k*delta)))
	return(dchisq(vt,df=d,ncp = tau,log=T)+log(val1))
}

aux_g<-function(v,w,x,delta,vt,vu,set){
    #print(v[1:10])
    n = length(v)
    I_c = rep(1,n)
    I = rep(0,n)
    for(i in 1:n)
        if(v[i]<set$max_v)
            if(v[i]>set$min_v)
                if(w[i]<set$max_w)
                    if(w[i]>set$min_w){
                        I[i]= 1
                        I_c[i]= 0
                    }
    #print('hi0')
    temp = phi(w,x,delta,vt,vu)
    #print('hi')
    g = cos(v*w)*Re(temp)+sin(v*w)*Im(temp)
    #print('hi1')
    p = abs(g)*I+abs(g)*exp(-v-w)*I_c
    #print('hi3')
    #print(p)
    return(log(p))
}

phi<-function(y,x,delta,vt,vu){
    d = 4*x$k*x$v/x$sigma_v^2
    nu = sqrt(x$k^2-2*x$sigma_v^2*y*complex(imaginary=1))
    val1 = nu*exp(-0.5*(nu-x$k)*delta)*(1-exp(-x$k*(delta)))/(x$k*(1-exp(-nu*delta)))
    val2 = exp((vu+vt)/x$sigma_v^2*(x$k*(1+exp(-x$k*delta))/(1-exp(-x$k*delta))-nu*(1+exp(-nu*delta))/(1-exp(-nu*delta))))
    v1 = sqrt(vt*vu)*4*nu*exp(-0.5*nu*delta)/(x$sigma_v^2*(1-exp(-nu*delta)))
    v2 = sqrt(vt*vu)*4*x$k*exp(-0.5*x$k*delta)/(x$sigma_v^2*(1-exp(-x$k*delta)))
    #print('start')
    #print(v1[10])
    #print(d)
    val3_1 = BesselI(v1,nu=0.5*d-1)
    #print(val3_1[2])
    val3_2 = BesselI(v2,nu=0.5*d-1)
    #print(val3_2[2])
    val3 = val3_1/val3_2
    #print(val3[2])
    return(val1*val2*val3)
}
