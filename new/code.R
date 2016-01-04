library(MASS)                                   ## library required to generate multivariate norm rv
library(MCMCpack)                               ## library required to generate inverse gamma rv

source('~/Desktop/FYP/Codes/new/posterior.R') # load the posterior functions
source('~/Desktop/FYP/Codes/sim_svvg.R') # load the price simulator
source('~/Desktop/FYP/Codes/new/mhrw.R') # load the MH-random walk
source('~/Desktop/FYP/Codes/new/update.R')

# need to fix delta for the sim

mcmc_para<-function(data,delta=1,iterations){
    set.seed(150)
    n = length(data)
    

    mu = rnorm(1,mean = 0,sd = 1)
    gamma = rnorm(1,mean = 0,sd=1)
    set = list(max_v = 15000,min_v = 0,max_w = 15000,min_w = 0)
    sigma_j = rgamma(1,1,1)
    k = 5
    theta = abs(rnorm(1))
    v_p =1/rgamma(1,10,1/10)
    rho = 0.5
    sigma_v = 0.25
    
    vt = rep(1,n*delta) # latent variable
    v_star = rep(1,n*delta) # latent variable
    z = rgamma(n*delta,1,1)    # latent variable
    w = rep(10000,n*delta)
    
    x = list(mu=mu,k=k,theta=theta,sigma_v=sigma_v,rho=rho,gamma=gamma,sigma_j=sigma_j,v=v_p)
    
    u1 = sigma_v^2/(2*v_p)
    u2 = sigma_v^2/(2*k)
    u3 = sqrt(2*k*v_p)
    z_list = t(as.matrix(z))
    w_list = t(as.matrix(w))
    v_list = t(as.matrix(vt))
    vs_list = t(as.matrix(v_star))
    for(i in 1:iterations){
        #print(z)
        print('mu')
        mu[i+1] = update_mu(mu,x,data[2:n],data[1:n-1],delta,z[2:n],v_star[2:n],vt[2:n],vt[1:n-1],w[2:n],set)$mu
        x$mu = mu[i+1]
        print('gam')
        gamma[i+1] = update_gamma(gamma,x,data[2:n],data[1:n-1],delta,z[2:n],v_star[2:n],vt[2:n],vt[1:n-1],w[2:n],set)$gamma
        x$gamma = gamma[i+1]
        
        print('sigma_j')
        sigma_j[i+1] = update_sigma_j(sigma_j,x,data[2:n],data[1:n-1],delta,z[2:n],v_star[2:n],vt[2:n],vt[1:n-1],w[2:n],set)$sigma_j
        x$sigma_j = sigma_j[i+1]
        
        print('rho')
        rho[i+1] = update_rho(rho,x,data[2:n],data[1:n-1],delta,z[2:n],v_star[2:n],vt[2:n],vt[1:n-1],w[2:n],set)$rho
        x$rho = rho[i+1]
        
        print('k')
        k[i+1] = update_k(u1,k,x,data[2:n],data[1:n-1],delta,z[2:n],v_star[2:n],vt[2:n],vt[1:n-1],w[2:n],set)$k
        x$k = k[i+1]
        
        print('v_p')
        v_p[i+1] = update_v(u2,v_p,x,data[2:n],data[1:n-1],delta,z[2:n],v_star[2:n],vt[2:n],vt[1:n-1],w[2:n],set)$v
        x$v = v_p[i+1]
        
        print('sigma_v')
        sigma_v[i+1] = update_sigma_v(u3,sigma_v,x,data[2:n],data[1:n-1],delta,z[2:n],v_star[2:n],vt[2:n],vt[1:n-1],w[2:n],set)$sigma_v
        x$sigma_v = sigma_v[i+1]
        
        print('z')
        z[2:n] = update_z(z_list,x,data[2:n],data[1:n-1],delta,z[2:n],v_star[2:n],vt[2:n],vt[1:n-1],w[2:n],set)
        z_list = rbind(z_list,z)
        print('v')
        vt = update_vt(v_list,x,data[2:n],data[1:n-1],delta,z[2:n],v_star[2:n],vt,w[2:n],set)
        v_list = rbind(v_list,vt)
        print('v_star')
        v_star[2:n] = update_vs(vs_list,x,data[2:n],data[1:n-1],delta,z[2:n],v_star[2:n],vt[2:n],vt[1:n-1],w[2:n],set)
        vs_list = rbind(vs_list,v_star)
        print('w')
        w[2:n] = update_w(w_list,x,data[2:n],data[1:n-1],delta,z[2:n],v_star[2:n],vt[2:n],vt[1:n-1],w[2:n],set)
        w_list = rbind(w_list,w)
        print('end latent')
        
        u1[i+1] = x$sigma_v^2/(2*x$v)
        u2[i+1] = x$sigma_v^2/(2*x$k)
        u3[i+1] = sqrt(2*x$k*x$v)
        print(i)
    }

    print(rejectionRate(mcmc(mu)))
    print(rejectionRate(mcmc(gamma)))
    print(rejectionRate(mcmc(sigma_j)))
    print(rejectionRate(mcmc(rho)))
    print(rejectionRate((mcmc(k))))
    print(rejectionRate(mcmc(v_p)))
    print(rejectionRate(mcmc(sigma_v)))
    #return(list(gamma=gamma,mu=mu,sigma_j=sigma_j,rho=rho,k=k,v=v_p,sigma_v= sigma_v))
    return(w_list)
}



