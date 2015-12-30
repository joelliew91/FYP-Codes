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
    set = list(max_v = 10000,min_v = 0,max_w = 10000,min_w = 0)
    sigma_j = rgamma(1,1,1)
    k = 5
    theta = abs(rnorm(1))
    v_p =1/rgamma(1,10,1/10)
    rho = 0.5
    sigma_v = 0.25
    
    vt = rep(1,n*delta) # latent variable
    v_star = rep(1,n*delta) # latent variable
    z = rgamma(n*delta,1,1)    # latent variable
    w = rep(1,n*delta)
    
    x = list(mu=mu,k=k,theta=theta,sigma_v=sigma_v,rho=rho,gamma=gamma,sigma_j=sigma_j,v=v_p)
    u1 = sigma_v^2/(2*v_p)
    u2 = sigma_v^2/(2*k)
    u3 = sqrt(2*k*v_p)

    for(i in 1:iterations){
        mu[i+1] = update_mu(mu,x,data[2:n],data[1:n-1],delta,z[2:n],v_star[2:n],vt[2:n],vt[1:n-1],w[2:n],set)$mu
        x$mu = mu[i+1]
        
        gamma[i+1] = update_gamma(gamma,x,data[2:n],data[1:n-1],delta,z[2:n],v_star[2:n],vt[2:n],vt[1:n-1],w[2:n],set)$gamma
        x$gamma = gamma[i+1]
        
        sigma_j[i+1] = update_sigma_j(sigma_j,x,data[2:n],data[1:n-1],delta,z[2:n],v_star[2:n],vt[2:n],vt[1:n-1],w[2:n],set)$sigma_j
        x$sigma_j = sigma_j[i+1]
        
        rho[i+1] = update_rho(rho,x,data[2:n],data[1:n-1],delta,z[2:n],v_star[2:n],vt[2:n],vt[1:n-1],w[2:n],set)$rho
        x$rho = rho[i+1]
        
        k[i+1] = update_k(u1,k,x,data[2:n],data[1:n-1],delta,z[2:n],v_star[2:n],vt[2:n],vt[1:n-1],w[2:n],set)$k
        x$k = k[i+1]
        
        v_p[i+1] = update_v(u2,v_p,x,data[2:n],data[1:n-1],delta,z[2:n],v_star[2:n],vt[2:n],vt[1:n-1],w[2:n],set)$v
        x$v = v_p[i+1]
        
        sigma_v[i+1] = update_sigma_v(u3,sigma_v,x,data[2:n],data[1:n-1],delta,z[2:n],v_star[2:n],vt[2:n],vt[1:n-1],w[2:n],set)$sigma_v
        x$sigma_v = sigma_v[i+1]
        
        u1[i+1] = x$sigma_v^2/(2*x$v)
        u2[i+1] = x$sigma_v^2/(2*x$k)
        u3[i+1] = sqrt(2*x$k*x$v)
    }

    print(rejectionRate(mcmc(mu)))
    print(rejectionRate(mcmc(gamma)))
    print(rejectionRate(mcmc(sigma_j)))
    print(rejectionRate(mcmc(rho)))
    print(rejectionRate((mcmc(k))))
    print(rejectionRate(mcmc(v_p)))
    print(rejectionRate(mcmc(sigma_v)))
    return(list(gamma=gamma,mu=mu,sigma_j=sigma_j,rho=rho,k=k,v=v_p,sigma_v= sigma_v))
}



