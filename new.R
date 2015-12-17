library(MASS)                                   ## library required to generate multivariate norm rv
library(MCMCpack)                               ## library required to generate inverse gamma rv
library(HI)                                     ## library required for Adaptive Rejection Metrop Sampling
library(truncnorm)                              ## library required for truncated normal

source('~/Desktop/FYP/Codes/posteriors.R') # load the posterior functions
source('~/Desktop/FYP/Codes/update.R')
source('~/Desktop/FYP/Codes/sim_svvg.R') # load the price simulator
source('~/Desktop/FYP/Codes/mhrw.R') # load the MH-random walk


# need to fix delta for the sim

mcmc_para<-function(data,delta=1,iterations){
    set.seed(150)
    n = length(data)
    
    Gt = rgamma(n*delta,1,1)    # latent variable
    Jt = rnorm(n*delta)     # latent variable
    vt = rep(1,n*delta) # latent variable
    mu = rnorm(1) # para 1
    #print(mu)
    k = abs(rnorm(1)) # para 2
    theta = abs(rnorm(1)) # para 3
    w_v = rinvgamma(1,2,200)
    phi = rnorm(1,mean=0,sd=sqrt(0.5*w_v))
    rho = 1/(1+w_v/phi^2) # para 5
    sigma.v = sqrt(w_v/(1-rho^2)) # para 4
    gamma = rnorm(1) # para 6
    sigma = sqrt(rinvgamma(1,2.5,5)) # 1/rgamma(1,2.5,5) # para 7
    v = 1/rgamma(1,10,1/10) # para 8
    para = list(mu=mu,k=k,theta=theta,sigma.v=sigma.v,rho=rho,gamma=gamma,sigma=sigma,v=v,vt=vt,w_v=w_v,phi=phi,Jt=Jt,Gt=Gt)
    iter = list(theta=theta,k=k)
    for(i in 1:iterations){
        
        para$mu = get_mu(data,para,delta)     ## Update mu g
        
        iter$theta[i+1] = get_theta(data,para,delta)     ## Update theta g
        iter$k[i+1] = get_k(data,para,delta)        ## Update k g
        
        para$theta = (para$theta*i+iter$theta[i+1])/(i+1)
        para$k = (para$k*i+iter$k[i+1])/(i+1)
        
        temp = get_w_v_phi(data,para,delta)
        para$w_v = temp$w_v
        para$phi = temp$phi
        para$sigma.v = sqrt(temp$w_v+temp$phi^2)
        para$rho = temp$phi/para$sigma.v
        
        para$gamma = get_gamma(data,para,delta)

        para$sigma = get_sigma(data,para,delta)   ## Update sigma ns
        para$v = get_v(para,delta)
        para$Jt[2:n*delta] = get_Jt(data,para,delta)               ## Update Jt ns
        #Gt = update_Gt(para,delta,Gt)
    }
    
    return(para)
}



