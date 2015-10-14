
sim_svvg<-function(years){
      y = c(1)                                        ## Y0 set at 0
      vol = c(1)                                      ## Volatility at t=0 is set at 1
      delta = 1/250
      mu = 0.05
      k = 0.015
      theta = 0.8
      sigma.v = 0.1
      rho = -0.4
      gam = -0.01
      sigma.y = 0.4
      v = 3.0
      sigma = 0
      E = matrix(c(1,rho*sigma.v,rho*sigma.v,sigma.v^2),2,2)
      epsilon = mvrnorm(years*250,mu=c(0,0),Sigma=E)  ## Generate 5000 moves ie 20yrs in the 
                                                      ## BM of log price and volatility
      
      epsilon.j = rnorm(years*250)                    ## Generate the Jump outcomes
      Gamma_path = rgamma(years*250,delta/v,v)
      J = gam*Gamma_path + sigma*sqrt(Gamma_path)*epsilon.j
      
      for(i in 1:(years*250)){
            Y = y[i] + mu*delta + sqrt(vol[i]*delta)*epsilon[i,1] + J[i]
            V = vol[i] + k*(theta-vol[i])*delta + sigma.v*sqrt(vol[i]*delta)*epsilon[i,2]
            y = rbind(y,Y)
            vol = rbind(vol,V)
      }
      ## plot(y,type='l',ylab='Log Stock Price',xlab = 'Time',main='SVVG')
      return(list(y=y,vol=vol))
}

est_svvg<-function(iter){

      return(para)
}

svvg_fn_v<-function(x){
      M = 1/10
      m = 10
      G = para$G
      n = length(G)
      val1 = ((cumprod(G)[n])^(delta/x))*(1/(gamma(delta/x)*x^(delta/x)))^T
      val2 = exp((cumsum(G)[n]+1/M)*-1/x)*(1/x)^(m+1)
      return(log(val1+val2))
}

svvg_fn_G<-function(x){
      va1 = x^(delta/v)*exp(para$J[i-1]^2/(x*2*para$sigma[i-1]^2))
      va2 = exp(-x*(para$gamma[i-1]^2/(2*para$sigma[i-1]^2)+1/v))
      return(log(va2*va1))
}

## Initiation
library(MASS)                                   ## library required to generate multivariate norm rv
library(MCMCpack)                               ## library required to generate inverse gamma rv
library(HI)                                     ## library required for Adaptive Rejection Metrop Sampling
set.seed(100)
delta = 1/250
para = as.data.frame(cbind(0,0,0.01,0.02,0.5,0,1,0.1,0,1)) ## Choose values randomly but within the domain of posterior
names(para) = c("mu","k","theta","sigma.v","rho","gamma","v","sigma","J","G")
T = 5000

data = sim_svvg(round(T/250))$y
for(i in 2:round(T/250)){
      ## posterior of v
      v = arms(runif(1,0,1),svvg_fn_v, function(x) (x>0)*(x<1),1)
      ## end of v
      print("v:")
      print(v)
      ## posterior of G
      G = arms(runif(1,0,50),svvg_fn_G, function(x) (x>0)*(x<50),1)
      ## end of G
      print("G:")
      print(G)
      ## posterior of J
      W = 1/((1-para$rho[i-1]^2)*para$v[i-1]*delta)+1/(G*para$sigma[i-1]^2)
      C = data[i]-data[i-1] - para$mu[i-1]*delta
      D = v - para$v[i-1] - para$k[i-1]*(para$theta[i-1]-para$v[i-1])*delta
      S = 1/((1-para$rho[i-1]^2)*para$v[i-1]*delta)*(C-para$rho[i-1]*D/para$sigma.v[i-1])+para$gamma[i-1]/para$sigma[i-1]^2
      J = rnorm(1,mean=S/W,sd = 1/sqrt(W))
      ## end of J
      
      ## posterior of sigma
      
      ## end for sigma
      ## Posterior of N
      ## End of N
      ## Posterior of mu
      ## End for mu
}