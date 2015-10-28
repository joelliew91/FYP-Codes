normal<-function(x,me,std){
      val = exp(-((x-me)/std)^2/2)/(sqrt(2*pi)*std)
      return(val)
}
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

svvg_fn_v<-function(x){
      M = 1/10
      m = 10
      G = para$G
      val1 = sum(log(G))*delta/x - T*log(x)*(delta/x)-T*lgamma(delta/x)
      val2 = -(m+1)*log(x)-1/x*(1/M+sum(G))
      return(val1+val2)
}

svvg_fn_G<-function(x){
    
      val2 = -para$J[i-1]^2/(x*2*para$sigma[i-1]^2)
      val1 = (delta/x-1.5)*log(x)
      val3 = -x*(1/v+para$gamma[i-1]^2/(2*para$sigma[i-1]^2))
      return(val1+val2+val3)
}

svvg_fn_vt<-function(x){
    val1 = 0
    val2 = 0

    if(i==2){
        e_y_2 = (data[i] - data[i-1] - para$mu[i-1]*delta)/sqrt(para$vt[i-1]*delta)
        e_v_2 = (x - para$vt[i-1] - delta*para$k[i-1]*(para$theta[i-1]-para$vt[i-1]))/(para$sigma.v[i-1]*sqrt(para$vt[i-1]*delta))
        val2 = -log(x)-0.5*(e_y_2^2-2*para$rho[i-1]*e_v_2*e_y_2+e_v_2^2)/(1-para$rho[i-1]^2)
    }
    else {
        e_y_1 = (data[i-1] - data[i-2] - para$mu[i-1]*delta)/sqrt(para$vt[i-2]*delta)
        e_y_2 = (data[i] - data[i-1] - para$mu[i-1]*delta)/sqrt(para$vt[i-1]*delta)
        e_v_1 = (para$vt[i-1] - para$vt[i-2] - delta*para$k[i-1]*(para$theta[i-1]-para$vt[i-1]))/(para$sigma.v[i-1]*sqrt(para$vt[i-2]*delta))
        e_v_2 = (x - para$vt[i-1] - delta*para$k[i-1]*(para$theta[i-1]-para$vt[i-1]))/(para$sigma.v[i-1]*sqrt(para$vt[i-1]*delta))
        val1 = -0.5*(-2*para$rho[i-1]*e_y_1+e_v_1^2)/(1-para$rho[i-1]^2)
        val2 = -log(x)-0.5*(e_y_2^2-2*para$rho[i-1]*e_v_2*e_y_2+e_v_2^2)/(1-para$rho[i-1]^2)
    }
    return(val1+val2)
}

## Initiation
library(MASS)                                   ## library required to generate multivariate norm rv
library(MCMCpack)                               ## library required to generate inverse gamma rv
library(HI)                                     ## library required for Adaptive Rejection Metrop Sampling
library(truncnorm)                              ## library required for truncated normal
set.seed(120)
delta = 1/250
para = as.data.frame(cbind(0.01,0.5,0.1,0.01,0.5,0,3,0.1,0,1,0.02)) ## Choose values randomly but within the domain of posterior
names(para) = c("mu","k","theta","sigma.v","rho","gamma","v","sigma","J","G","vt")
T = 5000

data = sim_svvg(round(T/250))$y
for(i in 2:T){
      para = rbind(para,rep(0,10))
      ## posterior of v
      v = arms(para$v[i-1],svvg_fn_v, function(x) (x>0)*(x<1000),50000)
      v = mean(v)
      para$v[i] = v
      ## end of v
      
      ## posterior of vt
      vt = arms(para$vt[i-1],svvg_fn_vt,function(x) (x>0)*(x=<1),50000)
      #hist(vt,prob=TRUE,main="vt")
      vt = mean(vt)
      para$vt[i] = vt
      ## end of vt
      
      diff_y = data[i] - data[i-1]
      print(diff_y)
      diff_v = para$vt[i] - para$vt[i-1]
      
      ## posterior of G
      G = arms(para$G[i-1],svvg_fn_G, function(x) (x>0)*(x<50),50000)
      G = mean(G)
      para$G[i] = G
      ## end of G

      ## posterior of J
      W = 1/((1-para$rho[i-1]^2)*para$vt[i-1]*delta)+1/(G*para$sigma[i-1]^2)
      C = diff_y - para$mu[i-1]*delta
      D = diff_v - para$k[i-1]*(para$theta[i-1]-para$vt[i-1])*delta
      S = 1/((1-para$rho[i-1]^2)*para$vt[i-1]*delta)*(C-para$rho[i-1]*D/para$sigma.v[i-1])+para$gamma[i-1]/para$sigma[i-1]^2
      J = S/W
      para$J[i] = J
      ## end of J
      
      ## posterior of sigma
      m = 2.5
      M = 5
      alpha = T/2 + m
      beta = (1/M + 0.5*sum((para$J - para$G*para$theta[i-1])^2/para$G))^-1
      sigma = beta/(alpha-1)
      para$sigma[i] = sigma
      ## end for sigma
      
      ## Posterior of gamma
      m = 0
      M = 1
      W = 1/sigma^2*sum(para$G)+1/M^2
      S = 1/sigma^2*sum(para$J)+m/M^2
      gam = S/W
      para$gamma[i] = gam
      ## end of gamma
      
      ## Posterior of mu
      m = 0
      M = 1
      C_n = data[2:i] - data[1:i-1] - delta*para$mu[i-1]
      D_n = para$vt[2:i] - para$vt[1:i-1] - para$k[i-1]*(para$theta[i-1]-para$v[1:i-1])*delta
      W = delta/(1-para$rho[i-1]^2)*sum(1/para$vt[1:i-1])+1/M^2
      S = sum((C_n-para$rho[i-1]*D_n/para$sigma.v[i-1])/para$vt[1:i-1])/(1-para$rho[i-1]^2)+m/M^2
      mu = S/W
      para$mu[i] = mu
      ## End for mu
      
      ## Posterior of Theta
      m = 0
      M = 1
      W = para$k[i-1]^2*delta/(para$sigma.v[i-1]^2*(1-para$rho[i-1]^2))*sum(1/para$vt[1:i-1])+1/M^2
      S = para$k[i-1]/(para$sigma.v[i-1]*(1-para$rho[i-1]^2))*sum((D_n/para$sigma.v[i-1]-para$rho[i-1]*C_n)/para$vt[1:i-1])+m/M^2
      u = runif(10000,min=0,max=1)
      me = S/W
      std = sqrt(1/W)
      print(me)
      print(std)
      sam_theta = qnorm(u*(1-pnorm(0,me,std))+pnorm(0,me,std),me,std)
      prob = dnorm(sam_theta,me,std)
      theta = sum(prob*sam_theta)/sum(prob)
      para$theta[i] = theta
      ## end for theta
      
      ## Posterior of rho and sigma.v
      m = 2
      M = 200
      C_n_adj = C_n/sqrt(delta*para$vt[1:i-1])
      D_n_adj = D_n/sqrt(delta*para$vt[1:i-1])
      W = sum(C_n_adj^2)+2
      S = sum(C_n_adj*D_n_adj)
      alpha = T/2 + m
      beta = (1/M - S^2/(2*W)+0.5*sum(D_n_adj^2))^-1
      w = beta/(alpha-1)
      phi = S/W
      sigma.v = w + phi^2
      rho = phi / sigma.v
      para$rho[i] = rho
      para$sigma.v[i] = sigma.v
      ## end of posterior of rho and sigma.v
      
      ## Posterior of k
      m = 0
      M = 1
      W = delta/((1-para$rho[i]^2)*para$sigma.v[i]^2)*sum((para$theta[i]-para$vt[1:i-1])^2/para$v[1:i-1])+1/M^2
      S = 1/(para$sigma.v[i-1]*(1-para$rho[i-1]^2))*sum(((para$theta[i]-para$v[1:i-1])*(D_n/para$sigma.v[i]-para$rho[i]*C_n))/para$vt[1:i-1])+m/M^2
      u = runif(10000,min=0,max=1)
      me1 = S/W
      std1 = 1/sqrt(W)
      print(me1)
      print(std1)
      sam_k = qnorm(u*(1-pnorm(0,me1,std1))+pnorm(0,me1,std1),me1,std1)
      prob1 = dnorm(sam_k,me1,std1)
      k = sum(prob1*sam_k)/sum(prob1)
      para$k[i] = k
      ## end of posterior of k
      print('HI')
}





## Posterior of N
#{A1 = (diff_y - para$mu[i-1]*delta - para$swig[i-1])/sqrt(delta*para$v[i-1])
#A2 = (diff_v - delta*para$mu[i-1])/sqrt(delta*para$v[i-1])
#B = (para$v[i] - para$v[i-1] - para$k[i-1]*(para$theta[i-1]-para$v[i-1])*delta)/(para$sigma[i]*sqrt(delta*para$v[i-1]))
#a1 = exp(-(A1^2-2*para$rho[i-1]*A1*B)/(2*(1-para$rho[i-1]^2)))
#a2 = exp(-(A2^2 - 2*para$rho[i-1]*A2*B)/(2*(1-para$rho[i-1]^2)))*(1-para$lambda[i-1])
#N = rbinom(1,1,a1/(a1+a2))
#para$N = N ## checked
## End of N

## Posterior of swig
#C = diff_y - para$mu[i-1]*delta
#D = diff_v - para$k[i-1]*(para$theta[i-1]-para$v[i-1])*delta
#S = para$N[i]*(C-para$rho[i-1]*D/para$sigma.v[i-1])/(delta*para$v[i-1]*(1-para$rho[i-1]^2)) + para$mu
#W =       
#swig = rnorm(1,mean=S/W,sd = 1/sqrt(W))
#para$swig[i] = swig
## end of swig