library(MASS)                                   ## library required to generate multivariate norm rv
library(MCMCpack)                               ## library required to generate inverse gamma rv
library(HI) 
generate_stable<-function(n,alpha,beta,gamma=0,delta=1){
      k<-function(a){
            1-abs(1-a)
      }
      phi = runif(n,min=-0.5*pi,max = 0.5*pi)
      phi0 = -0.5*pi*beta*k(alpha)/alpha
      w = rexp(n)
      ans = sin(alpha*(phi-phi0))/(cos(phi)^(1/alpha))
      ans = ans*(cos(phi-alpha*(phi-phi0))/w)^(1/alpha-1)
      
      return(ans*delta+gamma)
}
combine<-function(years){
      y = list(svmj=1,svcj=1,svls=1,svvg=1)                                        ## Y0 set at 0
      vol = list(svmj=1,svcj=1,svls=1,svvg=1)                                   ## Volatility at t=0 is set at 1
      delta = 1
      mu = 0.05
      k = 0.015
      theta = 0.8
      sigma.v = 0.1
      rho = -0.4
      mu.y = -3
      sigma.y = 3.5
      sigma.y.inf = 0.4
      lam.y = 0.015
      mu.v = 1
      rho.j = -0.4
      a = 1.8
      sigma = 0.3
      beta = -1
      gam = -0.01
      v = 3
      E = matrix(c(1,rho*sigma.v,rho*sigma.v,sigma.v^2),2,2)
      
      epsilon = mvrnorm(years*250,mu=c(0,0),Sigma=E)  ## Generate 5000 moves ie 20yrs in the 
                                                      ## BM of log price and volatility
      swig.v = rexp(years*250,rate = mu.v)
      swig = rnorm(years*250,mean = mu.y,sd = sigma.y)
      swig.y = swig+rho.j*swig.v
      jump = rbinom(years*250,size=1,prob = lam.y*delta)
      J = jump*swig
      J.y = jump*swig.y
      J.v = jump*swig.v
      
      for(i in 1:(years*250)){
            Y = y$svmj[i] + mu*delta + sqrt(vol$svmj[i]*delta)*epsilon[i,1] + J[i]
            V = vol$svmj[i] + k*(theta-vol$svmj[i])*delta + sigma.v*sqrt(vol$svmj[i]*delta)*epsilon[i,2]
            y$svmj[i+1] = Y
            vol$svmj[i+1] = V
      }
      
      for(i in 1:(years*250)){
            Y = y$svcj[i] + mu*delta + sqrt(vol$svcj[i]*delta)*epsilon[i,1] + J.y[i]
            V = vol$svcj[i] + k*(theta-vol$svcj[i])*delta + sigma.v*sqrt(vol$svcj[i]*delta)*epsilon[i,2]+J.v[i]
            y$svcj[i+1] = Y
            vol$svcj[i+1] = V
      }
      x = seq(1,length(y$svmj),by=1)
      
      epsilon.j = rnorm(years*250)                    ## Generate the Jump outcomes
      Gamma_path = rgamma(years*250,delta/v,v)
      J = gam*Gamma_path + sigma*sqrt(Gamma_path)*epsilon.j
      
      for(i in 1:(years*250)){
            Y = y$svvg[i] + mu*delta + sqrt(vol$svvg[i]*delta)*epsilon[i,1] + J[i]
            V = vol$svvg[i] + k*(theta-vol$svvg[i])*delta + sigma.v*sqrt(vol$svvg[i]*delta)*epsilon[i,2]
            y$svvg[i+1] = Y
            vol$svvg[i+1] = V
      }
      log_stable_path = generate_stable(years*250,alpha=a,beta=-1,gamma=0,delta=sigma*delta^(1/a))
      for(i in 1:(years*250)){
            Y = y$svls[i] + mu*delta + sqrt(vol$svls[i]*delta)*epsilon[i,1] + log_stable_path[i]
            V = vol$svls[i] + k*(theta-vol$svls[i])*delta + sigma.v*sqrt(vol$svls[i]*delta)*epsilon[i,2]
            y$svls[i+1] = Y
            vol$svls[i+1] = V
      }
      
      plot(x,y$svmj,type='l',ylim=range(c(y$svmj,y$svcj,y$svvg,y$svls)),col='black',ylab='Log Stock Price',xlab = 'Time',main='Plots of SVMJ,SVCJ,SVVG,SVLS')
      par(new=T)
      plot(x,y$svcj,ylim=range(c(y$svmj,y$svcj,y$svvg,y$svls)),xlab='',ylab='',type='l',col='red',axes=FALSE)
      par(new=T)
      plot(x,y$svvg,ylim=range(c(y$svmj,y$svcj,y$svvg,y$svls)),xlab='',ylab='',type='l',col='blue',axes=FALSE)
      par(new=T)
      plot(x,y$svls,ylim=range(c(y$svmj,y$svcj,y$svvg,y$svls)),xlab='',ylab='',type='l',col='purple',axes=FALSE)
      legend("topleft",legend = c("SVMJ","SVCJ","SVVG","SVLS"),col=c('black','red','blue','purple'),lwd=c(3,3,3,3))
      
      return()
}

