generate_stable<-function(n,alpha,beta,gamma,delta){
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
sim_svmj<-function(years){
      y = c(1)                                        ## Y0 set at 0
      vol = c(1)                                      ## Volatility at t=0 is set at 1
      delta = 1/250
      mu = 0.05
      k = 0.015
      theta = 0.8
      sigma.v = 0.1
      rho = -0.4
      a = 1.8
      sigma = 0.3
      beta = -1
      
      E = matrix(c(1,rho*sigma.v,rho*sigma.v,sigma.v^2),2,2)
      
      epsilon = mvrnorm(years*250,mu=c(0,0),Sigma=E)  ## Generate 5000 moves ie 20yrs in the 
                                                      ## BM of log price and volatility
      log_stable_path = generate_stable(years*250,alpha=a,beta=-1,gamma=0,delta=sigma*delta^(1/a))
      
      
      for(i in 1:(years*250)){
            Y = y[i] + mu*delta + sqrt(vol[i]*delta)*epsilon[i,1] + log_stable_path[i]
            V = vol[i] + k*(theta-vol[i])*delta + sigma.v*sqrt(vol[i]*delta)*epsilon[i,2]
            y = rbind(y,Y)
            vol = rbind(vol,V)
      }
      #par(mfrow=c(2,1))
      #plot(y,type='l',ylab='Log Stock Price',xlab = 'Time',main='SVMJ')
      #plot(exp(y),type='l',ylab='Stock Price',xlab = 'Time',main = 'Stock Price using SVMJ')
      return(list(y=y,vol=vol))
}
