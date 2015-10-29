sim_svmj<-function(years){
      y = c(1)                                        ## Y0 set at 0
      vol = c(1)                                      ## Volatility at t=0 is set at 1
      delta = 1/250
      mu = 0.05
      k = 0.015
      theta = 0.8
      sigma.v = 0.1
      rho = -0.4
      mu.y = -3
      sigma.y = 3.5
      lam.y = 0.015
      mu.v = 1
      rho.j = -0.4
      E = matrix(c(1,rho*sigma.v,rho*sigma.v,sigma.v^2),2,2)

      epsilon = mvrnorm(years*250,mu=c(0,0),Sigma=E)  ## Generate 5000 moves ie 20yrs in the 
      ## BM of log price and volatility
      swig.v = rexp(years*250,rate = mu.v)
      for(i in 1:length(swig.v)){
            swig.y = rnorm(1,mean = mu.y+rho.j*swig.v[i],sd = sigma.y)
      }
      
      jump = rbinom(years*250,size=1,prob = lam.y*delta)
      J.y = jump*swig.y
      J.v = jump*swig.v
      
      
      for(i in 1:(years*250)){
            Y = y[i] + mu*delta + sqrt(vol[i]*delta)*epsilon[i,1] + J.y[i]
            V = vol[i] + k*(theta-vol[i])*delta + sigma.v*sqrt(vol[i]*delta)*epsilon[i,2]+J.v[i]
            y = rbind(y,Y)
            vol = rbind(vol,V)
      }
      par(mfrow=c(2,1))
      plot(y,type='l',ylab='Log Stock Price',xlab = 'Time',main='SVCMJ')
      plot(exp(y),type='l',ylab='Stock Price',xlab = 'Time',main = 'Stock Price using SVCMJ')
      return(list(y=y,vol=vol))
}
