library(MASS)  
svmj_vt<-function(x){
      ep_y_1 = (diff_y[1] - mu*delta - N*swig) / sqrt(delta*vt[i-1])
      ep_v_1 = (x - vt[i-1] - k*delta*(theta-vt[i-1]))/(sigma.v*sqrt(vt[i-1]*delta))
      val1 = -(-2*rho*ep_y*ep_v+ep_v^2)/(2*(1-rho^2))
      val2 = -log(x)-()/(2*(1-rho^2))
      return(val1+val2)
}



test<-function(years){
      source("~/Desktop/Year 4 sem 1/FYP/Codes/svmj.R")
      data = sim_svmj(years)
      diff_y = data[2:length(data)] - data[1:length(data)-1]
      vt = c(1)
      delta = 1/250
      E = matrix(c(1,rho*sigma.v,rho*sigma.v,sigma.v^2),2,2)
      ep_v =  
      for(i in 2:length(data)){
            for(iter in 1:50000){
                  
            }
      }
          
      
}