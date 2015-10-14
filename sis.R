sis<-function(n,df,iter){
     sam1 = rt(n,df=df)
     t_hist = dt(sam1,df=2)
     n_hist = dnorm(sam1)
     w = dnorm(sam1)/dt(sam1,df=2)
     w = cbind(w)
     w = w/sum(w)
     for(i in 2:iter){
           temp = rt(n,df=df)
           sam1 = cbind(sam1,temp)
           t_hist = cbind(t_hist,dt(temp,df=df))
           n_hist = cbind(n_hist,dnorm(temp))
           w_t = w[,i-1]*n_hist[,i]/t_hist[,i]
           w_t = w_t/sum(w_t)
           w = cbind(w,w_t)
     }
     return(list(w=w,th=t_hist,nh=n_hist))
}