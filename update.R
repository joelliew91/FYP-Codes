get_mu<-function(data,para,delta){
    n = length(data)
    dY = data[2:n] - data[1:n-1]
    dVt = para$vt[2:n] - para$vt[1:n-1]
    dJt = para$Jt[2:n] - para$Jt[1:n-1]
    
    C = dY - para$mu*delta - dJt
    D = dVt - para$k*(para$theta-para$vt[1:n-1])*delta

    W = delta/(1-para$rho^2)*sum(para$vt)+1
    S = 1/(1-para$rho^2) * sum(1/para$vt[1:n-1]*(C-para$rho*D/para$sigma.v))

    return(S/W)
    
}

get_theta<-function(data,para,delta){
    n = length(data)
    dY = data[2:n] - data[1:(n-1)]
    dVt = para$vt[2:n] - para$vt[1:n-1]
    dJt = para$Jt[2:n] - para$Jt[1:n-1]
    
    C = dY - para$mu*delta - dJt
    D = dVt + para$vt[1:n-1]*para$k*delta
    W = para$k^2*delta/(para$sigma.v^2*(1-para$rho^2))*sum(1/para$vt[1:n-1])+1
    S = para$k/(para$sigma.v*(1-para$rho^2))*sum((D/para$sigma.v-para$rho*C)/para$vt[1:n-1])
    
    u = runif(1)
    mn = S/W
    std = sqrt(1/W)
    f0 = pnorm(0,mean=mn,sd=std)
    theta = qnorm((1-f0)*u+f0,mean=mn,sd=std)
    
    return(theta)
}

get_k<-function(data,para,delta){
    n = length(data)
    dY = data[2:n] - data[1:(n-1)]
    dVt = para$vt[2:n] - para$vt[1:n-1]
    dJt = para$Jt[2:n] - para$Jt[1:n-1]
    
    C = dY - para$mu*delta - dJt
    D = dVt
    W = delta/(para$sigma.v^2*(1-para$rho^2))*sum(1/para$vt[1:n-1]*(para$theta-para$vt[1:n-1])^2)+1
    S = 1/(para$sigma.v*(1-para$rho^2))*sum((para$theta-para$vt[1:n-1])*(D/para$sigma.v-para$rho*C)/para$vt[1:n-1])
    
    u = runif(1)
    mn = S/W
    std = sqrt(1/W)
    f0 = pnorm(0,mean=mn,sd=std)
    k = qnorm((1-f0)*u+f0,mean=mn,sd=std)
    return(k)
}

get_w_v_phi<-function(data,para,delta){
    n = length(data)
    dY = data[2:n] - data[1:(n-1)]
    dVt = para$vt[2:n] - para$vt[1:n-1]
    dJt = para$Jt[2:n] - para$Jt[1:n-1]
    
    C = dY - para$mu*delta - dJt
    D = (dVt - para$k*(para$theta-para$vt[1:n-1]))*delta/sqrt(delta*para$vt[1:n-1])
    
    S = sum(C*D)
    W = sum(C^2)+2
    alpha = n/2 + 2
    beta = 1/(0.5*sum(D^2)+1/200-S^2/(2*W))
    
    return(list(w_v=beta/(alpha-1),phi = S/W))
}

get_gamma<-function(data,para,delta){
    n = length(data)
    W = sum(para$Gt[1:n*delta-1])/para$sigma^2+1
    S = sum(para$Jt[1:n*delta-1])/para$sigma^2
    
    return(S/W)
}

get_sigma<-function(data,para,delta){
    n = length(data)*delta
    a = n/2 + 2.5
    b = (0.5*sum((para$Jt[1:n-1]-para$theta*para$Gt[1:n-1])^2/para$Gt[1:n-1])+1/5)^-1
    
    return(b/(a-1))
}

get_v<-function(para,delta){
    v = mean(arms(para$v,post_svvg_fn_v,function(x,...) (x>0)*(x<10),1000,para,delta))
    return(v)
}

get_Jt<-function(data,para,delta){
    n = length(data)
    dY = data[2:n] - data[1:n-1]
    dVt = para$vt[2:n*delta] - para$vt[1:n*delta-1]

    C = dY - para$mu*delta
    D = dVt - para$k*delta*(para$theta-para$vt[1:n*delta-1])
    
    W = 1/(para$sigma^2*para$Gt[2:n*delta])+1/((1-para$rho^2)*para$vt[1:n*delta-1]*delta)
    S = 1/((1-para$rho^2)*para$vt[1:n*delta-1]*delta)*(C-para$rho*D/para$sigma.v)+para$gamma/para$sigma^2

    return(S/W)
}


update_Gt<-function(para,delta,Gt){
    for(i in 1:(length(para$Gt)-1))
        Gt <- arms(para$Gt[i], svvg_fn_Gt, function(x,para,delta,i) (x>0)*(x<1000), 5000,para,delta,i)
    return(Gt)
}






