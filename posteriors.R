post_svvg_fn_v<-function(x,para,delta){
    M = 1/10
    m = 10
    n = length(para$Gt)*delta
    val1 = sum(log(para$Gt[1:n-1]))*delta/x - T*log(x)*(delta/x)-T*lgamma(delta/x)
    val2 = -(m+1)*log(x)-1/x*(1/M+sum(para$Gt[1:n-1]))
    return(val1+val2)
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


svvg_fn_Gt<-function(x,para,delta){
    n = delta*length(data)
    val2 = -para$Jt[1:n-1]^2/(x*2*para$sigma^2)
    val1 = (delta/para$v-1.5)*log(x)
    val3 = -x*(1/para$v+para$gamma^2/(2*para$sigma^2))
    return(val1+val2+val3)
}

svvg_fn_log_mu_post<-function(curr_val,para){
    val = -log(para$sigma)-0.5*((curr_val-para$mu)/para$sigma)^2
    return(val)
}

svvg_fn_log_k_post<-function(curr_val,para,u){
    k = exp(curr_val)+u
    val = -log(para$sigma)-0.5*((k-para$mu)/para$sigma)^2#+curr_val
    return(val)
}