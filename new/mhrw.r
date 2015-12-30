normrwmetrop_mu<-function(logf,proposal,x,...){
    val = x
    log_prev = logf(x,...)
    std = proposal$sd
    s = proposal$scale
    
    x$mu = x$mu + rnorm(1)*s*std
    R = exp(logf(x,...)-log_prev)
    
    if (is.na(R) == FALSE)
        if(R>runif(1))
            val = x
    return(val)

}

normrwmetrop_gamma<-function(logf,proposal,x,...){
    val = x
    log_prev = logf(x,...)
    std = proposal$sd
    s = proposal$scale
    
    x$gamma = x$gamma + rnorm(1)*s*std
    R = exp(logf(x,...)-log_prev)
    
    if (is.na(R) == FALSE)
    if(R>runif(1))
    val = x
    return(val)
    
}

normrwmetrop_sigma_j<-function(logf,proposal,x,...){
    val = x
    log_prev = logf(x,...)+log(x$sigma_j)
    std = proposal$sd
    s = proposal$scale
    x$sigma_j = exp(log(x$sigma_j) + rnorm(1)*s*std)
    R = exp(logf(x,...)-log_prev+log(x$sigma_j))
    
    if (is.na(R) == FALSE)
    if(R>runif(1))
    val = x
    return(val)
    
}

normrwmetrop_rho<-function(logf,proposal,x,...){
    val = x
    log_prev = logf(x,...)+log(x$rho+1)
    std = proposal$sd
    s = proposal$scale
    x$rho = exp(log((x$rho+1)/(1-x$rho)) + rnorm(1)*s*std)
    x$rho = (x$rho-1)/(x$rho+1)
    R = exp(logf(x,...)-log_prev+log(x$rho+1))
    
    if (is.na(R) == FALSE)
    if(R>runif(1))
    val = x
    return(val)
    
}

normrwmetrop_k<-function(logf,proposal,x,...){
    val = x
    log_prev = logf(x,...)+log(x$k-x$sigma_v^2/(2*x$v))
    std = proposal$sd
    s = proposal$scale
    x$k = exp(log(x$k-x$sigma_v^2/(2*x$v)) + rnorm(1)*s*std)+x$sigma_v^2/(2*x$v)
    R = exp(logf(x,...)-log_prev+log(x$k-x$sigma_v^2/(2*x$v)))

    if (is.na(R) == FALSE)
    if(R>runif(1))
    val = x
    
    return(val)
    
}

normrwmetrop_v<-function(logf,proposal,x,...){
    val = x
    log_prev = logf(x,...)+log(x$v-x$sigma_v^2/(2*x$k))
    std = proposal$sd
    s = proposal$scale
    x$v = exp(log(x$v-x$sigma_v^2/(2*x$k)) + rnorm(1)*s*std)+x$sigma_v^2/(2*x$k)
    R = exp(logf(x,...)-log_prev+log(x$v-x$sigma_v^2/(2*x$k)))
    
    if (is.na(R) == FALSE)
    if(R>runif(1))
    val = x
    
    return(val)
    
}

normrwmetrop_sigma_v<-function(logf,proposal,x,...){
    val = x
    u =sqrt(2*x$k*x$v)
    log_prev = logf(x,...)+log(x$sigma_v*(u-x$sigma_v)/u)
    std = proposal$sd
    s = proposal$scale
    x$sigma_v = exp(log(x$sigma_v/(u-x$sigma_v))+ rnorm(1)*s*std)
    x$sigma_v = x$sigma_v*u/(x$sigma_v+1)
    R = exp(logf(x,...)-log_prev+log(x$sigma_v*(u-x$sigma_v)/u))
    
    if (is.na(R) == FALSE)
    if(R>runif(1))
    val = x
    
    return(val)
    
}



