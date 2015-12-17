normrwmetrop<-function(logf,proposal,start,...){
    val = start
    log_prev = logf(start,...)
    std = proposal$sd
    s = proposal$scale
    #print(val)
    new = val + rnorm(1)*s*std
    #print(new)
    R = exp(logf(new,...)-log_prev)
    #print(log_prev)
    #print(R)
    if (is.na(R) == FALSE)
        if(R>runif(1))
            val = new
    
    
    return(val)

}




