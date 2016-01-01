update_mu<-function(mu,x,yt,yu,delta,z,v,vt,vu,w,set){
    std = sd(mu)
    if(is.na(std)||std==0)
        std = 1
    proposal = list(sd=std,scale=0.5)
    new = normrwmetrop_mu(posterior,proposal,x,yt,yu,delta,z,v,vt,vu,w,set)
    return(new)
}

update_gamma<-function(gamma,x,yt,yu,delta,z,v,vt,vu,w,set){
    std = sd(gamma)
    if(is.na(std)||std==0)
        std = 1
    proposal = list(sd=std,scale=0.7)
    new = normrwmetrop_gamma(posterior,proposal,x,yt,yu,delta,z,v,vt,vu,w,set)
    return(new)
}

update_sigma_j<-function(sigma_j,x,yt,yu,delta,z,v,vt,vu,w,set){
    std = sd(log(sigma_j))
    if(is.na(std)||std==0)
        std = 1
    proposal = list(sd=std,scale=0.5)
    new = normrwmetrop_sigma_j(posterior,proposal,x,yt,yu,delta,z,v,vt,vu,w,set)
    return(new)
}

update_rho<-function(rho,x,yt,yu,delta,z,v,vt,vu,w,set){
    std = sd(log((rho+1)/(1-rho)))
    if(is.na(std)||std==0)
        std = 1
    proposal = list(sd=std,scale=0.1)
    new = normrwmetrop_rho(posterior,proposal,x,yt,yu,delta,z,v,vt,vu,w,set)
    return(new)

}

update_k<-function(u1,k,x,yt,yu,delta,z,v,vt,vu,w,set){
    std = sd(log(k-u1))

    if(is.na(std)||std==0)
        std = 1
    proposal = list(sd=std,scale=1)
    new = normrwmetrop_k(posterior,proposal,x,yt,yu,delta,z,v,vt,vu,w,set)
    return(new)
}

update_v<-function(u2,v_p,x,yt,yu,delta,z,v,vt,vu,w,set){
    std = sd(log(v_p-u2))

    if(is.na(std)||std==0)
        std = 1
    proposal = list(sd=std,scale=1)
    new = normrwmetrop_v(posterior,proposal,x,yt,yu,delta,z,v,vt,vu,w,set)
    return(new)
}

update_sigma_v<-function(u3,sigma_v,x,yt,yu,delta,z,v,vt,vu,w,set){
    std = sd(log(sigma_v/(u3-sigma_v)))
    if(is.na(std)||std==0)
        std = 1
    proposal = list(sd=std,scale=1)
    new = normrwmetrop_sigma_v(posterior,proposal,x,yt,yu,delta,z,v,vt,vu,w,set)
    return(new)
}

update_z<-function(z_list,x,yt,yu,delta,z,v,vt,vu,w,set){
    s = 1
    n = dim(z_list)[2]

    for(i in 2:n){
        #print(i)
        std = sd(z_list[,i])
        if(is.na(std)||std==0)
            std = 1

        log_prev = posterior_latent_z(x,yt[i-1],yu[i-1],delta,z[i-1],v[i-1],vt[i-1],vu[i-1],w[i-1],set)
        old = z[i-1]
        new = z[i-1] + rnorm(1)*s*std
        z[i-1] = new
        R = exp(posterior_latent_z(x,yt[i-1],yu[i-1],delta,z[i-1],v[i-1],vt[i-1],vu[i-1],w[i-1],set)-log_prev)
        z[i-1] = old
        if (is.na(R) == FALSE)
            if(R>runif(1))
                z[i-1] = new
    }
    return(z)
}

update_w<-function(w_list,x,yt,yu,delta,z,v,vt,vu,w,set){
    s = 1
    n = dim(w_list)[2]
    
    for(i in 2:n){
        #print(i)
        std = sd(log(w_list[,i]))
        if(is.na(std)||std==0)
            std = 1
        #print('log_prev')
        log_prev = posterior_latent_w(x,yt[i-1],yu[i-1],delta,z[i-1],v[i-1],vt[i-1],vu[i-1],w[i-1],set)+log(w[i-1])
        old = w[i-1]
        new = exp(log(w[i-1]) + rnorm(1)*s*std)
        w[i-1] = new
        #print('new')
        R = exp(posterior_latent_w(x,yt[i-1],yu[i-1],delta,z[i-1],v[i-1],vt[i-1],vu[i-1],w[i-1],set)-log_prev+log(w[i-1]))
        w[i-1] = old
        if (is.na(R) == FALSE)
        if(R>runif(1))
        w[i-1] = new
    }
    return(w)
}

update_vt<-function(v_list,x,yt,yu,delta,z,v,vt,w,set){
    s = 1
    n = dim(v_list)[2]
    i=2
    while(i<n){
        #print(v[i-1:i])
        #print(i)
        std = sd(log(v_list[,i]))
        if(is.na(std)||std==0)
            std = 1
        log_prev = posterior_latent_vt(x,c(yt[i-1],yt[i]),c(yu[i-1],yu[i]),delta,c(z[i-1],z[i]),c(v[i-1],v[i]),c(vt[i],vt[i+1]),c(vt[i-1],vt[i]),c(w[i-1],w[i]),set)+log(vt[i])
        old = vt[i]
        new = exp(log(vt[i])+rnorm(1)*s*std)
        vt[i] = new
        #print('new')
        R = exp(posterior_latent_vt(x,c(yt[i-1],yt[i]),c(yu[i-1],yu[i]),delta,c(z[i-1],z[i]),c(v[i-1],v[i]),c(vt[i],vt[i+1]),c(vt[i-1],vt[i]),c(w[i-1],w[i]),set)-log_prev+log(vt[i]))
        vt[i] = old
        if (is.na(R) == FALSE)
            if(R>runif(1))
                vt[i] = new
        i = i+1
    }
    return(vt)
}

update_vs<-function(vs,x,yt,yu,delta,z,v,vt,vu,w,set){
    s = 1
    n = dim(vs)[2]
    for(i in 2:n){
        std = sd(log(vs[,i]))
        if(is.na(std)||std==0)
            std = 1
        log_prev = posterior_latent_vs(x,yt[i-1],yu[i-1],delta,z[i-1],v[i-1],vt[i-1],vu[i-1],w[i-1],set)+log(v[i-1])
        old = v[i-1]
        new = exp(log(v[i-1]) + rnorm(1)*s*std)
        v[i-1] = new
        #print('new')
        R = exp(posterior_latent_vs(x,yt[i-1],yu[i-1],delta,z[i-1],v[i-1],vt[i-1],vu[i-1],w[i-1],set)-log_prev+log(v[i-1]))
        #print(R)
        v[i-1] = old
        if (is.na(R) == FALSE)
            if(R>runif(1))
                v[i-1] = new

    }
    return(v)

}




