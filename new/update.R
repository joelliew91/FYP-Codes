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


