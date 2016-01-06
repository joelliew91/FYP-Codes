#include "include.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <random>

using namespace std;

double accept(std::vector<double> hist){
    int n = hist.size();
    int count = 0;
    for(int i=1;i<n;i++){
        if(hist[i] != hist[i-1])
            count++;
    }
    return 1.0*count/n;
}
double normal(double s1,double s2){
    return sqrt(-2*log(s1))*cos(2*M_PI*s2);
}


double uniform(unsigned s){
    return s%10000/10000.0;
}

double mean(vector<double> history){
    int n = history.size();
    double sum = 0;
    for(int i =0;i<n;i++)
        sum += history[i];

    return sum/n;
}

double sd(vector<double> history){
    int n = history.size();
    if(n==1) return 1.0;
    double m = mean(history);
    double sum = 0;
    for(int i=0;i<n;i++)
        sum += pow(history[i]-m,2);
    
    if(sum==0) return 1.0;
    
    return sqrt(sum/(n-1));
}



double update_sigma_v(){
    double std = sd(sigma_v_h);
    double s = 2.0;
    double old_post = 0;
    double new_post = 0;
    double s1 = uniform(rand());
    double s2 = uniform(rand());
    double fac = normal(s1,s2);
    double u_s = sqrt(2*k*v_p);
    double new_temp = exp(log(sigma_v/(u_s-sigma_v)) + fac*s*std);
    new_temp = new_temp*u_s/(new_temp+1);
    double temp = sigma_v;
    for(int i=1;i<total;i++){
        old_post += posterior(z[i],v[i],v[i-1],v_star[i],price[i],price[i-1],w[i]);
        sigma_v = new_temp;
        new_post += posterior(z[i],v[i],v[i-1],v_star[i],price[i],price[i-1],w[i]);
        sigma_v = temp;
    }
    double u = rand()%1000/(1.0*1000);
    double R = exp(new_post-old_post+log(new_temp*(u_s-new_temp)/u_s)-log(temp*(u_s-temp)/u_s));
    
    if(R>u)
        sigma_v = new_temp;
    
    return log(sigma_v/(u_s-sigma_v));
}

double update_v_p(){
    double std = sd(v_p_h);
    double s = 2.0;
    double old_post = 0;
    double new_post = 0;
    double s1 = uniform(rand());
    double s2 = uniform(rand());
    double fac = normal(s1,s2);
    double u_s = pow(sigma_v,2)/(2*k);
    double new_temp = exp(log(v_p-u_s) + fac*s*std)+u_s;
    double temp = v_p;
    for(int i=1;i<total;i++){
        old_post += posterior(z[i],v[i],v[i-1],v_star[i],price[i],price[i-1],w[i]);
        v_p = new_temp;
        new_post += posterior(z[i],v[i],v[i-1],v_star[i],price[i],price[i-1],w[i]);
        v_p = temp;
    }
    double u = rand()%1000/(1.0*1000);
    double R = exp(new_post-old_post+log(new_temp-u_s)-log(temp-u_s));
    
    if(R>u)
        v_p = new_temp;
    
    return log(v_p-u_s);
}

double update_k(){
    double std = sd(k_h);
    double s = 2.0;
    double old_post = 0;
    double new_post = 0;
    double s1 = uniform(rand());
    double s2 = uniform(rand());
    double fac = normal(s1,s2);
    double u_s = pow(sigma_v,2)/(2*v_p);
    double new_temp = exp(log(k-u_s) + fac*s*std)+u_s;
    double temp = k;
    for(int i=1;i<total;i++){
        old_post += posterior(z[i],v[i],v[i-1],v_star[i],price[i],price[i-1],w[i]);
        k = new_temp;
        new_post += posterior(z[i],v[i],v[i-1],v_star[i],price[i],price[i-1],w[i]);
        k = temp;
    }
    double u = rand()%1000/(1.0*1000);
    double R = exp(new_post-old_post+log(new_temp-u_s)-log(temp-u_s));
    
    if(R>u)
        k = new_temp;
    
    return log(k-u_s);
}

double update_rho(){
    double std = sd(rho_h);
    double s = 1.0;
    double old_post = 0;
    double new_post = 0;
    double s1 = uniform(rand());
    double s2 = uniform(rand());
    double fac = normal(s1,s2);
    double new_temp = exp(log((1+rho)/(1-rho)) + fac*s*std);
    new_temp = (new_temp-1)/(1+new_temp);
    double temp = rho;
    for(int i=1;i<total;i++){
        old_post += posterior(z[i],v[i],v[i-1],v_star[i],price[i],price[i-1],w[i]);
        rho = new_temp;
        new_post += posterior(z[i],v[i],v[i-1],v_star[i],price[i],price[i-1],w[i]);
        rho = temp;
    }
    double u = rand()%1000/(1.0*1000);
    double R = exp(new_post-old_post+log(new_temp+1)-log(temp+1));
    
    if(R>u)
        rho = new_temp;
    
    return log((1+rho)/(1-rho));
}

double update_sigma_j(){
    double std = sd(sigma_j_h);
    double s = 1.0;
    double old_post = 0;
    double new_post = 0;
    double s1 = uniform(rand());
    double s2 = uniform(rand());
    double fac = normal(s1,s2);
    double new_temp = exp(log(sigma_j) + fac*s*std);
    double temp = sigma_j;
    for(int i=1;i<total;i++){
        old_post += posterior(z[i],v[i],v[i-1],v_star[i],price[i],price[i-1],w[i]);
        sigma_j = new_temp;
        new_post += posterior(z[i],v[i],v[i-1],v_star[i],price[i],price[i-1],w[i]);
        sigma_j = temp;
    }
    double u = rand()%1000/(1.0*1000);
    double R = exp(new_post-old_post+log(new_temp)-log(temp));
    
    if(R>u)
        sigma_j = new_temp;
    
    return log(sigma_j);
}

double update_gam(){
    double std = sd(gam_h);
    double s = 1.0;
    double old_post = 0;
    double new_post = 0;
    double s1 = uniform(rand());
    double s2 = uniform(rand());
    double fac = normal(s1,s2);
    double new_temp = gam + fac*s*std;
    double temp = gam;
    for(int i=1;i<total;i++){
        old_post += posterior(z[i],v[i],v[i-1],v_star[i],price[i],price[i-1],w[i]);
        gam = new_temp;
        new_post += posterior(z[i],v[i],v[i-1],v_star[i],price[i],price[i-1],w[i]);
        gam = temp;
    }
    double u = rand()%1000/(1.0*1000);
    double R = exp(new_post-old_post);
    
    if(R>u)
        gam = new_temp;
    
    return gam;
}

double update_mu(){
    double std = sd(mu_h);
    double s = 1.0;
    double old_post = 0;
    double new_post = 0;
    double s1 = uniform(rand());
    double s2 = uniform(rand());
    double fac = normal(s1,s2);
    double new_temp = mu + fac*s*std;
    double temp = mu;
    for(int i=1;i<total;i++){
        old_post += posterior(z[i],v[i],v[i-1],v_star[i],price[i],price[i-1],w[i]);
        mu = new_temp;
        new_post += posterior(z[i],v[i],v[i-1],v_star[i],price[i],price[i-1],w[i]);
        mu = temp;
    }
    double u = rand()%1000/(1.0*1000);
    double R = exp(new_post-old_post);
    
    if(R>u)
        mu = new_temp;

    return mu;
}