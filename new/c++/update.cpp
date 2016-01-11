#include "include.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <random>

using namespace std;

double scaler(std::vector<double> history){
    double acc;
    double scale;
    if(history.size()==1)
        scale = 1.0;
    else
        acc = accept(history);
    if(acc>0.85) scale = 5.0;
    if(acc<0.15) scale = 1/5;
    
    return scale;
}
double final_mean(std::vector<double> history,int cutoff){
    int n = history.size();
    double sum = 0;
    for(int i =cutoff;i<n;i++)
        sum += history[i];
    
    return sum/(1.0*n);
}

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
    flag = 1;
    double std = sd(sigma_v_h);
    double s = 1;//scaler(sigma_v_h);
    double s1 = uniform(rand());
    double s2 = uniform(rand());
    double fac = normal(s1,s2);
    double u_s = sqrt(2*k*v_p);
    double new_temp = exp(log(sigma_v/(u_s-sigma_v)) + fac*s*std);
    new_temp = new_temp*u_s/(new_temp+1);
    double d = 0.5*4*k*v_p/pow(new_temp,2);
    if(1<d && d<dmax){
        double temp = sigma_v;
        double old_post = posterior();
        sigma_v = new_temp;
        double new_post = posterior();
        sigma_v = temp;
        double u = rand()%1000/(1.0*1000);
        double R = exp(new_post-old_post+log(new_temp*(u_s-new_temp)/u_s)-log(temp*(u_s-temp)/u_s));
        if(R>u)
            if(flag)
            sigma_v = new_temp;
    }
    
    flag = 1;
    return log(sigma_v/(u_s-sigma_v));
}

double update_v_p(){
    flag = 1;
    double std = sd(v_p_h);
    double s = 1;//scaler(v_p_h);
    double s1 = uniform(rand());
    double s2 = uniform(rand());
    double fac = normal(s1,s2);
    double u_s = pow(sigma_v,2)/(2*k);
    double new_temp = exp(log(v_p-u_s) + fac*s*std)+u_s;
    double temp = v_p;
    double d = 0.5*4*k*new_temp/pow(sigma_v,2);
    if(1<d && d<dmax){
        double old_post = posterior();
        v_p = new_temp;
        double new_post = posterior();
        v_p = temp;
        double u = rand()%1000/(1.0*1000);
        double R = exp(new_post-old_post+log(new_temp-u_s)-log(temp-u_s));
        if(R>u)
            if(flag)
            v_p = new_temp;
        
    }
    return log(v_p-u_s);
}

double update_k(){
    flag = 1;
    double std = sd(k_h);
    double s = 1;//scaler(k_h);
    double s1 = uniform(rand());
    double s2 = uniform(rand());
    double fac = normal(s1,s2);
    double u_s = pow(sigma_v,2)/(2*v_p);
    double new_temp = exp(log(k-u_s) + fac*s*std)+u_s;
    double temp = k;
    double d = 0.5*4*new_temp*v_p/pow(sigma_v,2);
    if(1<d && d<dmax){
        double old_post = posterior();
        k = new_temp;
        double new_post = posterior();
        k= temp;
        double u = rand()%1000/(1.0*1000);
        double R = exp(new_post-old_post+log(new_temp-u_s)-log(temp-u_s));
        if(R>u)
            if(flag)
                k = new_temp;

    }

    
    return log(k-u_s);
}

double update_rho(){
    flag = 1;
    double std = sd(rho_h);
    double s = 1.2;//scaler(rho_h);
    double s1 = uniform(rand());
    double s2 = uniform(rand());
    double fac = normal(s1,s2);
    double new_temp = exp(log((1+rho)/(1-rho)) + fac*s*std);
    new_temp = (new_temp-1)/(1+new_temp);
    double temp = rho;
    double old_post = posterior();
    rho = new_temp;
    double new_post = posterior();
    rho = temp;

    double u = rand()%1000/(1.0*1000);
    double R = exp(new_post-old_post+log(new_temp+1)-log(temp+1));
    
    if(R>u)
        if(flag)
            rho = new_temp;
    flag = 1;
    return log((1+rho)/(1-rho));
}

double update_sigma_j(){
    flag = 1;
    double std = sd(sigma_j_h);
    double s = 0.2;//scaler(sigma_j_h);
    double s1 = uniform(rand());
    double s2 = uniform(rand());
    double fac = normal(s1,s2);
    double new_temp = exp(log(sigma_j) + fac*s*std);
    double temp = sigma_j;

    double old_post = posterior();
    sigma_j = new_temp;
    double new_post = posterior();
    sigma_j = temp;

    double u = rand()%1000/(1.0*1000);
    double R = exp(new_post-old_post+log(new_temp)-log(temp));
    
    if(R>u)
        if(flag)
            sigma_j = new_temp;

    return log(sigma_j);
}

double update_gam(){
    flag=1;
    double std = sd(gam_h);
    double s = 0.5;//scaler(gam_h);
    double s1 = uniform(rand());
    double s2 = uniform(rand());
    double fac = normal(s1,s2);
    
    double new_temp = gam + fac*s*std;
    double temp = gam;
    
    double old_post = posterior();
    gam = new_temp;
    double new_post = posterior();
    gam = temp;

    double u = rand()%1000/(1.0*1000);
    double R = exp(new_post-old_post);
    
    if(R>u)
        if(flag)
            gam = new_temp;
    return gam;
}

double update_mu(){
    flag = 1;
    double std = sd(mu_h);
    double s = 1;//scaler(mu_h);
    double s1 = uniform(rand());
    double s2 = uniform(rand());
    double fac = normal(s1,s2);
    double new_temp = mu + fac*s*std;
    double temp = mu;
    double old_post = posterior();
    mu = new_temp;
    double new_post = posterior();
    mu = temp;
    double u = rand()%1000/(1.0*1000);
    double R = exp(new_post-old_post);
    if(R>u)
        if(flag)
            mu = new_temp;
    return mu;
}