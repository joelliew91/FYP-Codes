#include "include.h"
#include <math.h>
#include <vector>
#include <random>
#include <iostream>

using namespace std;

void accept_latent(double** hist){
    double max = 0;
    double min =1;
    for(int i=1;i<total;i++){
        int count = 0;
        for(int j=1;j<=iterations;j++)
            if(hist[i][j] != hist[i][j-1])
                count++;
        double rate = count/(1.0*iterations);
        if(rate>max)
            max = rate;
        if(rate<min)
            min = rate;
    }
    cout<<"min: "<<min<<" max: "<<max<<endl;
    return ;
}
double mean(int iter,double* history){

    double sum = 0;
    for(int i =0;i<iter;i++)
        sum += history[i];
    
    return sum/iter;
}

double sd(int iter,double* history){
    if (iter==1) return 1.0;
    
    double m = mean(iter,history);
    double sum = 0;
    for(int i=0;i<iter;i++)
        sum += pow(history[i]-m,2);
    
    if(sum==0) return 1.0;
    
    return sqrt(sum/(iter-1));
}

void update_latent_v(int iter,double** hist){
    hist[0][iter] = log(1);
    for(int i=1;i<total-1;i++){
        double std = sd(iter,hist[i]);
        double s = 1.0;
        double old_post = 0;
        double new_post = 0;
        double s1 = uniform(rand());
        double s2 = uniform(rand());
        double fac = normal(s1,s2);
        double new_temp = exp(log(v[i]) + fac*s*std);
        double temp = v[i];
        //cout<<new_temp<<endl;
        double upper = 2;
        if(i==total-1) upper = 1;
        if(i != total-1){
            old_post = posterior(z[i],v[i],v[i-1],v_star[i],price[i],price[i-1],w[i])+posterior(z[i+1],v[i+1],v[i],v_star[i+1],price[i+1],price[i],w[i+1]);
            v[i] = new_temp;
            new_post = posterior(z[i],v[i],v[i-1],v_star[i],price[i],price[i-1],w[i])+posterior(z[i+1],v[i+1],v[i],v_star[i+1],price[i+1],price[i],w[i+1]);
            v[i] = temp;
        }
        else{
            old_post = posterior(z[i],v[i],v[i-1],v_star[i],price[i],price[i-1],w[i]);
            v[i] = new_temp;
            new_post = posterior(z[i],v[i],v[i-1],v_star[i],price[i],price[i-1],w[i]);
            v[i] = temp;
        }
        //cout<<new_post<<endl;
        double u = rand()%1000/(1.0*1000);
        double log_R = new_post-old_post+log(new_temp)-log(temp);
        if(log_R>log(u))
            v[i] = new_temp;
        hist[i][iter] = log(v[i]);
    }
    
    
    return ;
}

void update_latent_v_s(int iter,double** hist){
    hist[0][iter] = log(1);
    for(int i=1;i<total;i++){
        double std = sd(iter,hist[i]);
        double s = 1.0;
        double old_post = 0;
        double new_post = 0;
        double s1 = uniform(rand());
        double s2 = uniform(rand());
        double fac = normal(s1,s2);
        double new_temp = exp(log(v_star[i]) + fac*s*std);
        double temp = v_star[i];
        //cout<<new_temp<<endl;
        old_post = posterior(z[i],v[i],v[i-1],v_star[i],price[i],price[i-1],w[i])+log(temp);
        v_star[i] = new_temp;
        new_post = posterior(z[i],v[i],v[i-1],v_star[i],price[i],price[i-1],w[i])+log(new_temp);
        v_star[i] = temp;
        
        double u = rand()%1000/(1.0*1000);
        double R = exp(new_post-old_post);
        if(R>u)
            v_star[i] = new_temp;
        hist[i][iter] = log(v_star[i]);
    }
    
    return ;
}

void update_latent_w(int iter,double** hist){
    hist[0][iter] = log(1);
    for(int i=1;i<total;i++){
        double std = sd(iter,hist[i]);
        double s = 1.0;
        double old_post = 0;
        double new_post = 0;
        double s1 = uniform(rand());
        double s2 = uniform(rand());
        double fac = normal(s1,s2);
        double new_temp = exp(log(w[i]) + fac*s*std);
        double temp = w[i];
        old_post = posterior(z[i],v[i],v[i-1],v_star[i],price[i],price[i-1],w[i])+log(temp);
        w[i] = new_temp;
        new_post = posterior(z[i],v[i],v[i-1],v_star[i],price[i],price[i-1],w[i])+log(new_temp);
        w[i] = temp;
        
        double u = rand()%1000/(1.0*1000);
        double R = exp(new_post-old_post);
        if(R>u)
            w[i] = new_temp;
        hist[i][iter] = log(w[i]);
    }
    
    return ;
}

void update_latent_z(int iter,double** hist){
    hist[0][iter] = 1;
    for(int i=1;i<total;i++){
        double std = sd(iter,hist[i]);
        double s = 1.0;
        double old_post = 0;
        double new_post = 0;
        double s1 = uniform(rand());
        double s2 = uniform(rand());
        double fac = normal(s1,s2);
        double new_temp = z[i] + fac*s*std;
        double temp = z[i];
        old_post = posterior(z[i],v[i],v[i-1],v_star[i],price[i],price[i-1],w[i]);
        z[i] = new_temp;
        new_post = posterior(z[i],v[i],v[i-1],v_star[i],price[i],price[i-1],w[i]);
        z[i] = temp;

        double u = rand()%1000/(1.0*1000);
        double R = exp(new_post-old_post);
        if(R>u)
            z[i] = new_temp;
        hist[i][iter]=z[i];
    }

    return ;
}