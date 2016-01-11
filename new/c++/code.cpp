#include <iostream>
#include <math.h>
#include "include.h"
#include <random>
#include <vector>
#include <chrono>
#include <complex_bessel.h>
#include <time.h>
using namespace std;

double mu;vector<double> mu_h;
double gam;vector<double> gam_h;
double sigma_j;vector<double> sigma_j_h;
double rho;vector<double> rho_h;
double k;vector<double> k_h;
double v_p;vector<double> v_p_h;
double sigma_v;vector<double> sigma_v_h;
double theta;vector<double> theta_h;
double *v;
double *v_star;
double *z;
double *w;
double max_v;
double max_w;
double min_v;
double min_w;
double dmax = 2;
int flag = 1;
int iterations=50000;
int delta=1;
vector<double> price;
int total;



int main(){
    srand(time(0));

    initialize();
    //complex<double> temp(-100000,-100000);
    //cout<<sp_bessel::besselI(10000000000,temp)<<endl;
    double ** latent_z = new double *[total];
    double ** latent_w = new double *[total];
    double ** latent_v = new double *[total];
    double ** latent_v_star = new double* [total];
    for(int i=0;i<total;i++){
        latent_z[i] = new double[iterations+1];
        latent_w[i] = new double[iterations+1];
        latent_v[i] = new double[iterations+1];
        latent_v_star[i] = new double[iterations+1];
    }

    for(int i=0;i<total;i++){
        latent_z[i][0] = 1;
        latent_w[i][0] = log(1);
        latent_v[i][0] = 1;
        latent_v_star[i][0] = 1;
    }

    for(int i=0;i<iterations;i++){
        cout<<"mu"<<endl;
        mu_h.push_back(update_mu());
        cout<<"gam"<<endl;
        gam_h.push_back(update_gam());
        cout<<"s_j"<<endl;
        sigma_j_h.push_back(update_sigma_j());
        cout<<"rho"<<endl;
        rho_h.push_back(update_rho());
        cout<<"k"<<endl;
        k_h.push_back(update_k());
        cout<<"v_p"<<endl;
        v_p_h.push_back(update_v_p());
        cout<<"s_v"<<endl;
        sigma_v_h.push_back(update_sigma_v());
        
        cout<<"lat_v"<<endl;
        update_latent_v(i+1,latent_v);
        cout<<"lat_vs"<<endl;
        update_latent_v_s(i+1,latent_v_star);
        cout<<"lat_z"<<endl;
        update_latent_z(i+1,latent_z);
        cout<<"lat_w"<<endl;
        update_latent_w(i+1,latent_w);
        
        cout<<i+1<<endl;
    }
    
    cout<<"mu: "<<accept(mu_h)<<" mean: "<<final_mean(mu_h,30000)<<endl;
    cout<<"gam: "<<accept(gam_h)<<" mean: "<<final_mean(gam_h,30000)<<endl;
    cout<<"sigma_j: "<<accept(sigma_j_h)<<" mean: "<<final_mean(sigma_j_h,30000)<<endl;
    cout<<"rho: "<<accept(rho_h)<<" mean: "<<final_mean(rho_h,30000)<<endl;
    cout<<"k: "<<accept(k_h)<<" mean: "<<final_mean(k_h,30000)<<endl;
    cout<<"v_p: "<<accept(v_p_h)<<" mean: "<<final_mean(v_p_h,30000)<<endl;
    cout<<"sigma_v: "<<accept(sigma_v_h)<<" mean: "<<final_mean(sigma_v_h,30000)<<endl;
    cout<<"lat_w ";accept_latent(latent_w);cout<<latent_mean(latent_w)<<endl;
    cout<<"lat_z ";accept_latent(latent_z);cout<<latent_mean(latent_z)<<endl;
    cout<<"lat_vs ";accept_latent(latent_v_star);cout<<latent_mean(latent_v_star)<<endl;
    cout<<"lat_v ";accept_latent(latent_v);cout<<latent_mean(latent_v)<<endl;

	return(0);
}
