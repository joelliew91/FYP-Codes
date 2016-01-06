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
int iterations=100;
int delta=1;
vector<double> price;
int total;



int main(){
    srand(time(0));

    initialize();
    
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
    
    cout<<"mu: "<<accept(mu_h)<<endl;
    cout<<"gam: "<<accept(gam_h)<<endl;
    cout<<"sigma_j: "<<accept(sigma_j_h)<<endl;
    cout<<"rho: "<<accept(rho_h)<<endl;
    cout<<"k: "<<accept(k_h)<<endl;
    cout<<"v_p: "<<accept(v_p_h)<<endl;
    cout<<"sigma_v: "<<accept(sigma_v_h)<<endl;
    accept_latent(latent_w);
    accept_latent(latent_z);
    accept_latent(latent_v_star);
    accept_latent(latent_v);

	return(0);
}
