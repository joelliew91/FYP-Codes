#include <iostream>
#include <math.h>
#include "include.h"
#include <random>
#include <vector>
#include "./libAmosBessel/libAmosBessel.h"

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
int main(){

    //cout<<"iterations"<<endl;
    //cout<<"delta"<<endl;
    //cin>> iterations;
    //cin>> delta;
    
    initialize();
    for(int i=0;i<iterations;i++){
        
    }
    double val = posterior(z[1],v[1],v[0],v_star[1],price[1],price[0],w[1]);
    cout<<val<<endl;
	return(0);
}
