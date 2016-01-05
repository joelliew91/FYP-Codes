#include "include.h"
#include <iostream>
#include <random>
#include <vector>
#include <fstream>
#include <string>


using namespace std;

void initialize(){
    max_v = 10000;
    min_v = 0;
    min_w = 0;
    max_w = 10000;
    
    default_random_engine generator;
    normal_distribution<double> distribution(0,10.0);
    mu = distribution(generator);
    gam = distribution(generator);
    gamma_distribution<double> distribution1(1.0,10.0);
    sigma_j = distribution1(generator);
    k = 5;
    theta = abs(distribution(generator));
    v_p = 1/distribution1(generator);
    rho = 0.5;
    sigma_v = 1.12;
    
    
    ifstream myfile("data.csv");
    int count = 0;
    string line;
    while(myfile.good()){
        getline(myfile,line);
        int pos = line.find(",");
        if(pos==3){
            price.push_back(stod(line.substr(pos+1)));
            //cout<<line<<endl;
        }
        
    }
    myfile.close();
    v = new double[delta*price.size()];
    v_star = new double[delta*price.size()];
    w = new double[delta*price.size()];
    z = new double[delta*price.size()];
    for(int i =0;i<price.size();i++){
        v[i] = 1;
        v_star[i] = 1;
        w[i] = 1;
        z[i] = 1;
    }
    mu_h.push_back(mu);
    sigma_j_h.push_back(sigma_j);
    theta_h.push_back(theta);
    k_h.push_back(k);
    v_p_h.push_back(v_p);
    sigma_v_h.push_back(sigma_v);
    rho_h.push_back(rho);
    gam_h.push_back(gam);
    return ;
}

