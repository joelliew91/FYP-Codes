#include "include.h"
#include <vector>
#include <iostream>
#include <math.h>
#include <boost/math/special_functions/bessel.hpp>
#include <complex>
#include <complex.h>
using namespace std;

double posterior(double z,double vt,double vu,double v_star,double yt,double yu,double w){
    //double post;
    double val1 = prior();
    double val2 = likelihood(z,vt,vu,v_star,yt,yu);
    double val3 = variance_gamma(z);
    double val4 = aux_g(v_star,w,vt,vu);
    return val1+val2+val3;
}

double prior(){
    double prior_mu,prior_gam,prior_sigj,prior_rho,prior_com;
    prior_mu = -pow(mu,2)/2;
    prior_gam = -gam;
    double rhob = (rho+1)/2;
    prior_rho = (1/7-1)*log(rhob);
    prior_sigj = -sigma_j;
    double I=0;
    double d = 4.0*k*v_p/pow(sigma_v,2);
    //cout<<k<<" "<<v_p<<" "<<1/pow(sigma_v,2)<<endl;
    //cout<<d<<endl;
    if(0.5*d<dmax)
        if(1<0.5*d)
            I=1;
    prior_com = log(I)-sigma_v-v_p;
    return(prior_mu+prior_gam+prior_rho+prior_sigj+prior_com);
}

double likelihood(double z,double vt,double vu,double v_star,double yt,double yu){
    double var = (1-pow(rho,2))*v_star;
    double v_bar = (vt-vu-k*v_p*delta+k*v_star)/sigma_v;
    double mn = yu + mu*delta + z + rho*v_bar;
    double val = -log(pow(var,0.5))-pow(yt-mn,2)/(2*var);
    return val;
}

double variance_gamma(double z){
    double alpha = delta*1.0;
    double sigma = sigma_j*sqrt(delta);
    double val1 = gam*pow(z,2)/sigma-log(alpha)*delta/alpha-log(tgamma(delta/alpha));
    double val2 = (delta/(2*alpha)-0.5)*log(pow(z,2)/(pow(gam,2)+2*pow(sigma,2)/alpha));
    double v = sqrt(pow(z,2)*(pow(gam,2)+2*pow(sigma,2)/alpha))/pow(sigma,2);
    double val3 = log(boost::math::cyl_bessel_k(delta/alpha-0.5,v));
    return val1+val2+val3;
}

double aux_g(double v,double w,double vt,double vu){
    int I = 0;
    if(v>min_v)
        if(v<max_v)
            if(w>min_w)
                if(w<max_w)
                    I=1;
    std::complex<double> val = phi(v,w,vt,vu);
    //cout<<val<<endl;
    return(I);
}

std::complex<double> phi(double v,double w,double vt,double vu){
    double de = 1.0*delta;
    double com = -2*pow(sigma_v,2)*w;
    complex<double> nu(pow(k,2),com);
    nu = sqrt(nu);
    double d = 0.5*4*k*v/pow(sigma_v,2)-1;
    complex<double> v1 = sqrt(vt*vu)*4.0*nu*exp(-0.5*nu*de)/(pow(sigma_v,2)*(1.0-exp(-nu*de)));
    complex<double> v2 = sqrt(vt*vu)*4.0*k*exp(-0.5*k*de)/(pow(sigma_v,2)*(1.0-exp(-k*de)));
    double re = real(v1);
    double im = imag(v1);
    //cdouble z = 1+1.5fi;
    //complex<double> val1 = boost::math::cyl_bessel_i(d,v1);///boost::math::cyl_bessel_i(d,v2);
    //cout<<val1<<endl;
    return(nu);
}





