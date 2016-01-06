#include "include.h"
#include <vector>
#include <iostream>
#include <math.h>
#include <boost/math/special_functions/bessel.hpp>
#include <complex>
#include <complex_bessel.h>
#include <boost/math/distributions/non_central_chi_squared.hpp>

using namespace std;

double posterior(double z,double vt,double vu,double v_star,double yt,double yu,double w){
    //double post;
    double val1 = prior();
    double val2 = likelihood(z,vt,vu,v_star,yt,yu);
    double val3 = variance_gamma(z);
    double val4 = aux_g(v_star,w,vt,vu);
    double val5 = transition(vt,vu);
    return val1+val2+val3+val4+val5;
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
    complex<double> a(v,0);
    double val3 = log(sp_bessel::besselK(delta/alpha-0.5,a).real());
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
    double g = cos(v*w)*val.real()+sin(v*w)*val.imag();
    double aux = abs(g)*I+abs(g)*exp(-v-w)*(1-I);
    return(aux);
}

std::complex<double> phi(double v,double w,double vt,double vu){
    double de = 1.0*delta;
    double com = -2*pow(sigma_v,2)*w;
    complex<double> nu(pow(k,2),com);
    nu = sqrt(nu);
    double d = 0.5*4*k*v/pow(sigma_v,2)-1;
    complex<double> v1 = sqrt(vt*vu)*4.0*nu*exp(-0.5*nu*de)/(pow(sigma_v,2)*(1.0-exp(-nu*de)));
    complex<double> v2 = sqrt(vt*vu)*4.0*k*exp(-0.5*k*de)/(pow(sigma_v,2)*(1.0-exp(-k*de)));
    
    //cout<<d<<endl;
    complex<double> val1 = sp_bessel::besselI(d,v1)/sp_bessel::besselI(d,v2);
    val1 = val1*nu*exp(-0.5*(nu-k)*de)*(1.0-exp(-k*de))/(k*(1.0-exp(-nu*de)));
    val1 = val1*exp((vu+vt)/pow(sigma_v,2)*(k*(1.0+exp(-k*de))/(1.0-exp(-k*de))-nu*(1.0+exp(-nu*de))/(1.0-exp(-nu*de))));
    //cout<<val1<<endl;
    return(val1);
}

double transition(double vt,double vu){
    double ncp = 4.0*k*exp(-k*delta)*vu/(pow(sigma_v,2)*(1-exp(-k*delta)));
    double  d = 4*k*k/pow(sigma_v,2);
    //complex<double> temp(sqrt(ncp*vt),0);
    //double b = sp_bessel::besselI(d/2-1.0,temp).real();
    //double nonchisq = 0.5*exp(-(vt+ncp)/2)*pow(vt/ncp,d/4-0.5)*b;
    //cout<<vt<<endl;
    //cout<<"vu: "<<vu<<" ncp: "<<ncp<<endl;
    double val;
    
    if(isinf(vt)||isinf(vu))
        val = 0;
    else{
        boost::math::non_central_chi_squared_distribution<long double> myNonCentralChiSquared(d, ncp);
        val = pdf(myNonCentralChiSquared,vt);
    }
    return log(val*pow(sigma_v,2)*(1-exp(-k*delta))/(4*k));
    //return log(nonchisq*pow(sigma_v,2)*(1-exp(-k*delta))/(4*k));
}




