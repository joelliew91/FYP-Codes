#include "include.h"
#include <vector>
#include <iostream>
#include <math.h>
#include <boost/math/special_functions/bessel.hpp>
#include <complex>
#include <complex_bessel.h>
#include <boost/math/distributions/non_central_chi_squared.hpp>

using namespace std;

double post_no_prior(double z,double vt,double vu,double v_star,double yt,double yu,double w){
    double val2 = likelihood(z,vt,vu,v_star,yt,yu);
    double val3 = variance_gamma(z);
    double val4 = aux_g(v_star,w,vt,vu);
    double val5 = transition(vt,vu);
    return val2+val3+val4+val5;
}
double posterior(){
    flag = 1;
    double val1 = prior();
    double sum= 0 ;
    for(int i =1;i<total;i++){
        sum += post_no_prior(z[i],v[i],v[i-1],v_star[i],price[i],price[i-1],w[i]);
        if(!flag){
            break;
        }
    }
    return sum+val1;
}

double prior(){
    double prior_mu,prior_gam,prior_sigj,prior_rho,prior_com;
    prior_mu = -pow(mu,2)/2;
    prior_gam = -gam;
    double rhob = (rho+1)/2;
    prior_rho = (1/7-1)*log(rhob);
    prior_sigj = -sigma_j;
    double d = 4.0*k*v_p/pow(sigma_v,2);
    if(0.5*d>dmax || 0.5*d<1)
        flag = 0;
    prior_com = -sigma_v-v_p;
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
    double va = sqrt(pow(z,2)*(pow(gam,2)+2*pow(sigma,2)/alpha))/pow(sigma,2);
    complex<double> a(va,0);
    complex<double> temp;
    temp = sp_bessel::besselK(delta/alpha-0.5,a);
    if(norm(temp)==0||isnan(va)||isinf(va))
        flag=0;
    double val3 = log(temp.real());
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
    if(norm(val) == 0)
        flag=0;
    double g = cos(v*w)*val.real()+sin(v*w)*val.imag();
    double aux = abs(g)*I+abs(g)*exp(-v-w)*(1-I);
    
    return(log(aux));
}

std::complex<double> phi(double v,double w,double vt,double vu){
    double de = 1.0*delta;
    double com = -2*pow(sigma_v,2)*w;
    complex<double> nu(pow(k,2),com);
    nu = sqrt(nu);
    double d = 0.5*4*k*v/pow(sigma_v,2)-1;
    complex<double> v1 = sqrt(vt*vu)*4.0*nu*exp(-0.5*nu*de)/(pow(sigma_v,2)*(1.0-exp(-nu*de)));
    complex<double> v2 = sqrt(vt*vu)*4.0*k*exp(-0.5*k*de)/(pow(sigma_v,2)*(1.0-exp(-k*de)));
    complex<double> den = sp_bessel::besselI(d,v2);
    complex<double> val1(0,0);
    if(norm(den)==0||isinf(d)||isnan(d)||isinf(v1.real())||isinf(v1.imag())||isinf(v2.real())||isinf(v2.imag())){
        flag = 0;
        return val1;
    }
    else{
        val1 = sp_bessel::besselI(d,v1)/den;
        val1 = val1*nu*exp(-0.5*(nu-k)*de)*(1.0-exp(-k*de))/(k*(1.0-exp(-nu*de)));
        val1 = val1*exp((vu+vt)/pow(sigma_v,2)*(k*(1.0+exp(-k*de))/(1.0-exp(-k*de))-nu*(1.0+exp(-nu*de))/(1.0-exp(-nu*de))));
    }
    return(val1);
}

double transition(double vt,double vu){
    double ncp = 4.0*k*exp(-k*delta)*vu/(pow(sigma_v,2)*(1-exp(-k*delta)));
    double  d = 4*k*k/pow(sigma_v,2);

    double val;
    
    if(isinf(vt)||isinf(vu)||isnan(vt)||isnan(ncp)){
        val = 0;
        flag = 0;
    }
    else{
        boost::math::non_central_chi_squared_distribution<long double> myNonCentralChiSquared(d, ncp);
        val = pdf(myNonCentralChiSquared,vt);
    }
    return log(val*pow(sigma_v,2)*(1-exp(-k*delta))/(4*k));

}




