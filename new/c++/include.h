#ifndef INIT
#define INIT

#include <vector>
#include <complex>

extern double mu;extern std::vector<double> mu_h;
extern double gam;extern std::vector<double> gam_h;
extern double sigma_j;extern std::vector<double> sigma_j_h;
extern double rho;extern std::vector<double> rho_h;
extern double k;extern std::vector<double> k_h;
extern double v_p;extern std::vector<double> v_p_h;
extern double sigma_v;extern std::vector<double> sigma_v_h;
extern double *v;
extern double *v_star;
extern double *w;
extern double *z;
extern double theta;extern std::vector<double> theta_h;
extern double max_v;
extern double max_w;
extern double min_v;
extern double min_w;
extern int iterations;
extern int delta;
extern std::vector<double> price;
extern double dmax;


void initialize();
double posterior(double z,double vt,double vu,double v_star,double yt,double yu,double w);
double prior();
double likelihood(double z,double vt,double vu,double v_star,double yt,double yu);
double variance_gamma(double z);
double aux_g(double v,double w,double vt,double vu);
std::complex<double> phi(double v,double w,double vt,double vu);



#endif