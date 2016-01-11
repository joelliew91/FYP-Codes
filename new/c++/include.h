#ifndef INIT
#define INIT

#include <vector>
#include <complex>
#include <random>

extern double mu;extern std::vector<double> mu_h;
extern double gam;extern std::vector<double> gam_h;
extern double sigma_j;extern std::vector<double> sigma_j_h;
extern double rho;extern std::vector<double> rho_h;
extern double k;extern std::vector<double> k_h;
extern double v_p;extern std::vector<double> v_p_h;
extern double sigma_v;extern std::vector<double> sigma_v_h;
extern double theta;extern std::vector<double> theta_h;

extern double *v;
extern double *v_star;
extern double *w;
extern double *z;

extern double max_v;
extern double max_w;
extern double min_v;
extern double min_w;

extern int iterations;
extern int delta;
extern std::vector<double> price;
extern double dmax;
extern int total;
extern int flag;

void initialize();
double post_no_prior(double z,double vt,double vu,double v_star,double yt,double yu,double w);
double posterior();
double prior();
double likelihood(double z,double vt,double vu,double v_star,double yt,double yu);
double variance_gamma(double z);
double aux_g(double v,double w,double vt,double vu);
std::complex<double> phi(double v,double w,double vt,double vu);
double transition(double vt,double vu);

double sd(std::vector<double> history);
double mean(std::vector<double> history);
double normal(double s1,double s2);
double uniform(unsigned s);
double accept(std::vector<double> hist);

double update_gam();
double update_mu();
double update_sigma_j();
double update_rho();
double update_k();
double update_v_p();
double update_sigma_v();

void update_latent_w(int iter,double** hist);
void update_latent_z(int iter,double** hist);
void update_latent_v_s(int iter,double** hist);
void update_latent_v(int iter,double** hist);

double mean(int iter,double* history);
double sd(int iter,double* history);
void accept_latent(double** hist);
double scaler(std::vector<double> history);

double final_mean(std::vector<double> history,int cutoff);
double accept_latent(int iter,double* hist);
double latent_mean(double ** hist);
#endif