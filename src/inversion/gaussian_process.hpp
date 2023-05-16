#ifndef GAUSSIAN_PROCESS_HPP
#define GAUSSIAN_PROCESS_HPP

typedef double (*KernelFunction_t)(double, double, double, double);

inline double SquaredExponential(double x1, double x2, double l, double s1s2 = 1.)
{
  return s1s2*exp(-(x1-x2)*(x1-x2)/(2.*l*l));
}

inline double OrnsteinUhlenbeck(double x1, double x2, double l, double s1s2 = 1.)
{
  return s1s2*exp(-std::fabs(x1-x2)/l);
}

void gp_covariance(KernelFunction_t kernel, double **cov,
  double const *x, double const *s, int n, double l);

void gp_covariance2(KernelFunction_t kernel, double **cov,
  double const *x1, double const *s1, int n1, 
  double const *x2, double const *s2, int n2, double l);

double gp_predict(KernelFunction_t kernel, double *arr2,
  double const *x2, double const *s2, int n2, double const *arr1,
  double const *x1, double const *s1, int n1, double len);

double gp_lnprior(KernelFunction_t kernel, double const *arr1, 
  double const *x1, double const *s1, int n1, double len);

#endif
