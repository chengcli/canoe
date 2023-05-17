// C/C++ headers
#include <iostream>
#include <cmath>

#include <climath/linalg.h>
#include <Eigen/Dense>
#include <utils/ndarrays.hpp>
#include "gaussian_process.hpp"

void gp_covariance(KernelFunction_t kernel, double **cov,
  double const *x, double const *s, int n, double l)
{
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      cov[i][j] = kernel(x[i], x[j], l, s[i]*s[j]);
    }
}

void gp_covariance2(KernelFunction_t kernel, double **cov,
  double const *x1, double const *s1, int n1, 
  double const *x2, double const *s2, int n2, double l)
{
  for (int i = 0; i < n1; ++i)
    for (int j = 0; j < n2; ++j) {
      cov[i][j] = kernel(x1[i], x2[j], l, s1[i]*s2[j]);
    }
}

double gp_predict(KernelFunction_t kernel, double *arr2,
  double const *x2, double const *s2, int n2, double const *arr1,
  double const *x1, double const *s1, int n1, double len)
{
  double **cov1, **cov2;
  NewCArray(cov1, n1, n1);
  NewCArray(cov2, n2, n1);

  gp_covariance(kernel, cov1, x1, s1, n1, len);
  gp_covariance2(kernel, cov2, x2, s2, n2, x1, s1, n1, len);

  // copy arr1 to b because the solution x=A^{-1}*b will be stored 
  // in b after calling lubksb
  int *indx = new int [n1];
  double *b = new double [n1];
  memcpy(b, arr1, n1*sizeof(double));

  ludcmp(cov1, n1, indx);
  lubksb(cov1, n1, indx, b);
  mvdot(arr2, cov2, b, n2, n1);

  double d = -0.5*vvdot(arr1, b, n1);

  delete[] b;
  delete[] indx;
  FreeCArray(cov1);
  FreeCArray(cov2);

  // returns prior probability
  return d;
}

double gp_lnprior(KernelFunction_t kernel, double const *arr1, 
  double const *x1, double const *s1, int n1, double len)
{
  double **cov1;
  NewCArray(cov1, n1, n1);

  gp_covariance(kernel, cov1, x1, s1, n1, len);

  int *indx = new int [n1];
  double d, *b = new double [n1];
  memcpy(b, arr1, n1*sizeof(double));

  ludcmp(cov1, n1, indx);
  lubksb(cov1, n1, indx, b);

  d = -0.5*vvdot(arr1, b, n1);

  delete[] b;
  delete[] indx;
  FreeCArray(cov1);

  // returns prior probability
  return d;
}
