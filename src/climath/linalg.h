#ifndef LINALG_H_
#define LINALG_H_

double vvdot(double const *a, double const *b, int n);
void mvdot(double *r, double **m, double const *v, int n1, int n2);
int ludcmp(double **a, int n, int *indx);
void lubksb(double **a, int n, int *indx, double *b);
void leastsq(double **A, double *b, int n1, int n2);
void tridiag(int n, double *a, double *b, double *c, double *r, double *u, int pivot_type);

void band_decomp(int n, int m1, int m2, double *a, double *al, int *index, double *d);
void band_back_sub(int n, int m1, int m2, double *a, double *al, int *index, double *b);
void band_multiply(int n, int m1, int m2, double *a, double *x, double *b);
void band_improve(int n, int m1, int m2, double *aorig, double *a, double *al, int *index,
                  double *b, double *x);

#endif
