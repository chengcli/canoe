#ifndef ROOT_H_
#define ROOT_H_

typedef double (*RootFunction_t)(double, void*);
int root(double x1, double x2, double xacc, double *x_root, RootFunction_t func, void *aux);

#endif
