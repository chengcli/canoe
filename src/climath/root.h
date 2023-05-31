#ifndef SRC_CLIMATH_ROOT_H_
#define SRC_CLIMATH_ROOT_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef double (*RootFunction_t)(double, void *);
int root(double x1, double x2, double xacc, double *x_root, RootFunction_t func,
         void *aux);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // SRC_CLIMATH_ROOT_H_
