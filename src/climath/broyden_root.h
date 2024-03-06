#ifndef SRC_CLIMATH_BROYDEN_ROOT_H_
#define SRC_CLIMATH_BROYDEN_ROOT_H_

#ifdef __cplusplus
extern "C" {
#endif

int broyden_root(int n, double *x,
                 void (*vecfunc)(int, double *, double *, void *), double tol_f,
                 int max_it, void *arg);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // SRC_CLIMATH_BROYDEN_ROOT_H_
