#ifndef SRC_CLIMATH_BROYDEN_ROOT_UTILS_H_
#define SRC_CLIMATH_BROYDEN_ROOT_UTILS_H_

#ifdef __cplusplus
extern "C" {
#endif

/*
 * DS stands for Dennis and Schnabel (1996).
 */
#define DS_EPS 1e-8
#define DS_SQRT_EPS sqrt(DS_EPS)
#define DS_TOL_X pow(DS_EPS, 2. / 3.)
#define DS_TOL_F pow(DS_EPS, 1. / 3.)
#define DS_TOL_MIN pow(DS_EPS, 2. / 3.)
#define DS_ALPHA 1.e-4
#define DS_MAX_STEP 100.

#define DS_LINE_STEP 0
#define DS_HOOK_STEP 1 /* Not yet implemented. */
#define DS_DOGLEG_STEP 2
/*
 * Values for status (retcode in DS96):
 */
#define DS_X_ACCEPTED 0
#define DS_X_NO_PROGRESS 1
#define DS_REDUCE_DELTA 2
#define DS_INCREASE_DELTA 3
#define DS_MAX_IT_EXCEEDED 4
#define DS_MAX_TAKEN_5 5
#define DS_INITIAL 6
#define DS_SINGULAR_JACOBIAN 7

double *fvector(int nl, int nh, char *calling_func);

void free_fvector(double *m, int nl, int nh, char *calling_func);

int global_step(int n, double *x_old, double f_old, double *g, double *r,
                double *sn, double max_step, double *delta, int step_type,
                int *status, double *x, double *f, double *fvec, void *arg,
                void (*vecfunc)(int, double *, double *, void *));

int line_search(int n, double *x_old, double f_old, double *g, double *sn,
                double max_step, int *status, double *x, double *f,
                double *fvec, void *arg,
                void (*vecfunc)(int, double *, double *, void *));

int dogleg_driver(int n, double *x_old, double f_old, double *g, double *r,
                  double *sn, double max_step, double *delta, int *status,
                  double *x, double *f, double *fvec, void *arg,
                  void (*vecfunc)(int, double *, double *, void *));

int dogleg_step(int n, double *g, double *r, double *sn, double newt_length,
                double max_step, double *delta, int *first_dog, double *s_hat,
                double *nu_hat, double *s);

int trust_region(int n, double *x_old, double f_old, double *g, double *s,
                 int newt_taken, double max_step, int step_type, double *r,
                 double *delta, int *status, double *x_prev, double *f_prev,
                 double *x, double *f, double *fvec, void *arg,
                 void (*vecfunc)(int, double *, double *, void *));

int qr_decompose(int n, double *r, double *c, double *d);

void qr_update(int n, double *r, double *qt, double *u, double *v);

void qr_rotate(int n, double *r, double *qt, int i, double a, double b);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // SRC_CLIMATH_BROYDEN_ROOT_UTILS_H_
