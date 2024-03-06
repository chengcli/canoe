// C
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// climath
#include "broyden_root_utils.h"

/*
 * Shift macros:
 */
#define R(i, j) r[j + (i) * n]
#define QT(i, j) qt[j + (i) * n]

#define TRUE 1
#define FALSE 0

#define MIN(x, y)                  \
  ({                               \
    const double _x = (double)(x); \
    const double _y = (double)(y); \
    _x < _y ? _x : _y;             \
  })

#define MAX(x, y)                  \
  ({                               \
    const double _x = (double)(x); \
    const double _y = (double)(y); \
    _x > _y ? _x : _y;             \
  })

#define NR_SIGN(a, b) ((b) > 0. ? fabs(a) : -fabs(a))

/*======================= fvector()
 * ============================================*/

/*
 *  Allocates memory for a 1D double array
 *  with range [nl..nh].
 */

#undef DEBUG

double *fvector(int nl, int nh, char *calling_func) {
  unsigned int len_safe;
  int nl_safe, nh_safe;
  double *m;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int idbms = 0;
  static char dbmsname[] = "fvector";

#if defined(DEBUG)
  fprintf(stderr, "fvector() called by %s \n", calling_func);
  fflush(stderr);
#endif

  if (nh < nl) {
    fprintf(stderr, "**error:%s, called by %s, range (%d,%d)\n", dbmsname,
            calling_func, nl, nh);
    exit(1);
  }

  nl_safe = (nl < 0) ? nl : 0;
  nh_safe = (nh > 0) ? nh : 0;
  len_safe = (unsigned)(nh_safe - nl_safe + 1);

  m = (double *)calloc(len_safe, sizeof(double));

  if (!m) {
    fprintf(stderr, "**error:%s, called by %s, nl=%d,nh=%d,len_safe=%d\n",
            dbmsname, calling_func, nl, nh, len_safe);
    exit(1);
  }
  m -= nl_safe;

  return m;
}

/*======================= end of fvector()
 * ====================================*/

/*======================= free_fvector()
 * ======================================*/

#undef DEBUG

/*
 *  Frees memory allocated by fvector().
 */

void free_fvector(double *m, int nl, int nh, char *calling_func) {
  int nl_safe;

  nl_safe = (nl < 0) ? nl : 0;
  m += nl_safe;
  free(m);

#if defined(DEBUG)
  fprintf(stderr, "free_fvector() called by %s \n", calling_func);
#endif

  return;
}

/*======================= global_step() ======================================*/

/*
 * Take a globally convergent step towards a vector root.
 * Returns max_taken.
 *
 * The function line_search() backtracks along the quasi-Newton direction.
 *
 * The function dogleg_driver() uses the model trust-region approach, where
 * delta is the radius of the trust region. A step is taken
 * in the steepest descent direction for small delta, in the quasi_Newton
 * direction for large delta, and on the connecting line segment for
 * intermediate delta.
 * See Dennis and Schnabel (1996, DS96).
 *
 * The value of step_type can be DS_LINE_STEP, DS_HOOK_STEP, or DS_DOGLEG_STEP.
 *
 * NOTE: Unlike in DS96, here R is R of QR, not the transpose of R.
 */

int global_step(int n, double *x_old, double f_old, double *g, double *r,
                double *sn, double max_step, double *delta, int step_type,
                int *status, double *x, double *f, double *fvec, void *arg,
                void (*vecfunc)(int, double *, double *, void *)) {
  int max_taken = FALSE;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int idbms = 0;
  static char dbmsname[] = "global_step";

  if (step_type == DS_LINE_STEP) {
    max_taken = line_search(n, x_old, f_old, g, sn, max_step, status, x, f,
                            fvec, arg, vecfunc);
  } else if (step_type == DS_HOOK_STEP) {
    fprintf(stderr, "**error:%s, DS_HOOK_STEP not yet implemented\n", dbmsname);
    exit(1);
  } else if (step_type == DS_DOGLEG_STEP) {
    max_taken = dogleg_driver(n, x_old, f_old, g, r, sn, max_step, delta,
                              status, x, f, fvec, arg, vecfunc);
  } else {
    fprintf(stderr, "**error:%s, unrecognized step_type=%d\n", dbmsname,
            step_type);
    exit(1);
  }

  return max_taken;
}

/*======================= end of global_step() ===============================*/

/*======================= line_search() ======================================*/

/*
 * Adapted from Numerical Recipes in C, 2nd ed., p. 385-386, which is
 * an implementation of Algorithm A6.3.1 of Dennis and Schnabel (1996).
 * Returns max_taken.
 * Assumes zero-based indexing.
 */

int line_search(int n, double *x_old, double f_old, double *g, double *sn,
                double max_step, int *status, double *x, double *f,
                double *fvec, void *arg,
                void (*vecfunc)(int, double *, double *, void *)) {
  int i, max_taken = FALSE;
  double a, b, lambda, lambda_prev, lambda_min, disc, f_prev, rhs1, rhs2,
      initial_slope, newt_length, rel_step_length, tmp, tmp_lambda;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int idbms = 0;
  static char dbmsname[] = "line_search";

  *status = DS_X_ACCEPTED;

  tmp = 0.;
  for (i = 0; i < n; i++) {
    tmp += sn[i] * sn[i];
  }
  newt_length = sqrt(tmp);

  if (newt_length > max_step) {
    tmp = max_step / newt_length;
    for (i = 0; i < n; i++) {
      sn[i] *= tmp;
    }
  }

  initial_slope = 0.;
  for (i = 0; i < n; i++) {
    initial_slope += g[i] * sn[i];
  }

  rel_step_length = 0.;
  for (i = 0; i < n; i++) {
    rel_step_length =
        MAX(rel_step_length, fabs(sn[i]) / MAX(fabs(x_old[i]), 1.));
  }

  lambda_min = DS_TOL_X / rel_step_length;
  lambda = 1.;

  /*
   * Iteration loop.
   */
  while (TRUE) {
    for (i = 0; i < n; i++) {
      x[i] = x_old[i] + lambda * sn[i];
    }

    /*
     * Calculate fvec[].
     */
    (*vecfunc)(n, x, fvec, arg);

    *f = 0.;
    for (i = 0; i < n; i++) {
      *f += fvec[i] * fvec[i];
    }
    *f *= .5;

    if (lambda < lambda_min) {
      /*
       * Convergence on dx.
       * For zero finding, calling program should verify the convergence.
       */
      for (i = 0; i < n; i++) {
        x[i] = x_old[i];
      }
      *status = DS_X_NO_PROGRESS;
      return max_taken;
    } else if (*f <= f_old + DS_ALPHA * lambda * initial_slope) {
      /* Sufficient function decrease. */
      if (lambda == 1. && newt_length > 0.99 * max_step) {
        max_taken = TRUE;
      }
      return max_taken;
    } else {
      if (lambda == 1.) {
        /* First backtrack uses a quadratic model. */
        tmp_lambda = -initial_slope / (2. * (*f - f_old - initial_slope));
      } else {
        /* Subsequent backtracks use a cubic model. */
        rhs1 = *f - initial_slope * lambda - f_old;
        rhs2 = f_prev - initial_slope * lambda_prev - f_old;
        a = (rhs1 / (lambda * lambda) - rhs2 / (lambda_prev * lambda_prev)) /
            (lambda - lambda_prev);
        b = (-lambda_prev * rhs1 / (lambda * lambda) +
             lambda * rhs2 / (lambda_prev * lambda_prev)) /
            (lambda - lambda_prev);
        if (a == 0.) {
          tmp_lambda = -initial_slope / (2. * b);
        } else {
          disc = b * b - 3. * a * initial_slope;
          if (disc < 0.) {
            fprintf(stderr, "**error:%s,roundoff problem\n", dbmsname);
            exit(1);
          } else {
            tmp_lambda = (-b + sqrt(disc)) / (3. * a);
          }
        }
        if (tmp_lambda > .5 * lambda) {
          tmp_lambda = .5 * lambda;
        }
      }
    }
    lambda_prev = lambda;
    f_prev = *f;
    lambda = MAX(tmp_lambda, .1 * lambda);
  }

  /* Never get here. */
}

/*======================= end of line_search() ===============================*/

/*======================= dogleg_driver() ====================================*/

/*
 * Adapted from Dennis and Schnabel (1996), Appendix A, Algorithm A6.4.3.
 * Returns max_taken;
 * Assumes zero-based indexing.
 */
int dogleg_driver(int n, double *x_old, double f_old, double *g, double *r,
                  double *sn, double max_step, double *delta, int *status,
                  double *x, double *f, double *fvec, void *arg,
                  void (*vecfunc)(int, double *, double *, void *)) {
  int i, max_taken, newt_taken, first_dog = TRUE;
  double newt_length, f_prev, tmp, *s, *s_hat, *nu_hat, *x_prev;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int idbms = 0;
  static char dbmsname[] = "dogleg_driver";

  /*
   * Allocate memory.
   */
  s = fvector(0, n - 1, dbmsname);
  s_hat = fvector(0, n - 1, dbmsname);
  nu_hat = fvector(0, n - 1, dbmsname);
  x_prev = fvector(0, n - 1, dbmsname);

  *status = DS_INITIAL;

  tmp = 0.;
  for (i = 0; i < n; i++) {
    tmp += sn[i] * sn[i];
  }
  newt_length = sqrt(tmp);

  while (*status >= DS_REDUCE_DELTA) {
    /*
     * Find new step.
     */
    newt_taken = dogleg_step(n, g, r, sn, newt_length, max_step, delta,
                             &first_dog, s_hat, nu_hat, s);
    /*
     * Check new point and update trust region.
     */
    max_taken = trust_region(n, x_old, f_old, g, s, newt_taken, max_step,
                             DS_DOGLEG_STEP, r, delta, status, x_prev, &f_prev,
                             x, f, fvec, arg, vecfunc);
  }

  /*
   * Free allocated memory.
   */
  free_fvector(x_prev, 0, n - 1, dbmsname);
  free_fvector(nu_hat, 0, n - 1, dbmsname);
  free_fvector(s_hat, 0, n - 1, dbmsname);
  free_fvector(s, 0, n - 1, dbmsname);

  return max_taken;
}

/*======================= end of dogleg_driver() =============================*/

/*======================= dogleg_step() ======================================*/

/*
 * Adapted from Dennis and Schnabel (1996), Appendix A, Algorithm A6.4.4.
 * Returns newt_taken.
 * Assumes zero-based indexing.
 */
int dogleg_step(int n, double *g, double *r, double *sn, double newt_length,
                double max_step, double *delta, int *first_dog, double *s_hat,
                double *nu_hat, double *s) {
  int i, j, newt_taken;
  static int eta_warned = FALSE;
  double alpha, beta, al_be, eta, lambda, tmp, tmp_nu, tmp_cauchy;
  static double cauchy_length;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int idbms = 0;
  static char dbmsname[] = "dogleg_step";

  if (newt_length <= *delta) {
    /*
     * s is Newton step.
     */
    newt_taken = TRUE;
    for (i = 0; i < n; i++) {
      s[i] = sn[i];
    }
    *delta = newt_length;
  } else {
    /*
     * Newton step is too long, find s on double-dogleg curve.
     */
    newt_taken = FALSE;
    if (*first_dog == TRUE) {
      /*
       * Calculate double-dogleg curve.
       */
      *first_dog = FALSE;

      alpha = 0.;
      for (i = 0; i < n; i++) {
        alpha += g[i] * g[i];
      }

      beta = 0.;
      for (i = 0; i < n; i++) {
        tmp = 0.;
        for (j = i; j < n; j++) {
          /*
           * NOTE: Unlike DS96, here R is R of QR, not the transpose of R.
           */
          tmp += R(i, j) * g[j];
        }
        beta += tmp * tmp;
      }
      al_be = alpha / beta;
      cauchy_length = al_be * sqrt(alpha);

      tmp = 0.;
      for (i = 0; i < n; i++) {
        s_hat[i] = -g[i] * al_be;
        tmp += g[i] * sn[i];
      }
      eta = 0.2 + 0.8 * alpha * al_be / fabs(tmp);

      /*
       * Check range of eta, which should be [0,1]:
       */
      if (eta > 1.) {
        /*
         * If eta > 1.0 print a one-time warning.
         * This behavior is known to be caused by optimization (-O)
         * for gcc (at least versions 2.7.2.3.f.1 and 2.96).
         */
        if (!eta_warned) {
          fprintf(
              stderr,
              "eta=%g > 1, setting eta=1. This warning will not be repeated.\n"
              "Check for optimization error (-O vs -g).",
              eta);
          eta_warned = TRUE;
        }
        eta = 1.;
      }

      for (i = 0; i < n; i++) {
        nu_hat[i] = eta * sn[i] - s_hat[i];
      }

      if (*delta == -1) {
        /*
         * First iteration, and no initial trust region was
         * provided by the user.
         */
        *delta = MIN(cauchy_length, max_step);
      }
    }

    if (eta * newt_length <= *delta) {
      /*
       * Take partial step in Newton direction.
       */
      for (i = 0; i < n; i++) {
        s[i] = ((*delta) / newt_length) * sn[i];
      }
    } else if ((cauchy_length) >= (*delta)) {
      /*
       * Take step in steepest descent direction.
       */
      for (i = 0; i < n; i++) {
        s[i] = ((*delta) / (cauchy_length)) * s_hat[i];
        /*
         * Screen for NaN.
         */
        if (!isfinite(s[i])) {
          fprintf(stderr,
                  "**error:%s, s[%d]=%g,*delta=%g,cauchy_length=%g,s_hat=%g",
                  dbmsname, i, s[i], *delta, cauchy_length, s_hat[i]);
          exit(1);
        }
      }
    } else {
      /*
       * Take convex-combination step.
       */
      tmp = 0.;
      tmp_nu = 0.;
      for (i = 0; i < n; i++) {
        tmp += nu_hat[i] * s_hat[i];
        tmp_nu += nu_hat[i] * nu_hat[i];
      }
      tmp_cauchy = cauchy_length * cauchy_length - (*delta) * (*delta);
      lambda = (-tmp + sqrt(tmp * tmp - tmp_nu * tmp_cauchy)) / tmp_nu;
      for (i = 0; i < n; i++) {
        s[i] = s_hat[i] + lambda * nu_hat[i];
        /*
         * Screen for NaN.
         */
        if (!isfinite(s[i])) {
          fprintf(stderr, "**error:%s, s[%d]=%g,lambda=%g,nu_hat=%g,s_hat=%g",
                  dbmsname, i, s[i], lambda, nu_hat[i], s_hat[i]);
          exit(1);
        }
      }
    }
  }

  return newt_taken;
}

/*======================= end of dogleg_step() ===============================*/

/*======================= trust_region() =====================================*/

/*
 * Adapted from Dennis and Schnabel (1996), Appendix A, Algorithm A6.4.5.
 * Returns max_taken.
 * Assumes zero-based indexing.
 */

int trust_region(int n, double *x_old, double f_old, double *g, double *s,
                 int newt_taken, double max_step, int step_type, double *r,
                 double *delta, int *status, double *x_prev, double *f_prev,
                 double *x, double *f, double *fvec, void *arg,
                 void (*vecfunc)(int, double *, double *, void *)) {
  int i, j, max_taken = FALSE;
  double initial_slope, step_length, rel_step_length, delta_tmp, tmp;
  double df, df_tol, df_pred;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int idbms = 0;
  static char dbmsname[] = "trust_region";

  /*
   * Screen for NaN.
   */
  for (i = 0; i < n; i++) {
    if (!isfinite(s[i])) {
      fprintf(stderr, "**error:%s,s[%d]=%g\n", dbmsname, i, s[i]);
      exit(1);
    }
  }

  tmp = 0.;
  for (i = 0; i < n; i++) {
    tmp += s[i] * s[i];
  }
  step_length = sqrt(tmp);

  /* Take step. */
  for (i = 0; i < n; i++) {
    x[i] = x_old[i] + s[i];
  }

  /* Calculate fvec[]. */
  (*vecfunc)(n, x, fvec, arg);

  /* Compute f. */
  *f = 0.;
  for (i = 0; i < n; i++) {
    *f += fvec[i] * fvec[i];
  }
  *f *= .5;

  df = *f - f_old;

  initial_slope = 0.;
  for (i = 0; i < n; i++) {
    initial_slope += g[i] * s[i];
  }

  if (*status != DS_INCREASE_DELTA) {
    *f_prev = 0.;
  }

  df_tol = DS_ALPHA * initial_slope;

  if (*status == DS_INCREASE_DELTA && (*f >= *f_prev || df > df_tol)) {
    /*
     * Retreat.
     */
    *status = DS_X_ACCEPTED;
    for (i = 0; i < n; i++) {
      x[i] = x_prev[i];
    }
    *f = *f_prev;
    *delta *= .5;
  } else if (df >= df_tol) {
    /*
     * The value of f is too large.
     */
    rel_step_length = 0.;
    for (i = 0; i < n; i++) {
      rel_step_length = MAX(rel_step_length, fabs(s[i]) / MAX(fabs(x[i]), 1.));
    }

    if (rel_step_length < DS_TOL_X) {
      /*
       * The step is too small.
       */
      *status = DS_X_NO_PROGRESS;
      for (i = 0; i < n; i++) {
        x[i] = x_old[i];
      }
    } else {
      /*
       * Reduce delta.
       */
      *status = DS_REDUCE_DELTA;
      delta_tmp = (-initial_slope * step_length) / (2. * (df - initial_slope));
      if (delta_tmp < (*delta) * 0.1) {
        *delta *= 0.1;
      } else if (delta_tmp > (*delta) * 0.5) {
        *delta *= 0.5;
      } else {
        *delta = delta_tmp;
      }
    }
  } else {
    /*
     * The value of f is sufficiently small.
     */
    df_pred = initial_slope;

    if (step_type == DS_HOOK_STEP) {
      fprintf(stderr, "**error:%s, DS_HOOK_STEP not yet implemented\n",
              dbmsname);
      exit(1);
    } else if (step_type == DS_DOGLEG_STEP) {
      for (i = 0; i < n; i++) {
        tmp = 0.;
        for (j = i; j < n; j++) {
          /*
           * NOTE: Unlike DS96, here R is R of QR, not the transpose of R.
           */
          tmp += R(i, j) * s[j];
        }
        df_pred += tmp * tmp * .5;
      }
    } else {
      fprintf(stderr, "**error:%s, unrecognized step_type=%d\n", dbmsname,
              step_type);
      exit(1);
    }

    if (((*status) != DS_REDUCE_DELTA &&
         fabs(df_pred - df) <= 0.1 * fabs(df)) ||
        (df <= initial_slope && newt_taken == FALSE &&
         (*delta) <= 0.99 * max_step)) {
      /*
       * Double delta.
       */
      *status = DS_INCREASE_DELTA;
      for (i = 0; i < n; i++) {
        x_prev[i] = x[i];
      }
      *f_prev = *f;
      *delta = MIN((*delta) * 2., max_step);
    } else {
      /*
       * Accept x, choose delta for next iteration.
       */
      *status = DS_X_ACCEPTED;
      if (step_length > 0.99 * max_step) {
        max_taken = TRUE;
      }
      if (df >= 0.1 * df_pred) {
        /*
         * Decrease delta.
         */
        *delta *= 0.5;
      } else if (df <= 0.75 * df_pred) {
        /*
         * Increase delta.
         */
        *delta = MIN((*delta) * 2., max_step);
      } else {
        /*
         * Leave delta unchanged.
         */
        ;
      }
    }
  }

  return max_taken;
}

/*======================= end of trust_region() ==============================*/

/*======================= qr_decompose() =====================================*/

/*
 * Adapted from Numerical Recipes in C, 2nd ed., p.99.
 * Returns FALSE if normal and TRUE if singular.
 * Assumes zero-based indexing.
 */

int qr_decompose(int n, double *r, double *c, double *d) {
  int i, j, k, singular;
  double sigma, sum, tau, scale;

  singular = FALSE;

  for (k = 0; k < n - 1; k++) {
    /*
     * Put scale=0. inside k loop as in Dennis & Schnabel (1996),
     * Algorithm A3.2.1, p.305, rather than outside the loop as in
     * Numerical Recipes' qrdcmp(), p. 99.
     */
    scale = 0.;
    for (i = k; i < n; i++) {
      scale = MAX(scale, fabs(R(i, k)));
    }
    if (scale == 0.) {
      /*
       * Singular case.
       */
      c[k] = 0.;
      d[k] = 0.;
      singular = TRUE;
    } else {
      for (i = k; i < n; i++) {
        R(i, k) /= scale;
      }
      sum = 0.;
      for (i = k; i < n; i++) {
        sum += R(i, k) * R(i, k);
      }
      sigma = NR_SIGN(sqrt(sum), R(k, k));
      R(k, k) += sigma;
      c[k] = sigma * R(k, k);
      d[k] = -scale * sigma;
      for (j = k + 1; j < n; j++) {
        sum = 0.;
        for (i = k; i < n; i++) {
          sum += R(i, k) * R(i, j);
        }
        tau = sum / c[k];
        for (i = k; i < n; i++) {
          R(i, j) -= tau * R(i, k);
        }
      }
    }
  }
  d[n - 1] = R(n - 1, n - 1);

  if (d[n - 1] == 0.) {
    singular = TRUE;
  }

  return singular;
}

/*======================= end of qr_decompose() ==============================*/

/*======================= qr_update() ========================================*/

/*
 * Adapted from Numerical Recipes in C, 2nd ed., p. 101.
 * Assumes zero-based indexing.
 */

void qr_update(int n, double *r, double *qt, double *u, double *v) {
  int i, j, k;
  double tmp;

  /* Find largest k such that u[k] != 0. */
  for (k = n - 1; k > 0; k--) {
    if (u[k]) {
      break;
    }
  }

  for (i = k - 1; i >= 0; i--) {
    qr_rotate(n, r, qt, i, u[i], -u[i + 1]);
    if (u[i] == 0.) {
      u[i] = fabs(u[i + 1]);
    } else if (fabs(u[i]) > fabs(u[i + 1])) {
      tmp = u[i + 1] / u[i];
      u[i] = fabs(u[i]) * sqrt(1. + tmp * tmp);
    } else {
      tmp = u[i] / u[i + 1];
      u[i] = fabs(u[i + 1]) * sqrt(1. + tmp * tmp);
    }
  }

  for (j = 0; j < n; j++) {
    R(0, j) += u[0] * v[j];
  }
  for (i = 0; i < k; i++) {
    qr_rotate(n, r, qt, i, R(i, i), -R(i + 1, i));
  }

  return;
}

/*======================= end of qr_update() =================================*/

/*======================= qr_rotate() ========================================*/

/*
 * Adapted from Numerical Recipes in C, 2nd ed., p. 101.
 * Assumes zero-based indexing.
 */

void qr_rotate(int n, double *r, double *qt, int i, double a, double b) {
  int j;
  double c, factor, s, w, y;

  if (a == 0.) {
    c = 0.;
    s = b > 0. ? 1. : -1.;
  } else if (fabs(a) > fabs(b)) {
    factor = b / a;
    c = NR_SIGN(1. / sqrt(1. + factor * factor), a);
    s = factor * c;
  } else {
    factor = a / b;
    s = NR_SIGN(1. / sqrt(1. + factor * factor), b);
    c = factor * s;
  }

  for (j = 0; j < n; j++) {
    y = R(i, j);
    w = R(i + 1, j);
    R(i, j) = c * y - s * w;
    R(i + 1, j) = s * y + c * w;
  }

  for (j = 0; j < n; j++) {
    y = QT(i, j);
    w = QT(i + 1, j);
    QT(i, j) = c * y - s * w;
    QT(i + 1, j) = s * y + c * w;
  }

  return;
}

/*======================= end of qr_rotate() =================================*/
