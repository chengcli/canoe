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

/*======================= broyden_root() =====================================*/

/*
 * Finds a vector root, x[], of the vector function vecfunc().
 * This globally convergent algorithm is adapted from
 * Numerical Recipes in C, 2nd ed., p. 389-392, which is based on
 * Dennis and Schnabel (1996).
 *
 * Call with an initial guess for x[].
 * Assumes zero-based indexing.
 *
 * The function global_step() performs the globally-convergent step.
 * The argument step_type specifies the type of step taken. Currently, the
 * valid choices are DS_LINE_STEP and DS_DOGLEG_STEP, the latter being more
 * sophisticated and reliable (the "DS" stands for Dennis and Schnabel 1996).
 *
 * Returns 0 on normal execution and an error code if
 * the routine has failed, has converged to a local minimum, or can make
 * no further progress, in which case one should retry with a
 * different initial guess.
 */

int broyden_root(int n, double *x,
                 void (*vecfunc)(int, double *, double *, void *), double tol_f,
                 int max_it, void *arg) {
  int k, j, i, restart, singular, skip, num_bytes, old_max_taken,
      it = 0, max_taken = FALSE, count_max_taken = 0, status = DS_X_ACCEPTED;
  static int nold = 0;
  double sum, denom, f, f_old, max_step, test, h, tmp, delta = -1;
  static double *c, *d, *x_old, *fvec, *fvec2, *fvec_old, *g, *sn, *s, *t, *w,
      *qt, *r;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int idbms = 0;
  static char dbmsname[] = "broyden_root";

  if (n > nold) {
    if (nold != 0) {
      /*
       * Free previously allocated memory:
       */
      free_fvector(r, 0, nold * nold - 1, dbmsname);
      free_fvector(qt, 0, nold * nold - 1, dbmsname);
      free_fvector(w, 0, nold - 1, dbmsname);
      free_fvector(t, 0, nold - 1, dbmsname);
      free_fvector(s, 0, nold - 1, dbmsname);
      free_fvector(sn, 0, nold - 1, dbmsname);
      free_fvector(g, 0, nold - 1, dbmsname);
      free_fvector(fvec_old, 0, nold - 1, dbmsname);
      free_fvector(fvec2, 0, nold - 1, dbmsname);
      free_fvector(fvec, 0, nold - 1, dbmsname);
      free_fvector(x_old, 0, nold - 1, dbmsname);
      free_fvector(d, 0, nold - 1, dbmsname);
      free_fvector(c, 0, nold - 1, dbmsname);
    }
    /*
     * Allocate memory:
     */
    c = fvector(0, n - 1, dbmsname);
    d = fvector(0, n - 1, dbmsname);
    x_old = fvector(0, n - 1, dbmsname);
    fvec = fvector(0, n - 1, dbmsname);
    fvec2 = fvector(0, n - 1, dbmsname);
    fvec_old = fvector(0, n - 1, dbmsname);
    g = fvector(0, n - 1, dbmsname);
    sn = fvector(0, n - 1, dbmsname);
    s = fvector(0, n - 1, dbmsname);
    t = fvector(0, n - 1, dbmsname);
    w = fvector(0, n - 1, dbmsname);
    qt = fvector(0, n * n - 1, dbmsname);
    r = fvector(0, n * n - 1, dbmsname);

    nold = n;
  } else {
    /* Clear working memory: */
    num_bytes = n * sizeof(double);
    memset(c, 0, num_bytes);
    memset(d, 0, num_bytes);
    memset(x_old, 0, num_bytes);
    memset(fvec, 0, num_bytes);
    memset(fvec2, 0, num_bytes);
    memset(fvec_old, 0, num_bytes);
    memset(g, 0, num_bytes);
    memset(sn, 0, num_bytes);
    memset(s, 0, num_bytes);
    memset(t, 0, num_bytes);
    memset(w, 0, num_bytes);
    num_bytes *= n;
    memset(qt, 0, num_bytes);
    memset(r, 0, num_bytes);
  }

  /*
   * Calculate fvec[].
   */
  (*vecfunc)(n, x, fvec, arg);

  f = 0.;
  for (i = 0; i < n; i++) {
    f += fvec[i] * fvec[i];
  }
  f *= 0.5;

  /*
   * Test if initial guess is a root.
   * NOTE: We do not compare to the more stringent 0.01*tol_f used in
   *       NR and DS96 at iteration zero, so that our solutions are
   *       consistent.
   */
  test = 0.;
  for (i = 0; i < n; i++) {
    test = MAX(test, fabs(fvec[i]));
  }
  if (test < tol_f) {
    /* initial x[] is a root. */
    return status;
  }

  /*
   * Calculate max_step for globally convergent step.
   */
  sum = 0.;
  for (i = 0; i < n; i++) {
    sum += x[i] * x[i];
  }
  max_step = DS_MAX_STEP * MAX(sqrt(sum), (double)n);

  /*
   * Main iteration loop.
   */
  restart = TRUE;
  for (it = 0; it < max_it; it++) {
    if (restart == TRUE) {
      /*
       * Compute forward-difference approximation to Jacobian.
       */
      for (j = 0; j < n; j++) {
        tmp = x[j];
        h = DS_SQRT_EPS * fabs(tmp);
        if (h == 0.) {
          h = DS_SQRT_EPS;
        }
        x[j] = tmp + h;

        (*vecfunc)(n, x, fvec2, arg);

        h = x[j] - tmp;
        x[j] = tmp;
        for (i = 0; i < n; i++) {
          R(i, j) = (fvec2[i] - fvec[i]) / h;
        }
      }

      /*
       * Calculate QR decomposition.
       */
      singular = qr_decompose(n, r, c, d);
      if (singular) {
        status = DS_SINGULAR_JACOBIAN;
        return status;
      }

      /* Compute transpose, QT. */
      memset(qt, 0, n * n * sizeof(double));
      for (i = 0; i < n; i++) {
        QT(i, i) = 1.;
      }
      for (k = 0; k < n; k++) {
        if (c[k]) {
          for (j = 0; j < n; j++) {
            sum = 0.;
            for (i = k; i < n; i++) {
              sum += R(i, k) * QT(i, j);
            }
            sum /= c[k];
            for (i = k; i < n; i++) {
              QT(i, j) -= sum * R(i, k);
            }
          }
        }
      }
      /* Form R explicitly. */
      for (i = 0; i < n; i++) {
        R(i, i) = d[i];
        for (j = 0; j < i; j++) {
          R(i, j) = 0.;
        }
      }
    } else {
      for (i = 0; i < n; i++) {
        s[i] = x[i] - x_old[i];
      }
      for (i = 0; i < n; i++) {
        sum = 0.;
        for (j = i; j < n; j++) {
          sum += R(i, j) * s[j];
        }
        t[i] = sum;
      }
      skip = TRUE;
      for (i = 0; i < n; i++) {
        sum = 0.;
        for (j = 0; j < n; j++) {
          sum += QT(j, i) * t[j];
        }
        w[i] = fvec[i] - fvec_old[i] - sum;
        if (fabs(w[i]) >= DS_EPS * (fabs(fvec[i]) + fabs(fvec_old[i]))) {
          skip = FALSE;
        } else {
          w[i] = 0.;
        }
      }
      if (skip == FALSE) {
        for (i = 0; i < n; i++) {
          sum = 0.;
          for (j = 0; j < n; j++) {
            sum += QT(i, j) * w[j];
          }
          t[i] = sum;
        }

        denom = 0.;
        for (i = 0; i < n; i++) {
          denom += s[i] * s[i];
        }
        /* Store s/(s.s) in s. */
        for (i = 0; i < n; i++) {
          s[i] /= denom;
        }

        /*
         * Update r and qt.
         */
        qr_update(n, r, qt, t, s);

        for (i = 0; i < n; i++) {
          if (R(i, i) == 0.) {
            fprintf(stderr, "**error:%s, R(%d,%d) singular\n", dbmsname, i, i);
            exit(1);
          }
          d[i] = R(i, i);
        }
      }
    }

    for (i = 0; i < n; i++) {
      sum = 0.;
      for (j = 0; j < n; j++) {
        sum += QT(i, j) * fvec[j];
      }
      g[i] = sum;
    }
    for (i = n - 1; i >= 0; i--) {
      sum = 0.;
      for (j = 0; j <= i; j++) {
        sum += R(j, i) * g[j];
      }
      g[i] = sum;
    }

    /*
     * Store old x,fvec,f.
     */
    for (i = 0; i < n; i++) {
      /*
       * Screen for NaN:
       */
      if (!isfinite(x[i]) || !isfinite(fvec[i])) {
        fprintf(stderr, "**error:%s, x[%d]=%g fvec[%d]=%g\n", dbmsname, i, x[i],
                i, fvec[i]);
        exit(1);
      }
      x_old[i] = x[i];
      fvec_old[i] = fvec[i];
    }
    f_old = f;

    /*
     * Compute right-hand side of linear equations, sn[].
     */
    for (i = 0; i < n; i++) {
      sum = 0.;
      for (j = 0; j < n; j++) {
        sum += QT(i, j) * fvec[j];
      }
      sn[i] = -sum;
    }

    /*
     * Solve R.x = sn.
     * See Numerical Recipes in C, 2nd ed., p. 100.
     */
    sn[n - 1] /= d[n - 1];
    for (i = n - 2; i >= 0; i--) {
      sum = 0.;
      for (j = i + 1; j < n; j++) {
        sum += R(i, j) * sn[j];
      }
      sn[i] = (sn[i] - sum) / d[i];
    }

    /*
     * Calculate new x,f,fvec[].
     */
    old_max_taken = max_taken;
    max_taken = global_step(n, x_old, f_old, g, r, sn, max_step, &delta,
                            DS_DOGLEG_STEP, &status, x, &f, fvec, arg, vecfunc);

    /*
     * Screen for NaN:
     */
    for (i = 0; i < n; i++) {
      if (!isfinite(x[i]) || !isfinite(fvec[i])) {
        fprintf(stderr,
                "**error:%s, after global_step(): x[%d]=%g fvec[%d]=%g\n",
                dbmsname, i, x[i], i, fvec[i]);
        exit(1);
      }
    }

    /*
     * Screen for repeated maximum steps.
     */
    if (max_taken == TRUE) {
      if (old_max_taken == TRUE) {
        count_max_taken++;
      } else {
        count_max_taken = 1;
      }
    } else {
      count_max_taken = 0;
    }

    /*
     * Test for convergence.
     */
    test = 0.;
    for (i = 0; i < n; i++) {
      test = MAX(test, fabs(fvec[i]));
    }
    if (test < tol_f) {
      status = DS_X_ACCEPTED;
      return status;
    }

    if (count_max_taken >= 5) {
      status = DS_MAX_TAKEN_5;
      return status;
    }

    if (status == DS_X_NO_PROGRESS) {
      if (restart == TRUE) {
        return status;
      } else {
        test = 0.;
        denom = MAX(f, 0.5 * (double)n);
        for (i = 0; i < n; i++) {
          tmp = fabs(g[i]) * MAX(fabs(x[i]), 1.) / denom;
          test = MAX(test, tmp);
        }
        if (test < DS_TOL_MIN) {
          return status;
        } else {
          /*
           * Try reinitializing the Jacobian.
           */
          restart = TRUE;
        }
      }
    } else {
      restart = FALSE;
      test = 0.;
      for (i = 0; i < n; i++) {
        tmp = (fabs(x[i] - x_old[i])) / MAX(fabs(x[i]), 1.);
        test = MAX(test, tmp);
        if (test < DS_TOL_X) {
          /*
           * Convergence.
           */
          return status;
        }
      }
    }
  }

  status = DS_MAX_IT_EXCEEDED;
  return status;
}

/*======================= end of broyden_root() ==============================*/
