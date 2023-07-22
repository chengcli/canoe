#ifndef SRC_CLIMATH_ROOT_HPP_
#define SRC_CLIMATH_ROOT_HPP_

#include <cmath>
#include <cstdio>
#include <cstdlib>

#undef MAX_IT
#define MAX_IT 100

#undef UNLIKELY_VAL
#define UNLIKELY_VAL -1.11111e+30

template <typename FLOAT, typename FUNC>
int root(FLOAT x1, FLOAT x2, FLOAT xacc, FLOAT *x_root, FUNC func) {
  int iter, compare;
  FLOAT
  fh, fl, fm, fnew, s, xh, xl, xm, xnew;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
  int
    idbms=0;
   */
  static char dbmsname[] = "_FindRoot";

  fl = func(x1);
  fh = func(x2);
  if ((fl > 0. && fh < 0.) || (fl < 0. && fh > 0.)) {
    xl = x1;
    xh = x2;
    /* Set *x_root to an unlikely value: */
    *x_root = UNLIKELY_VAL;

    for (iter = 0; iter < MAX_IT; iter++) {
      xm = 0.5 * (xl + xh);
      fm = func(xm);
      s = sqrt(fm * fm - fl * fh);
      if (s == 0.) {
        return 0;
      }
      xnew = xm + (xm - xl) * ((fl > fh ? 1. : -1.) * fm / s);

      if (fabs(xnew - *x_root) <= xacc) {
        return 0;
      }
      *x_root = xnew;

      fnew = func(*x_root);
      if (fnew == 0.) {
        return 0;
      }

      if ((fnew > 0. ? fabs(fm) : -fabs(fm)) != fm) {
        xl = xm;
        fl = fm;
        xh = *x_root;
        fh = fnew;
      } else if ((fnew > 0. ? fabs(fl) : -fabs(fl)) != fl) {
        xh = *x_root;
        fh = fnew;
      } else if ((fnew > 0. ? fabs(fh) : -fabs(fh)) != fh) {
        xl = *x_root;
        fl = fnew;
      } else {
        fprintf(stderr, "**error:%s, should never get here\n", dbmsname);
        exit(1);
      }
      if (fabs(xh - xl) <= xacc) {
        return 0;
      }
    }
    fprintf(stderr,
            "**error:%s, exceeded MAX_IT = %d, current root calc = %e\n",
            dbmsname, MAX_IT, *x_root);
    exit(1);
  } else {
    if (fl == 0.) {
      *x_root = x1;
      return 0;
    }
    if (fh == 0.) {
      *x_root = x2;
      return 0;
    }

    compare = fabs(fl) >= fabs(fh) ? 1 : -1;
    if (compare < 0) {
      return -1;
    } else {
      return 1;
    }
  }

  /* Should never get here. */
  fprintf(stderr, "**error:%s, should never reach here\n", dbmsname);
  exit(1);
}

#undef MAX_IT
#undef UNLIKELY_VAL

#endif  // SRC_CLIMATH_ROOT_HPP_
