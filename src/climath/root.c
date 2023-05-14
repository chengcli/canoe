#define MAX_IT 100
#define UNLIKELY_VAL -1.11111e+30

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "root.h"

int root(double x1, double x2, double xacc, double *x_root, RootFunction_t func, void *aux)
{
  int
    iter,
    compare;
  double 
    fh,fl,fm,fnew,
    s,xh,xl,xm,xnew;
  static char
    name[]="_root";

  fl = func(x1, aux);
  fh = func(x2, aux);
  if ((fl > 0. && fh < 0.) || (fl < 0. && fh > 0.)) {
    xl      = x1;
    xh      = x2;
    /* Set *x_root to an unlikely value: */
    *x_root = UNLIKELY_VAL;

    for (iter = 0; iter < MAX_IT; iter++) {
      xm = 0.5*(xl+xh);
      fm = func(xm, aux);
      s  = sqrt(fm*fm-fl*fh);
      if (s == 0.) {
        return 0;
      }
      xnew = xm+(xm-xl)*((fl > fh ? 1. : -1.)*fm/s);

      if (fabs(xnew-*x_root) <= xacc) {
        return 0;
      }
      *x_root = xnew;

      fnew    = func(*x_root, aux);
      if (fnew == 0.) {
        return 0;
      }

      if ((fnew > 0. ? fabs(fm) : -fabs(fm)) != fm) {
        xl = xm;
        fl = fm;
        xh = *x_root;
        fh = fnew;
      }
      else if ((fnew > 0. ? fabs(fl) : -fabs(fl)) != fl) {
        xh = *x_root;
        fh = fnew;
      }
      else if ((fnew > 0. ? fabs(fh) : -fabs(fh)) != fh) {
        xl = *x_root;
        fl = fnew;
      }
      else {
        fprintf(stderr, "### FATAL ERROR in _root function: should never get here\n");
        exit(1);
      }
      if (fabs(xh-xl) <= xacc) {
        return 0;
      }
    }
    fprintf(stderr, "### FATAL ERROR in _root function: exceeded MAX_IT = ");
    fprintf(stderr, "%d\n", MAX_IT);
    fprintf(stderr, "current root calc = ");
    fprintf(stderr, "%12.4g\n", *x_root);
    exit(1);
  }
  else {
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
    }
    else {
      return 1;
    }
  }

  /* Should never get here. */
  fprintf(stderr, "### FATAL ERROR in _root function: should never get here\n");
  exit(1);
}

#undef MAX_IT
#undef UNLIKELY_VAL
