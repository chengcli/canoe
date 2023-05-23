#include <stdio.h>
#include <stdlib.h>

#include "core.h"
#include "linalg.h"

#define WITHOUT_PIVOTING 0
#define WITH_PIVOTING 1

/*
 * Cubic spline routine.
 * Adapted from Numerical Recipes in C, p. 96.
 * The tridiagonal system is solved with pivoting
 * to avoid numerical instability.
 * Assumes zero-offset arrays.
 *
 * NOTE: We stripe the data into one array with float_triplet to get
 *       a cache-aware memory layout.
 */

void spline(int n, struct float_triplet *table, double y1_bot, double y1_top) {
  int j, jm1, jp1;
  double dx_a, dx_b, dx_c, *a, *b, *c, *r, *y2;
  /*
   * The following are part of DEBUG_MILESTONE(.) statements:
   */
  int idbms = 0;
  static char dbmsname[] = "spline";

  /*
   * Check validity of n.
   */
  if (n < 3) {
    fprintf(stderr, "**error:%s, n=%d < 3\n", dbmsname, n);
    exit(1);
  }

  /*
   * Allocate memory.
   */
  a = (double *)malloc(n * sizeof(double));
  b = (double *)malloc(n * sizeof(double));
  c = (double *)malloc(n * sizeof(double));
  r = (double *)malloc(n * sizeof(double));
  y2 = (double *)malloc(n * sizeof(double));

  for (j = 0; j < n; j++) {
    y2[j] = (table + j)->z;
  }

  for (j = 0; j < n; j++) {
    jm1 = j - 1;
    jp1 = j + 1;
    if (jm1 >= 0) {
      dx_a = (table + j)->x - (table + jm1)->x;
      if (dx_a <= 0.) {
        fprintf(stderr, "**error:%s, x[%d]=%g x[%d]=%g\n", dbmsname, jm1,
                (table + jm1)->x, j, (table + j)->x);
        exit(1);
      }
      if (jp1 < n) {
        dx_b = (table + jp1)->x - (table + jm1)->x;
        if (dx_b <= 0.) {
          fprintf(stderr, "**error:%s, x[%d]=%g x[%d]=%g\n", dbmsname, jm1,
                  (table + jm1)->x, jp1, (table + jp1)->x);
          exit(1);
        }
      }
    }
    if (jp1 < n) {
      dx_c = (table + jp1)->x - (table + j)->x;
      if (dx_c <= 0.) {
        fprintf(stderr, "**error:%s, x[%d]=%g x[%d]= %g\n", dbmsname, j,
                (table + j)->x, jp1, (table + jp1)->x);
        exit(1);
      }
    }

    if (j == 0) {
      if (y1_bot > 0.99e+30) {
        y2[j] = 0.;
      } else {
        b[j] = dx_c / 3.;
        c[j] = dx_c / 6.;
        r[j] = ((table + jp1)->y - (table + j)->y) / dx_c - y1_bot;
      }
    } else if (j == n - 1) {
      if (y1_top > 0.99e+30) {
        y2[j] = 0.;
      } else {
        a[j] = dx_a / 6.;
        b[j] = dx_a / 3.;
        r[j] = y1_top - ((table + j)->y - (table + jm1)->y) / dx_a;
      }
    } else {
      a[j] = dx_a / 6.;
      b[j] = dx_b / 3.;
      c[j] = dx_c / 6.;
      r[j] = ((table + jp1)->y - (table + j)->y) / dx_c -
             ((table + j)->y - (table + jm1)->y) / dx_a;
    }
  }

  if (y1_bot > 0.99e+30) {
    if (y1_top > 0.99e+30) {
      /* y2 = 0 on both ends. */
      tridiag(n - 2, a + 1, b + 1, c + 1, r + 1, y2 + 1, WITH_PIVOTING);
    } else {
      /* y2 = 0 at start. */
      tridiag(n - 1, a + 1, b + 1, c + 1, r + 1, y2 + 1, WITH_PIVOTING);
    }
  } else {
    if (y1_top > 0.99e+30) {
      /* y2 = 0 at end. */
      tridiag(n - 1, a, b, c, r, y2, WITH_PIVOTING);
    } else {
      /* y2 needed at both ends */
      tridiag(n, a, b, c, r, y2, WITH_PIVOTING);
    }
  }

  for (j = 0; j < n; j++) {
    (table + j)->z = y2[j];
  }

  /*
   * Free allocated memory.
   */
  free(y2);
  free(r);
  free(c);
  free(b);
  free(a);

  return;
}

/*======================= end of spline() ===================================*/

/*======================= splint() ==========================================*/

/*
 *  Evaluates cubic-spline interpolations.
 *  This version assumes you have already found the correct position
 *  in the tables, unlike the Numerical Recipes version.
 *  The function find_place_in_table() may be used to find the position.
 *
 *  NOTE: We stripe the data into one array with float_triplet to get
 *        a cache-aware memory layout.
 */

double splint(double xx, struct float_triplet *table, double dx) {
  double a, b, ans;

  a = ((table + 1)->x - xx) / dx;
  b = (xx - table->x) / dx;

  ans = a * table->y + b * (table + 1)->y +
        ((a * a * a - a) * (table)->z + (b * b * b - b) * (table + 1)->z) * dx *
            dx / 6;

  return ans;
}

/*======================= end of splint() ===================================*/
