#include <math.h>

#include "core.h"

/*
 * Find place in table using bisection.
 * Adapted from EPIC model 430
 *
 * NOTE: Finds place relative to the x component of the float_triplet table.
 */
int find_place_in_table(int n, struct float_triplet *table, double x,
                        double *dx, int il) {
  int im, iu, ascend, cmplo, cmphi;

  iu = n;
  ascend = (table[n - 1].x >= table[0].x);
  cmplo = fcmp(x, table[0].x);
  cmphi = fcmp(x, table[n - 1].x);

  /*
   * Check x against table endpoints.
   */
  if (ascend) {
    if (cmplo <= 0) {
      il = 0;
      *dx = table[il + 1].x - table[il].x;
      return il;
    }
    if (cmphi >= 0) {
      il = n - 2;
      *dx = table[il + 1].x - table[il].x;
      return il;
    }
  } else {
    if (cmplo >= 0) {
      il = 0;
      *dx = table[il + 1].x - table[il].x;
      return il;
    }
    if (cmphi <= 0) {
      il = n - 2;
      *dx = table[il + 1].x - table[il].x;
      return il;
    }
  }

  /*
   * Use bisection to search table.
   */
  while (iu - il > 1) {
    im = (iu + il) >> 1;
    if (fcmp(x, table[im].x) >= 0 == ascend) {
      il = im;
    } else {
      iu = im;
    }
  }

  *dx = table[il + 1].x - table[il].x;
  return il;
}
