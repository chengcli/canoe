#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*! Given a matrix a[0..n-1][0..n-1], this routine replaces it by the LU
 * decomposition of a rowwise permutation of itself. a and n are input. a is
 * output, arranged as in NRIC equation (2.3.14) ; indx[0..n-1] is an output
 * vector that records the row permutation effected by the partial pivoting; d
 * is output as +/- 1 depending on whether the number of row interchanges was
 * evenor odd, respectively. This routine is used in combination with lubksbto
 * solve linear equationsor invert a matrix. adapted from Numerical Recipes in
 * C, 2nd Ed., p. 46.
 */
int ludcmp(double **a, int n, int *indx) {
  int i, imax, j, k, d;
  double big, dum, sum, temp;
  double *vv = (double *)malloc(n * sizeof(double));

  d = 1;
  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if ((temp = fabs(a[i][j])) > big) big = temp;
    if (big == 0.0) {
      fprintf(stderr, "Singular matrix in routine ludcmp");
      exit(1);
    }
    vv[i] = 1.0 / big;
  }
  for (j = 0; j < n; j++) {
    for (i = 0; i < j; i++) {
      sum = a[i][j];
      for (k = 0; k < i; k++) sum -= a[i][k] * a[k][j];
      a[i][j] = sum;
    }
    big = 0.0;
    for (i = j; i < n; i++) {
      sum = a[i][j];
      for (k = 0; k < j; k++) sum -= a[i][k] * a[k][j];
      a[i][j] = sum;
      if ((dum = vv[i] * fabs(sum)) >= big) {
        big = dum;
        imax = i;
      }
    }
    if (j != imax) {
      for (k = 0; k < n; k++) {
        dum = a[imax][k];
        a[imax][k] = a[j][k];
        a[j][k] = dum;
      }
      d = -d;
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if (j != n - 1) {
      dum = (1.0 / a[j][j]);
      for (i = j + 1; i < n; i++) a[i][j] *= dum;
    }
  }
  free(vv);

  return d;
}
