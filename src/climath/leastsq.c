#include <stdlib.h>
#include <string.h>

#include "linalg.h"

/*! solve least square problem A.x = b
 * A[0..n1-1][0..n2-1] is input
 * b[0..n1-1] is input/output. Input dimension is n1, output dimension is n2
 * n1 >= n2
 */
void leastsq(double **A, double *b, int n1, int n2) {
  double *c = (double *)malloc(n1 * sizeof(double));
  memcpy(c, b, n1 * sizeof(double));

  double **B;
  B = (double **)malloc(n2 * sizeof(double *));
  B[0] = (double *)malloc(n2 * n2 * sizeof(double));
  for (int i = 0; i < n2; ++i) B[i] = B[0] + i * n2;

  for (int i = 0; i < n2; ++i) {
    // calculate A^T.A
    for (int j = 0; j < n2; ++j) {
      B[i][j] = 0.;
#pragma GCC ivdep
      for (int k = 0; k < n1; ++k) B[i][j] += A[k][i] * A[k][j];
    }

    // calculate A^T.b
    b[i] = 0.;
#pragma GCC ivdep
    for (int j = 0; j < n1; ++j) b[i] += A[j][i] * c[j];
  }

  // calculate (A^T.A)^{-1}.(A^T.b)
  int *indx = (int *)malloc(n2 * sizeof(int));
  ludcmp(B, n2, indx);
  lubksb(B, n2, indx, b);

  free(c);
  free(indx);
  free(B[0]);
  free(B);
}
