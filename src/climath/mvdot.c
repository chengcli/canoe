#include <assert.h>
/*! Matrix vector dot product: m.v
 * m[0..n1-1][0..n2-1] is input
 * v[0..n2-1] in input
 * r[0..n1-1] is output
 */
void mvdot(double *r, double **m, double const *v, int n1, int n2)
{
  assert(r != v); // r and v cannot be the same
  for (int i = 0; i < n1; ++i) {
    r[i] = 0.;
    #pragma GCC ivdep
    for (int j = 0; j < n2; ++j)
      r[i] += m[i][j]*v[j];
  }
}
