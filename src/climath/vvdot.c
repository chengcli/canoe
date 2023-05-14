/*! vector vector dot product: a.b
 * a[0..n-1] is input
 * b[0..n-1] is input
 */
double vvdot(double const *a, double const *b, int n)
{
  double result = 0.;
  #pragma GCC ivdep
  for (int i = 0; i < n; ++i)
    result += a[i]*b[i];
  return result;
}
