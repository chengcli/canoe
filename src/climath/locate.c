/*! Given an array xx[0..n-1] , and given a value x , 
 * returns a value j such that x is between xx[j] and xx[j+1]. 
 * xx must be monotonic, either increasing or decreasing. 
 * j=0 or j=n is returned to indicate that x is out of range.
 * adapted from Numerical Recipes in C, 2nd Ed., p. 117.
 */
int locate(double const *xx, double x, int n)
{
  xx -= 1;  // zero-offset to unit-offset

  int j;
  int ju,jm,jl;
  int ascnd;

  jl = 0;
  ju = n+1;
  ascnd = (xx[n] >= xx[1]);
  while (ju-jl > 1) {
    jm = (ju+jl) >> 1;
    if (x >= xx[jm] == ascnd)
      jl = jm;
    else
      ju = jm;
  }
  if (x == xx[1]) j = 1;
  else if(x == xx[n]) j = n-1;
  else j = jl;

  j -= 1; // unit-offset to zero-offset
  return j;
}
